import os
import sys
import tempfile
from pathlib import Path

import yaml
from django.conf import settings
from django.test import SimpleTestCase

from celery_blast.entrez_query import (
    EntrezQueryValidationError,
    build_remote_blast_command,
    run_remote_blast_command,
    update_snakemake_config_entrez_query,
    validate_entrez_query,
)


HARMLESS_PAYLOADS = [
    'txid2[Organism] AND "ribosomal protein"',
    "Escherichia coli[organism] AND 'DNA polymerase'",
    'bacteria[organism]; touch SHOULD_NOT_EXIST',
    'bacteria[organism] $(touch SHOULD_NOT_EXIST)',
    'bacteria[organism] `touch SHOULD_NOT_EXIST`',
    'bacteria[organism] | touch SHOULD_NOT_EXIST',
    'bacteria[organism] > SHOULD_NOT_EXIST < input',
    'β-proteobacteria[organism] AND kinase[Title]',
    '((txid2[Organism:exp]) AND ("outer membrane"[Title] OR porin[Title])) NOT plasmid[Title]',
]


class EntrezQuerySafetyTests(SimpleTestCase):
    def test_valid_queries_round_trip_through_yaml_without_value_changes(self):
        for query in HARMLESS_PAYLOADS:
            with self.subTest(query=query):
                with tempfile.TemporaryDirectory() as temp_dir:
                    config_path = os.path.join(temp_dir, "snakefile_config")
                    with open(config_path, "w", encoding="utf-8") as config_file:
                        config_file.write("project_id: 1\nentrez_query:\nquery_sequence: query.faa\n")

                    update_snakemake_config_entrez_query(config_path, query)

                    with open(config_path, encoding="utf-8") as config_file:
                        config = yaml.safe_load(config_file)
                    self.assertEqual(query, config["entrez_query"])

    def test_invalid_control_characters_raise_validation_error(self):
        for query in ("abc\x00def", "abc\ndef", "abc\rdef", "abc\tdef"):
            with self.subTest(query=repr(query)):
                with self.assertRaises(EntrezQueryValidationError):
                    validate_entrez_query(query)

    def test_query_is_one_logical_blast_argument(self):
        query = 'bacteria[organism]"; touch SHOULD_NOT_EXIST; echo "$(id)" | cat > file'

        command = build_remote_blast_command(
            search_strategy="blastp",
            database="nr",
            outfmt="6 qseqid sseqid",
            output_path="blast.table",
            word_size=3,
            e_value=0.001,
            num_alignments=100,
            query_path="query.faa",
            entrez_query=query,
        )

        self.assertIn("-entrez_query", command)
        query_index = command.index("-entrez_query") + 1
        self.assertEqual(query, command[query_index])
        self.assertEqual(1, command.count(query))

    def test_runner_does_not_execute_secondary_command_from_argument_text(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            marker_path = os.path.join(temp_dir, "marker")
            stderr_path = os.path.join(temp_dir, "stderr.log")
            payload = 'query"; touch "{}"; echo "'.format(marker_path)
            command = [
                sys.executable,
                "-c",
                "import sys; assert len(sys.argv) == 2",
                payload,
            ]

            run_remote_blast_command(command, stderr_path)

            self.assertFalse(os.path.exists(marker_path))

    def test_remote_snakefiles_do_not_shell_execute_entrez_command_strings(self):
        snakefiles = [
            Path(settings.BASE_DIR) / "static/snakefiles/reciprocal_blast/remote_searches/Snakefile",
            Path(settings.BASE_DIR) / "static/snakefiles/one_way_blast/remote_searches/Snakefile",
        ]

        for snakefile in snakefiles:
            with self.subTest(snakefile=snakefile):
                text = snakefile.read_text(encoding="utf-8")
                self.assertNotIn("shell(cmd_string)", text)
                self.assertNotIn(' -entrez_query \\"', text)
                self.assertIn("build_remote_blast_command", text)
                self.assertIn("run_remote_blast_command", text)
