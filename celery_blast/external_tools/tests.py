from unittest.mock import ANY, patch

from django.test import SimpleTestCase

from external_tools.genbank_download_clinker_synteny import (
    build_download_genbank_command,
    run_genbank_download_command,
)
from external_tools.tasks import (
    build_clinker_command,
    build_cdd_download_command,
    build_cdd_extract_command,
    build_fasttree_command,
    build_mafft_command,
    build_rpsblast_command,
    build_shiptv_command,
    run_external_tool_command,
)
from external_tools.entrez_search_service import (
    build_efetch_fasta_from_input_command,
    build_entrez_search_command,
    build_protein_query_to_fasta_command,
    build_pubmed_to_fasta_command,
    run_edirect_command,
)


class EntrezSearchServiceCommandTests(SimpleTestCase):
    def test_build_entrez_search_command_passes_query_as_argument(self):
        command = build_entrez_search_command(
            'protein',
            50,
            'kinase "quoted"',
            '/tmp/search.tsv',
        )

        self.assertEqual('bash', command[0])
        self.assertEqual('-o', command[1])
        self.assertEqual('pipefail', command[2])
        self.assertIn('esearch -db "$1" -query "$2"', command[4])
        self.assertEqual('protein', command[6])
        self.assertEqual('kinase "quoted"', command[7])
        self.assertEqual('50', command[8])
        self.assertEqual('Id Caption Title Organism', command[9])
        self.assertEqual('/tmp/search.tsv', command[10])

    def test_build_efetch_fasta_from_input_command_passes_paths_as_arguments(self):
        command = build_efetch_fasta_from_input_command('/tmp/ids.txt', '/tmp/out.faa')

        self.assertEqual(
            [
                'bash',
                '-c',
                'efetch -db protein -format fasta -input "$1" > "$2"',
                'efetch_fasta',
                '/tmp/ids.txt',
                '/tmp/out.faa',
            ],
            command,
        )

    def test_build_pubmed_to_fasta_command_preserves_active_download_shape(self):
        command = build_pubmed_to_fasta_command('hydra symbiont', '/tmp/out.faa')

        self.assertIn('elink -target protein | efetch -format fasta -start 1 -stop 100', command[4])
        self.assertEqual('hydra symbiont', command[6])
        self.assertEqual('/tmp/out.faa', command[7])

    def test_build_pubmed_to_fasta_command_can_use_refseq_filter(self):
        command = build_pubmed_to_fasta_command('hydra symbiont', '/tmp/out.faa', refseq_only=True)

        self.assertIn('elink -target protein | efilter -source refseq | efetch -format fasta', command[4])
        self.assertNotIn('-start 1 -stop 100', command[4])

    def test_build_protein_query_to_fasta_command(self):
        command = build_protein_query_to_fasta_command('WP_12345', '/tmp/out.faa')

        self.assertEqual('WP_12345', command[6])
        self.assertEqual('/tmp/out.faa', command[7])
        self.assertIn('esearch -db protein -query "$1"', command[4])

    def test_run_edirect_command_delegates_to_runner_without_checking_returncode(self):
        class Result:
            returncode = 2

        command = ['bash', '-c', 'exit 2']
        with patch('external_tools.entrez_search_service.run_external_command', return_value=Result()) as run_command:
            returncode = run_edirect_command(command, timeout=30)

        self.assertEqual(2, returncode)
        run_command.assert_called_once_with(command, timeout=30, shell=False, check=False)


class ExternalToolCommandTests(SimpleTestCase):
    def test_build_mafft_command_passes_paths_as_arguments(self):
        command = build_mafft_command('/tmp/input with spaces.faa', '/tmp/output.msa')

        self.assertEqual(
            [
                'bash',
                '-c',
                'mafft "$1" > "$2"',
                'mafft_alignment',
                '/tmp/input with spaces.faa',
                '/tmp/output.msa',
            ],
            command,
        )

    def test_build_fasttree_command_passes_paths_as_arguments(self):
        command = build_fasttree_command('/tmp/input.msa', '/tmp/output tree.nwk')

        self.assertEqual(
            [
                'bash',
                '-c',
                'fasttree -lg "$1" > "$2"',
                'fasttree_phylogeny',
                '/tmp/input.msa',
                '/tmp/output tree.nwk',
            ],
            command,
        )

    def test_build_shiptv_command_uses_argv_without_shell(self):
        command = build_shiptv_command('/tmp/tree.nwk', '/tmp/metadata.tsv', '/tmp/tree.html')

        self.assertEqual(
            [
                'shiptv',
                '--newick', '/tmp/tree.nwk',
                '--metadata', '/tmp/metadata.tsv',
                '--output-html', '/tmp/tree.html',
            ],
            command,
        )

    def test_build_clinker_command_expands_to_argv(self):
        command = build_clinker_command(['/tmp/a.gbk', '/tmp/b.gbk'], '/tmp/clinker.html')

        self.assertEqual(
            [
                'clinker',
                '/tmp/a.gbk',
                '/tmp/b.gbk',
                '-p',
                '/tmp/clinker.html',
                '-i',
                '0.25',
            ],
            command,
        )

    def test_build_cdd_download_command(self):
        command = build_cdd_download_command('https://example.test/Cdd_LE.tar.gz', '/tmp/Cdd_LE.tar.gz')

        self.assertEqual(
            ['wget', 'https://example.test/Cdd_LE.tar.gz', '-q', '-O', '/tmp/Cdd_LE.tar.gz'],
            command,
        )

    def test_build_cdd_extract_command(self):
        command = build_cdd_extract_command('/tmp/Cdd_LE.tar.gz', '/tmp/CDD')

        self.assertEqual(['tar', '-zxvf', '/tmp/Cdd_LE.tar.gz', '-C', '/tmp/CDD/'], command)

    def test_build_rpsblast_command(self):
        command = build_rpsblast_command(
            '/tmp/query.faa',
            '/tmp/cdd',
            '/tmp/domains.tsv',
            {
                'rps_e_value': 0.001,
                'rps_num_threads': 4,
                'rps_max_hsps': 2,
                'rps_num_alignments': 100,
            },
        )

        self.assertEqual('rpsblast', command[0])
        self.assertEqual('/tmp/query.faa', command[2])
        self.assertEqual('/tmp/cdd', command[4])
        self.assertEqual('/tmp/domains.tsv', command[8])
        self.assertEqual('0.001', command[10])
        self.assertEqual('4', command[12])
        self.assertEqual('2', command[14])
        self.assertEqual('100', command[16])

    def test_run_external_tool_command_delegates_to_runner_with_check_enabled(self):
        class Result:
            returncode = 0

        command = ['shiptv', '--version']
        with patch('external_tools.tasks.run_external_command', return_value=Result()) as run_command:
            returncode = run_external_tool_command(command, timeout=45)

        self.assertEqual(0, returncode)
        run_command.assert_called_once_with(
            command,
            timeout=45,
            shell=False,
            logger=ANY,
            check=True,
            cleanup_exceptions=(),
        )

    def test_run_external_tool_command_delegates_cleanup_exceptions(self):
        class Result:
            returncode = 0

        cleanup_exceptions = (RuntimeError,)
        command = ['rpsblast', '-query', '/tmp/query.faa']
        with patch('external_tools.tasks.run_external_command', return_value=Result()) as run_command:
            returncode = run_external_tool_command(command, timeout=45, cleanup_exceptions=cleanup_exceptions)

        self.assertEqual(0, returncode)
        self.assertEqual(cleanup_exceptions, run_command.call_args.kwargs['cleanup_exceptions'])


class GenbankDownloadCommandTests(SimpleTestCase):
    def test_build_download_genbank_command_passes_paths_as_arguments(self):
        command = build_download_genbank_command(
            'https://example.test/file with spaces.gbk.gz',
            '/tmp/output file.gbk',
        )

        self.assertEqual('bash', command[0])
        self.assertEqual('-o', command[1])
        self.assertEqual('pipefail', command[2])
        self.assertEqual('wget -qO- "$1" | gzip -d > "$2"', command[4])
        self.assertEqual('https://example.test/file with spaces.gbk.gz', command[6])
        self.assertEqual('/tmp/output file.gbk', command[7])

    def test_run_genbank_download_command_delegates_to_runner_with_check_enabled(self):
        class Result:
            returncode = 0

        command = ['bash', '-c', 'true']
        with patch(
            'external_tools.genbank_download_clinker_synteny.run_external_command',
            return_value=Result(),
        ) as run_command:
            returncode = run_genbank_download_command(command, timeout=25)

        self.assertEqual(0, returncode)
        run_command.assert_called_once_with(command, timeout=25, shell=False, check=True)
