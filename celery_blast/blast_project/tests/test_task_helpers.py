from subprocess import SubprocessError, TimeoutExpired
from unittest.mock import ANY, patch

from django.test import SimpleTestCase
from django.test.utils import override_settings

from celery_blast.processes import ExternalCommandError, ExternalCommandTimeout
from blast_project.tasks import (
    build_makeblastdb_command,
    build_species_taxids_command,
    execute_reciprocal_snakemake_workflow,
    run_database_setup_command,
    run_makeblastdb_command,
    run_species_taxids_command,
)


class FakeProgressRecorder:
    def __init__(self):
        self.progress_updates = []

    def set_progress(self, current, total, description):
        self.progress_updates.append((current, total, description))


class FakeExternalTools:
    def __init__(self):
        self.msa_task_id = None
        self.phylo_task_id = None

    def update_for_all_query_sequences_msa_task(self, task_id):
        self.msa_task_id = task_id

    def update_for_all_query_sequences_phylo_task(self, task_id):
        self.phylo_task_id = task_id


class BlastProjectTaskHelperTests(SimpleTestCase):
    def test_run_database_setup_command_delegates_to_runner(self):
        class Result:
            returncode = 0

        with patch('blast_project.tasks.run_external_command', return_value=Result()) as run_command:
            returncode = run_database_setup_command(['wget', 'url'], 600)

        self.assertEqual(0, returncode)
        run_command.assert_called_once_with(
            ['wget', 'url'],
            timeout=600,
            shell=False,
            logger=ANY,
            check=True,
        )

    def test_run_database_setup_command_maps_timeout(self):
        with patch('blast_project.tasks.run_external_command',
                   side_effect=ExternalCommandTimeout(['wget'], message='timeout')):
            with self.assertRaises(TimeoutExpired):
                run_database_setup_command(['wget'], 600)

    def test_run_database_setup_command_maps_command_error(self):
        with patch('blast_project.tasks.run_external_command',
                   side_effect=ExternalCommandError(['wget'], returncode=1)):
            with self.assertRaises(SubprocessError):
                run_database_setup_command(['wget'], 600)

    def test_build_makeblastdb_command_with_taxmap_file(self):
        with patch('blast_project.tasks.BLAST_DATABASE_DIR', 'media/databases/'):
            command = build_makeblastdb_command(4, 'media/databases/4/TEST.database', taxmap_file=True)

        self.assertEqual(
            [
                'makeblastdb',
                '-in', 'media/databases/4/TEST.database',
                '-out', 'media/databases/4/TEST.database',
                '-taxid_map', 'media/databases/4/acc_taxmap.table',
                '-dbtype', 'prot',
                '-input_type', 'fasta',
                '-parse_seqids',
            ],
            command,
        )

    def test_build_makeblastdb_command_with_taxonomic_node(self):
        with patch('blast_project.tasks.BLAST_DATABASE_DIR', 'media/databases/'):
            command = build_makeblastdb_command(4, 'media/databases/4/TEST.database', taxonomic_node=1140)

        self.assertEqual(
            [
                'makeblastdb',
                '-in', 'media/databases/4/TEST.database',
                '-out', 'media/databases/4/TEST.database',
                '-taxid', '1140',
                '-dbtype', 'prot',
                '-input_type', 'fasta',
                '-parse_seqids',
            ],
            command,
        )

    def test_build_makeblastdb_command_requires_taxmap_or_taxid(self):
        with self.assertRaises(SubprocessError):
            build_makeblastdb_command(4, 'media/databases/4/TEST.database')

    def test_build_species_taxids_command_passes_args_safely(self):
        command = build_species_taxids_command(1140, '/data/taxids.txt')

        self.assertEqual(
            [
                'bash',
                '-c',
                'get_species_taxids.sh -t "$1" >> "$2" 2>&1',
                'get_species_taxids',
                '1140',
                '/data/taxids.txt',
            ],
            command,
        )

    def test_run_species_taxids_command_delegates_to_runner_without_checking_returncode(self):
        class Result:
            returncode = 2

        command = build_species_taxids_command(1140, '/data/taxids.txt')
        with patch('blast_project.tasks.run_external_command', return_value=Result()) as run_command:
            returncode = run_species_taxids_command(command)

        self.assertEqual(2, returncode)
        run_command.assert_called_once_with(
            command,
            timeout=200,
            shell=False,
            logger=ANY,
            check=False,
        )

    @override_settings(SUBPROCESS_TIME_LIMIT=120)
    def test_run_makeblastdb_command_delegates_to_runner(self):
        class Result:
            returncode = 0

        command = ['makeblastdb', '-in', 'db.faa']
        with patch('blast_project.tasks.run_external_command', return_value=Result()) as run_command:
            returncode = run_makeblastdb_command(command)

        self.assertEqual(0, returncode)
        run_command.assert_called_once_with(
            command,
            timeout=120,
            shell=False,
            logger=ANY,
            check=True,
            cleanup_exceptions=ANY,
        )

    @override_settings(SUBPROCESS_TIME_LIMIT=90)
    def test_execute_reciprocal_snakemake_workflow_builds_command_and_updates_models(self):
        class Result:
            returncode = 0

        updated = []
        progress_recorder = FakeProgressRecorder()
        external_tools = FakeExternalTools()

        def update_project(project_id, task_id):
            updated.append((project_id, task_id))

        with patch('blast_project.tasks.run_external_command', return_value=Result()) as run_command, \
                patch('blast_project.tasks.create_external_tools_after_snakemake_workflow_finishes') as create_tools, \
                patch('blast_project.tasks.ExternalTools') as external_tools_model:
            external_tools_model.objects.get_external_tools_based_on_project_id.return_value = external_tools

            returncode = execute_reciprocal_snakemake_workflow(
                project_id=3,
                working_dir='media/blast_projects/3/',
                config_file='media/blast_projects/3/snakefile_config',
                snakefile_dir='static/snakefiles/reciprocal_blast/Snakefile',
                task_result_updater=update_project,
                task_id='task-456',
                project_type='local',
                progress_recorder=progress_recorder,
            )

        self.assertEqual(0, returncode)
        self.assertEqual([(3, 'task-456')], updated)
        self.assertEqual('task-456', external_tools.msa_task_id)
        self.assertEqual('task-456', external_tools.phylo_task_id)
        self.assertEqual([(25, 100, 'PROGRESS'), (100, 100, 'SUCCESS')], progress_recorder.progress_updates)
        create_tools.assert_called_once_with(3, 'local')
        run_command.assert_called_once_with(
            [
                'snakemake',
                '--snakefile', 'static/snakefiles/reciprocal_blast/Snakefile',
                '--cores', '1',
                '--configfile', 'media/blast_projects/3/snakefile_config',
                '--directory', 'media/blast_projects/3/',
            ],
            timeout=90,
            shell=False,
            logger=ANY,
            check=True,
            cleanup_exceptions=ANY,
        )
