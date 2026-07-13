from unittest.mock import ANY, patch

from django.test import SimpleTestCase
from django.test.utils import override_settings

from one_way_blast.tasks import execute_snakemake_workflow


class FakeProgressRecorder:
    def __init__(self):
        self.progress_updates = []

    def set_progress(self, current, total, description):
        self.progress_updates.append((current, total, description))


class OneWayTaskHelperTests(SimpleTestCase):
    @override_settings(SUBPROCESS_TIME_LIMIT=45)
    def test_execute_snakemake_workflow_builds_command_and_updates_task(self):
        updated = []
        progress_recorder = FakeProgressRecorder()

        def update_project(project_id, task_id):
            updated.append((project_id, task_id))

        class Result:
            returncode = 0

        with patch('one_way_blast.tasks.run_external_command', return_value=Result()) as run_command:
            returncode = execute_snakemake_workflow(
                project_id=7,
                working_dir='media/one_way_blast/7/',
                config_file='media/one_way_blast/7/snakefile_config',
                snakefile_dir='static/snakefiles/one_way_blast/Snakefile',
                task_result_updater=update_project,
                task_id='task-123',
                progress_recorder=progress_recorder,
            )

        self.assertEqual(0, returncode)
        self.assertEqual([(7, 'task-123')], updated)
        self.assertEqual([(100, 100, 'SUCCESS')], progress_recorder.progress_updates)
        run_command.assert_called_once_with(
            [
                'snakemake',
                '--snakefile', 'static/snakefiles/one_way_blast/Snakefile',
                '--cores', '1',
                '--configfile', 'media/one_way_blast/7/snakefile_config',
                '--directory', 'media/one_way_blast/7/',
                '--keep-incomplete',
            ],
            timeout=45,
            shell=False,
            logger=ANY,
            check=True,
        )
