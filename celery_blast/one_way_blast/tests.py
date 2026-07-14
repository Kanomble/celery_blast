import os
from unittest.mock import ANY, patch

from django.conf import settings
from django.contrib.auth.models import User
from django.test import SimpleTestCase, TestCase
from django.test.utils import override_settings
from django.urls import reverse
from django_celery_results.models import TaskResult

from blast_project.models import BlastSettings
from one_way_blast.models import OneWayRemoteBlastProject
from one_way_blast.tasks import execute_snakemake_workflow


class FakeProgressRecorder:
    def __init__(self):
        self.progress_updates = []

    def set_progress(self, current, total, description):
        self.progress_updates.append((current, total, description))


def expected_snakemake_environment():
    env = os.environ.copy()
    pythonpath = env.get('PYTHONPATH')
    if pythonpath:
        env['PYTHONPATH'] = settings.BASE_DIR + os.pathsep + pythonpath
    else:
        env['PYTHONPATH'] = settings.BASE_DIR
    return env


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
                '--snakefile', os.path.abspath(
                    os.path.join(settings.BASE_DIR, 'static/snakefiles/one_way_blast/Snakefile')
                ),
                '--cores', '1',
                '--configfile', os.path.abspath(
                    os.path.join(settings.BASE_DIR, 'media/one_way_blast/7/snakefile_config')
                ),
                '--directory', os.path.abspath(os.path.join(settings.BASE_DIR, 'media/one_way_blast/7/')),
                '--keep-incomplete',
            ],
            timeout=45,
            shell=False,
            logger=ANY,
            check=True,
            env=expected_snakemake_environment(),
        )


class OneWayRemotePollingViewTests(TestCase):
    def setUp(self):
        self.user = User.objects.create_user(
            'testuser',
            'test_email@email.com',
            password='test',
        )
        self.settings = BlastSettings.objects.create(
            e_value=0.001,
            word_size=3,
            num_alignments=10000,
            max_target_seqs=10000,
            num_threads=1,
            max_hsps=500,
        )
        self.project = OneWayRemoteBlastProject.objects.create(
            r_project_title='remote one-way polling project',
            r_project_query_sequences='query.faa',
            r_project_user=self.user,
            r_project_settings=self.settings,
            r_project_database='nr',
            r_search_strategy='blastp',
            r_project_execution_task_result=TaskResult.objects.create(
                task_id='remote-one-way-progress',
                status='PROGRESS',
            ),
        )

    def test_remote_detail_page_renders_expected_polling_configuration(self):
        self.client.login(username='testuser', password='test')

        response = self.client.get(
            reverse('one_way_remote_project_details', kwargs={'project_id': self.project.id})
        )

        self.assertEqual(200, response.status_code)
        self.assertContains(response, 'id="progress_bar"', html=False)
        self.assertContains(
            response,
            reverse(
                'ajax_call_to_snakemake_logfiles',
                kwargs={'remote': 1, 'project_id': self.project.id},
            ),
            html=False,
        )
        self.assertContains(response, "type: 'get'", html=False)
        self.assertContains(response, '}, 20000);', html=False)
        self.assertContains(response, 'error: function(data)', html=False)
        self.assertContains(response, 'data.responseJSON.error', html=False)
        self.assertContains(response, 'error polling progress data', html=False)
        self.assertContains(response, 'clearInterval(pollInterval);', html=False)

    def test_remote_polling_endpoint_returns_backend_json_error(self):
        response = self.client.get(
            reverse(
                'ajax_call_to_snakemake_logfiles',
                kwargs={'remote': 1, 'project_id': 999999},
            ),
            HTTP_X_REQUESTED_WITH='XMLHttpRequest',
        )

        self.assertEqual(400, response.status_code)
        self.assertIn('error', response.json())

    @override_settings(SUBPROCESS_TIME_LIMIT=45)
    def test_execute_one_way_remote_snakemake_workflow_builds_absolute_command(self):
        progress_recorder = FakeProgressRecorder()

        class Result:
            returncode = 0

        with patch('one_way_blast.tasks.run_external_command', return_value=Result()) as run_command:
            returncode = execute_snakemake_workflow(
                project_id=9,
                working_dir='media/one_way_blast/remote_searches/9/',
                config_file='media/one_way_blast/remote_searches/9/snakefile_config',
                snakefile_dir='static/snakefiles/one_way_blast/remote_searches/Snakefile',
                task_result_updater=lambda project_id, task_id: None,
                task_id='task-remote',
                progress_recorder=progress_recorder,
            )

        self.assertEqual(0, returncode)
        self.assertEqual([(100, 100, 'SUCCESS')], progress_recorder.progress_updates)
        run_command.assert_called_once_with(
            [
                'snakemake',
                '--snakefile', os.path.abspath(
                    os.path.join(settings.BASE_DIR, 'static/snakefiles/one_way_blast/remote_searches/Snakefile')
                ),
                '--cores', '1',
                '--configfile', os.path.abspath(
                    os.path.join(settings.BASE_DIR, 'media/one_way_blast/remote_searches/9/snakefile_config')
                ),
                '--directory', os.path.abspath(
                    os.path.join(settings.BASE_DIR, 'media/one_way_blast/remote_searches/9/')
                ),
                '--keep-incomplete',
            ],
            timeout=45,
            shell=False,
            logger=ANY,
            check=True,
            env=expected_snakemake_environment(),
        )
