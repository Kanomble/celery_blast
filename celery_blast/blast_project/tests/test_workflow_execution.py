from concurrent.futures import ThreadPoolExecutor

from django.contrib.auth.models import User
from django.db import close_old_connections
from django.test import TransactionTestCase
from django_celery_results.models import TaskResult

from blast_project.models import BlastProject, BlastSettings, RemoteBlastProject, WorkflowLifecycle
from blast_project.workflow_execution import (
    claim_workflow_execution,
    finish_workflow_execution,
    lifecycle_from_task_status,
    start_reciprocal_workflow,
)
from celery_blast.resource_governance import WorkflowQuotaExceeded
from refseq_transactions.models import BlastDatabase


def create_settings():
    return BlastSettings.objects.create(
        e_value=0.001,
        word_size=3,
        num_alignments=10,
        max_target_seqs=10,
        num_threads=1,
        max_hsps=10,
    )


class ReciprocalWorkflowExecutionTests(TransactionTestCase):
    reset_sequences = True

    def setUp(self):
        self.owner = User.objects.create_user('workflow-owner', password='test')
        self.database = BlastDatabase.objects.create(
            database_name='workflow database',
            database_description='test database',
            assembly_entries=1,
            path_to_database_file='testfiles/databases/workflow',
        )
        self.project = BlastProject.objects.create(
            project_title='workflow local project',
            search_strategy='blastp',
            project_query_sequences='query.faa',
            project_user=self.owner,
            project_forward_settings=create_settings(),
            project_backward_settings=create_settings(),
            project_forward_database=self.database,
            project_backward_database=self.database,
            species_name_for_backward_blast='Species',
        )
        self.remote_project = RemoteBlastProject.objects.create(
            r_project_title='workflow remote project',
            r_search_strategy='blastp',
            r_project_query_sequences='query.faa',
            r_project_user=self.owner,
            r_project_forward_settings=create_settings(),
            r_project_backward_settings=create_settings(),
            r_project_forward_database='nr',
            r_project_backward_database=self.database,
            r_species_name_for_backward_blast='Species',
        )

    def task_sender(self, enqueued):
        def sender(project_id, execution_token):
            enqueued.append((project_id, execution_token))

        return sender

    def test_two_sequential_local_execution_requests_enqueue_only_one_task(self):
        enqueued = []

        first = start_reciprocal_workflow(self.project.id, task_sender=self.task_sender(enqueued))
        second = start_reciprocal_workflow(self.project.id, task_sender=self.task_sender(enqueued))

        self.assertTrue(first.started)
        self.assertFalse(second.started)
        self.assertEqual(1, len(enqueued))
        self.assertEqual(first.task_id, second.task_id)
        self.project.refresh_from_db()
        self.assertEqual(WorkflowLifecycle.QUEUED, self.project.project_workflow_state)
        self.assertEqual(first.task_id, self.project.project_workflow_execution_token)

    def test_concurrent_local_execution_requests_enqueue_only_one_task(self):
        enqueued = []

        def run_start():
            close_old_connections()
            try:
                return start_reciprocal_workflow(self.project.id, task_sender=self.task_sender(enqueued))
            finally:
                close_old_connections()

        with ThreadPoolExecutor(max_workers=2) as executor:
            results = list(executor.map(lambda _: run_start(), range(2)))

        self.assertEqual(1, sum(1 for result in results if result.started))
        self.assertEqual(1, len(enqueued))
        self.project.refresh_from_db()
        self.assertEqual(WorkflowLifecycle.QUEUED, self.project.project_workflow_state)

    def test_active_local_task_statuses_cannot_restart(self):
        active_statuses = ['PENDING', 'STARTED', 'PROGRESS', 'RETRY']
        for status in active_statuses:
            project = BlastProject.objects.create(
                project_title='active {}'.format(status),
                search_strategy='blastp',
                project_query_sequences='query.faa',
                project_user=self.owner,
                project_forward_settings=create_settings(),
                project_backward_settings=create_settings(),
                project_forward_database=self.database,
                project_backward_database=self.database,
                species_name_for_backward_blast='Species',
            )
            task_result = TaskResult.objects.create(task_id='active-{}'.format(status), status=status)
            project.project_execution_snakemake_task = task_result
            project.project_workflow_state = lifecycle_from_task_status(status)
            project.project_workflow_execution_token = task_result.task_id
            project.save()
            enqueued = []

            result = start_reciprocal_workflow(project.id, task_sender=self.task_sender(enqueued))

            self.assertFalse(result.started, status)
            self.assertEqual([], enqueued, status)
            project.refresh_from_db()
            self.assertEqual(lifecycle_from_task_status(status), project.project_workflow_state)

    def test_failed_local_project_can_restart(self):
        failed_task = TaskResult.objects.create(task_id='failed-local', status='FAILURE')
        self.project.project_execution_snakemake_task = failed_task
        self.project.project_workflow_state = WorkflowLifecycle.FAILED
        self.project.project_workflow_execution_token = failed_task.task_id
        self.project.save()
        enqueued = []

        result = start_reciprocal_workflow(self.project.id, task_sender=self.task_sender(enqueued))

        self.assertTrue(result.started)
        self.assertEqual(1, len(enqueued))
        self.project.refresh_from_db()
        self.assertEqual(WorkflowLifecycle.QUEUED, self.project.project_workflow_state)
        self.assertNotEqual('failed-local', self.project.project_workflow_execution_token)

    def test_user_exceeding_active_workflow_quota_cannot_enqueue_another_workflow(self):
        active_task = TaskResult.objects.create(task_id='quota-active', status='STARTED')
        self.project.project_execution_snakemake_task = active_task
        self.project.project_workflow_state = WorkflowLifecycle.RUNNING
        self.project.project_workflow_execution_token = active_task.task_id
        self.project.save()
        second_project = BlastProject.objects.create(
            project_title='quota second project',
            search_strategy='blastp',
            project_query_sequences='query.faa',
            project_user=self.owner,
            project_forward_settings=create_settings(),
            project_backward_settings=create_settings(),
            project_forward_database=self.database,
            project_backward_database=self.database,
            species_name_for_backward_blast='Species',
        )
        enqueued = []

        with self.assertRaises(WorkflowQuotaExceeded):
            start_reciprocal_workflow(second_project.id, task_sender=self.task_sender(enqueued))

        self.assertEqual([], enqueued)
        second_project.refresh_from_db()
        self.assertEqual(WorkflowLifecycle.STARTABLE, second_project.project_workflow_state)

    def test_duplicate_delivery_claims_only_once(self):
        enqueued = []
        result = start_reciprocal_workflow(self.project.id, task_sender=self.task_sender(enqueued))

        first_claim = claim_workflow_execution(self.project.id, result.task_id, 'local')
        second_claim = claim_workflow_execution(self.project.id, result.task_id, 'local')

        self.assertTrue(first_claim)
        self.assertFalse(second_claim)
        self.project.refresh_from_db()
        self.assertEqual(WorkflowLifecycle.RUNNING, self.project.project_workflow_state)

    def test_stale_task_cannot_mark_newer_local_execution_successful_or_failed(self):
        old_task = TaskResult.objects.create(task_id='old-local', status='STARTED')
        new_task = TaskResult.objects.create(task_id='new-local', status='STARTED')
        self.project.project_execution_snakemake_task = new_task
        self.project.project_workflow_state = WorkflowLifecycle.RUNNING
        self.project.project_workflow_execution_token = new_task.task_id
        self.project.save()

        self.assertFalse(finish_workflow_execution(self.project.id, old_task.task_id, 'local', successful=True))
        self.assertFalse(finish_workflow_execution(self.project.id, old_task.task_id, 'local', successful=False))

        self.project.refresh_from_db()
        self.assertEqual(WorkflowLifecycle.RUNNING, self.project.project_workflow_state)
        self.assertEqual(new_task.task_id, self.project.project_workflow_execution_token)

    def test_remote_workflows_receive_equivalent_start_protection(self):
        enqueued = []

        first = start_reciprocal_workflow(
            self.remote_project.id,
            remote=True,
            task_sender=self.task_sender(enqueued),
        )
        second = start_reciprocal_workflow(
            self.remote_project.id,
            remote=True,
            task_sender=self.task_sender(enqueued),
        )

        self.assertTrue(first.started)
        self.assertFalse(second.started)
        self.assertEqual(1, len(enqueued))
        self.remote_project.refresh_from_db()
        self.assertEqual(WorkflowLifecycle.QUEUED, self.remote_project.r_project_workflow_state)

    def test_failed_remote_project_can_restart(self):
        failed_task = TaskResult.objects.create(task_id='failed-remote', status='FAILURE')
        self.remote_project.r_project_execution_snakemake_task = failed_task
        self.remote_project.r_project_workflow_state = WorkflowLifecycle.FAILED
        self.remote_project.r_project_workflow_execution_token = failed_task.task_id
        self.remote_project.save()
        enqueued = []

        result = start_reciprocal_workflow(
            self.remote_project.id,
            remote=True,
            task_sender=self.task_sender(enqueued),
        )

        self.assertTrue(result.started)
        self.assertEqual(1, len(enqueued))
        self.remote_project.refresh_from_db()
        self.assertEqual(WorkflowLifecycle.QUEUED, self.remote_project.r_project_workflow_state)
        self.assertNotEqual('failed-remote', self.remote_project.r_project_workflow_execution_token)
