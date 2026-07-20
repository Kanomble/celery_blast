from unittest.mock import patch

from django.contrib.auth.models import User
from django.test import TestCase
from django.urls import reverse
from django_celery_results.models import TaskResult

from blast_project.setup_state import cathi_setup_is_running
from external_tools.models import DomainDatabase


class CathiSetupStateTests(TestCase):
    def test_setup_is_running_for_active_setup_task(self):
        DomainDatabase.objects.create(
            domain_database_loaded=False,
            domain_database_download_task_result=TaskResult.objects.create(
                task_id='setup-progress',
                status='PROGRESS',
            ),
        )

        self.assertTrue(cathi_setup_is_running())

    def test_setup_is_not_running_after_successful_load(self):
        DomainDatabase.objects.create(
            domain_database_loaded=True,
            domain_database_download_task_result=TaskResult.objects.create(
                task_id='setup-success',
                status='SUCCESS',
            ),
        )

        self.assertFalse(cathi_setup_is_running())


class CathiSetupAccessGateTests(TestCase):
    def setUp(self):
        self.user = User.objects.create_user(
            username='setup-gate-user',
            email='setup-gate@example.test',
            password='test',
        )
        self.client.login(username='setup-gate-user', password='test')
        DomainDatabase.objects.create(
            domain_database_loaded=False,
            domain_database_download_task_result=TaskResult.objects.create(
                task_id='setup-active',
                status='PROGRESS',
            ),
        )

    def test_one_way_project_creation_redirects_while_setup_is_running(self):
        response = self.client.get(reverse('one_way_project_creation'))

        self.assertRedirects(response, reverse('blast_project_dashboard'))

    def test_refseq_summary_download_does_not_enqueue_while_setup_is_running(self):
        with patch('refseq_transactions.views.download_refseq_assembly_summary.delay') as delay:
            response = self.client.post(
                reverse(
                    'download_refseq_assembly_summary',
                    kwargs={'summary_file': 'RefSeq'},
                )
            )

        self.assertRedirects(response, reverse('refseq_transactions_dashboard'))
        delay.assert_not_called()
