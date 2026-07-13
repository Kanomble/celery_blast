from django.contrib.auth.models import User
from django.http import Http404, HttpResponse
from django.test import Client, RequestFactory, TestCase

from blast_project.access import (
    get_owned_local_project_or_404,
    get_owned_remote_project_or_404,
    project_owner_required,
)
from blast_project.models import BlastProject, BlastSettings, RemoteBlastProject
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


class ProjectAccessTests(TestCase):
    def setUp(self):
        self.owner = User.objects.create_user('owner', password='test')
        self.other_user = User.objects.create_user('other', password='test')
        self.database = BlastDatabase.objects.create(
            database_name='access test database',
            database_description='test database',
            assembly_entries=1,
            path_to_database_file='testfiles/databases/access',
        )
        self.local_project = BlastProject.objects.create(
            project_title='access local project',
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
            r_project_title='access remote project',
            r_search_strategy='blastp',
            r_project_query_sequences='query.faa',
            r_project_user=self.owner,
            r_project_forward_settings=create_settings(),
            r_project_backward_settings=create_settings(),
            r_project_forward_database='nr',
            r_project_backward_database=self.database,
            r_species_name_for_backward_blast='Species',
        )
        self.client = Client()

    def test_local_project_owner_can_fetch_project(self):
        project = get_owned_local_project_or_404(self.owner, self.local_project.id)

        self.assertEqual(self.local_project, project)

    def test_local_project_non_owner_gets_404(self):
        with self.assertRaises(Http404):
            get_owned_local_project_or_404(self.other_user, self.local_project.id)

    def test_remote_project_owner_can_fetch_project(self):
        project = get_owned_remote_project_or_404(self.owner, self.remote_project.id)

        self.assertEqual(self.remote_project, project)

    def test_remote_project_non_owner_gets_404(self):
        with self.assertRaises(Http404):
            get_owned_remote_project_or_404(self.other_user, self.remote_project.id)

    def test_project_owner_decorator_blocks_non_owner_before_view_runs(self):
        called = []

        @project_owner_required()
        def protected_view(request, project_id):
            called.append(project_id)
            return HttpResponse('ok')

        request = RequestFactory().get('/protected')
        request.user = self.other_user

        with self.assertRaises(Http404):
            protected_view(request, project_id=self.local_project.id)

        self.assertEqual([], called)

    def test_project_owner_decorator_attaches_owned_project(self):
        @project_owner_required()
        def protected_view(request, project_id):
            self.assertEqual(self.local_project, request.owned_blast_project)
            return HttpResponse('ok')

        request = RequestFactory().get('/protected')
        request.user = self.owner

        response = protected_view(request, project_id=self.local_project.id)

        self.assertEqual(200, response.status_code)

    def test_project_details_returns_404_for_other_users_project(self):
        self.client.force_login(self.other_user)

        response = self.client.get('/blast_project/{}/project_details'.format(self.local_project.id))

        self.assertEqual(404, response.status_code)

    def test_remote_project_execution_requires_login(self):
        response = self.client.post('/blast_project/{}/remote_project_execution'.format(self.remote_project.id))

        self.assertEqual(302, response.status_code)
        self.assertIn('/blast_project/login', response['Location'])
