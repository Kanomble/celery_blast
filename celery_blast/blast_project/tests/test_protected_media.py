import os
import tempfile
from pathlib import Path
from unittest.mock import patch

from django.conf import settings
from django.contrib.auth.models import User
from django.http import Http404
from django.test import SimpleTestCase, TestCase, override_settings
from django.urls import reverse

from blast_project.models import BlastProject, BlastSettings
from blast_project.protected_media import resolve_project_file, serve_protected_project_file
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


class ProtectedProjectArtifactTests(TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tempdir.cleanup)
        self.media_root = Path(self.tempdir.name) / 'media'
        self.project_root = self.media_root / 'blast_projects' / '1'
        self.project_root.mkdir(parents=True)
        (self.project_root / 'result.txt').write_text('owned result', encoding='utf-8')
        (self.project_root / 'plot_amount_hits_of_target_taxon.png').write_bytes(b'png-data')

        self.owner = User.objects.create_user('owner', password='test')
        self.other_user = User.objects.create_user('other', password='test')
        self.database = BlastDatabase.objects.create(
            database_name='protected media database',
            database_description='test database',
            assembly_entries=1,
            path_to_database_file='testfiles/databases/access',
        )
        self.project = BlastProject.objects.create(
            id=1,
            project_title='protected media project',
            search_strategy='blastp',
            project_query_sequences='query.faa',
            project_user=self.owner,
            project_forward_settings=create_settings(),
            project_backward_settings=create_settings(),
            project_forward_database=self.database,
            project_backward_database=self.database,
            species_name_for_backward_blast='Species',
        )

    def artifact_url(self, relative_path):
        return reverse(
            'project_artifact',
            kwargs={
                'project_id': self.project.id,
                'remote_or_local': 'local',
                'relative_path': relative_path,
            },
        )

    def get_with_project_root(self, url):
        project_dir = str(self.media_root / 'blast_projects') + os.sep
        with patch('blast_project.views.BLAST_PROJECT_DIR', project_dir):
            return self.client.get(url)

    def test_anonymous_user_cannot_retrieve_project_artifact(self):
        response = self.get_with_project_root(self.artifact_url('result.txt'))

        self.assertNotEqual(200, response.status_code)

    def test_owner_can_retrieve_project_artifact(self):
        self.client.force_login(self.owner)

        response = self.get_with_project_root(self.artifact_url('result.txt'))

        self.assertEqual(200, response.status_code)
        self.assertEqual(b'owned result', b''.join(response.streaming_content))

    def test_other_user_cannot_retrieve_project_artifact(self):
        self.client.force_login(self.other_user)

        response = self.get_with_project_root(self.artifact_url('result.txt'))

        self.assertEqual(404, response.status_code)

    def test_traversal_and_absolute_paths_are_rejected(self):
        self.client.force_login(self.owner)

        unsafe_paths = [
            '../secret.txt',
            '%2e%2e/secret.txt',
            '/etc/passwd',
            'C%3A%5CWindows%5Cwin.ini',
            'nested/%2e%2e/secret.txt',
        ]
        for unsafe_path in unsafe_paths:
            with self.subTest(unsafe_path=unsafe_path):
                response = self.get_with_project_root(self.artifact_url(unsafe_path))
                self.assertEqual(404, response.status_code)

    def test_symlink_escape_is_rejected(self):
        outside_file = Path(self.tempdir.name) / 'outside.txt'
        outside_file.write_text('outside', encoding='utf-8')
        link_path = self.project_root / 'escape.txt'
        try:
            os.symlink(outside_file, link_path)
        except (OSError, NotImplementedError):
            self.skipTest('symlink creation is not available in this environment')

        with self.assertRaises(Http404):
            resolve_project_file(self.project_root, 'escape.txt')

    def test_missing_files_return_safe_404(self):
        self.client.force_login(self.owner)

        response = self.get_with_project_root(self.artifact_url('missing.txt'))

        self.assertEqual(404, response.status_code)
        self.assertNotContains(response, str(self.project_root), status_code=404)

    def test_direct_media_url_is_not_served_by_django(self):
        response = self.client.get('/media/blast_projects/{}/result.txt'.format(self.project.id))

        self.assertEqual(404, response.status_code)

    def test_existing_authorized_result_view_still_works(self):
        (self.project_root / 'reciprocal_results.html').write_text('<p>results</p>', encoding='utf-8')
        self.client.force_login(self.owner)
        project_dir = str(self.media_root / 'blast_projects') + os.sep

        with patch('blast_project.views.BLAST_PROJECT_DIR', project_dir):
            response = self.client.get(reverse('reciprocal_results', kwargs={'project_id': self.project.id}))

        self.assertEqual(200, response.status_code)
        self.assertContains(response, '<p>results</p>', html=False)

    @override_settings(
        PROTECTED_MEDIA_USE_X_ACCEL=True,
        PROTECTED_MEDIA_X_ACCEL_SIZE_THRESHOLD=0,
    )
    def test_x_accel_redirect_uses_internal_relative_path(self):
        with override_settings(PROTECTED_MEDIA_ROOT=str(self.media_root)):
            response = serve_protected_project_file(self.project_root, 'result.txt')

        self.assertEqual('/_protected_media/blast_projects/1/result.txt', response['X-Accel-Redirect'])
        self.assertNotIn(str(self.media_root), response['X-Accel-Redirect'])


class NginxProtectedMediaConfigTests(SimpleTestCase):
    def test_nginx_denies_public_media_and_uses_internal_location(self):
        nginx_conf = Path(settings.BASE_DIR).parent / 'nginx' / 'nginx.conf'
        if not nginx_conf.exists():
            self.skipTest('nginx configuration is not mounted in this runtime')
        content = nginx_conf.read_text(encoding='utf-8')

        self.assertIn('location /media', content)
        self.assertIn('return 404;', content)
        self.assertIn('location /_protected_media/', content)
        self.assertIn('internal;', content)
