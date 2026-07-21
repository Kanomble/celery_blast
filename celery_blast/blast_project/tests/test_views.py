#cmd for running these tests: python manage.py test blast_project.tests.test_views
#different tags: long fast view view_combination
from django.test import Client, TestCase, tag
from django.contrib.auth.models import User
from django.urls import reverse
from pathlib import Path

from blast_project.models import BlastProject, BlastSettings, RemoteBlastProject
from external_tools.models import DomainDatabase
from refseq_transactions.models import AssemblyLevels, BlastDatabase
from django.utils import timezone
from django_celery_results.models import TaskResult
#from blast_project.py_django_db_services import get_all_succeeded_databases
from django.core.files.uploadedfile import SimpleUploadedFile
from blast_project.py_services import check_if_taxdb_exists, delete_project_and_associated_directories_by_id
from shutil import rmtree
from os.path import isdir
from time import sleep
from unittest.mock import patch

@tag('blastprojectviews')
class BlastProjectViewsTestCase(TestCase):
    c = Client()
    def setUp(self) -> None:
        print("[+] Setting up model data for testing")
        user = User.objects.create_user(
            'testuser',
            'test_email@email.com',
            password='test',
            last_login=timezone.now()
        )
        AssemblyLevels.objects.create(
            assembly_level='Complete Genome'
        )
        AssemblyLevels.objects.create(
            assembly_level='Chromosome'
        )
        AssemblyLevels.objects.create(
            assembly_level='Contig'
        )
        AssemblyLevels.objects.create(
            assembly_level='Scaffold'
        )

        fw_settings = BlastSettings.objects.create(
            e_value=0.001,
            word_size=3,
            num_alignments=10000,
            max_target_seqs=10000,
            num_threads=1,
            max_hsps=500
        )

        bw_settings = BlastSettings.objects.create(
            e_value=0.001,
            word_size=3,
            num_alignments=1,
            max_target_seqs=1,
            num_threads=1,
            max_hsps=500
        )

        celery_task = TaskResult.objects.create(
            status="SUCCESS"
        )

        fw_and_bw_db = BlastDatabase.objects.create(
            database_name = 'Curvibacter sp. aep1-3 database',
            database_description = 'Full Assembly Content',
            assembly_entries = 1,
            path_to_database_file="testfiles/databases/curvibacter_test_db",
            database_download_and_format_task=celery_task
        )

        assembly_level = ['Complete Genome', 'Chromosome', 'Contig', 'Scaffold']
        assembly_levels = AssemblyLevels.objects.filter(assembly_level__in=assembly_level)
        for assembly_lvl in assembly_levels:
            fw_and_bw_db.assembly_levels.add(assembly_lvl)

        fw_and_bw_db.save()

        BlastProject.objects.create(
            project_title='test project 1',
            search_strategy='blastp',
            project_query_sequences='lps_transport.faa',
            project_user=user,
            project_forward_settings=fw_settings,
            project_backward_settings=bw_settings,
            project_forward_database=fw_and_bw_db,
            project_backward_database=fw_and_bw_db,
            species_name_for_backward_blast='Curvibacter sp. AEP1-3'
        )

    def _set_domain_database_state(self, loaded=False, status=None):
        DomainDatabase.objects.all().delete()
        task_result = None
        if status is not None:
            task_result = TaskResult.objects.create(
                task_id='domain-database-{}'.format(status.lower()),
                status=status,
            )
        return DomainDatabase.objects.create(
            domain_database_loaded=loaded,
            domain_database_download_task_result=task_result,
        )

    def _get_logged_in_dashboard(self):
        self.c.login(username='testuser', password='test')
        return self.c.get('/blast_project/')

    def _create_remote_project(self, status=None):
        task_result = None
        if status is not None:
            task_result = TaskResult.objects.create(
                task_id='remote-project-{}'.format(status.lower()),
                status=status,
            )

        user = User.objects.get(username='testuser')
        database = BlastDatabase.objects.get(database_name='Curvibacter sp. aep1-3 database')
        return RemoteBlastProject.objects.create(
            r_project_title='test remote project {}'.format(status or 'not-started'),
            r_search_strategy='blastp',
            r_project_query_sequences='lps_transport.faa',
            r_project_user=user,
            r_project_forward_settings=BlastSettings.objects.create(
                e_value=0.001,
                word_size=3,
                num_alignments=10000,
                max_target_seqs=10000,
                num_threads=1,
                max_hsps=500,
            ),
            r_project_backward_settings=BlastSettings.objects.create(
                e_value=0.001,
                word_size=3,
                num_alignments=1,
                max_target_seqs=1,
                num_threads=1,
                max_hsps=500,
            ),
            r_project_forward_database='nr',
            r_project_backward_database=database,
            r_species_name_for_backward_blast='Curvibacter sp. AEP1-3',
            r_project_execution_snakemake_task=task_result,
        )

    def _create_remote_bokeh_template(self, project_id):
        bokeh_dir = Path('media/blast_projects/remote_projects') / str(project_id)
        bokeh_dir.mkdir(parents=True, exist_ok=True)
        self.addCleanup(rmtree, str(bokeh_dir), True)
        bokeh_file = bokeh_dir / 'interactive_bokeh_plot.html'
        bokeh_file.write_text('<div id="remote-bokeh-placeholder">Remote plot</div>', encoding='utf-8')
        return bokeh_file

    def tearDown(self) -> None:
        print("[*] Deleting all BlastProjects from test database")
        for project in BlastProject.objects.all():
            print("\t[*] Trying to delete project: {} with id {}".format(project.project_title, project.id))
            delete_project_and_associated_directories_by_id(project.id)

    @tag('fast','view')
    def test_register_get(self):
        response = self.c.get('/blast_project/register')
        self.assertEqual(200,response.status_code)
        self.assertTemplateUsed(response, template_name='blast_project/register.html')

    @tag('fast','view')
    def test_register_post(self):
        response = self.c.post('/blast_project/register',
                               {'username':'testuser2',
                                'password':'nvm',
                                'email':'testuser2@email.de'})
        self.assertEqual(200,response.status_code)
        self.assertTemplateUsed(response, template_name='blast_project/register.html')

    @tag('fast','view')
    def test_login_get(self):
        response = self.c.get('/blast_project/login')
        self.assertEqual(200,response.status_code)
        self.assertTemplateUsed(response, template_name='blast_project/login.html')

    @tag('fast','view')
    def test_login_post(self):
        response = self.c.post('/blast_project/login',{'username':'testuser','password':'test'})
        self.assertEqual(302,response.status_code)
        self.assertRedirects(response, '/blast_project/',status_code=302,
                             target_status_code=200,fetch_redirect_response=True)

    @tag('fast','view')
    def test_dashboard_view_get_not_logged_in(self):
        response = self.c.get('/blast_project/')
        self.assertEqual(302,response.status_code)
        self.assertRedirects(response, '/blast_project/login?next=/blast_project/',status_code=302,
                             target_status_code=200,fetch_redirect_response=True)

    @tag('fast','view')
    def test_dashboard_view_get(self):
        self.c.login(username='testuser',password='test')
        project = BlastProject.objects.get(project_title='test project 1')
        remote_project = self._create_remote_project()
        response = self.c.get('/blast_project/')
        self.assertEqual(200,response.status_code)
        self.assertTemplateUsed(response,'blast_project/blast_project_dashboard.html')
        #response.context['blast_projects'] is a django.db.models.query.QuerySet class
        self.assertEqual(response.context['blast_projects'][0],
                         project)
        self.assertEqual(response.context['RemoteBlastProjects'][0],
                         remote_project)
        self.assertEqual(response.context['ActiveBlastDatabases'][0],
                         project.project_forward_database)

    @tag('fast', 'view')
    def test_dashboard_complete_state_has_refactored_information_hierarchy(self):
        self._set_domain_database_state(loaded=True, status='SUCCESS')
        self._create_remote_project()

        response = self._get_logged_in_dashboard()

        self.assertEqual(200, response.status_code)
        self.assertContains(response, '<h1', count=1, html=False)
        self.assertContains(
            response,
            'href="{}"'.format(reverse('project_creation')),
            html=False,
        )
        self.assertContains(response, 'Create a reciprocal BLAST project')
        self.assertContains(response, 'id="dashboard_workflows"', html=False)
        self.assertContains(response, 'id="dashboard_project_summary"', html=False)
        self.assertContains(response, 'Remote reciprocal projects')
        self.assertContains(response, 'id="dashboard_help"', html=False)

        html = response.content.decode()
        self.assertNotIn('onclick="displayDivElement', html)

        help_sections = {
            'overview': 'overview_heading',
            'target_search': 'target_search_heading',
            'database_management': 'database_management_heading',
            'project_creation': 'project_creation_heading',
        }
        for section_id, heading_id in help_sections.items():
            self.assertContains(response, 'id="{}"'.format(section_id), html=False)
            self.assertRegex(
                html,
                r'<section[^>]+id="{}"[^>]+aria-labelledby="{}"'.format(
                    section_id,
                    heading_id,
                ),
            )

        help_controls = {
            'dashboard-help-overview-button': 'overview',
            'dashboard-help-target-search-button': 'target_search',
            'dashboard-help-database-management-button': 'database_management',
            'dashboard-help-project-creation-button': 'project_creation',
        }
        for button_id, section_id in help_controls.items():
            self.assertRegex(
                html,
                r'<button[^>]+id="{}"[^>]+type="button"[^>]+aria-controls="{}"[^>]+aria-expanded="true"'.format(
                    button_id,
                    section_id,
                ),
            )

    @tag('fast', 'view')
    def test_dashboard_cdd_database_missing_renders_setup_control(self):
        self._set_domain_database_state(loaded=False)

        response = self._get_logged_in_dashboard()

        self.assertEqual(200, response.status_code)
        self.assertContains(response, 'id="project_dashboard"', html=False)
        self.assertContains(
            response,
            'action="{}"'.format(reverse('setup_cathi_view')),
            html=False,
        )
        self.assertContains(response, 'Start CATHI Tool Set-Up Procedure')

    @tag('fast', 'view')
    def test_dashboard_cdd_download_progress_renders_progress_control(self):
        self._set_domain_database_state(loaded=False, status='PROGRESS')

        response = self._get_logged_in_dashboard()

        self.assertEqual(200, response.status_code)
        self.assertContains(response, 'id="cdd_download_progress"', html=False)
        self.assertContains(response, 'id="download_button"', html=False)
        self.assertContains(response, 'domain_database_download_progress')

    @tag('fast', 'view')
    def test_dashboard_cdd_download_failure_renders_restart_control(self):
        self._set_domain_database_state(loaded=False, status='FAILURE')

        response = self._get_logged_in_dashboard()

        self.assertEqual(200, response.status_code)
        self.assertContains(response, 'id="cdd_download_failure"', html=False)
        self.assertContains(
            response,
            'action="{}"'.format(reverse('delete_domain_database')),
            html=False,
        )
        self.assertContains(response, 'Set-Up Failed: Delete Unfinished Task And Restart')

    @tag('fast', 'view')
    def test_project_creation_get_preserves_project_type_selection_and_forms(self):
        self._set_domain_database_state(loaded=True, status='SUCCESS')
        self.c.login(username='testuser', password='test')

        with patch('blast_project.views.check_if_taxdb_exists', return_value=True):
            response = self.c.get(reverse('project_creation'))

        self.assertEqual(200, response.status_code)
        self.assertTemplateUsed(response, 'blast_project/project_creation_dashboard.html')
        self.assertContains(response, '<h1', count=1, html=False)
        self.assertContains(response, 'id="local"', html=False)
        self.assertContains(response, 'id="remote"', html=False)
        self.assertContains(response, 'for="local"', html=False)
        self.assertContains(response, 'for="remote"', html=False)
        self.assertContains(response, 'id="project_creation_local"', html=False)
        self.assertContains(response, 'id="project_creation_remote"', html=False)
        self.assertContains(response, 'id="advanced_settings"', html=False)
        self.assertContains(response, 'id="advanced_remote_settings"', html=False)

        expected_local_field_names = (
            'project_title',
            'query_sequence_file',
            'query_sequence_text',
            'project_forward_database',
            'project_backward_database',
            'species_name_for_backward_blast',
        )
        for field_name in expected_local_field_names:
            self.assertContains(response, 'name="{}"'.format(field_name), html=False)

        expected_remote_field_names = (
            'r_project_title',
            'r_query_sequence_file',
            'r_query_sequence_text',
            'r_project_forward_database',
            'r_project_backward_database',
            'r_entrez_query',
            'r_species_name_for_backward_blast',
        )
        for field_name in expected_remote_field_names:
            self.assertContains(response, 'name="{}"'.format(field_name), html=False)

        expected_shared_settings_names = (
            'fw_e_value',
            'fw_word_size',
            'fw_num_alignments',
            'fw_max_target_seqs',
            'fw_num_threads',
            'fw_max_hsps',
            'bw_e_value',
            'bw_word_size',
            'bw_num_alignments',
            'bw_max_target_seqs',
            'bw_num_threads',
            'bw_max_hsps',
            'bitscore_filter',
            'max_amount_of_rbh_for_msa_and_phylogeny',
            'trimal_gt',
            'trimal_st',
            'trimal_cons',
            'mview_sort',
            'mview_coloring',
        )
        for field_name in expected_shared_settings_names:
            self.assertContains(
                response,
                'name="{}"'.format(field_name),
                count=2,
                html=False,
            )

        project_creation_url = reverse('project_creation')
        self.assertContains(
            response,
            'action="{}"'.format(project_creation_url),
            count=2,
            html=False,
        )

        html = response.content.decode()
        self.assertRegex(
            html,
            r'<input[^>]+type="hidden"[^>]+value="local"[^>]+id="project_type"[^>]+name="project_type"',
        )
        self.assertRegex(
            html,
            r'<input[^>]+type="hidden"[^>]+value="remote"[^>]+id="project_type"[^>]+name="project_type"',
        )

    @tag('fast', 'view')
    def test_project_creation_local_invalid_submission_renders_accessible_errors(self):
        self._set_domain_database_state(loaded=True, status='SUCCESS')
        self.c.login(username='testuser', password='test')

        post_data = {
            'project_type': 'local',
            'project_title': 'Retained local project title',
            'query_sequence_text': '',
            'fw_e_value': 'not-a-number',
            'fw_word_size': '3',
            'fw_num_alignments': '10000',
            'fw_max_target_seqs': '10000',
            'fw_num_threads': '1',
            'fw_max_hsps': '500',
            'bw_e_value': '0.001',
            'bw_word_size': '3',
            'bw_num_alignments': '1',
            'bw_max_target_seqs': '1',
            'bw_num_threads': '1',
            'bw_max_hsps': '500',
            'bitscore_filter': '50',
            'max_amount_of_rbh_for_msa_and_phylogeny': '500',
            'trimal_gt': '0.8',
            'trimal_st': '0.001',
            'trimal_cons': '60',
            'mview_sort': 'coverage',
            'mview_coloring': 'any',
        }

        with patch('blast_project.views.check_if_taxdb_exists', return_value=True):
            response = self.c.post(reverse('project_creation'), post_data)

        self.assertEqual(200, response.status_code)
        self.assertTemplateUsed(response, 'blast_project/project_creation_dashboard.html')
        self.assertContains(response, 'action="{}"'.format(reverse('project_creation')), html=False)
        self.assertContains(response, 'id="local_project_error_summary"', html=False)
        self.assertContains(response, 'Retained local project title')
        self.assertContains(response, 'value="Retained local project title"', html=False)
        self.assertContains(response, 'id="fw_e_value_errors"', html=False)
        self.assertContains(response, 'Enter a number.')
        self.assertContains(response, 'aria-invalid="true"', html=False)
        self.assertContains(response, 'for="project_title"', html=False)
        self.assertContains(response, 'for="query_sequence_file"', html=False)
        self.assertContains(response, 'for="query_sequence_text"', html=False)
        self.assertContains(response, 'Upload a FASTA file or enter protein identifiers separated by commas')
        self.assertContains(response, 'value="local" id="project_type" name="project_type"', html=False)

        for field_name in (
                'project_title',
                'query_sequence_file',
                'query_sequence_text',
                'project_forward_database',
                'project_backward_database',
                'species_name_for_backward_blast',
                'project_type',
        ):
            self.assertContains(response, 'name="{}"'.format(field_name), html=False)

    @tag('fast', 'view')
    def test_project_creation_remote_invalid_submission_renders_accessible_errors(self):
        self._set_domain_database_state(loaded=True, status='SUCCESS')
        self.c.login(username='testuser', password='test')

        post_data = {
            'project_type': 'remote',
            'r_project_title': 'Retained remote project title',
            'r_query_sequence_text': '',
            'r_project_forward_database': 'nr',
            'r_entrez_query': '',
            'fw_e_value': '0.001',
            'fw_word_size': '3',
            'fw_num_alignments': '10000',
            'fw_max_target_seqs': '10000',
            'fw_num_threads': '1',
            'fw_max_hsps': '500',
            'bw_e_value': 'bad-number',
            'bw_word_size': '3',
            'bw_num_alignments': '1',
            'bw_max_target_seqs': '1',
            'bw_num_threads': '1',
            'bw_max_hsps': '500',
            'bitscore_filter': '50',
            'max_amount_of_rbh_for_msa_and_phylogeny': '500',
            'trimal_gt': '0.8',
            'trimal_st': '0.001',
            'trimal_cons': '60',
            'mview_sort': 'coverage',
            'mview_coloring': 'any',
        }

        with patch('blast_project.views.check_if_taxdb_exists', return_value=True):
            response = self.c.post(reverse('project_creation'), post_data)

        self.assertEqual(200, response.status_code)
        self.assertTemplateUsed(response, 'blast_project/project_creation_dashboard.html')
        self.assertContains(response, 'action="{}"'.format(reverse('project_creation')), html=False)
        self.assertContains(response, 'id="remote_project_error_summary"', html=False)
        self.assertContains(response, 'Retained remote project title')
        self.assertContains(response, 'value="Retained remote project title"', html=False)
        self.assertContains(response, 'id="r_entrez_query_errors"', html=False)
        self.assertContains(response, 'id="bw_e_value_errors"', html=False)
        self.assertContains(response, 'Enter a number.')
        self.assertContains(response, 'aria-invalid="true"', html=False)
        self.assertContains(response, 'for="r_project_title"', html=False)
        self.assertContains(response, 'for="r_query_sequence_file"', html=False)
        self.assertContains(response, 'for="r_query_sequence_text"', html=False)
        self.assertContains(response, 'Upload a FASTA file or enter protein identifiers separated by commas')
        self.assertContains(response, 'value="remote" id="project_type" name="project_type"', html=False)

        for field_name in (
                'r_project_title',
                'r_query_sequence_file',
                'r_query_sequence_text',
                'r_project_forward_database',
                'r_project_backward_database',
                'r_entrez_query',
                'r_species_name_for_backward_blast',
                'project_type',
        ):
            self.assertContains(response, 'name="{}"'.format(field_name), html=False)

    @tag('fast','view')
    def test_project_details_get(self):
        project = BlastProject.objects.get(project_title='test project 1')
        self.c.login(username='testuser', password='test')
        response = self.c.get('/blast_project/'+str(project.id)+'/project_details')
        self.assertEqual(200,response.status_code)
        self.assertEqual(response.context['BlastProject'],project)
        self.assertNotContains(response, 'id="myTable"', html=False)
        self.assertContains(response, "initializeReciprocalResultTable('myTable');", html=False)
        self.assertContains(response, 'javascript/load_ajax_table.js', html=False)

    @tag('fast', 'view')
    def test_project_details_get_with_success_task_and_missing_reciprocal_summary(self):
        project = BlastProject.objects.get(project_title='test project 1')
        project.project_execution_snakemake_task = TaskResult.objects.create(
            status='SUCCESS',
            task_id='detail-success-missing-summary',
        )
        project.save()

        self.c.login(username='testuser', password='test')
        response = self.c.get('/blast_project/'+str(project.id)+'/project_details')

        self.assertEqual(200, response.status_code)
        self.assertTemplateUsed(response, 'blast_project/project_details_dashboard.html')

    @tag('fast', 'view')
    def test_project_details_progress_state_renders_status_and_preserves_polling_target(self):
        project = BlastProject.objects.get(project_title='test project 1')
        project.project_execution_snakemake_task = TaskResult.objects.create(
            status='PROGRESS',
            task_id='detail-progress',
        )
        project.save()

        self.c.login(username='testuser', password='test')
        response = self.c.get('/blast_project/'+str(project.id)+'/project_details')

        self.assertEqual(200, response.status_code)
        self.assertContains(response, 'Status: In progress')
        self.assertContains(response, 'id="progress_bar"', html=False)
        self.assertContains(
            response,
            'data-progress-poll-url="{}"'.format(
                reverse('ajax_call_to_logfiles', kwargs={'project_id': project.id}),
            ),
            html=False,
        )
        self.assertContains(response, 'data-progress-bar-id="progress_bar"', html=False)
        self.assertContains(response, 'data-progress-poll-interval="20000"', html=False)
        self.assertContains(response, 'Loading ...')
        self.assertNotContains(response, 'Delete Blast Project')
        self.assertNotContains(response, "jQuery.ajax({")
        self.assertNotContains(response, 'id="myTable"', html=False)
        self.assertContains(response, "initializeReciprocalResultTable('myTable');", html=False)
        self.assertContains(response, 'javascript/load_ajax_table.js', html=False)
        self.assertNotContains(
            response,
            reverse('ajax_call_to_remote_logfiles', kwargs={'project_id': project.id}),
            html=False,
        )

    @tag('fast', 'view')
    def test_ajax_call_to_logfiles_returns_progress_for_ajax_request(self):
        project = BlastProject.objects.get(project_title='test project 1')
        project.project_execution_snakemake_task = TaskResult.objects.create(
            status='SUCCESS',
            task_id='detail-progress-endpoint',
        )
        project.save()

        self.c.login(username='testuser', password='test')
        response = self.c.get(
            reverse('ajax_call_to_logfiles', kwargs={'project_id': project.id}),
            HTTP_X_REQUESTED_WITH='XMLHttpRequest',
        )

        self.assertEqual(200, response.status_code)
        self.assertEqual('application/json', response['Content-Type'].split(';')[0])
        self.assertEqual({'progress': 100}, response.json())

    @tag('fast', 'view')
    def test_ajax_call_to_logfiles_returns_json_errors_for_invalid_requests(self):
        project = BlastProject.objects.get(project_title='test project 1')
        project.project_execution_snakemake_task = TaskResult.objects.create(
            status='SUCCESS',
            task_id='detail-progress-invalid-requests',
        )
        project.save()
        url = reverse('ajax_call_to_logfiles', kwargs={'project_id': project.id})

        self.c.login(username='testuser', password='test')

        missing_ajax_response = self.c.get(url)
        self.assertEqual(400, missing_ajax_response.status_code)
        self.assertEqual('application/json', missing_ajax_response['Content-Type'].split(';')[0])
        self.assertEqual({'error': 'AJAX request required'}, missing_ajax_response.json())

        unsupported_method_response = self.c.post(
            url,
            HTTP_X_REQUESTED_WITH='XMLHttpRequest',
        )
        self.assertEqual(405, unsupported_method_response.status_code)
        self.assertEqual('GET', unsupported_method_response['Allow'])
        self.assertEqual('application/json', unsupported_method_response['Content-Type'].split(';')[0])
        self.assertEqual({'error': 'GET method required'}, unsupported_method_response.json())

        missing_project_response = self.c.get(
            reverse('ajax_call_to_logfiles', kwargs={'project_id': project.id + 1000}),
            HTTP_X_REQUESTED_WITH='XMLHttpRequest',
        )
        self.assertEqual(404, missing_project_response.status_code)
        self.assertEqual('application/json', missing_project_response['Content-Type'].split(';')[0])
        self.assertEqual({'error': 'Project not found'}, missing_project_response.json())

    @tag('fast', 'view')
    def test_ajax_call_to_logfiles_returns_json_for_unavailable_or_failed_progress(self):
        project = BlastProject.objects.get(project_title='test project 1')
        url = reverse('ajax_call_to_logfiles', kwargs={'project_id': project.id})

        self.c.login(username='testuser', password='test')

        unavailable_response = self.c.get(
            url,
            HTTP_X_REQUESTED_WITH='XMLHttpRequest',
        )
        self.assertEqual(409, unavailable_response.status_code)
        self.assertEqual('application/json', unavailable_response['Content-Type'].split(';')[0])
        self.assertEqual(
            {'error': 'Project execution status is unavailable'},
            unavailable_response.json(),
        )

        project.project_execution_snakemake_task = TaskResult.objects.create(
            status='PROGRESS',
            task_id='detail-progress-processing-failure',
        )
        project.save()
        with patch('blast_project.views.read_task_logs_summary_table', side_effect=Exception):
            failure_response = self.c.get(
                url,
                HTTP_X_REQUESTED_WITH='XMLHttpRequest',
            )

        self.assertEqual(500, failure_response.status_code)
        self.assertEqual('application/json', failure_response['Content-Type'].split(';')[0])
        self.assertEqual({'error': 'Unable to read progress'}, failure_response.json())

    @tag('fast', 'view')
    def test_ajax_call_to_logfiles_redirects_unauthenticated_requests(self):
        project = BlastProject.objects.get(project_title='test project 1')
        self.c.logout()

        response = self.c.get(
            reverse('ajax_call_to_logfiles', kwargs={'project_id': project.id}),
            HTTP_X_REQUESTED_WITH='XMLHttpRequest',
        )

        self.assertEqual(302, response.status_code)
        self.assertIn(reverse('login'), response['Location'])

    @tag('fast', 'view')
    def test_project_details_completed_state_renders_local_result_link(self):
        project = BlastProject.objects.get(project_title='test project 1')
        project.project_execution_snakemake_task = TaskResult.objects.create(
            status='SUCCESS',
            task_id='detail-success-result-link',
        )
        project.save()

        self.c.login(username='testuser', password='test')
        response = self.c.get('/blast_project/'+str(project.id)+'/project_details')

        self.assertEqual(200, response.status_code)
        self.assertContains(response, 'Status: Success')
        self.assertContains(
            response,
            'href="{}"'.format(reverse('reciprocal_results', kwargs={'project_id': project.id})),
            html=False,
        )
        self.assertContains(response, "initializeReciprocalResultTable('myTable');", html=False)
        self.assertContains(response, 'javascript/load_ajax_table.js', html=False)
        self.assertContains(response, 'Query Sequence Information')
        self.assertContains(response, 'Reciprocal Result Summary')

    @tag('fast', 'view')
    def test_project_details_failure_state_renders_failure_information(self):
        project = BlastProject.objects.get(project_title='test project 1')
        project.project_execution_snakemake_task = TaskResult.objects.create(
            status='FAILURE',
            task_id='detail-failure',
        )
        project.save()

        self.c.login(username='testuser', password='test')
        response = self.c.get('/blast_project/'+str(project.id)+'/project_details')

        self.assertEqual(200, response.status_code)
        self.assertContains(response, 'Status: Failure')
        self.assertContains(response, 'An error occurred during execution of the pipeline.')
        self.assertContains(response, 'Delete Blast Project')
        self.assertContains(response, 'id="progress_bar"', html=False)
        self.assertNotContains(response, 'id="myTable"', html=False)
        self.assertContains(response, "initializeReciprocalResultTable('myTable');", html=False)
        self.assertContains(response, 'javascript/load_ajax_table.js', html=False)

    @tag('fast', 'view')
    def test_remote_project_details_not_started_does_not_render_result_table(self):
        remote_project = self._create_remote_project()

        self.c.login(username='testuser', password='test')
        response = self.c.get('/blast_project/'+str(remote_project.id)+'/remote_project_details')

        self.assertEqual(200, response.status_code)
        self.assertNotContains(response, 'id="myTable"', html=False)
        self.assertContains(response, "initializeReciprocalResultTable('myTable');", html=False)
        self.assertContains(response, 'javascript/load_ajax_table.js', html=False)

    @tag('fast', 'view')
    def test_remote_project_details_completed_state_renders_remote_result_link(self):
        remote_project = self._create_remote_project(status='SUCCESS')
        self._create_remote_bokeh_template(remote_project.id)

        self.c.login(username='testuser', password='test')
        response = self.c.get('/blast_project/'+str(remote_project.id)+'/remote_project_details')

        self.assertEqual(200, response.status_code)
        self.assertTemplateUsed(response, 'blast_project/remote_project_details_dashboard.html')
        self.assertContains(response, 'Status: Success')
        self.assertContains(
            response,
            'href="{}"'.format(reverse('remote_reciprocal_results', kwargs={'project_id': remote_project.id})),
            html=False,
        )
        self.assertContains(response, "initializeReciprocalResultTable('myTable');", html=False)
        self.assertContains(response, 'javascript/load_ajax_table.js', html=False)
        self.assertContains(response, 'id="interactive_bokeh_plot"', html=False)
        self.assertContains(response, 'Remote plot')

    @tag('fast', 'view')
    def test_remote_project_details_progress_and_failure_states_render_status(self):
        progress_project = self._create_remote_project(status='PROGRESS')
        failure_project = self._create_remote_project(status='FAILURE')

        self.c.login(username='testuser', password='test')
        progress_response = self.c.get('/blast_project/'+str(progress_project.id)+'/remote_project_details')
        failure_response = self.c.get('/blast_project/'+str(failure_project.id)+'/remote_project_details')

        self.assertEqual(200, progress_response.status_code)
        self.assertContains(progress_response, 'Status: In progress')
        self.assertContains(progress_response, 'id="progress_bar"', html=False)
        self.assertContains(
            progress_response,
            'data-progress-poll-url="{}"'.format(
                reverse('ajax_call_to_remote_logfiles', kwargs={'project_id': progress_project.id}),
            ),
            html=False,
        )
        self.assertContains(progress_response, 'data-progress-bar-id="progress_bar"', html=False)
        self.assertContains(progress_response, 'data-progress-poll-interval="20000"', html=False)
        self.assertContains(progress_response, 'Loading ...')
        self.assertNotContains(progress_response, 'Delete Blast Project')
        self.assertNotContains(progress_response, 'id="myTable"', html=False)
        self.assertContains(progress_response, "initializeReciprocalResultTable('myTable');", html=False)
        self.assertContains(progress_response, 'javascript/load_ajax_table.js', html=False)
        self.assertNotContains(
            progress_response,
            reverse('ajax_call_to_logfiles', kwargs={'project_id': progress_project.id}),
            html=False,
        )

        self.assertEqual(200, failure_response.status_code)
        self.assertContains(failure_response, 'Status: Failure')
        self.assertContains(failure_response, 'An error occurred during execution of the pipeline.')
        self.assertContains(failure_response, 'Delete Blast Project')
        self.assertNotContains(failure_response, 'id="myTable"', html=False)
        self.assertContains(failure_response, "initializeReciprocalResultTable('myTable');", html=False)
        self.assertContains(failure_response, 'javascript/load_ajax_table.js', html=False)

    @tag('fast', 'view')
    def test_ajax_call_to_remote_logfiles_returns_progress_for_ajax_request(self):
        remote_project = self._create_remote_project(status='SUCCESS')

        self.c.login(username='testuser', password='test')
        response = self.c.get(
            reverse('ajax_call_to_remote_logfiles', kwargs={'project_id': remote_project.id}),
            HTTP_X_REQUESTED_WITH='XMLHttpRequest',
        )

        self.assertEqual(200, response.status_code)
        self.assertEqual('application/json', response['Content-Type'].split(';')[0])
        self.assertEqual({'progress': 100}, response.json())

    @tag('fast', 'view')
    def test_ajax_call_to_remote_logfiles_returns_json_errors_for_invalid_requests(self):
        remote_project = self._create_remote_project(status='SUCCESS')
        url = reverse('ajax_call_to_remote_logfiles', kwargs={'project_id': remote_project.id})

        self.c.login(username='testuser', password='test')

        missing_ajax_response = self.c.get(url)
        self.assertEqual(400, missing_ajax_response.status_code)
        self.assertEqual('application/json', missing_ajax_response['Content-Type'].split(';')[0])
        self.assertEqual({'error': 'AJAX request required'}, missing_ajax_response.json())

        unsupported_method_response = self.c.post(
            url,
            HTTP_X_REQUESTED_WITH='XMLHttpRequest',
        )
        self.assertEqual(405, unsupported_method_response.status_code)
        self.assertEqual('GET', unsupported_method_response['Allow'])
        self.assertEqual('application/json', unsupported_method_response['Content-Type'].split(';')[0])
        self.assertEqual({'error': 'GET method required'}, unsupported_method_response.json())

        missing_project_response = self.c.get(
            reverse('ajax_call_to_remote_logfiles', kwargs={'project_id': remote_project.id + 1000}),
            HTTP_X_REQUESTED_WITH='XMLHttpRequest',
        )
        self.assertEqual(404, missing_project_response.status_code)
        self.assertEqual('application/json', missing_project_response['Content-Type'].split(';')[0])
        self.assertEqual({'error': 'Project not found'}, missing_project_response.json())

    @tag('fast', 'view')
    def test_ajax_call_to_remote_logfiles_returns_json_for_unavailable_or_failed_progress(self):
        remote_project = self._create_remote_project()
        url = reverse('ajax_call_to_remote_logfiles', kwargs={'project_id': remote_project.id})

        self.c.login(username='testuser', password='test')

        unavailable_response = self.c.get(
            url,
            HTTP_X_REQUESTED_WITH='XMLHttpRequest',
        )
        self.assertEqual(409, unavailable_response.status_code)
        self.assertEqual('application/json', unavailable_response['Content-Type'].split(';')[0])
        self.assertEqual(
            {'error': 'Project execution status is unavailable'},
            unavailable_response.json(),
        )

        remote_project.r_project_execution_snakemake_task = TaskResult.objects.create(
            status='PROGRESS',
            task_id='remote-progress-processing-failure',
        )
        remote_project.save()
        with patch('blast_project.views.read_task_logs_summary_table', side_effect=Exception):
            failure_response = self.c.get(
                url,
                HTTP_X_REQUESTED_WITH='XMLHttpRequest',
            )

        self.assertEqual(500, failure_response.status_code)
        self.assertEqual('application/json', failure_response['Content-Type'].split(';')[0])
        self.assertEqual({'error': 'Unable to read progress'}, failure_response.json())

    @tag('fast', 'view')
    def test_ajax_call_to_remote_logfiles_redirects_unauthenticated_requests(self):
        remote_project = self._create_remote_project(status='SUCCESS')
        self.c.logout()

        response = self.c.get(
            reverse('ajax_call_to_remote_logfiles', kwargs={'project_id': remote_project.id}),
            HTTP_X_REQUESTED_WITH='XMLHttpRequest',
        )

        self.assertEqual(302, response.status_code)
        self.assertIn(reverse('login'), response['Location'])

    @tag('long','view','project_creation_view')
    #two test stages 1. post request testing 2. database testing
    #throws exception if the project directory does exist! --> IntegrityError
    def test_project_creation_view_post(self):
        blast_database = BlastDatabase.objects.get(database_name='Curvibacter sp. aep1-3 database')
        upload_file = open('testfiles/testsequences/lps_transport.faa', 'rb')
        post_dict = {'project_title': 'Test Project Title 1',
                     'species_name_for_backward_blast':'Curvibacter sp. AEP1-3',
                     'user_email':'lukas.becker@hhu.de',
                     'project_forward_database':blast_database.pk,
                     'project_backward_database':blast_database.pk,
                     'query_sequence_file': SimpleUploadedFile(upload_file.name, upload_file.read()),
                     'fw_e_value':0.001,
                     'fw_word_size':3,
                     'fw_num_alignments':10000,
                     'fw_max_target_seqs':10000,
                     'fw_num_threads':1,
                     'fw_max_hsps':500,
                     'bw_e_value': 0.001,
                     'bw_word_size': 3,
                     'bw_num_alignments': 10000,
                     'bw_max_target_seqs': 10000,
                     'bw_num_threads': 1,
                     'bw_max_hsps': 500,
                     }

        #file_dict = {'query_sequence_file': SimpleUploadedFile(upload_file.name, upload_file.read())}
        if check_if_taxdb_exists():
            print("[+] Logging in as testuser")
            self.c.login(username='testuser', password='test')
            print("[+] Send POST request to sever")
            response = self.c.post('/blast_project/project_creation',post_dict)
            sleep(1)
            print(response)
            print(BlastProject.objects.all())
            print("[+] Load new blast project from database and delete static and project directories")
            new_project = BlastProject.objects.get(project_title='Test Project Title 1')
            if isdir(new_project.get_project_dir() + '/'):
                rmtree(new_project.get_project_dir() + '/')
            if isdir('static/images/result_images/'+str(new_project.id)):
                rmtree('static/images/result_images/'+str(new_project.id))

            self.assertEqual(302,response.status_code)
            self.assertEqual(new_project.project_title, 'Test Project Title 1')
            self.assertRedirects(response, '/blast_project/'+str(new_project.id)+'/project_details', status_code=302,
                                 target_status_code=200, fetch_redirect_response=True)
        else:
            #taxid does not exist -> the test function would trigger a download process of the
            #taxonomy db please consider to restart the server application or to download the taxonomy db manually (into the database directory)
            print("[-] There is no taxonomy database loaded, skipping this test ...")
            self.assertTrue(False)

    @tag('fast','view','delete_project_view')
    def test_delete_project_view_get(self):
        project = BlastProject.objects.get(project_title='test project 1')
        self.c.login(username='testuser', password='test')
        response = self.c.post('/blast_project/'+str(project.id)+'/project_deletion')
        self.assertEqual(200,response.status_code)
        self.assertTemplateUsed(response,'blast_project/success.html')

    @tag('long','view','view_combination')
    def test_create_new_project_and_delete_it(self):
        blast_database = BlastDatabase.objects.get(
            database_name='Curvibacter sp. aep1-3 database')
        upload_file = open('testfiles/testsequences/lps_transport.faa', 'rb')
        post_dict = {'project_title': 'Test Project Title 2',
                     'species_name_for_backward_blast':'Curvibacter sp. AEP1-3',
                     'user_email':'lukas.becker@hhu.de',
                     'project_forward_database':blast_database.pk,
                     'project_backward_database':blast_database.pk,
                     'query_sequence_file':
                         SimpleUploadedFile(upload_file.name, upload_file.read()),
                     'fw_e_value':0.001,
                     'fw_word_size':3,
                     'fw_num_alignments':10000,
                     'fw_max_target_seqs':10000,
                     'fw_num_threads':1,
                     'fw_max_hsps':500,
                     'bw_e_value': 0.001,
                     'bw_word_size': 3,
                     'bw_num_alignments': 10000,
                     'bw_max_target_seqs': 10000,
                     'bw_num_threads': 1,
                     'bw_max_hsps': 500,
                     }
        #file_dict = {'query_sequence_file':
        # SimpleUploadedFile(upload_file.name, upload_file.read())}
        if check_if_taxdb_exists():
            number_project_before = len(BlastProject.objects.all())
            self.c.login(username='testuser', password='test')
            response_post = self.c.post(
                '/blast_project/project_creation',post_dict)
            sleep(2)
            for project in BlastProject.objects.all():
                print("[+] Project Name: {}".format(project.project_title))

            new_project = BlastProject.objects.get(
                project_title='Test Project Title 2')

            number_project_before_deletion = len(BlastProject.objects.all())
            response_delete = self.c.post(
                '/blast_project/'+str(new_project.id)+'/project_deletion',
                {'project_id':new_project.id})
            number_project_after_deletion = len(BlastProject.objects.all())

            for project in BlastProject.objects.all():
                print("[+] Project Name: {}".format(project.project_title))

            self.assertEqual(200,response_delete.status_code)
            self.assertEqual(302,response_post.status_code)
            self.assertTrue(number_project_before < number_project_before_deletion)
            self.assertTrue(number_project_after_deletion == number_project_before)
            self.assertEqual(new_project.project_title, 'Test Project Title 2')
            self.assertRedirects(response_post,
                                 '/blast_project/'+str(new_project.id)+'/project_details',
                                 status_code=302,target_status_code=200,
                                 fetch_redirect_response=True)
        else:
            #taxid does not exist and would be downloaded by test
            print("[-] There is no taxonomy database loaded, skipping this test ...")
            self.assertTrue(False)


