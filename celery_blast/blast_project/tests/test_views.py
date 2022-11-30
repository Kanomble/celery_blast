#cmd for running these tests: python manage.py test blast_project.tests.test_views
from django.test import Client
from django.test import TestCase
from django.contrib.auth.models import User
from blast_project.models import BlastProject, BlastSettings
from refseq_transactions.models import AssemblyLevels, BlastDatabase
from django.utils import timezone
from django_celery_results.models import TaskResult

class BlastProjectViewsTestCase(TestCase):
    c = Client()
    def setUp(self) -> None:
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

    def test_register_get(self):
        response = self.c.get('/blast_project/register')
        self.assertEqual(200,response.status_code)
        self.assertTemplateUsed(response, template_name='blast_project/register.html')

    def test_register_post(self):
        response = self.c.post('/blast_project/register',
                               {'username':'testuser2',
                                'password':'nvm',
                                'email':'testuser2@email.de'})
        self.assertEqual(200,response.status_code)
        self.assertTemplateUsed(response, template_name='blast_project/register.html')

    def test_login_get(self):
        response = self.c.get('/blast_project/login')
        self.assertEqual(200,response.status_code)
        self.assertTemplateUsed(response, template_name='blast_project/login.html')

    def test_login_post(self):
        response = self.c.post('/blast_project/login',{'username':'testuser','password':'test'})
        self.assertEqual(302,response.status_code)
        self.assertRedirects(response, '/blast_project/',status_code=302,
                             target_status_code=200,fetch_redirect_response=True)

    def test_dashboard_view_get_not_logged_in(self):
        response = self.c.get('/blast_project/')
        self.assertEqual(302,response.status_code)
        self.assertRedirects(response, '/blast_project/login?next=/blast_project/',status_code=302,
                             target_status_code=200,fetch_redirect_response=True)

    def test_dashboard_view_get(self):
        self.c.login(username='testuser',password='test')
        project = BlastProject.objects.get(project_title='test project 1')
        response = self.c.get('/blast_project/')
        self.assertEqual(200,response.status_code)
        self.assertTemplateUsed(response,'blast_project/blast_project_dashboard.html')
        #response.context['blast_projects'] is a django.db.models.query.QuerySet class
        self.assertEqual(response.context['blast_projects'][0],
                         project)
        self.assertEqual(response.context['ActiveBlastDatabases'][0],
                         project.project_forward_database)

    def test_project_details_get(self):
        project = BlastProject.objects.get(project_title='test project 1')
        self.c.login(username='testuser', password='test')
        response = self.c.get('/blast_project/'+str(project.id)+'/project_details')
        self.assertEqual(200,response.status_code)
        self.assertEqual(response.context['BlastProject'],project)