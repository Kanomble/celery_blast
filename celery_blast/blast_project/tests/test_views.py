#cmd for running these tests: python manage.py test blast_project.tests.test_views
#different tags: long fast view view_combination
from django.test import Client, TestCase, tag
from django.contrib.auth.models import User
from blast_project.models import BlastProject, BlastSettings
from refseq_transactions.models import AssemblyLevels, BlastDatabase
from django.utils import timezone
from django_celery_results.models import TaskResult
#from blast_project.py_django_db_services import get_all_succeeded_databases
from django.core.files.uploadedfile import SimpleUploadedFile
from blast_project.py_services import check_if_taxdb_exists
from shutil import rmtree
from os.path import isdir
from time import sleep

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
        response = self.c.get('/blast_project/')
        self.assertEqual(200,response.status_code)
        self.assertTemplateUsed(response,'blast_project/blast_project_dashboard.html')
        #response.context['blast_projects'] is a django.db.models.query.QuerySet class
        self.assertEqual(response.context['blast_projects'][0],
                         project)
        self.assertEqual(response.context['ActiveBlastDatabases'][0],
                         project.project_forward_database)

    @tag('fast','view')
    def test_project_details_get(self):
        project = BlastProject.objects.get(project_title='test project 1')
        self.c.login(username='testuser', password='test')
        response = self.c.get('/blast_project/'+str(project.id)+'/project_details')
        self.assertEqual(200,response.status_code)
        self.assertEqual(response.context['BlastProject'],project)

    @tag('long','view')
    #two test stages 1. post request testing 2. database testing
    #throws exception if the project directory does exist! --> IntegrityError
    def test_project_creation_view_post(self):
        blast_database = BlastDatabase.objects.get(database_name='Curvibacter sp. aep1-3 database')
        upload_file = open('testfiles/testsequences/lps_transport.faa', 'rb')
        post_dict = {'project_title': 'Test Project Title',
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
            self.c.login(username='testuser', password='test')
            response = self.c.post('/blast_project/project_creation',post_dict)
            new_project = BlastProject.objects.get(project_title='Test Project Title')
            if isdir(new_project.get_project_dir() + '/'):
                rmtree(new_project.get_project_dir() + '/')
            if isdir('static/images/result_images/'+str(new_project.id)):
                rmtree('static/images/result_images/'+str(new_project.id))

            self.assertEqual(302,response.status_code)
            self.assertEqual(new_project.project_title, 'Test Project Title')
            self.assertRedirects(response, '/blast_project/'+str(new_project.id)+'/project_details', status_code=302,
                                 target_status_code=200, fetch_redirect_response=True)
        else:
            #taxid does not exist -> the test function would trigger a download process of the
            #taxonomy db please consider to restart the server application or to download the taxonomy db manually (into the database directory)
            print("[-] There is no taxonomy database loaded, skipping this test ...")
            self.assertTrue(False)

    @tag('fast','view')
    def test_delete_project_view_get(self):
        project = BlastProject.objects.get(project_title='test project 1')
        self.c.login(username='testuser', password='test')
        response = self.c.get('/blast_project/'+str(project.id)+'/project_deletion')
        self.assertEqual(200,response.status_code)
        self.assertTemplateUsed(response,'blast_project/project_details_dashboard.html')

    @tag('long','view','view_combination')
    def test_create_new_project_and_delete_it(self):
        blast_database = BlastDatabase.objects.get(
            database_name='Curvibacter sp. aep1-3 database')
        upload_file = open('testfiles/testsequences/lps_transport.faa', 'rb')
        post_dict = {'project_title': 'Test Project Title',
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
            new_project = BlastProject.objects.get(
                project_title='Test Project Title')
            number_project_before_deletion = len(BlastProject.objects.all())
            response_delete = self.c.delete(
                '/blast_project/'+str(new_project.id)+'/project_deletion',
                {'project_id':new_project.id})
            number_project_after_deletion = len(BlastProject.objects.all())
            self.assertEqual(200,response_delete.status_code)
            self.assertEqual(302,response_post.status_code)
            self.assertTrue(number_project_before < number_project_before_deletion)
            self.assertTrue(number_project_after_deletion == number_project_before)
            self.assertEqual(new_project.project_title, 'Test Project Title')
            self.assertRedirects(response_post,
                                 '/blast_project/'+str(new_project.id)+'/project_details',
                                 status_code=302,target_status_code=200,
                                 fetch_redirect_response=True)
        else:
            #taxid does not exist and would be downloaded by test
            print("[-] There is no taxonomy database loaded, skipping this test ...")
            self.assertTrue(False)

    @tag('fast','view')
    def test_uploadgenome_view_get(self):
        self.c.login(username='testuser', password='test')
        response = self.c.get('/blast_project/upload_genomes/')
        self.assertEqual(200, response.status_code)

    #TODO fix celery database handling for tests ...
    @tag('upload_genomes','special')
    def test_uploadgenome_view_post(self):
        self.c.login(username='testuser', password='test')
        genome_file = open(
            'testfiles/upload_genomes/concatenated_genome_files/other_symbionts.faa', 'rb')
        organism_names = open(
            'testfiles/upload_genomes/concatenated_genome_files/organisms.txt', 'rb')
        assembly_accessions = open(
            'testfiles/upload_genomes/concatenated_genome_files/accessions.txt', 'rb')
        assembly_levels = open(
            'testfiles/upload_genomes/concatenated_genome_files/levels.txt', 'rb')
        taxmap_file = open(
            'testfiles/upload_genomes/concatenated_genome_files/acc_map.tab', 'rb')

        post_dict={'database_title':"All Hydra Symbionts Test",
                   'database_description':"Collection Of Symbiont Genomes",
                   'assembly_entries':5,
                   'user_email':'lukas.becker@hhu.de',
                   'genome_fasta_file':
                       SimpleUploadedFile(genome_file.name, genome_file.read()),
                   'organism_name_file':
                       SimpleUploadedFile(organism_names.name, organism_names.read()),
                   'assembly_accessions_file':
                       SimpleUploadedFile(assembly_accessions.name, assembly_accessions.read()),
                   'assembly_level_file':
                       SimpleUploadedFile(assembly_levels.name, assembly_levels.read()),
                   'taxmap_file':
                       SimpleUploadedFile(taxmap_file.name, taxmap_file.read())
                   }
        response = self.c.post('/blast_project/upload_genomes/',data=post_dict)

        print("[+] ... Starting to sleep for 20 seconds to wait for makeblastdb ...")
        #sleep(20)

        self.assertTrue(200,response.status_code)
        self.assertTemplateUsed(response,'blast_project/success.html')