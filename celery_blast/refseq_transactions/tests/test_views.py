from time import sleep
from unittest.mock import patch

from blast_project.models import BlastProject, BlastSettings
# from blast_project.py_django_db_services import get_all_succeeded_databases
# from django.core.files.uploadedfile import SimpleUploadedFile
# from blast_project.py_services import check_if_taxdb_exists
# from shutil import rmtree
# from os.path import isdir
from blast_project.py_services import delete_blastdb_and_associated_directories_by_id
from django.contrib.auth.models import User
from django.test import Client, TestCase, tag
from django.utils import timezone
from django_celery_results.models import TaskResult
from refseq_transactions.models import AssemblyLevels, BlastDatabase


class RefseqTransactionsViewsTestCase(TestCase):
    c = Client()

    def tearDown(self) -> None:
        super().tearDown()
        print("[+] Tear down refseqtrasactions celery test data")
        print("\t[*] Trying to delete test directories and models ...")
        for database in BlastDatabase.objects.all():
            print("\t[*] Deleting database: {}".format(database.database_name))
            delete_blastdb_and_associated_directories_by_id(database.id)

    def setUp(self) -> None:
        super().setUp()
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
            database_name='Curvibacter sp. aep1-3 database',
            database_description='Full Assembly Content',
            assembly_entries=1,
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

    @tag("refseq_transactions_views")
    def test_dashboard_view(self):
        self.c.login(username='testuser', password='test')
        response = self.c.get("/refseq_transactions/")
        self.assertEqual(200, response.status_code)
        self.assertTemplateUsed(response, 'refseq_transactions/refseq_transactions_dashboard.html')

    @tag("refseq_transactions_views")
    def test_invalid_metadata_post_keeps_metadata_form_when_summary_files_exist(self):
        self.c.login(username='testuser', password='test')
        post_dict = {
            "database_name": "invalid metadata post",
            "database_description": "missing database summary file",
            "taxid_uploaded_file": "",
            "taxid_text_field": "",
        }

        with patch('refseq_transactions.views.refseq_file_exists', return_value=True), \
                patch('refseq_transactions.views.genbank_file_exists', return_value=True), \
                patch('refseq_transactions.views.sleep', return_value=None):
            response = self.c.post('/refseq_transactions/create_refseq_database_metadata', data=post_dict)

        self.assertEqual(200, response.status_code)
        self.assertTemplateUsed(response, 'refseq_transactions/refseq_transactions_dashboard.html')
        self.assertTrue(response.context['refseq_exists'])
        self.assertTrue(response.context['genbank_exists'])
        self.assertIn('InactiveBlastDatabases', response.context)
        self.assertIn('FailedDatabases', response.context)
        self.assertContains(response, 'id="post_form_container"', html=False)
        self.assertContains(response, 'BLAST database Metadata')
        self.assertContains(response, 'This field is required')

    @tag("refseq_transactions_views")
    def test_summary_download_redirects_to_progress_page(self):
        self.c.login(username='testuser', password='test')

        class Result:
            id = 'summary-download-task'

        with patch('refseq_transactions.views.download_refseq_assembly_summary.delay', return_value=Result()) as delay:
            response = self.c.post('/refseq_transactions/RefSeq/download_refseq_assembly_summary')

        delay.assert_called_once_with('RefSeq')
        self.assertRedirects(
            response,
            '/refseq_transactions/RefSeq/summary-download-task/assembly_summary_download_progress',
            status_code=302,
            target_status_code=200,
            fetch_redirect_response=True,
        )

    @tag("refseq_transactions_views")
    def test_summary_download_progress_ajax_reports_success(self):
        self.c.login(username='testuser', password='test')
        TaskResult.objects.create(
            task_id='summary-success',
            status='SUCCESS',
            result='{"progress": {"percent": 100, "description": "SUCCESS"}}',
        )

        response = self.c.get('/refseq_transactions/summary-success/assembly_summary_download_progress_ajax')

        self.assertEqual(200, response.status_code)
        self.assertEqual(100, response.json()['progress'])
        self.assertTrue(response.json()['complete'])

    @tag("refseq_transactions_views")
    def test_database_preview_creation_redirects_to_progress_page(self):
        self.c.login(username='testuser', password='test')

        class Result:
            id = 'database-preview-task'

        post_dict = {
            "assembly_levels": ["Complete Genome", "Chromosome"],
            "database_summary_file": "RefSeq",
            "database_name": "database preview async",
            "database_description": "database preview async test",
            "taxid_uploaded_file": "",
            "taxid_text_field": "",
            "user_email": "test_email@email.com",
        }

        with patch('refseq_transactions.views.create_blast_database_preview.delay', return_value=Result()) as delay, \
                patch('refseq_transactions.forms.check_if_database_title_exists', return_value=False):
            response = self.c.post('/refseq_transactions/create_refseq_database_metadata', data=post_dict)

        delay.assert_called_once()
        self.assertEqual('database preview async', delay.call_args.args[0]['database_name'])
        self.assertRedirects(
            response,
            '/refseq_transactions/database-preview-task/database_preview_creation_progress',
            status_code=302,
            target_status_code=200,
            fetch_redirect_response=True,
        )

    @tag("refseq_transactions_views")
    def test_database_preview_creation_progress_ajax_reports_success(self):
        self.c.login(username='testuser', password='test')
        TaskResult.objects.create(
            task_id='database-preview-success',
            status='SUCCESS',
            result='7',
        )

        response = self.c.get('/refseq_transactions/database-preview-success/database_preview_creation_progress_ajax')

        self.assertEqual(200, response.status_code)
        self.assertEqual(100, response.json()['progress'])
        self.assertTrue(response.json()['complete'])

    @tag("refseq_transactions_views")
    def test_create_and_delete_blast_database_model_and_directory(self):
        self.c.login(username='testuser', password='test')
        post_dict = {'assembly_levels': ["Complete Genome", "Chromosome"],
                     "database_name": "test database curvibacter 2",
                     "database_description": "database for testing purposes",
                     "taxid_uploaded_file": "curvibacter_hydra.taxids"
                     }

        # this should trigger the "else" branche of the py_refseq_transactions.py function: create_blast_database_model_and_directory
        response = self.c.post('/refseq_transactions/create_refseq_database_metadata', data=post_dict)
        print("\t[*] Sleeping now ...")
        sleep(2)
        new_database = BlastDatabase.objects.get(database_name='test database curvibacter 2')
        print("\t[*] New Database: ", new_database.database_name)
        response_delete = self.c.post('/refseq_transactions/{}/delete_database'.format(new_database.id))
        print("\t[*] Deleting database ...")
        self.assertTrue(200, response.status_code)
        self.assertTrue(302, response_delete.status_code)
        self.assertRedirects(response,
                             '/refseq_transactions/',
                             status_code=302, target_status_code=200,
                             fetch_redirect_response=True)
        self.assertRedirects(response_delete, '/refseq_transactions/', status_code=302, target_status_code=200)

    '''test_create_blast_database_model_and_directory

    functions that gets executed during this POST request:

    filtered_table = read_current_assembly_summary_with_pandas(assembly_levels)
    filtered_table = filter_duplicates_by_ftp_path(filtered_table)
    new_blastdb = create_and_save_refseq_database_model(
        database_name=database_name,
        database_description=database_description,
        assembly_levels=assembly_levels,
        assembly_entries=len(filtered_table))

    refseq_database_table_path = create_blastdatabase_directory(new_blastdb.id)

    write_pandas_table_to_project_dir(refseq_database_table_path,
                                      filtered_table,
                                      database_name)

    '''

    @tag("refseq_transactions_views")
    def test_create_blast_database_model_and_directory(self):
        print(BlastDatabase.objects.all())
        self.c.login(username='testuser', password='test')
        post_dict = {'assembly_levels': ["Complete Genome", "Chromosome"],
                     "database_name": "test database curvibacter 1",
                     "database_description": "database for testing purposes 1",
                     "taxid_uploaded_file": "curvibacter_hydra.taxids"
                     }

        # this should trigger the "else" branche of the py_refseq_transactions.py function: create_blast_database_model_and_directory
        response = self.c.post('/refseq_transactions/create_refseq_database_metadata', data=post_dict)

        print("\t[*] Sleeping now ...")
        sleep(2)
        new_database = BlastDatabase.objects.get(database_name='test database curvibacter 1')
        print("\t[*] New Database: ", new_database.database_name)

        self.assertTrue(new_database.database_description == "database for testing purposes 1")
        self.assertTrue(200, response.status_code)
        self.assertRedirects(response,
                             '/refseq_transactions/',
                             status_code=302, target_status_code=200,
                             fetch_redirect_response=True)
