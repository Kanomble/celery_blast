from time import sleep

from celery.contrib.testing.worker import start_worker
from celery_blast.celery import app
from django.contrib.auth.models import User
from django.test import TransactionTestCase, tag, Client
from django.utils import timezone
from refseq_transactions.models import AssemblyLevels, BlastDatabase


@tag('refseqtransactions_celery_tests')
class RefseqTransactionsTestCase(TransactionTestCase):
    c = Client()

    @classmethod
    def setUpClass(cls):
        print("[+] Starting celery worker process")
        super().setUpClass()
        cls.celery_worker = start_worker(app, perform_ping_check=False)
        cls.celery_worker.__enter__()

    def setUp(self) -> None:
        super().setUp()
        print("[+] Setting up model data for refseqtransactions celery tests")
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

    @classmethod
    def tearDownClass(cls):
        print("[+] Tear down refseqtransactions test celery worker")
        super().tearDownClass()
        cls.celery_worker.__exit__(None, None, None)

    def tearDown(self) -> None:
        super().tearDown()
        print("[+] Tear down refseqtrasactions celery test data")
        print("\t[*] Trying to delete test directories and models ...")
        for database in BlastDatabase.objects.all():
            print("\t[*] Deleting database: {}".format(database.database_name))
            # delete_blastdb_and_associated_directories_by_id(database.id)

    @tag("refseq_transaction_celery_tasks")
    def test_download_and_format_refseq_database_task(self):
        self.c.login(username='testuser', password='test')
        post_dict = {'assembly_levels': ["Complete Genome", "Chromosome"],
                     "database_name": "test database curvibacter 2",
                     "database_description": "database for testing purposes",
                     "taxid_uploaded_file": "curvibacter_hydra.taxids"
                     }

        # this should trigger the "else" branche of the py_refseq_transactions.py function: create_blast_database_model_and_directory
        response_post = self.c.post('/refseq_transactions/create_refseq_database_metadata', data=post_dict)
        print("\t[*] Sleeping now ...")
        sleep(2)
        new_database = BlastDatabase.objects.get(database_name='test database curvibacter 2')
        print("\t[*] New Database: ", new_database.database_name)

        print("[*] Trying to download and format BLAST Database")
        response_download_task = self.c.post(
            "/refseq_transactions/{}/download_and_format_blast_database".format(new_database.id))
        print("\t[*] Sleeping 15s ...")
        sleep(15)
        # path_to_database = new_database.get_pandas_table_name() + '.database.pal'
        # out = check_output(['blastdbcheck', '-db', path_to_database])
        # print(out)
        self.assertTrue(new_database.database_description == "database for testing purposes 1")
        self.assertTrue(200, response_post.status_code)
        self.assertRedirects(response_post,
                             '/refseq_transactions/',
                             status_code=302, target_status_code=200,
                             fetch_redirect_response=True)
        self.assertTrue(200, response_download_task.status_code)
        self.assertRedirects(response_download_task.status_code,
                             '/refseq_transactions/',
                             status_code=302, target_status_code=200,
                             fetch_redirect_response=True)
