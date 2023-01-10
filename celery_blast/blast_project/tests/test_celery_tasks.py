from django.test import TransactionTestCase, tag, Client
from celery.contrib.testing.worker import start_worker
from django.core.files.uploadedfile import SimpleUploadedFile
from django.contrib.auth.models import User
from blast_project.models import BlastProject, BlastSettings
from refseq_transactions.models import AssemblyLevels, BlastDatabase
from django.utils import timezone
from django_celery_results.models import TaskResult
from celery_blast.celery import app
from time import sleep
from django.conf import settings

@tag('special')
class GenomeUploadViewsTestCase(TransactionTestCase):
    c = Client()

    @classmethod
    def setUpClass(cls):
        print("[+] Starting celery worker process")
        super().setUpClass()
        cls.celery_worker = start_worker(app,perform_ping_check=False)
        cls.celery_worker.__enter__()

    @classmethod
    def tearDownClass(cls):
        print("[+] Tear down genomeuploadview celery test data")
        super().tearDownClass()
        cls.celery_worker.__exit__(None,None,None)

    def setUp(self) -> None:
        super().setUp()
        print("[+] Setting up genomeuploadview celery test data")
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
    def test_uploadgenome_view_get(self):
        self.c.login(username='testuser', password='test')
        response = self.c.get('/blast_project/upload_genomes/')
        self.assertEqual(200, response.status_code)

    @tag('upload_genomes','special')
    def test_uploadgenome_view_post(self):
        print("[+] Setting up test for the uploadgenome view")
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
        print("[+] Sending POST request to server ...")
        response = self.c.post('/blast_project/upload_genomes/',data=post_dict)
        print("[+] response: {}".format(response.status_code))
        sleep(5)
        print("[+] Loading new database object from database")
        database = BlastDatabase.objects.get(database_name='All Hydra Symbionts Test')
        print("[+] Celery Task For Database Formatting: ", database.database_download_and_format_task.status)
        sleep(10)
        print("[+] Reloading object ...")
        database = BlastDatabase.objects.get(database_name='All Hydra Symbionts Test')
        print("[+] Celery Task For Database Formatting: ", database.database_download_and_format_task.status)
        self.assertTrue(200,response.status_code)
        self.assertTemplateUsed(response,'blast_project/success.html')

@tag('database_statistics_task')
class DatabaseStatisticsTransactionTestCase(TransactionTestCase):
    c = Client()

    @classmethod
    def setUpClass(cls):
        print("[+] Starting celery worker process")
        super().setUpClass()
        cls.celery_worker = start_worker(app, perform_ping_check=False)
        cls.celery_worker.__enter__()

    @classmethod
    def tearDownClass(cls):
        print("[+] Tear down genomeuploadview celery test data")
        super().tearDownClass()
        cls.celery_worker.__exit__(None, None, None)

    def setUp(self) -> None:
        super().setUp()
        print("[+] Setting up genomeuploadview celery test data")
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
            database_name='high quality cyanobacteria database',
            database_description='-',
            assembly_entries=251,
            path_to_database_file="testfiles/databases/2",
            database_download_and_format_task=celery_task
        )

        assembly_level = ['Complete Genome', 'Chromosome', 'Contig', 'Scaffold']
        assembly_levels = AssemblyLevels.objects.filter(assembly_level__in=assembly_level)
        for assembly_lvl in assembly_levels:
            fw_and_bw_db.assembly_levels.add(assembly_lvl)

        fw_and_bw_db.save()

        BlastProject.objects.create(
            project_title='kaiABC operon synechocystis sp. PCC6803 vs cyanobacteria',
            search_strategy='blastp',
            project_query_sequences='kaiABC.faa',
            project_user=user,
            project_forward_settings=fw_settings,
            project_backward_settings=bw_settings,
            project_forward_database=fw_and_bw_db,
            project_backward_database=fw_and_bw_db,
            species_name_for_backward_blast='synechocystis sp. PCC6803'
        )

    def test_celery_database_statistics_task(self):
        setattr(settings, 'BLAST_PROJECT_DIR', 'testfiles/blast_project/')
        setattr(settings, 'BLAST_DATABASE_DIR', 'testfiles/databases/')
        print(settings.BLAST_PROJECT_DIR)
        print(settings.BLAST_DATABASE_DIR)
