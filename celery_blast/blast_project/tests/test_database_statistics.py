from django.test import TestCase, tag
import pandas as pd
from django.contrib.auth.models import User
from blast_project.models import BlastProject, BlastSettings
from refseq_transactions.models import AssemblyLevels, BlastDatabase
from django.utils import timezone
from django_celery_results.models import TaskResult

from blast_project.py_database_statistics import add_taxonomic_information_to_db, \
                                                extract_taxonomic_information, calculate_database_statistics
from os.path import isfile
from django.conf import settings

class DatabaseStatisticsTest(TestCase):
    def setUp(self) -> None:
        super().setUp()
        print("[+] Setting up DatabaseStatistics test data")
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

    path_to_test_files = 'testfiles/database_statistics/'
    path_to_log_files = 'testfiles/database_statistics/log/'

    def test_add_taxonomic_information_to_db(self):
        logfile_path = self.path_to_log_files + 'add_taxonomic_information_to_db.log'
        database_dataframe = pd.read_csv(self.path_to_test_files + 'HIGH_QUALITY_CYANOBACTERIA_DATABASE', index_col=0)
        taxonomy_df = add_taxonomic_information_to_db('lukas.becker@hhu.de',
                                                      logfile_path,
                                                      list(database_dataframe['taxid'].unique()))

        self.assertTrue(isfile(logfile_path))
        self.assertTrue(len(taxonomy_df) == len(list(database_dataframe['taxid'].unique())))

    def test_extract_taxonomic_information(self):
        logfile_path = self.path_to_log_files + 'extract_taxonomic_information.log'
        database_dataframe = pd.read_csv(self.path_to_test_files + 'HIGH_QUALITY_CYANOBACTERIA_DATABASE', index_col=0)
        reciprocal_results_dataframe = pd.read_csv(self.path_to_test_files + 'reciprocal_results_with_taxonomy.csv', index_col=0)
        tax_counts = extract_taxonomic_information(logfile_path, False,
                                                   reciprocal_results_dataframe,database_dataframe,'order')

        for query_list in tax_counts:
            if query_list.name == "WP_010874243":
                synechococcus = query_list['Synechococcales']
                nostocales = query_list['Nostocales']
                chroococcales = query_list['Chroococcales']
        self.assertTrue(len(tax_counts) == 3)
        self.assertTrue(type(tax_counts) == list)
        self.assertTrue(type(tax_counts[0]) == pd.core.series.Series)
        self.assertTrue(synechococcus == 92)
        self.assertTrue(nostocales == 66)
        self.assertTrue(chroococcales == 24)

    def test_calculate_database_statistics(self):
        setattr(settings, 'BLAST_PROJECT_DIR', 'testfiles/blast_project/')
        setattr(settings, 'BLAST_DATABASE_DIR', 'testfiles/databases/')
        print(settings.BLAST_PROJECT_DIR)
        print(settings.BLAST_DATABASE_DIR)
        project = BlastProject.objects.get(project_title="kaiABC operon synechocystis sp. PCC6803 vs cyanobacteria")
        database = BlastDatabase.objects.get(database_name='high quality cyanobacteria database')
        print()
        print("+++++++++++++++++++++++++++++++++++++++")
        print(project.get_project_query_sequence_filepath())
        print(database.path_to_database_file)
        print()
        #(project_id: int,logfile:str,user_email:str, taxonomic_units:list)
        return_val = calculate_database_statistics(project.id,'lukas.becker@hhu.de','testfiles/blast_project/2/log/database_statistcis_test.log',['phylum','class','order'])
        self.assertTrue(database.path_to_database_file == 'testfiles/databases/2')
        self.assertTrue(project.get_project_query_sequence_filepath(),'testfiles/blast_project/2/kaiABC.faa')
        self.assertTrue(return_val == 0)
