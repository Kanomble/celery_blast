from django.test import TestCase
#this is the correct import
from blast_project.models import BlastProject, BlastSettings
from refseq_transactions.models import AssemblyLevels, BlastDatabase
from django.contrib.auth.models import User
from django.utils import timezone
import os

class BlastProjectTestCase(TestCase):
    def setUp(self):

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

        fw_and_bw_db = BlastDatabase.objects.create(
            database_name = 'Test Database 1',
            database_description = 'Full Assembly Content',
            assembly_entries = 2
        )

        BlastProject.objects.create(
            project_title='test project 1',
            search_strategy='blastp',
            project_query_sequences='lps_transport.faa',
            project_user = user,
            project_forward_settings = fw_settings,
            project_backward_settings = bw_settings,
            project_forward_database = fw_and_bw_db,
            project_backward_database = fw_and_bw_db,
            species_name_for_backward_blast = 'Curvibacter sp. AEP1-3'
        )

    def test_blast_database_model(self):
        blast_database = BlastDatabase.objects.get(database_name='Test Database 1')
        self.assertEqual('Test Database 1',blast_database.database_name)

    def test_blast_database_pandas_table_name(self):
        blast_database = BlastDatabase.objects.get(database_name='Test Database 1')
        self.assertEqual(blast_database.get_pandas_table_name(),'TEST_DATABASE_1')

    def test_blast_database_assembly_level_creation(self):
        blast_database = BlastDatabase.objects.get(database_name='Test Database 1')
        assembly_level = ['Complete Genome', 'Chromosome', 'Contig', 'Scaffold']
        assembly_levels = AssemblyLevels.objects.filter(assembly_level__in=assembly_level)
        for assembly_lvl in assembly_levels:
            blast_database.assembly_levels.add(assembly_lvl)

        level1 = blast_database.assembly_levels.get(assembly_level='Complete Genome')

        self.assertEqual(assembly_levels.count(),4)
        self.assertEqual(level1.assembly_level,'Complete Genome')

    def test_blast_project_model(self):
        blast_project = BlastProject.objects.get(project_title='test project 1')
        self.assertEqual(blast_project.project_title,'test project 1')

    def test_blast_project_query_sequence_file(self):
        blast_project = BlastProject.objects.get(project_title='test project 1')
        self.assertEqual(blast_project.get_list_of_query_sequences(filepath='testfiles/blast_project/'),
                         ['WP_087495344',
                          'WP_087495343',
                          'WP_087497333',
                          'WP_087495837',
                          'WP_087496015',
                          'WP_087493791',
                          'WP_087493790',
                          'WP_087493929'])

    def test_blast_project_write_sanekmake_config_file(self):
        blast_project = BlastProject.objects.get(project_title='test project 1')

        if os.path.isdir('testfiles/blast_project/'+str(blast_project.id)) == False:
            os.mkdir('testfiles/blast_project/'+str(blast_project.id))
        blast_project.write_snakemake_configuration_file(filepath='testfiles/blast_project/')
        with open('testfiles/blast_project/'+str(blast_project.id) +'/snakefile_config') as sconfig:
            lines = sconfig.readlines()
            config_dict = {}
            for line in lines:
                line = line.split(":")
                config_dict[line[0]] = line[1]
        self.assertEqual(len(config_dict.keys()),18)