from django.test import TestCase
#this is the correct import
from blast_project.models import BlastProject, BlastDatabase, BlastSettings, AssemblyLevels
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
            project_query_sequences='testfiles/query_sequences/lps_transport.faa',
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