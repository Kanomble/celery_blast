from django.test import TestCase, tag
#this is the correct import
from blast_project.models import BlastProject, BlastSettings, RemoteBlastProject
from refseq_transactions.models import AssemblyLevels, BlastDatabase
from django.contrib.auth.models import User
from django.utils import timezone
import os
from pathlib import Path
import tempfile

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

        remote_fw_settings = BlastSettings.objects.create(
            e_value=0.001,
            word_size=3,
            num_alignments=10000,
            max_target_seqs=10000,
            num_threads=1,
            max_hsps=500
        )

        remote_bw_settings = BlastSettings.objects.create(
            e_value=0.001,
            word_size=3,
            num_alignments=1,
            max_target_seqs=1,
            num_threads=1,
            max_hsps=500
        )

        RemoteBlastProject.objects.create(
            r_project_title='test remote project 1',
            r_search_strategy='blastp',
            r_project_query_sequences='lps_transport.faa',
            r_project_user=user,
            r_project_forward_settings=remote_fw_settings,
            r_project_backward_settings=remote_bw_settings,
            r_project_forward_database='nr',
            r_project_backward_database=fw_and_bw_db,
            r_species_name_for_backward_blast='Curvibacter sp. AEP1-3'
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

    def test_blast_project_missing_reciprocal_information_table_returns_empty_dict(self):
        blast_project = BlastProject.objects.get(project_title='test project 1')
        self.assertEqual(blast_project.read_reciprocal_information_table(filepath='testfiles/missing_blast_project/'), {})

    def test_remote_blast_project_missing_reciprocal_information_table_returns_empty_dict(self):
        blast_project = RemoteBlastProject.objects.get(r_project_title='test remote project 1')
        self.assertEqual(blast_project.read_reciprocal_information_table(filepath='testfiles/missing_blast_project/'), {})

    def test_query_information_table_escapes_fasta_header_like_html(self):
        blast_project = BlastProject.objects.get(project_title='test project 1')

        with tempfile.TemporaryDirectory() as tempdir:
            project_dir = Path(tempdir) / str(blast_project.id)
            project_dir.mkdir()
            (project_dir / 'query_sequence_information.csv').write_text(
                "Query\tFeatures\tDescription\n"
                "WP_TEST\t['<b>domain</b>', '<script>alert(1)</script>']\t"
                "header <img src=x onerror=alert(1)>\n",
                encoding='utf-8',
            )

            html = blast_project.read_query_information_table(filepath=str(Path(tempdir)) + os.sep)

        self.assertIn('&lt;b&gt;domain&lt;/b&gt;', html)
        self.assertIn('&lt;script&gt;alert(1)&lt;/script&gt;', html)
        self.assertIn('&lt;img src=x onerror=alert(1)&gt;', html)
        self.assertNotIn('<script>alert(1)</script>', html)
        self.assertNotIn('<img src=x onerror=alert(1)>', html)

    @tag('biological')
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

    @tag('biological')
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
