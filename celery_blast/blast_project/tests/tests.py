from django.test import TestCase
#this is the correct import
from blast_project.models import BlastProject, BlastDatabase, BlastSettings, AssemblyLevels

class BlastProjectTestCase(TestCase):
    def setUp(self):

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

        BlastSettings.objects.create(
            e_value=0.001,
            word_size=3,
            num_alignments=10000,
            max_target_seqs=10000,
            num_threads=1,
            max_hsps=500
        )

        BlastSettings.objects.create(
            e_value=0.001,
            word_size=3,
            num_alignments=1,
            max_target_seqs=1,
            num_threads=1,
            max_hsps=500
        )

        BlastDatabase.objects.create(
            database_name = 'Test Database 1',
            database_description = 'Full Assembly Content',
            assembly_entries = 2
        )


    def test_blast_database_model(self):
        blast_database = BlastDatabase.objects.get(database_name='Test Database 1')
        self.assertEqual('Test Database 1',blast_database.database_name)

    def test_blast_database_assembly_level_creation(self):
        blast_database = BlastDatabase.objects.get(database_name='Test Database 1')
        assembly_level = ['Complete Genome', 'Chromosome', 'Contig', 'Scaffold']
        assembly_levels = AssemblyLevels.objects.filter(assembly_level__in=assembly_level)
        print(assembly_levels)
        for assembly_lvl in assembly_levels:
            blast_database.assembly_levels.add(assembly_lvl)
        level1 = blast_database.assembly_levels.get(assembly_level='Complete Genome')
        self.assertEqual(assembly_levels.count(),4)
        self.assertEqual(level1.assembly_level,'Complete Genome')