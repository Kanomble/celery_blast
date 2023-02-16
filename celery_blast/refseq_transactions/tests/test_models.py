from django.test import TestCase
from refseq_transactions.models import AssemblyLevels


class AssemblyLevelsTest(TestCase):
    def setUp(self) -> None:
        if len(AssemblyLevels.objects.all()) < 4:
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

    def test_assembly_levels(self):
        assembly_levels = AssemblyLevels.objects.all()
        complete = AssemblyLevels.objects.get(assembly_level="Complete Genome")
        self.assertTrue(complete in assembly_levels)
