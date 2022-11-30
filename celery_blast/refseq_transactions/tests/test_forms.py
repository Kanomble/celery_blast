from django.test import TestCase
from django.utils import timezone
from django.core.files.uploadedfile import SimpleUploadedFile
from django_celery_results.models import TaskResult
from refseq_transactions.models import BlastDatabase, AssemblyLevels
from refseq_transactions.forms import RefseqDatabaseForm

class RefseqDatabaseFormTest(TestCase):
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

    def test_refseqdatabaseform_is_valid(self):
        upload_file = open('testfiles/taxonomic_nodes/curvibacter_hydra_test.taxids', 'rb')
        post_dict = {'assembly_levels':["Complete Genome","Chromosome"],
                     "database_name":"test database curvibacter 1",
                     "database_description":"database for testing purposes",
                     }
        file_dict = {'query_sequence_file': SimpleUploadedFile(upload_file.name, upload_file.read())}

        form = RefseqDatabaseForm(data=post_dict,files=file_dict)
        upload_file.close()
        self.assertTrue(form.is_valid())