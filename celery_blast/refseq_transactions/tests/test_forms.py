from django.core.files.uploadedfile import SimpleUploadedFile
from django.test import TestCase, tag
from refseq_transactions.forms import RefseqDatabaseForm
from refseq_transactions.models import AssemblyLevels


@tag("refseqdatabaseformtests")
class RefseqDatabaseFormTest(TestCase):

    @classmethod
    def setUpClass(cls):
        print("[+] Setting up model data for refseqdatabaseform tests")
        super().setUpClass()
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

    @classmethod
    def tearDownClass(cls):
        super().tearDownClass()

    def test_refseqdatabaseform_is_valid(self):
        upload_file = open('testfiles/taxonomic_nodes/curvibacter_hydra_test.taxids', 'rb')
        post_dict = {'assembly_levels': ["Complete Genome", "Chromosome"],
                     "database_name": "test database curvibacter 1",
                     "database_description": "database for testing purposes",
                     }
        file_dict = {'taxid_file': SimpleUploadedFile(upload_file.name, upload_file.read())}

        form = RefseqDatabaseForm(data=post_dict, files=file_dict)
        upload_file.close()
        self.assertTrue(form.is_valid())

    def test_refseqdatabaseform_with_uploadedtaxidfile_is_valid(self):
        upload_file = open('testfiles/taxonomic_nodes/curvibacter_hydra_test.taxids', 'rb')
        post_dict = {'assembly_levels': ["Complete Genome", "Chromosome"],
                     "database_name": "test database curvibacter 1",
                     "database_description": "database for testing purposes",
                     "taxid_uploaded_file": "curvibacter_hydra.taxids"
                     }
        form = RefseqDatabaseForm(data=post_dict)
        print(form.errors)
        self.assertTrue(form.is_valid())

    def test_refseqdatabaseform_is_not_valid_wrong_assembly_levels(self):
        post_dict = {'assembly_levels': ["U", "Chromosome"],
                     "database_name": "test database curvibacter 1",
                     "database_description": "database for testing purposes",
                     }

        form = RefseqDatabaseForm(data=post_dict)
        self.assertFalse(form.is_valid())

    def test_refseqdatabaseform_is_not_valid_no_database_name(self):
        post_dict = {'assembly_levels': ["Complete Genome", "Chromosome"],
                     "database_name": "",
                     "database_description": "database for testing purposes",
                     }

        form = RefseqDatabaseForm(data=post_dict)
        self.assertFalse(form.is_valid())
