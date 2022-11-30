from django.test import TestCase
from django.contrib.auth.models import User
from django.utils import timezone
from blast_project.forms import CreateTaxonomicFileForMultipleScientificNames, \
    ProjectCreationForm, UploadGenomeForm, UploadMultipleFilesGenomeForm
from refseq_transactions.models import BlastDatabase, AssemblyLevels
from blast_project.py_biopython import get_species_taxid_by_name
from blast_project.py_django_db_services import check_if_taxid_is_in_database, \
    check_if_sequences_are_in_database, get_all_succeeded_databases
from django.core.files.uploadedfile import SimpleUploadedFile
from django_celery_results.models import TaskResult

class UploadMultipleFilesGenomeFormTestCase(TestCase):
    def setUp(self) -> None:
        User.objects.create_user(
            'testuser',
            'test_email@email.com',
            password='test',
            last_login=timezone.now()
        )

    def test_uploadmultiplefilesgenomeform_is_valid(self):
        upload_file = open('testfiles/testsequences/lps_transport.faa', 'rb')
        file_dict = {'genome_file_field_0': SimpleUploadedFile(upload_file.name, upload_file.read())}

        data_dict={
            "database_title":"Test Database Upload Multiple Files",
            "database_description":"Curvibacter Test Upload Multiple Files Database",
            "organism_name_0":"Curvibacter sp. AEP1-3",
            "extra_field_count":0,
            "user_email":"lukas.becker@hhu.de"
        }
        user = User.objects.get(username="testuser")
        form = UploadMultipleFilesGenomeForm(user=user,data=data_dict,files=file_dict)
        self.assertTrue(form.is_valid())

    def test_uploadmultiplefilesgenomeform_is_not_valid(self):
        upload_file = open('testfiles/testsequences/lps_transport.faa', 'rb')
        file_dict = {'genome_file_field_0': SimpleUploadedFile(upload_file.name, upload_file.read())}

        data_dict={
            "database_title":"Test Database Upload Multiple Files",
            "database_description":"Curvibacter Test Upload Multiple Files Database",
            "organism_name_0":"XYZ",
            "extra_field_count":0,
            "user_email":"lukas.becker@hhu.de"
        }
        user = User.objects.get(username="testuser")
        form = UploadMultipleFilesGenomeForm(user=user,data=data_dict,files=file_dict)
        self.assertFalse(form.is_valid())

    def test_uploadmultiplefilesgenomeform_is_not_valid_wrong_filename(self):
        upload_file = open('testfiles/taxonomic_nodes/curvibacter_hydra_test.taxids', 'rb')
        file_dict = {'genome_file_field_0': SimpleUploadedFile(upload_file.name, upload_file.read())}

        data_dict={
            "database_title":"Test Database Upload Multiple Files",
            "database_description":"Curvibacter Test Upload Multiple Files Database",
            "organism_name_0":"XYZ",
            "extra_field_count":0,
            "user_email":"lukas.becker@hhu.de"
        }
        user = User.objects.get(username="testuser")
        form = UploadMultipleFilesGenomeForm(user=user,data=data_dict,files=file_dict)
        self.assertFalse(form.is_valid())

class CreateTaxonomicFileForMultipleScientificNamesTestCase(TestCase):
    def setUp(self):
        User.objects.create_user(
            'testuser',
            'test_email@email.com',
            password='test',
            last_login=timezone.now()
        )

    def test_form_is_valid(self):
        form = CreateTaxonomicFileForMultipleScientificNames(
            user=User.objects.get(username="testuser"),
            data={"filename":"test file",
                  "species_names":"Curvibacter sp. AEP1-3","user_email":"test_email@email.com"})
        self.assertTrue(form.is_valid())

class ProjectCreationFormTestCase(TestCase):
    def setUp(self) -> None:
        user = User.objects.create_user(
            'testuser',
            'test_email@email.com',
            password='test',
            last_login=timezone.now()
        )

        celery_task = TaskResult.objects.create(
            status="SUCCESS"
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

        blast_database = BlastDatabase.objects.create(
            database_name = 'Curvibacter sp. aep1-3 database',
            database_description = 'Full Assembly Content',
            assembly_entries = 1,
            path_to_database_file="testfiles/databases/curvibacter_test_db",
            database_download_and_format_task=celery_task
        )

        assembly_level = ['Complete Genome', 'Chromosome', 'Contig', 'Scaffold']
        assembly_levels = AssemblyLevels.objects.filter(assembly_level__in=assembly_level)
        for assembly_lvl in assembly_levels:
            blast_database.assembly_levels.add(assembly_lvl)

        blast_database.save()

    def test_get_species_taxid_by_name(self):
        taxonomic_nodes = get_species_taxid_by_name("lukas.becker@hhu.de", "Curvibacter sp. AEP1-3")
        self.assertEqual(taxonomic_nodes,['1844971', '1531298'])

    def test_true_check_if_taxid_is_in_database(self):
        taxonomic_nodes = get_species_taxid_by_name("lukas.becker@hhu.de", "Curvibacter sp. AEP1-3")
        blast_database = BlastDatabase.objects.get(
            database_name="Curvibacter sp. aep1-3 database")
        val_check = check_if_taxid_is_in_database(blast_database.id,taxonomic_nodes)
        self.assertEqual(True,val_check)

    def test_false_check_if_taxid_is_in_database(self):
        taxonomic_nodes = get_species_taxid_by_name("lukas.becker@hhu.de", "Curvibacter delicatus")
        blast_database = BlastDatabase.objects.get(
            database_name="Curvibacter sp. aep1-3 database")
        val_check = check_if_taxid_is_in_database(blast_database.id,taxonomic_nodes)
        #['80879'] not in db
        self.assertEqual(False,val_check)

    def test_true_check_if_sequences_are_in_database(self):
        #lps transporter
        protein_identifier = ['WP_087495344',
                          'WP_087495343',
                          'WP_087497333',
                          'WP_087495837',
                          'WP_087496015',
                          'WP_087493791',
                          'WP_087493790',
                          'WP_087493929']
        blast_database = BlastDatabase.objects.get(
            database_name="Curvibacter sp. aep1-3 database")
        val_check = check_if_sequences_are_in_database(blast_database.id,protein_identifier)
        self.assertEqual(True,val_check)

    def test_false_check_if_sequences_are_in_database(self):
        #lps transporter
        protein_identifier = ['WP_087495344',
                          'WP_087495343',
                          'WP_087497333',
                          'WP_087495837',
                          'WP_087496015',
                          'WP_087493791',
                          'WP_087493790',
                          'QPASDKJ']
        blast_database = BlastDatabase.objects.get(
            database_name="Curvibacter sp. aep1-3 database")
        val_check = check_if_sequences_are_in_database(blast_database.id,protein_identifier)
        self.assertEqual(['QPASDKJ'],val_check)


    def test_project_form_is_valid(self):
        #primary key of the referencing modelchoicefield object
        blast_database = get_all_succeeded_databases()[0]
        user = User.objects.get(username="testuser")
        upload_file = open('testfiles/testsequences/lps_transport.faa', 'rb')
        post_dict = {'project_title': 'Test Project Title',
                     'species_name_for_backward_blast':'Curvibacter sp. AEP1-3',
                     'user_email':'lukas.becker@hhu.de',
                     'project_forward_database':blast_database.pk,
                     'project_backward_database':blast_database.pk,
                     }
        file_dict = {'query_sequence_file': SimpleUploadedFile(upload_file.name, upload_file.read())}
        form = ProjectCreationForm(user=user,data=post_dict,files=file_dict)
        self.assertTrue(form.is_valid())
        self.assertEqual(2,len(form.cleaned_data['species_name_for_backward_blast']))

    def test_project_form_is_not_valid_wrong_species_name(self):
        #primary key of the referencing modelchoicefield object
        blast_database = get_all_succeeded_databases()[0]
        user = User.objects.get(username="testuser")
        upload_file = open('testfiles/testsequences/lps_transport.faa', 'rb')
        #species curvibacter delicatus does not reside in backward database
        post_dict = {'project_title': 'Test Project Title',
                     'species_name_for_backward_blast':'Curvibacter delicatus',
                     'user_email':'lukas.becker@hhu.de',
                     'project_forward_database':blast_database.pk,
                     'project_backward_database':blast_database.pk,
                     }
        file_dict = {'query_sequence_file': SimpleUploadedFile(upload_file.name, upload_file.read())}
        form = ProjectCreationForm(user=user,data=post_dict,files=file_dict)
        self.assertFalse(form.is_valid())
        self.assertEqual(form.errors['species_name_for_backward_blast'],
        ["specified taxonomic node: ['80879'] does not reside in the selected BACKWARD database: Curvibacter sp. aep1-3 database"])


    def test_project_form_is_not_valid_sequence_identifier_not_in_backward_db(self):
        #primary key of the referencing modelchoicefield object
        blast_database = get_all_succeeded_databases()[0]
        user = User.objects.get(username="testuser")
        upload_file = open('testfiles/testsequences/lps_transporter_corrupted.faa', 'rb')
        #species curvibacter delicatus does not reside in backward database
        post_dict = {'project_title': 'Test Project Title',
                     'species_name_for_backward_blast':'Curvibacter sp. AEP1-3',
                     'user_email':'lukas.becker@hhu.de',
                     'project_forward_database':blast_database.pk,
                     'project_backward_database':blast_database.pk,
                     }
        file_dict = {'query_sequence_file': SimpleUploadedFile(upload_file.name, upload_file.read())}
        form = ProjectCreationForm(user=user,data=post_dict,files=file_dict)
        self.assertFalse(form.is_valid())
        self.assertEqual(form.errors['query_sequence_file'],
                         ["following sequences do not reside in your backward database: ['PXYR_6321']"])

