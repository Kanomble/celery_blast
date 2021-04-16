from django.db import models

# Create your models here.
from os.path import isdir
from os import mkdir
from django.db import models
from django.contrib.auth.models import User
from django_celery_results.models import TaskResult
from django.core.exceptions import ValidationError
from .managers import BlastProjectManager

class BlastSettings(models.Model):
    e_value = models.DecimalField(
        max_digits=30,
        decimal_places=15,
        default=0.0001)
    word_size = models.IntegerField(
        default=3)
    num_threads = models.IntegerField(
        default=1)

    num_alignments = models.IntegerField()
    max_target_seqs = models.IntegerField()
    max_hsps = models.IntegerField()


class AssemblyLevels(models.Model):
    assembly_level = models.CharField(max_length=50)

class BlastDatabase(models.Model):
    # attribute fields
    database_name = models.CharField(
        max_length=200,
        blank=False, unique=True,
        verbose_name="database name")
    database_description = models.CharField(
        max_length=200,
        verbose_name="short description of database purpose")

    assembly_entries = models.IntegerField(
        verbose_name="number of assembly entries that should get downloaded")
    timestamp = models.DateTimeField(
        auto_now=True,
        verbose_name="date of database creation")

    # nullable fields
    # possibility to add a taxonomic file
    attached_taxonomic_node_file = models.CharField(
        max_length=300,
        blank=True, null=True,
        verbose_name="associated taxonomic file, which was used to limit assembly entries in db creation by taxids")
    path_to_database_file = models.CharField(
        max_length=300,
        blank=True, null=True,
        verbose_name="after makeblastdb task has finished this field is set automatically with the path to the BLAST database")

    # relationships
    database_download_and_format_task = models.OneToOneField(
        TaskResult,
        on_delete=models.CASCADE,
        blank=True, null=True,
        verbose_name="django_celery_results taskresult model for download and formatting procedure")

    # use the assembly_levels.SQL script for uploading the four existing assembly levels into the database
    assembly_levels = models.ManyToManyField(
        to=AssemblyLevels,
        verbose_name="possible assembly levels within this BLAST database")

    # functions
    def __str__(self):
        return "BLAST database: {}, created {} with {} entries.\n\t Database description: {}".format(
            self.database_name,
            self.timestamp,
            self.assembly_entries,
            self.database_description)


class BlastProject(models.Model):
    BLAST_SEARCH_PROGRAMS = [('blastp', 'blastp'), ('blastn', 'blastn')]

    # attribute fields
    project_title = models.CharField(
        max_length=200, blank=False, unique=True,
        verbose_name="title of this blast project")
    search_strategy = models.CharField(
        max_length=20,
        choices=BLAST_SEARCH_PROGRAMS,
        default='blastp',
        verbose_name="BLAST program that is used inside snakemake execution")
    project_query_sequences = models.CharField(
        max_length=200,
        verbose_name="query sequence filepath for the forward BLAST")
    timestamp = models.DateTimeField(auto_now=True)

    # relationships
    project_user = models.ForeignKey(
        User,
        on_delete=models.CASCADE,
        verbose_name="user who created the project")

    project_forward_settings = models.OneToOneField(
        BlastSettings,
        on_delete=models.CASCADE,
        verbose_name="settings for the forward BLAST execution",
        related_name='project_forward_settings')

    project_backward_settings = models.OneToOneField(
        BlastSettings,
        on_delete=models.CASCADE,
        verbose_name="settings for the backward BLAST execution",
        related_name='project_backward_settings')

    # forward database and backward database
    # each project can have ONE forward and ONE backward BlastDatabase
    # --> two foreign keys
    project_forward_database = models.ForeignKey(
        BlastDatabase,
        on_delete=models.CASCADE,
        verbose_name="associated forward BLAST database",
        related_name="project_forward_database")

    project_backward_database = models.ForeignKey(
        BlastDatabase,
        on_delete=models.CASCADE,
        verbose_name="associated backward BLAST database",
        related_name="project_backward_database")
    # one to one relationship
    project_execution_snakemake_task = models.OneToOneField(
        TaskResult,
        on_delete=models.SET_NULL,
        blank=True, null=True,
        verbose_name="django_celery_results taskresult model for this project")

    # customized initialization can be added in BlastProjectManager (e.g. direct creation of project directory
    objects = BlastProjectManager()

    # overwritten functions
    def __str__(self):
        return "Reciprocal BLAST Project, created {} by {} with fw db {} and bw db {}".format(
            self.timestamp, self.project_user.name,
            self.project_forward_database.database_name,
            self.project_backward_database.database_name)

    def get_project_username(self):
        return self.project_user.name

    def get_project_useremail(self):
        return self.project_user.email

    def get_project_dir(self):
        return 'media/blast_projects/' + str(self.id)

    def if_executed_return_associated_taskresult_model(self):
        # executed
        if self.project_execution_snakemake_task != None:
            return self.project_execution_snakemake_task
        # not executed
        else:
            return None

    # invoke like: project.objects.initialize_project_directory()
    def initialize_project_directory(self):
        # check if blast_project was previously created / check if media/blast_project directory exists
        if (isdir('media/blast_projects/' + str(self.id)) or isdir('media/blast_projects/') == False):
            raise ValidationError("project directory exists")
        else:
            try:
                mkdir('media/blast_projects/' + str(self.id))
            except Exception as e:
                raise ValidationError("couldnt create project directory : {}".format(e))


