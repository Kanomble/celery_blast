from django.db import models
from .managers import BlastDatabaseManager
from django_celery_results.models import TaskResult

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

    uploaded_files = models.BooleanField(
        default=False,blank=True,null=True
    )

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

    objects = BlastDatabaseManager()
    # functions
    def __str__(self):
        return "BLAST database: {}, created {} with {} entries.\n\t Database description: {}".format(
            self.database_name,
            self.timestamp,
            self.assembly_entries,
            self.database_description)

    def get_pandas_table_name(self):
        return self.database_name.replace(' ','_').upper()

    #can be deleted?
    #TODO deprecated
    def get_database_palfile_for_snakemake_config(self):
        return self.path_to_database_file + '/' + self.database_name.replace(' ','_').upper() + '.database.pal'
