from os.path import isdir, isfile
import pandas as pd
from os import mkdir
from django.db import models
from django.contrib.auth.models import User
from django_celery_results.models import TaskResult
from blast_project.models import BlastSettings, BlastDatabase
from .managers import OneWayBlastProjectManager, OneWayRemoteBlastProjectManager
from django.db import IntegrityError
# Create your models here.

#TODO documentation
class OneWayBlastProject(models.Model):
    project_title = models.CharField(
        max_length=200, blank=False, unique=True,
        verbose_name="title of this one way blast project"
    )

    project_query_sequences = models.CharField(
        max_length=200,
        verbose_name="query sequence filepath"
    )

    project_user = models.ForeignKey(
        User,
        on_delete=models.CASCADE,
        verbose_name="user who created this project"
    )

    project_settings = models.OneToOneField(
        BlastSettings,
        on_delete=models.CASCADE,
        verbose_name="BLAST settings"
    )

    project_database = models.ForeignKey(
        BlastDatabase,
        on_delete=models.CASCADE,
        verbose_name="associated BLAST database"
    )

    project_execution_task_result = models.OneToOneField(
        TaskResult,
        on_delete=models.SET_NULL,
        blank=True, null=True,
        verbose_name="django_celery_results taskresult model for this project"
    )

    timestamp = models.DateTimeField(auto_now=True)

    objects = OneWayBlastProjectManager()

    def __str__(self):
        return "One Way BLAST Project, created {} by {} with database {}".format(
            self.timestamp,self.project_user.username, self.project_database.database_name
        )

    def get_project_username(self):
        return self.project_user.name

    def get_project_useremail(self):
        return self.project_user.email

    def get_project_dir(self):
        return 'media/one_way_blast/' + str(self.id)

    def initialize_project_directory(self):
        # check if blast_project was previously created / check if media/blast_project directory exists
        if (isdir('media/one_way_blast/' + str(self.id)) or isdir('media/one_way_blast/') == False):
            raise IntegrityError("project directory exists")
        else:
            try:
                mkdir('media/one_way_blast/' + str(self.id))
                if (isdir('static/images/result_images/one_way_blast/' + str(self.id)) == False):
                    mkdir('static/images/result_images/one_way_blast/' + str(self.id))
            except Exception as e:
                raise IntegrityError("couldnt create project directory : {}".format(e))

    def write_snakemake_configuration_file(self):
        try:
            snk_config_file = open('media/one_way_blast/' + str(self.id) + '/snakefile_config', 'w')
            # database path from media/blast_projects/project_id as working directory for snakemake
            snk_config_file.write('project_id: ' + str(self.id) + "\n")
            snk_config_file.write('blastdb: ' + "\"" + "../../databases/" + str(
                self.project_database.id) + "/" + self.project_database.get_pandas_table_name() + ".database\"\n")
            snk_config_file.write('query_sequence: ' + "\"" + self.project_query_sequences + "\"\n")
            snk_config_file.write('user_email: '+str(self.project_user.email)+"\n")


            settings_dict = self.project_settings.get_values_as_dict()

            # print(bw_dict)
            for key in settings_dict.keys():
                snk_config_file.write(key + ': ' + settings_dict[key] + "\n")

            snk_config_file.close()

        except Exception as e:
            raise IntegrityError(
                "couldnt write snakemake configuration file in directory with exception : {}".format(e))

    def read_query_information_table(self, filepath='media/one_way_blast/'):
        try:
            path_to_information_table = filepath+str(self.id)+"/"+"query_sequence_information.csv"
            if isfile(path_to_information_table):
                table = pd.read_table(path_to_information_table,header=0,sep="\t",index_col=0)
                table = table.fillna(value='')
                table = table.to_html(classes='my_class" id="myTable')
                return table
            else:
                return "there is no query_sequence_information.csv in the project directory"
        except Exception as e:
            raise Exception("[-] ERROR during pandas parsing of query_sequence_information csv file with exception: {}".format(e))

#TODO documentation - on_delete=models.CASCADE!?
class OneWayRemoteBlastProject(models.Model):

    BLAST_SEARCH_PROGRAMS = [('blastp', 'blastp'), ('blastn', 'blastn')]
    BLAST_REMOTE_DATABASES = [('nr','nr'),('nt','nt')]

    r_project_title = models.CharField(
        max_length=200, blank=False, unique=True,
        verbose_name="title of this one way blast project"
    )

    r_project_query_sequences = models.CharField(
        max_length=200,
        verbose_name="query sequence filepath"
    )

    r_project_user = models.ForeignKey(
        User,
        on_delete=models.CASCADE,
        verbose_name="user who created this project"
    )

    r_project_settings = models.OneToOneField(
        BlastSettings,
        on_delete=models.CASCADE,
        verbose_name="BLAST settings"
    )

    r_search_strategy = models.CharField(
        max_length=20,
        choices=BLAST_SEARCH_PROGRAMS,
        default='blastp',
        verbose_name="BLAST program that is used inside snakemake execution")

    r_project_database = models.CharField(
        max_length=20,
        choices=BLAST_REMOTE_DATABASES,
        default='nr',
        verbose_name="BLAST program that is used inside snakemake execution")

    r_project_execution_task_result = models.OneToOneField(
        TaskResult,
        on_delete=models.SET_NULL,
        blank=True, null=True,
        verbose_name="django_celery_results taskresult model for this project"
    )

    #path to taxonomic node file ...
    r_entrez_query = models.CharField(
        max_length=300,
        blank=True, null=True,
        verbose_name="associated taxonomic file, which was used to limit assembly entries in db creation by taxids")

    r_timestamp = models.DateTimeField(auto_now=True)

    objects = OneWayRemoteBlastProjectManager()

    def __str__(self):
        return "One Way BLAST Project, created {} by {} with database {}".format(
            self.r_timestamp,self.r_project_user.username, self.r_project_database
        )

    def get_project_username(self):
        return self.r_project_user.name

    def get_project_useremail(self):
        return self.r_project_user.email

    def get_project_dir(self):
        return 'media/one_way_blast/remote_searches/' + str(self.id)

    def initialize_project_directory(self):
        # check if blast_project was previously created / check if media/blast_project directory exists
        if (isdir('media/one_way_blast/remote_searches/' + str(self.id)) or isdir('media/one_way_blast/remote_searches/') == False):
            raise IntegrityError("project directory exists")
        else:
            try:
                mkdir('media/one_way_blast/remote_searches/' + str(self.id))
                if (isdir('static/images/result_images/one_way_blast/remote_searches/' + str(self.id)) == False):
                    mkdir('static/images/result_images/one_way_blast/remote_searches/' + str(self.id))
            except Exception as e:
                raise IntegrityError("couldnt create project directory : {}".format(e))

    def write_snakemake_configuration_file(self):
        try:
            snk_config_file = open('media/one_way_blast/remote_searches/' + str(self.id) + '/snakefile_config', 'w')
            # database path from media/blast_projects/project_id as working directory for snakemake
            snk_config_file.write('project_id: ' + str(self.id) + "\n")
            snk_config_file.write('blastdb: '  + str(self.r_project_database) + "\n")
            snk_config_file.write('query_sequence: ' + "\"" + self.r_project_query_sequences + "\"\n")
            snk_config_file.write('search_strategy: ' + str(self.r_search_strategy) +"\n")
            snk_config_file.write('entrez_query:'+"\n")
            snk_config_file.write('user_email: '+str(self.r_project_user.email)+"\n")

            settings_dict = self.r_project_settings.get_values_as_dict()

            # print(bw_dict)
            for key in settings_dict.keys():
                snk_config_file.write(key + ': ' + settings_dict[key] + "\n")

            snk_config_file.close()

        except Exception as e:
            raise IntegrityError(
                "couldnt write snakemake configuration file in directory with exception : {}".format(e))

    def read_query_information_table(self, filepath='media/one_way_blast/remote_searches/'):
        try:
            path_to_information_table = filepath+str(self.id)+"/"+"query_sequence_information.csv"
            print("[*] {}".format(path_to_information_table))
            if isfile(path_to_information_table):
                table = pd.read_table(path_to_information_table,header=0,sep="\t",index_col=0)
                table = table.fillna(value='')
                table = table.to_html(classes='my_class" id="myTable')
                return table
            else:
                return "there is no query_sequence_information.csv in the project directory"
        except Exception as e:
            raise Exception("[-] ERROR during pandas parsing of query_sequence_information csv file with exception: {}".format(e))