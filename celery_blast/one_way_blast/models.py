from os.path import isdir
from os import mkdir
from django.db import models
from django.contrib.auth.models import User
from django_celery_results.models import TaskResult
from blast_project.models import BlastSettings, BlastDatabase
from .managers import OneWayBlastProjectManager
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
        return 'media/blast_projects/' + str(self.id)

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


            settings_dict = self.project_settings.get_values_as_dict()

            # print(bw_dict)
            for key in settings_dict.keys():
                snk_config_file.write(key + ': ' + settings_dict[key] + "\n")

            snk_config_file.close()

        except Exception as e:
            raise IntegrityError(
                "couldnt write snakemake configuration file in directory with exception : {}".format(e))
