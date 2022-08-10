# Create your models here.
from os.path import isdir, isfile
from os import mkdir
from django.db import models
from django.contrib.auth.models import User
from django_celery_results.models import TaskResult
from django.db import IntegrityError
from .managers import BlastProjectManager, BlastDatabaseManager
import pandas as pd

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

    #TODO documentation
    def values_as_fw_or_bw_dict(self,fwOrBw):
        settings_dict = {
            fwOrBw+'_e_value' : str(self.e_value),
            fwOrBw+'_word_size' : str(self.word_size),
            fwOrBw+'_num_threads' : str(self.num_threads),
            fwOrBw+'_num_alignments' : str(self.num_alignments),
            fwOrBw+'_max_target_seqs' : str(self.max_target_seqs),
            fwOrBw+'_max_hsps' : str(self.max_hsps)
        }
        return settings_dict

    def get_values_as_dict(self):
        settings_dict = {
            'e_value' : str(self.e_value),
            'word_size' : str(self.word_size),
            'num_threads' : str(self.num_threads),
            'num_alignments' : str(self.num_alignments),
            'max_hsps' : str(self.max_hsps)
        }
        return settings_dict

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
        verbose_name="user who created this project")

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

    #each project can have one forward BlastDatabase
    project_forward_database = models.ForeignKey(
        BlastDatabase,
        on_delete=models.CASCADE,
        verbose_name="associated forward BLAST database",
        related_name="project_database"
    )

    project_backward_database = models.ForeignKey(
        BlastDatabase,
        on_delete=models.CASCADE,
        verbose_name="associated backward BLAST database",
        related_name="project_backward_database"
    )

    species_name_for_backward_blast = models.CharField(
        max_length=200, blank=False,
        verbose_name="species name for the backward database and query sequences"
    )

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
        return "Reciprocal BLAST Project, created {} by {} with fw db {}".format(
            self.timestamp, self.project_user.username,
            self.project_forward_database.database_name)

    def get_list_of_query_sequences(self, filepath='media/blast_projects/'):
        try:
            query_sequence_file_path = self.get_project_query_sequence_filepath(filepath)
            query_file = open(query_sequence_file_path, 'r')
            qseqids = []
            for line in query_file.readlines():
                if ">" in line:
                    qseqid = line.split(" ")[0].split(">")[1].split(".")[0]
                    qseqids.append(qseqid)
            return qseqids
        except Exception as e:
            raise IntegrityError(
                "[-] couldnt extract query sequence ids from query sequence file : {} with exception : {}".format(
                    query_sequence_file_path, e))

    def get_project_username(self):
        return self.project_user.name

    def get_project_useremail(self):
        return self.project_user.email

    def get_project_dir(self):
        return 'media/blast_projects/' + str(self.id)

    def get_project_query_sequence_filepath(self, filepath='media/blast_projects/'):
        return filepath + str(self.id) + '/' + self.project_query_sequences

    def if_executed_return_associated_taskresult_model(self):
        # executed
        if self.project_execution_snakemake_task != None:
            return self.project_execution_snakemake_task
        # not executed
        else:
            return None

    ''' initialize_project_directory
    
        invokation is done with BlastProject.objects.initialize_project_directory() or 
        within the BlastProjectManager with blast_project.iniialize_project_directory()
        
        this function is executed within the create_blast_project function of the manager, 
        it directly creates a directory for the project and a directory for result images in the static folder
    '''
    def initialize_project_directory(self,filepath='media/blast_projects/'):
        # check if blast_project was previously created / check if media/blast_project directory exists
        if (isdir(filepath + str(self.id)) or isdir(filepath) == False):
            raise IntegrityError("project directory exists")
        else:
            try:
                mkdir(filepath + str(self.id))
                if(isdir('static/images/result_images/'+str(self.id)) == False):
                    mkdir('static/images/result_images/'+str(self.id))
            except Exception as e:
                raise IntegrityError("couldnt create project directory : {}".format(e))


    #TODO documentation
    def write_snakemake_configuration_file(self, filepath='media/blast_projects/'):
        try:
            snk_config_file = open(filepath + str(self.id)+'/snakefile_config','w')
            #database path from media/blast_projects/project_id as working directory for snakemake
            snk_config_file.write('project_id: '+str(self.id)+"\n")
            snk_config_file.write('blastdb: ' +"\"" +"../../databases/" + str(self.project_forward_database.id) + "/" + self.project_forward_database.get_pandas_table_name() + ".database\"\n")
            snk_config_file.write('backwarddb: '+"\""+"../../databases/"+str(self.project_backward_database.id) + "/" + self.project_backward_database.get_pandas_table_name() + ".database\"\n")
            snk_config_file.write('query_sequence: '+"\""+self.project_query_sequences+"\"\n")
            snk_config_file.write('bw_taxid: '+str(self.species_name_for_backward_blast[1])+"\n")
            snk_config_file.write('user_email: '+str(self.project_user.email)+"\n")
            bw_dict=self.project_backward_settings.values_as_fw_or_bw_dict('bw')
            fw_dict=self.project_forward_settings.values_as_fw_or_bw_dict('fw')

            for key_bw in bw_dict.keys():
                snk_config_file.write(key_bw+': '+bw_dict[key_bw]+"\n")

            for key_fw in fw_dict.keys():
                snk_config_file.write(key_fw+': '+fw_dict[key_fw]+"\n")

            snk_config_file.close()

        except Exception as e:
            raise IntegrityError("couldnt write snakemake configuration file in directory with exception : {}".format(e))

    def read_query_information_table(self, filepath='media/blast_projects/'):
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
