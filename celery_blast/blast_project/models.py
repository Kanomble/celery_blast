# Create your models here.
from os.path import isdir, isfile
from ast import literal_eval
from django.db import models
from django.contrib.auth.models import User
from django_celery_results.models import TaskResult
from django.db import IntegrityError
from .managers import BlastProjectManager, RemoteBlastProjectManager
from refseq_transactions.models import BlastDatabase
from os.path import isdir, isfile
from os import mkdir, listdir

import pandas as pd
from celery_blast.settings import BLAST_PROJECT_DIR, BLAST_DATABASE_DIR, REMOTE_BLAST_PROJECT_DIR


# TODO USE THOSE paths for setting up the project directory and snakemake config file
# from celery_blast.settings import BLAST_DATABASE_DIR, BLAST_PROJECT_DIR

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

    '''values_as_fw_or_bw_dict
    
        This function is used in the BlastProject model function write_snakemake_configuration_file.
        It returns a dictionary composed of the provided project settings.
        
        :params fwOrBw - "bw", "fw"
            :type str
        :returns settings_dict
            :type dict[str]=type(setting)
            
    '''

    def values_as_fw_or_bw_dict(self, fwOrBw: str) -> dict:
        settings_dict = {
            fwOrBw + '_e_value': str(self.e_value),
            fwOrBw + '_word_size': str(self.word_size),
            fwOrBw + '_num_threads': str(self.num_threads),
            fwOrBw + '_num_alignments': str(self.num_alignments),
            fwOrBw + '_max_target_seqs': str(self.max_target_seqs),
            fwOrBw + '_max_hsps': str(self.max_hsps)
        }
        return settings_dict

    '''get_values_as_dict
        
        This function returns a dictionary of the provided settings without max_target_seqs.
        
        :returns settings_dict
            :type dict[str]=type(settings)
            
    '''

    def get_values_as_dict(self):
        settings_dict = {
            'e_value': str(self.e_value),
            'word_size': str(self.word_size),
            'num_threads': str(self.num_threads),
            'num_alignments': str(self.num_alignments),
            'max_hsps': str(self.max_hsps)
        }
        return settings_dict


class BlastProject(models.Model):
    # currently just blastp is enabled
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

    # each project can have one forward BlastDatabase
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
        related_name="project_execution_snakemake_task",
        verbose_name="django_celery_results taskresult model for this projects snakemake pipeline")

    project_database_statistics_task = models.OneToOneField(
        TaskResult,
        on_delete=models.SET_NULL,
        blank=True, null=True,
        related_name="project_database_statistics_task",
        verbose_name="django_celery_results taskresult model for this projects database statistics")

    project_database_statistics_task_selection = models.OneToOneField(
        TaskResult,
        on_delete=models.SET_NULL,
        blank=True, null=True,
        related_name="project_database_statistics_task_selection",
        verbose_name="django_celery_results taskresult model for this "
                     "projects database statistics phylogenetic inference based on user selection task")

    # customized initialization can be added in BlastProjectManager (e.g. direct creation of project directory
    objects = BlastProjectManager()

    # overwritten functions
    def __str__(self):
        return "Reciprocal BLAST Project, created {} by {} with fw db {}".format(
            self.timestamp, self.project_user.username,
            self.project_forward_database.database_name)

    '''get_list_of_query_sequences
        
        Returns a list of query sequences without the additional .version of the query.
        
        :param self
            :type int
        :param filepath
            :type str
        
        :returns qseqids
            :type list[str]
    '''

    def get_list_of_query_sequences(self, filepath=BLAST_PROJECT_DIR):
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

    '''get_fasta_header_of_query_sequences
        
        This function extract information of the protein fasta header.
        
        :param self
            :type int
        :param filepath
            :type str
    
    '''

    def get_fasta_header_of_query_sequences(self, filepath=BLAST_PROJECT_DIR):
        try:
            query_sequence_file_path = self.get_project_query_sequence_filepath(filepath)
            query_file = open(query_sequence_file_path, 'r')
            qseqids_info = []
            for line in query_file.readlines():
                if ">" in line:
                    qseqid_info = ' '.join(line.split(">")[1].split(" ")[1:])
                    qseqid_info = qseqid_info.rstrip()
                    if len(qseqid_info) <= 1:
                        qseqid_info = "no information provided"
                    qseqids_info.append(qseqid_info)
            return qseqids_info
        except Exception as e:
            raise IntegrityError(
                "[-] couldnt extract query sequence ids from query sequence file : {} with exception : {}".format(
                    query_sequence_file_path, e))

    def get_project_username(self):
        return self.project_user.name

    def get_project_useremail(self):
        return self.project_user.email

    def get_project_dir(self):
        return BLAST_PROJECT_DIR + str(self.id)

    def get_project_query_sequence_filepath(self, filepath=BLAST_PROJECT_DIR):
        return filepath + str(self.id) + '/' + self.project_query_sequences

    def if_executed_return_associated_taskresult_model(self):
        # executed
        if self.project_execution_snakemake_task != None:
            return self.project_execution_snakemake_task
        # not executed
        else:
            return None

    ''' initialize_project_directory
    
        Invocation is done with BlastProject.objects.initialize_project_directory() or 
        within the BlastProjectManager with blast_project.iniialize_project_directory()
        
        Throws an exception if the specified filepath does not exists or if the 
        filepath + project_id directory exists.
        This function is executed within the create_blast_project function of the manager class, 
        it directly creates a directory for the project and a directory for result images in the static folder
            
    '''

    def initialize_project_directory(self, filepath=BLAST_PROJECT_DIR):
        # check if blast_project was previously created / check if media/blast_project directory exists
        if (isdir(filepath + str(self.id)) == True):
            raise IntegrityError("project directory exists")
        elif isdir(filepath) == False:
            raise IntegrityError("{} is not a directory".format(filepath))
        else:
            try:
                mkdir(filepath + str(self.id))
            except Exception as e:
                raise IntegrityError("couldnt create project directory : {}".format(e))

    '''write_snakemake_configuration_file
        
        This function is used during project creation. It writes all necessary settings for snakemake into
        the file snakefile_config in BLAST_PROJECT_DIR. Additional parameters that want to be used during
        the snakemake run have to be added here. The default filepath is media/blast_projects/.
        
        :params filepath=BLAST_PROJECT_DIR
            :type str
    '''

    def write_snakemake_configuration_file(self,project_settings, filepath=BLAST_PROJECT_DIR):
        try:
            with open(filepath + str(self.id) + '/snakefile_config', 'w') as snk_config_file:
                # database path from media/blast_projects/project_id as working directory for snakemake
                snk_config_file.write('project_id: ' + str(self.id) + "\n")
                snk_config_file.write('blastdb: ' + "\"" + "../../databases/" + str(
                    self.project_forward_database.id) + "/" + self.project_forward_database.get_pandas_table_name() + ".database\"\n")
                snk_config_file.write('backwarddb: ' + "\"" + "../../databases/" + str(
                    self.project_backward_database.id) + "/" + self.project_backward_database.get_pandas_table_name() + ".database\"\n")
                snk_config_file.write('query_sequence: ' + "\"" + self.project_query_sequences + "\"\n")
                snk_config_file.write('bw_taxid: ' + str(self.species_name_for_backward_blast[1]) + "\n")
                snk_config_file.write('user_email: ' + str(self.project_user.email) + "\n")
                bw_dict = self.project_backward_settings.values_as_fw_or_bw_dict('bw')
                fw_dict = self.project_forward_settings.values_as_fw_or_bw_dict('fw')

                for key_bw in bw_dict.keys():
                    snk_config_file.write(key_bw + ': ' + bw_dict[key_bw] + "\n")

                for key_fw in fw_dict.keys():
                    snk_config_file.write(key_fw + ': ' + fw_dict[key_fw] + "\n")

                snk_config_file.write("bitscore_filter: "+str(project_settings.cleaned_data['bitscore_filter']) + "\n")
                snk_config_file.write("max_rbhs_for_phylo: "+str(project_settings.cleaned_data['max_amount_of_rbh_for_msa_and_phylogeny']) + "\n")
                snk_config_file.write("trimal_gt: "+str(project_settings.cleaned_data['trimal_gt']) + "\n")
                snk_config_file.write("trimal_st: "+str(project_settings.cleaned_data['trimal_st']) + "\n")
                snk_config_file.write("trimal_cons: "+str(project_settings.cleaned_data['trimal_cons']) + "\n")
                snk_config_file.write("mview_sort: "+str(project_settings.cleaned_data['mview_sort']) + "\n")
                snk_config_file.write("mview_coloring: "+str(project_settings.cleaned_data['mview_coloring']) + "\n")

        except Exception as e:
            raise IntegrityError(
                "couldnt write snakemake configuration file in directory with exception : {}".format(e))

    '''read_query_information_table
        
        This function is getting executed within the project_details_dashboard.html website.
        
    '''

    def read_query_information_table(self, filepath=BLAST_PROJECT_DIR):
        def clean_feature_column(features):
            new_feature_column = []
            try:
                result_string = ''
                for idx, feature in enumerate(literal_eval(features)):
                    if idx % 2 == 0:
                        result_string += feature + " "
                    # TODO add functional linebreak
                    elif idx % 2 == 1:
                        result_string += feature + "     "
                new_feature_column.append(result_string)
                return pd.Series(new_feature_column)
            except Exception as e:
                raise Exception("ERROR:exception: {}".format(e))

        try:
            path_to_information_table = filepath + str(self.id) + '/query_sequence_information.csv'
            if isfile(path_to_information_table):
                table = pd.read_table(path_to_information_table, header=0, sep="\t", index_col=0)
                table.Features = table.Features.apply(clean_feature_column)
                table = table.fillna(value='')
                table = table.to_html(classes='my_class" id="myTable')
                return table
            else:
                return "there is no query_sequence_information.csv in the project directory"
        except Exception as e:
            raise Exception(
                "[-] ERROR during pandas parsing of query_sequence_information csv file with exception: {}".format(e))

    '''check_for_reciprocal_result_table
        
        This function checks if the reciprocal_result.html file is in the project directory or not.
    '''

    def check_for_reciprocal_result_table(self):
        if isfile(self.get_project_dir() + '/' + 'reciprocal_results.html'):
            return True
        else:
            return False

    '''return_list_of_all_logfiles

        Returns a list with all filepaths to logfiles of the snakemake run. 
        If there is no "log" directory, in the BLAST_PROJECT_DIR the function returns an empty list.

        :returns filelist
            :type list
    '''

    def return_list_of_all_logfiles(self) -> list:
        try:
            path_to_project_dir = BLAST_PROJECT_DIR + str(self.id) + "/log"
            if isdir(path_to_project_dir) == False:
                return []
            else:
                files = listdir(path_to_project_dir)
                filelist = []
                for file in files:
                    # check if "file" is a query sequence sub-directory
                    if isdir(path_to_project_dir + '/' + file):
                        for query_log in listdir(path_to_project_dir + '/' + file):
                            filelist.append(file + '/' + query_log)
                    else:
                        filelist.append(file)
                return filelist
        except Exception as e:
            raise Exception(
                "[-] ERROR creating list of all logfiles for project: {} with exception: {}".format(self.id, e))


class RemoteBlastProject(models.Model):
    # currently just blastp is enabled
    BLAST_SEARCH_PROGRAMS = [('blastp', 'blastp'), ('blastn', 'blastn')]

    # attribute fields
    r_project_title = models.CharField(
        max_length=200, blank=False, unique=True,
        verbose_name="title of this blast project")
    r_search_strategy = models.CharField(
        max_length=20,
        choices=BLAST_SEARCH_PROGRAMS,
        default='blastp',
        verbose_name="BLAST program that is used inside snakemake execution")
    r_project_query_sequences = models.CharField(
        max_length=200,
        verbose_name="query sequence filepath for the forward BLAST")
    r_timestamp = models.DateTimeField(auto_now=True)

    # relationships
    r_project_user = models.ForeignKey(
        User,
        on_delete=models.CASCADE,
        verbose_name="user who created this project")

    r_project_forward_settings = models.OneToOneField(
        BlastSettings,
        on_delete=models.CASCADE,
        verbose_name="settings for the forward BLAST execution",
        related_name='r_project_forward_settings')

    r_project_backward_settings = models.OneToOneField(
        BlastSettings,
        on_delete=models.CASCADE,
        verbose_name="settings for the backward BLAST execution",
        related_name='r_project_backward_settings')

    # each project can have one forward BlastDatabase
    r_project_forward_database = models.CharField(
        max_length=20,
        default='nr',
        verbose_name="BLAST program that is used inside snakemake execution",
    )

    r_entrez_query = models.CharField(
        max_length=500,
        default='',
        blank=True, null=True,
        verbose_name="Optional ENTREZ query for filtering the associated remote database",
    )

    r_project_backward_database = models.ForeignKey(
        BlastDatabase,
        on_delete=models.CASCADE,
        verbose_name="associated backward BLAST database",
        related_name="r_project_backward_database"
    )

    r_species_name_for_backward_blast = models.CharField(
        max_length=200, blank=False,
        verbose_name="species name for the backward database and query sequences"
    )

    # one to one relationship
    r_project_execution_snakemake_task = models.OneToOneField(
        TaskResult,
        on_delete=models.SET_NULL,
        blank=True, null=True,
        related_name="r_project_execution_snakemake_task",
        verbose_name="django_celery_results taskresult model for this projects snakemake pipeline")

    r_project_database_statistics_task = models.OneToOneField(
        TaskResult,
        on_delete=models.SET_NULL,
        blank=True, null=True,
        related_name="r_project_database_statistics_task",
        verbose_name="django_celery_results taskresult model for this projects database statistics")

    r_project_database_statistics_task_selection = models.OneToOneField(
        TaskResult,
        on_delete=models.SET_NULL,
        blank=True, null=True,
        related_name="r_project_database_statistics_task_selection",
        verbose_name="django_celery_results taskresult model for this "
                     "projects database statistics phylogenetic inference based on user selection task")

    # customized initialization can be added in BlastProjectManager (e.g. direct creation of project directory
    objects = RemoteBlastProjectManager()

    # overwritten functions
    def __str__(self):
        return "Reciprocal BLAST Project, created {} by {} with fw db {}".format(
            self.r_timestamp, self.r_project_user.username, self.r_project_forward_database)

    '''get_list_of_query_sequences

        Returns a list of query sequences without the additional .version of the query.

        :param self
            :type int
        :param filepath
            :type str

        :returns qseqids
            :type list[str]
    '''

    def get_list_of_query_sequences(self, filepath=REMOTE_BLAST_PROJECT_DIR):
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

    '''get_fasta_header_of_query_sequences

        This function extract information of the protein fasta header.

        :param self
            :type int
        :param filepath
            :type str

    '''

    def get_fasta_header_of_query_sequences(self, filepath=REMOTE_BLAST_PROJECT_DIR):
        try:
            query_sequence_file_path = self.get_project_query_sequence_filepath(filepath)
            query_file = open(query_sequence_file_path, 'r')
            qseqids_info = []
            for line in query_file.readlines():
                if ">" in line:
                    qseqid_info = ' '.join(line.split(">")[1].split(" ")[1:])
                    qseqid_info = qseqid_info.rstrip()
                    if len(qseqid_info) <= 1:
                        qseqid_info = "no information provided"
                    qseqids_info.append(qseqid_info)
            return qseqids_info
        except Exception as e:
            raise IntegrityError(
                "[-] couldnt extract query sequence ids from query sequence file : {} with exception : {}".format(
                    query_sequence_file_path, e))

    def get_project_username(self):
        return self.r_project_user.name

    def get_project_useremail(self):
        return self.r_project_user.email

    def get_project_dir(self):
        return REMOTE_BLAST_PROJECT_DIR + str(self.id)

    def get_project_query_sequence_filepath(self, filepath=REMOTE_BLAST_PROJECT_DIR):
        return filepath + str(self.id) + '/' + self.r_project_query_sequences

    def if_executed_return_associated_taskresult_model(self):
        # executed
        if self.r_project_execution_snakemake_task != None:
            return self.r_project_execution_snakemake_task
        # not executed
        else:
            return None

    ''' initialize_project_directory

        Invocation is done with BlastProject.objects.initialize_project_directory() or 
        within the BlastProjectManager with blast_project.iniialize_project_directory()

        Throws an exception if the specified filepath does not exists or if the 
        filepath + project_id directory exists.
        This function is executed within the create_blast_project function of the manager class, 
        it directly creates a directory for the project and a directory for result images in the static folder

    '''

    def initialize_project_directory(self, filepath=REMOTE_BLAST_PROJECT_DIR):
        # check if blast_project was previously created / check if media/blast_project directory exists
        if (isdir(filepath + str(self.id)) == True):
            raise IntegrityError("project directory exists")
        elif isdir(filepath) == False:
            raise IntegrityError("{} is not a directory".format(filepath))
        else:
            try:
                mkdir(filepath + str(self.id))
            except Exception as e:
                raise IntegrityError("couldnt create project directory : {}".format(e))

    '''write_snakemake_configuration_file

        This function is used during project creation. It writes all necessary settings for snakemake into
        the file snakefile_config in BLAST_PROJECT_DIR. Additional parameters that want to be used during
        the snakemake run have to be added here. The default filepath is media/blast_projects/.

        :params filepath=BLAST_PROJECT_DIR
            :type str
    '''

    def write_snakemake_configuration_file(self, project_settings, filepath=REMOTE_BLAST_PROJECT_DIR):
        try:
            with open(filepath + str(self.id) + '/snakefile_config', 'w') as snk_config_file:
                # database path from media/blast_projects/project_id as working directory for snakemake
                snk_config_file.write('project_id: ' + str(self.id) + "\n")
                snk_config_file.write('blastdb: ' + str(
                    self.r_project_forward_database) + "\n")
                snk_config_file.write('entrez_query:' + "\n")
                snk_config_file.write('backwarddb: ' + "\"" + "../../../databases/" + str(
                    self.r_project_backward_database.id) + "/" + self.r_project_backward_database.get_pandas_table_name() + ".database\"\n")
                snk_config_file.write('query_sequence: ' + "\"" + self.r_project_query_sequences + "\"\n")
                snk_config_file.write('bw_taxid: ' + str(self.r_species_name_for_backward_blast[1]) + "\n")
                snk_config_file.write('user_email: ' + str(self.r_project_user.email) + "\n")
                bw_dict = self.r_project_backward_settings.values_as_fw_or_bw_dict('bw')
                fw_dict = self.r_project_forward_settings.values_as_fw_or_bw_dict('fw')

                for key_bw in bw_dict.keys():
                    snk_config_file.write(key_bw + ': ' + bw_dict[key_bw] + "\n")

                for key_fw in fw_dict.keys():
                    snk_config_file.write(key_fw + ': ' + fw_dict[key_fw] + "\n")

                snk_config_file.write(
                    "bitscore_filter: " + str(project_settings.cleaned_data['bitscore_filter']) + "\n")
                snk_config_file.write("max_rbhs_for_phylo: " + str(
                    project_settings.cleaned_data['max_amount_of_rbh_for_msa_and_phylogeny']) + "\n")
                snk_config_file.write("trimal_gt: " + str(project_settings.cleaned_data['trimal_gt']) + "\n")
                snk_config_file.write("trimal_st: " + str(project_settings.cleaned_data['trimal_st']) + "\n")
                snk_config_file.write("trimal_cons: " + str(project_settings.cleaned_data['trimal_cons']) + "\n")
                snk_config_file.write("mview_sort: " + str(project_settings.cleaned_data['mview_sort']) + "\n")
                snk_config_file.write("mview_coloring: " + str(project_settings.cleaned_data['mview_coloring']) + "\n")

        except Exception as e:
            raise IntegrityError(
                "couldnt write snakemake configuration file in directory with exception : {}".format(e))

    '''read_query_information_table

        This function is getting executed within the project_details_dashboard.html website.

    '''

    def read_query_information_table(self, filepath=REMOTE_BLAST_PROJECT_DIR):
        def clean_feature_column(features):
            new_feature_column = []
            try:
                result_string = ''
                for idx, feature in enumerate(literal_eval(features)):
                    if idx % 2 == 0:
                        result_string += feature + " "
                    # TODO add functional linebreak
                    elif idx % 2 == 1:
                        result_string += feature + "     "
                new_feature_column.append(result_string)
                return pd.Series(new_feature_column)
            except Exception as e:
                raise Exception("ERROR:exception: {}".format(e))

        try:
            path_to_information_table = filepath + str(self.id) + '/query_sequence_information.csv'
            if isfile(path_to_information_table):
                table = pd.read_table(path_to_information_table, header=0, sep="\t", index_col=0)
                table.Features = table.Features.apply(clean_feature_column)
                table = table.fillna(value='')
                table = table.to_html(classes='my_class" id="myTable')
                return table
            else:
                return "there is no query_sequence_information.csv in the project directory"
        except Exception as e:
            raise Exception(
                "[-] ERROR during pandas parsing of query_sequence_information csv file with exception: {}".format(e))

    '''check_for_reciprocal_result_table

        This function checks if the reciprocal_result.html file is in the project directory or not.
    '''

    def check_for_reciprocal_result_table(self):
        if isfile(self.get_project_dir() + '/' + 'reciprocal_results.html'):
            return True
        else:
            return False

    '''return_list_of_all_logfiles

        Returns a list with all filepaths to logfiles of the snakemake run. 
        If there is no "log" directory, in the BLAST_PROJECT_DIR the function returns an empty list.

        :returns filelist
            :type list
    '''

    def return_list_of_all_logfiles(self) -> list:
        try:
            path_to_project_dir = REMOTE_BLAST_PROJECT_DIR + str(self.id) + "/log"
            if isdir(path_to_project_dir) == False:
                return []
            else:
                files = listdir(path_to_project_dir)
                filelist = []
                for file in files:
                    # check if "file" is a query sequence sub-directory
                    if isdir(path_to_project_dir + '/' + file):
                        for query_log in listdir(path_to_project_dir + '/' + file):
                            filelist.append(file + '/' + query_log)
                    else:
                        filelist.append(file)
                return filelist
        except Exception as e:
            raise Exception(
                "[-] ERROR creating list of all logfiles for project: {} with exception: {}".format(self.id, e))
