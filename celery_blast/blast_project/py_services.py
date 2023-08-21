import pandas as pd

from .models import BlastProject
from refseq_transactions.models import BlastDatabase
from external_tools.models import DomainDatabase

from os.path import isdir, isfile
from os import mkdir, listdir, remove, rmdir
from subprocess import check_output
from shutil import rmtree, make_archive
from wsgiref.util import FileWrapper
from django.db import IntegrityError, transaction
from celery_blast.settings import BLAST_PROJECT_DIR, BLAST_DATABASE_DIR, CDD_DATABASE_URL, TAXDB_URL

'''download_project_directory
    
    This function compresses the specified directory into a .zip file and creates the 
    class FileWrapper.
    
    :param directory
        :type str
    
    :returns file_wrapper_archive
        :type FileWrapper
'''


def download_project_directory(directory:str)->FileWrapper:
    try:
        path_to_zip = make_archive(directory,"zip",directory)
        return FileWrapper(open(path_to_zip, 'rb'))
    except Exception as e:
        raise Exception("[-] ERROR creating zip directory: {} with exception: {}".format(directory, e))

'''read_task_logs_summary_table
    
    This function loads the task_logfile.txt file into a pandas dataframe.
    It is used within the project_details view to track the progress of the snakemake pipeline.
    
    :returns logfiles_table
        :type pd.DataFrame
'''


def read_task_logs_summary_table() -> pd.DataFrame:
    try:
        data_path = BLAST_PROJECT_DIR + 'task_logfiles'
        logfiles_table = pd.read_table(data_path, sep="\t", header=0)
        return logfiles_table
    except Exception as e:
        raise Exception("[-] ERROR reading task_logfile.txt file in {} with exception: {}".format(data_path, e))


'''check_if_taxdb_exists
    
    This function is used before project creation to check if the taxonomy database exists in the database directory.
    
    :returns True or False
        :type boolean
'''


def check_if_taxdb_exists() -> bool:
    if isfile(BLAST_DATABASE_DIR + 'taxdb.btd') and isfile(BLAST_DATABASE_DIR + 'taxdb.bti'):
        return True
    else:
        return False


def check_if_file_exists(file_path) -> bool:
    return isfile(file_path)


'''list_taxonomic_files

utilization in create_taxonomic_file_view and refseqdatabaseform 
returns a list of all files and their corresponding total line length in the media/taxonomic_node_files folder that end with .taxids

'''


def list_taxonomic_files():
    try:
        files_in_taxonomic_node_files = listdir('media/taxonomic_node_files/')
        files = []
        length = []
        for file in files_in_taxonomic_node_files:
            lines = 0
            with open('media/taxonomic_node_files/' + file) as f:
                for line in f:
                    lines = lines + 1
            if file.endswith('.taxids'):
                files.append(file)
                length.append(lines)
        # [file for file in files_in_taxonomic_node_files if file.endswith('.taxids')]
        return files, length
    except Exception as e:
        raise Exception('exception ocurred in blast_project/py_services.list_taxonomic_files : {}'.format(e))


''' delete_refseqgenome_and_associated_directories_by_id

    deletes the blastdatabase entry by id
    removes the blastdatabase directory
    
    raises an integrity error if exception occurs
    
    :param database_id
        :type int
'''


def delete_blastdb_and_associated_directories_by_id(database_id):
    try:
        with transaction.atomic():
            blastdatabase = BlastDatabase.objects.get(id=database_id)
            if isdir(BLAST_DATABASE_DIR + str(database_id)):
                rmtree(BLAST_DATABASE_DIR + str(database_id))
            blastdatabase.delete()
    except Exception as e:
        raise IntegrityError("couldnt delete blast database entry : {}".format(e))


'''delete_project_and_associated_directories_by_id

    This function deletes the BlastProject with id project_id from the database.
    It also removes all associated directories with the rmtree shutils function.
    It is used in the delete_project_view in blast_project.views.py.
    
    :param porject_id
        :type int
    
'''


def delete_project_and_associated_directories_by_id(project_id: int) -> None:
    try:
        with transaction.atomic():
            project = BlastProject.objects.get(id=project_id)
            if isdir(BLAST_PROJECT_DIR + str(project_id)):
                rmtree(BLAST_PROJECT_DIR + str(project_id))
            project.delete()
    except Exception as e:
        raise IntegrityError("couldnt delete blast project entry : {}".format(e))


''' create_blastdatabase_directory
    
    creates an directory with the name database_id in the /media/databases folder 
    
    :param database_id
        :type int
'''


def create_blastdatabase_directory(database_id, database_filepath=BLAST_DATABASE_DIR):
    try:
        mkdir(database_filepath + str(database_id))
        return database_filepath + str(database_id)
    except Exception as e:
        raise IntegrityError(
            'something went wrong during database directory creation: {}'.format(e))


'''upload_file
    simple function for uploading a file to a specific server side location, defined by destination
'''


def upload_file(project_file, destination: str):
    try:
        with open(destination, 'wb+') as dest:
            for chunk in project_file.chunks():
                dest.write(chunk)
    except Exception as e:
        raise IntegrityError(
            'exception during file upload of : {} : exception : {}'.format(project_file.name, e))


'''get_html_results
    This function is used in blast_project.views.load_reciprocal_result_html_table_view, it returns
    the plain HTML as a string. 
    
    :param project_id
        :type int
    :param filename - html result filename
        :type str
    :param html_result_path
        :type str
    
    :returns data - string representation of a pandas html table
        :type list[str]
'''


def get_html_results(project_id: int, filename: str, html_result_path=BLAST_PROJECT_DIR) -> list:
    try:
        with open(html_result_path + str(project_id) + "/" + filename) as res:
            data = res.readlines()
        return data
    except Exception as e:
        raise FileNotFoundError("[-] ERROR: Couldn't read file {} with Exception: {}".format(filename, e))


'''html_table_exists
    Checks if the specified file exists.
    
    :param project_id
        :type int
    :param filename
        :type str
    :param html_result_path
        :type str
        
    :returns True/False
        :type bool
'''


def html_table_exists(project_id:int, filename:str, html_result_path=BLAST_PROJECT_DIR)->bool:
    if (isfile(html_result_path + str(project_id) + "/" + filename)):
        return True
    else:
        return False


'''concatenate_genome_fasta_files_in_db_dir
    This function is used during uploading of custom genome files.
    It concatenates all genome files into a database file, that is used by makeblastdb.
    
    :param path_to_database
        :type str
    :param database_title
        :type str
    :param genome_files
        :type list
'''
def concatenate_genome_fasta_files_in_db_dir(path_to_database:str, database_title:str, genome_files:list):
    try:
        database_name = database_title.replace(' ', '_').upper() + '.database'
        with open(path_to_database + database_name, 'w') as dbfile:
            for file in genome_files:
                with open(path_to_database + file, 'r') as gfile:
                    for line in gfile.readlines():
                        dbfile.write(line)
    except Exception as e:
        raise IntegrityError('couldnt concatenate database files : {}'.format(e))


'''check_domain_database_status
    This function checks if the domain database is loaded by utilizing the subprocess
    check_output function in conjunction with the blastdbcheck tool from the BLAST+ software suite.
    Sets the domain_database_loaded variable of the DomainDatabase class to True if the database directory
    is healthy.
    
    :returns returncode
        :type bool
'''
def check_domain_database_status():
    try:
        cdd_db_path = BLAST_DATABASE_DIR + "CDD/"
        # testing domain database integrity
        out = check_output(['blastdbcheck', '-db', cdd_db_path + 'Cdd'])
        returncode = False
        # integrity test failed
        domain_database_query_set = DomainDatabase.objects.all()

        if 'No errors reported for 1 alias(es)' in str(out):
            if len(domain_database_query_set) != 1:
                DomainDatabase.objects.all().delete()
                domain_database_model = DomainDatabase(domain_database_loaded=True)
                domain_database_model.save()
            else:
                domain_database_query_set[0].domain_database_loaded = True
                domain_database_query_set[0].save()
                returncode = True
        else:
            if len(domain_database_query_set) != 1:
                DomainDatabase.objects.all().delete()
                domain_database_model = DomainDatabase(domain_database_loaded=False)
                domain_database_model.save()
            else:
                domain_database_query_set[0].domain_database_loaded = False
                domain_database_query_set[0].save()
                returncode = False

        return returncode
    except Exception as e:
        raise Exception("check_domain_database_status has thrown an error: {}".format(e))

def delete_domain_database():
    try:
        cdd_db_path = BLAST_DATABASE_DIR + "CDD/"

        for file in listdir(cdd_db_path):
            file_to_remove = cdd_db_path + file
            remove(file_to_remove)
        rmdir(cdd_db_path)
        domain_database_query_set = DomainDatabase.objects.all()
        if len(domain_database_query_set) != 1:
            raise Exception("DomainDatabase model in PostgreSQL is corrupted, there are multiple entries.")
        else:
            domain_database_query_set[0].domain_database_loaded = False
            domain_database_query_set[0].save()

    except Exception as e:
        raise Exception("ERROR during deletion of domain database in: {} with exception: {}".format(cdd_db_path, e))