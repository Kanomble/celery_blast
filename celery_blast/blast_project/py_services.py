import os

import pandas as pd

from .models import BlastProject, RemoteBlastProject
from refseq_transactions.models import BlastDatabase
from external_tools.models import DomainDatabase

from os.path import isdir, isfile, join, relpath
from os import mkdir, listdir, remove, rmdir, walk
from subprocess import check_output
from shutil import rmtree, make_archive
from wsgiref.util import FileWrapper
from django.db import IntegrityError, transaction
from celery_blast.settings import BLAST_PROJECT_DIR, BLAST_DATABASE_DIR, CDD_DATABASE_URL, TAXDB_URL, REMOTE_BLAST_PROJECT_DIR
import zipfile
import io

'''list_all_available_logfiles
    
    This functions uses the os module to get a list of all available logfiles for a certain reciprocal BLAST
    project.
    
    :param path_to_logfiles
        :type str
    :param path_to_task_table
        :type str
    
    :returns direct_logfiles, query_specific_logfiles
        :type tuple(list[list, list], dict)
'''
def list_all_available_logfiles(path_to_logfiles: str, path_to_task_table) -> tuple:
    try:
        direct_logfiles = {}
        query_specific_logfiles_dir = []
        if isdir(path_to_logfiles):
            for file in listdir(path_to_logfiles):
                if isdir(path_to_logfiles + file) == False:
                    direct_logfiles[file] = path_to_logfiles + file
                else:
                    query_specific_logfiles_dir.append(path_to_logfiles + file + '/')
            query_specific_logfiles = {}
            for dir in query_specific_logfiles_dir:
                key = dir.split("/")[-2]
                query_specific_logfiles[key] = []
                for file in os.listdir(dir):
                    query_specific_logfiles[key].append(file)
        else:
            raise Exception("[-] There is no such path available on the server: {}".format(path_to_logfiles))

        logfile_table = pd.read_csv(path_to_task_table)
        direct_logfiles = logfile_table[logfile_table.logfile.isin(list(direct_logfiles.keys()))]
        database_logfiles = direct_logfiles[direct_logfiles.progress == 0.0]
        direct_logfiles = direct_logfiles[direct_logfiles.progress > 0.0]

        direct_logfiles = zip(list(direct_logfiles.logfile), list(direct_logfiles.progress))
        database_logfiles = zip(list(database_logfiles.logfile), list(database_logfiles.progress))
        return direct_logfiles, query_specific_logfiles, database_logfiles
    except Exception as e:
        raise Exception("[-] ERROR during logfile listing with exception: {}".format(e))


'''download_project_directory
    
    This function compresses the specified directory into a .zip file and creates the 
    class FileWrapper.
    
    :param base_directory
        :type str
    
    :returns zip_buffer
        :type io.BytesIO
'''


def download_project_directory(base_directory:str)->io.BytesIO:
    try:
        def add_files_to_zip(zip_file, directory):
            for root, _, filenames in walk(directory):
                for filename in filenames:
                    file_path = join(root, filename)
                    relative_path = relpath(file_path, directory)
                    zip_file.write(file_path, relative_path)
        def generate_zip_archive(directory):
            zip_buffer = io.BytesIO()
            with zipfile.ZipFile(zip_buffer, 'a', zipfile.ZIP_DEFLATED) as zip_file:
                add_files_to_zip(zip_file, directory)
            zip_buffer.seek(0)
            return zip_buffer
        # Create an in-memory zip archive
        zip_buffer = generate_zip_archive(base_directory)
        return zip_buffer
    except Exception as e:
        raise Exception("[-] ERROR creating zip directory: {} with exception: {}".format(base_directory, e))

'''read_task_logs_summary_table
    
    This function loads the task_logfile.txt file into a pandas dataframe.
    It is used within the project_details view to track the progress of the snakemake pipeline.
    
    :param remote_or_local
        :type str
        
    :returns logfiles_table
        :type pd.DataFrame
'''


def read_task_logs_summary_table(remote_or_local:str) -> pd.DataFrame:
    try:
        if remote_or_local == "local":
            data_path = BLAST_PROJECT_DIR + 'task_logfiles'
        elif remote_or_local == "remote":
            data_path = REMOTE_BLAST_PROJECT_DIR + 'task_logfiles'
        else:
            raise Exception("[-] project is neither local nor remote ...")
        logfiles_table = pd.read_table(data_path, sep=",", header=0)
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

def delete_remote_project_and_associated_directories_by_id(project_id: int) -> None:
    try:
        with transaction.atomic():
            project = RemoteBlastProject.objects.get(id=project_id)
            if isdir(REMOTE_BLAST_PROJECT_DIR + str(project_id)):
                rmtree(REMOTE_BLAST_PROJECT_DIR + str(project_id))
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

def get_remote_html_results(project_id: int, filename: str, html_result_path=REMOTE_BLAST_PROJECT_DIR) -> list:
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


'''check_blast_database_integrity
    
    This function checks if the specified database is correctly formatted.
    Similar to the check_domain_database_status() function.
    
    :param database_id
        :type int
    :returns returncode
        :type bool

'''
def check_blast_database_integrity(database_id:int)->bool:
    try:
        returncode = False
        blastdb = BlastDatabase.objects.get(id=database_id)

        if isdir(BLAST_DATABASE_DIR + str(database_id)):
            blastdb_path = BLAST_DATABASE_DIR + str(database_id) + "/" + blastdb.get_pandas_table_name() + ".database"
            out = check_output(['blastdbcheck', '-db', blastdb_path])
            if "Result=SUCCESS. No errors reported for 1 volume(s)" in str(out) or "Result=SUCCESS. No errors reported for 1 alias(es)" in str(out):
                returncode = True
            else:
                returncode = False
        else:
            returncode = False

        return returncode
    except Exception as e:
        raise Exception("[-] ERROR during check up of BLAST database integrity with ID: {} and exception: {}".format(database_id, e))

'''check_domain_database_status
    This function checks if the domain database is loaded by utilizing the subprocess
    check_output function in conjunction with the blastdbcheck tool from the BLAST+ software suite.
    Sets the "domain_database_loaded" variable of the DomainDatabase class to True if the database directory
    is functional. If not it sets the value to True.
    
    :returns returncode
        :type bool
'''
def check_domain_database_status()->bool:
    try:
        cdd_db_path = BLAST_DATABASE_DIR + "CDD/"
        # testing domain database integrity

        if isdir(cdd_db_path + 'Cdd/') :
                try:
                    out = check_output(['blastdbcheck', '-db', cdd_db_path + 'Cdd'])
                except Exception as e:
                    out = str(e)

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
        else:
            if isdir(cdd_db_path) == False:
                mkdir(cdd_db_path)
            domain_database_query_set = DomainDatabase.objects.all()
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

'''delete_domain_database

    This function is executed within the delete_domain_database_view function.
    It deletes the domain database sitting in BLAST_DATABASE_DIR + "CDD/" and the associated 
    DomainDatabase model. It then creates a new model with the attribute domain_database_loaded
    set to false.

'''
def delete_domain_database():
    try:
        cdd_db_path = BLAST_DATABASE_DIR + "CDD/"

        if isdir(cdd_db_path):
            for file in listdir(cdd_db_path):
                file_to_remove = cdd_db_path + file
                remove(file_to_remove)
            rmdir(cdd_db_path)

        DomainDatabase.objects.all().delete()
        domain_database_model = DomainDatabase(domain_database_loaded=False)
        domain_database_model.save()
        return 0
    except Exception as e:
        raise Exception("ERROR during deletion of domain database in: {} with exception: {}".format(cdd_db_path, e))