from os import getcwd, chdir, mkdir, remove
from os.path import isdir, isfile
from subprocess import Popen

from blast_project.py_django_db_services import get_database_by_id
from .py_services import write_blastdatabase_snakemake_configfile

from celery import shared_task
from celery.utils.log import get_task_logger

#logger for celery worker instances
logger = get_task_logger(__name__)

@shared_task
def download_blast_databases(database_id):
    blast_database_summary_table = get_database_by_id(database_id).get_pandas_table_name()
    logger.info('load blastdatabase: {} from database'.format(blast_database_summary_table))

    #check if blast summary file is available
    bdb_summary_table_path = 'media/databases/'+str(database_id)+'/'+blast_database_summary_table
    if(isfile(bdb_summary_table_path ) == False):
        logger.warning('there is no database table with path : {}',format(bdb_summary_table_path))
        raise Exception("couldnt download assembly files, there is no summary table with path : {}".format(bdb_summary_table_path ))

    write_blastdatabase_snakemake_configfile(database_id)


''' download_refseq_assembly_summary_file
    
    this function downloads the current assembly_summary_refseq.txt file into 
    the media directory /databases/refseq_summary_file/ if this directory does not exists the function
    creates the directory
    
    it uses subprocess.Popen to invoke wget, stdout (e.g. percentage is printed out inside celery worker console)
    
    :returns download_file
        :type str
'''
@shared_task
def download_refseq_assembly_summary_file():
    try:
        refseq_url = "ftp://ftp.ncbi.nih.gov/genomes/refseq/assembly_summary_refseq.txt"
        current_working_directory = getcwd() #/blast/reciprocal_blast
        path_to_assembly_file_location = current_working_directory + '/media/databases/refseq_summary_file'
        timeout=300

        logger.info('setup filepath parameter:\n\t cwd : {} \n\t path_to_assembly_file_location : {}'
                    .format(current_working_directory,path_to_assembly_file_location))

        if(isdir(path_to_assembly_file_location) == False):
            logger.warning('path_to_assembly_file_location : {} does not exists, trying to create it with mkdir ...')
            mkdir(path_to_assembly_file_location)

        chdir(path_to_assembly_file_location)

        if(isfile('assembly_summary_refseq.txt')):
            logger.warning('assembly_summary_refseq.txt exists deleting it in order to download a newer version')
            remove('assembly_summary_refseq.txt')

        # invoke wget program
        wget_process = Popen(['wget', refseq_url])
        # communicate with subprocess : https://docs.python.org/3/library/subprocess.html#subprocess.Popen.communicate
        # wait for process to terminate and set returncode attribute
        logger.info(
            'waiting for popen instance {} to finish with timeout set to {}'
                .format(wget_process.pid, timeout))
        returncode = wget_process.wait(timeout=timeout)
        logger.info('returncode : {}'.format(returncode))


        chdir(current_working_directory)

        logger.info('download completed')

        return refseq_url.split('/')[-1]
    except Exception as e:
        raise Exception("couldn't download assembly_summary_refseq.txt file : {}".format(e))
