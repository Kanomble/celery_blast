from os import getcwd, chdir, mkdir, remove
from os.path import isdir, isfile
from subprocess import Popen, PIPE as subPIPE, STDOUT as subSTDOUT

from blast_project.py_django_db_services import get_database_by_id, update_blast_database_with_task_result_model
from .py_services import write_blastdatabase_snakemake_configfile
from django_celery_results.models import TaskResult
from time import sleep
from celery import shared_task, task, current_task
from celery.utils.log import get_task_logger
from celery_progress.backend import ProgressRecorder
#logger for celery worker instances
logger = get_task_logger(__name__)

#TODO documentation
#task_track_started --> no need for celery-progress?
@shared_task(bind=True)
def download_blast_databases(self, database_id):


    progress_recorder = ProgressRecorder(self)
    progress_recorder.set_progress(0, 100,"started process")

    try:
        update_blast_database_with_task_result_model(database_id, str(self.request.id))
    except Exception as e:
        logger.warning("couldnt save taskresult into blastdatabase ... ")
        raise Exception("coulndt save taskresult into blastdatabase ...")

    #check if blast summary file is available
    blast_database_summary_table = get_database_by_id(database_id).get_pandas_table_name()
    logger.info('load blastdatabase: {} from database'.format(blast_database_summary_table))
    bdb_summary_table_path = 'media/databases/' + str(database_id) + '/' + blast_database_summary_table
    if (isfile(bdb_summary_table_path) == False):
        logger.warning('there is no database table with path : {}', format(bdb_summary_table_path))
        raise Exception(
            "couldnt download assembly files, there is no summary table with path : {}".format(bdb_summary_table_path))

    # summary file and ProgressRecord
    try:
        write_blastdatabase_snakemake_configfile(database_id, str(self.request.id))
    except Exception as e:
        logger.warning("couldnt write snakemake configuration file")
        raise Exception("couldnt write snakemake configuration file")


    progress_recorder.set_progress(20, 100, "starting snakemake process")

    try:
        snakemake_working_dir = 'media/databases/' + str(database_id) + '/'
        snakemake_config_file = 'media/databases/' + str(database_id) + '/snakefile_config'
        snakefile_dir = 'static/snakefiles/assembly_download/Snakefile'
        #cmd to copy
        #snakemake --snakefile 'static/snakefiles/assembly_download/Snakefile' --cores 2 --configfile 'media/databases/12/snakefile_config' --directory 'media/databases/12/' --wms-monitor 'http://172.23.0.5:5000' --latency-wait 10 --dry-run
        snakemake_process = Popen(['snakemake','--snakefile',snakefile_dir,'--wms-monitor','http://172.23.0.5:5000','--cores','1','--configfile',snakemake_config_file,'--directory',snakemake_working_dir,'--latency-wait','10'], shell=False, stdout=subPIPE, stderr=subSTDOUT)
        logger.info('waiting for popen snakemake instance {} to finish'.format(snakemake_process.pid))
        progress_recorder.set_progress(40, 100, "waiting for snakemake to finish")
        returncode = snakemake_process.wait(timeout=4000) # 66 Minutes

        if (returncode != 0):
            logger.warning('subprocess Popen snakemake instance resulted in an error!')
            raise Exception('Popen hasnt succeeded ...')

        logger.info('returncode : {}'.format(returncode))
        logger.info('blast database download and formatting completed')
        progress_recorder.set_progress(100, 100, "complete")
        return returncode
    except Exception as e:
        raise Exception("couldnt download assembly files, error during execution of snakemake, with exception : ".format(e))

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
