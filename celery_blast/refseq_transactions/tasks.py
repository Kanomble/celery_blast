from os import getcwd, chdir, mkdir, remove
from os.path import isdir, isfile
from subprocess import Popen, PIPE as subPIPE, STDOUT as subSTDOUT

from blast_project.py_django_db_services import get_database_by_id, update_blast_database_with_task_result_model
from .py_services import write_blastdatabase_snakemake_configfile,\
    get_ftp_paths_and_taxids_from_summary_file, write_progress_database_transactions,\
    get_current_progress_database_transactions, get_bdb_summary_table_name

from time import sleep
from celery import shared_task
from celery.utils.log import get_task_logger
from celery_progress.backend import ProgressRecorder
#logger for celery worker instances
logger = get_task_logger(__name__)

#TODO documentation task 3
@shared_task()
def format_available_databases(path_to_database, database_files_taxids_dict, progress_recorder):
    try:
        progress_steps = 100 / ((len(database_files_taxids_dict.keys())) * 2)
        progress = 50
        progress_recorder.set_progress(progress, 100, "starting makeblastdb processes")
        for assembly_file in database_files_taxids_dict.keys():
            database = path_to_database + assembly_file
            taxid = database_files_taxids_dict[assembly_file]
            cmd = "makeblastdb -in {} -dbtype prot -out {} -taxid {}".format(database,database,taxid)

            if (isfile(database+'.pdb') == True):
                logger.info('file {} exists, scipping formatting'.format(database))
            else:
                proc = Popen("makeblastdb -in {} -dbtype prot -out {} -taxid {} -parse_seqids"
                             .format(database,database,taxid), shell=True)
                returncode = proc.wait(timeout=600)#10*60 Sec.
                logger.info('processed database with makeblastdb : {}'.format(cmd))
                if(returncode != 0):
                    logger.warning("something went wrong during creation of database: {}".format(database))

            progress += progress_steps
            progress_recorder.set_progress(progress, 100, "starting makeblastdb processes")

        progress_recorder.set_progress(99, 100, "ended database formatting")
        return database_files_taxids_dict.keys()
    except Exception as e:
        logger.warning('couldnt format downloaded assembly files to blast databases with exception : {}'.format(e))
        raise Exception('couldnt format downloaded assembly files to blast databases with exception : {}'.format(e))

#TODO documentation task 2
@shared_task()
def download_wget_ftp_paths(path_to_database,dictionary_ftp_paths_taxids,progress_recorder):
    try:
        error_log = open(path_to_database+'download_error.log', 'w')
        transform_ftp_path = lambda file: file.split('/')[-1].rstrip(file[-3:])

        #progress_steps = round(100 / (len(dictionary_ftp_paths_taxids.keys()) * 2), 3)

        downloaded_files = {}
        #database_id for format_available_databases

        progress_steps = 100/((len(dictionary_ftp_paths_taxids.keys()))*2)
        progress = 0
        progress_recorder.set_progress(progress, 100, "starting wget processes")

        for file in dictionary_ftp_paths_taxids.keys():

            gunzip_output = transform_ftp_path(file)
            if(isfile(path_to_database + gunzip_output) == True):
                downloaded_files[gunzip_output] = dictionary_ftp_paths_taxids[file]
                logger.info('file {} exists, scipping download'.format(gunzip_output))
                progress += progress_steps
                progress_recorder.set_progress(progress, 100, "downloaded {}".format(gunzip_output))
            elif(isfile(path_to_database + gunzip_output) == False):
                for attempt in range(10):
                    try:


                        proc = Popen('wget -qO- {} | gzip -d > {}'.format(file, path_to_database+gunzip_output), shell=True)
                        returncode = proc.wait(timeout=300)  # 66 Minutes
                        if(returncode != 0):
                            raise Exception
                        #downloaded_files[gunzip_output] = dictionary_ftp_paths_taxids[file]
                        #logger.info('downloaded : {} returncode : {}'.format(path_to_database + gunzip_output,returncode))


                    except Exception as e:
                        logger.warning("download exception : {}".format(e))
                        error_log.write("{} {}\n".format(file, attempt))
                        logger.warning('next download attempt of file : {} with attempt : {}'.format(file,attempt))
                    else:
                        logger.info('downloaded : {} returncode : {}'.format(path_to_database + gunzip_output,returncode))
                        downloaded_files[gunzip_output] = dictionary_ftp_paths_taxids[file]
                        progress += progress_steps
                        progress_recorder.set_progress(progress, 100, "downloaded {}".format(gunzip_output))
                        break
            else:
                error_log.write('couldnt download: {} '.format(file))
                logger.warning('couldnt download: {} after 10 attempts'.format(file))
                progress += progress_steps
                progress_recorder.set_progress(progress, 100, "failed trying to download {}".format(file))
                #error_log.close()
                #raise Exception
        error_log.close()
        return downloaded_files
    except Exception as e:
        logger.warning('couldnt download assemblies with exception : {}'.format(e))
        raise Exception('couldnt download assemblies with exception : {}'.format(e))

#TODO documentation task 1
@shared_task(bind=True)
def download_blast_databases_based_on_summary_file(self, database_id):
    try:
        progress_recorder = ProgressRecorder(self)
        logger.info('starting downloading task')
        working_dir = getcwd()
        logger.info('working dir : {}'.format(working_dir))
        path_to_database = 'media/databases/' + str(database_id) + '/'




        progress_recorder.set_progress(0,100,"started downloading")

        dictionary_ftp_paths_taxids = get_ftp_paths_and_taxids_from_summary_file(database_id)

        logger.info('trying to update blast database with current task')
        update_blast_database_with_task_result_model(database_id, str(self.request.id))

        dictionary_ftp_paths_taxids = download_wget_ftp_paths(path_to_database,dictionary_ftp_paths_taxids,progress_recorder)

        database_files = format_available_databases(path_to_database,dictionary_ftp_paths_taxids,progress_recorder)
        progress_recorder.set_progress(99, 100, "writing alias file")
        database_pandas_table_name = get_bdb_summary_table_name(database_id)
        alias_filename = 'media/databases/' + str(database_id) + '/' + database_pandas_table_name+'.database.pal'
        logger.info('starting to write database alias file : {}'.format(alias_filename))
        alias_file = open(alias_filename,'w')
        alias_file.write("TITLE {}\n".format(database_pandas_table_name+'.database'))
        alias_file.write("DBLIST")
        for database_file in database_files:
            alias_file.write(" \""+database_file+"\"")
        alias_file.write("\n")
        alias_file.close()
        progress_recorder.set_progress(100, 100, "writing alias file")
        return 0
    except Exception as e:
        logger.warning('couldnt perform task because of exception : {}'.format(e))
        raise Exception('couldnt perform compute ftp path task with exception {}'.format(e))


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
        snakemake_process = Popen(['snakemake',
                                   '--snakefile',snakefile_dir,
                                   '--wms-monitor','http://172.23.0.6:5000',
                                   '--cores','1',
                                   '--configfile',snakemake_config_file,
                                   '--directory',snakemake_working_dir,
                                   '--latency-wait','10',
                                   '--quiet',
                                   '--keep-incomplete'], shell=False, stdout=subPIPE, stderr=subSTDOUT)
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
        path_to_assembly_file_location = current_working_directory + '/media/databases/refseq_summary_file/'
        timeout=300

        logger.info('setup filepath parameter:\n\t cwd : {} \n\t path_to_assembly_file_location : {}'
                    .format(current_working_directory,path_to_assembly_file_location))

        if(isdir(path_to_assembly_file_location) == False):
            logger.warning('path_to_assembly_file_location : {} does not exists, trying to create it with mkdir ...')
            mkdir(path_to_assembly_file_location)

        path_to_assembly_file = path_to_assembly_file_location+'assembly_summary_refseq.txt'
        if(isfile(path_to_assembly_file)):
            logger.warning('assembly_summary_refseq.txt exists deleting it in order to download a newer version')
            remove(path_to_assembly_file)

        # invoke wget program

        wget_process = Popen(['wget', refseq_url,'-O',path_to_assembly_file])
        # communicate with subprocess : https://docs.python.org/3/library/subprocess.html#subprocess.Popen.communicate
        # wait for process to terminate and set returncode attribute
        logger.info(
            'waiting for popen instance {} to finish with timeout set to {}'
                .format(wget_process.pid, timeout))
        returncode = wget_process.wait(timeout=timeout)
        logger.info('returncode : {}'.format(returncode))


        logger.info('download completed')

        return refseq_url.split('/')[-1]
    except Exception as e:
        raise Exception("couldn't download assembly_summary_refseq.txt file : {}".format(e))
