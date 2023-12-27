# The best practice is to create a common logger for all of your tasks at the top of your module:

import os
import psutil
from subprocess import Popen, STDOUT as subSTDOUT, SubprocessError, TimeoutExpired
from django.conf import settings
from celery import shared_task
from external_tools.models import ExternalTools
from celery.utils.log import get_task_logger
from celery_progress.backend import ProgressRecorder
from celery.exceptions import SoftTimeLimitExceeded
from .py_django_db_services import update_blast_project_with_task_result_model, \
    update_blast_database_with_task_result_model, create_external_tools_after_snakemake_workflow_finishes, \
    update_blast_project_with_database_statistics_task_result_model, get_all_blast_databases, \
    update_remote_blast_project_with_task_result_model
from .py_database_statistics import calculate_database_statistics
from celery_blast.settings import BLAST_DATABASE_DIR, BLAST_PROJECT_DIR, TAXDB_URL, TAXONOMIC_NODES, \
    CDD_DATABASE_URL, REMOTE_BLAST_PROJECT_DIR

# logger for celery worker instances
logger = get_task_logger(__name__)

'''download_and_format_taxdb

    This task downloads the taxonomy database from the NCBI FTP site.
    WARNING: This function uses the default BLAST_DATABASE_DIR filepath, overwriting BLAST_DATABASE_DIR does not
    affect the function.
    
    :param self
        :type TaskResult object?
    :returns returncode
        :type int

'''


@shared_task(bind=True)
def download_and_format_taxdb(self):
    logger.info("INFO:NO TAXONOMY DATABASE")
    if os.path.isfile(BLAST_DATABASE_DIR + "taxdb.tar.gz"):
        os.remove(BLAST_DATABASE_DIR + "taxdb.tar.gz")
    if os.path.isfile(BLAST_DATABASE_DIR + "taxdb.tar"):
        os.remove(BLAST_DATABASE_DIR + "taxdb.tar")
    logger.info("INFO:STARTING TO DOWNLOAD TAXONOMY DATABASE")

    try:
        taxdb_ftp_path = TAXDB_URL
        path_to_taxdb_location = BLAST_DATABASE_DIR
        path_to_taxdb_location = path_to_taxdb_location + 'taxdb.tar.gz'

        proc = Popen(["wget", taxdb_ftp_path, "-q", "-O", path_to_taxdb_location], shell=False)
        returncode = proc.wait(timeout=600)
        if returncode != 0:
            raise SubprocessError
        logger.info("INFO:TRYING TO DECOMPRESS THE TAXONOMY DATABASE")

        proc = Popen(["tar", "-zxvf", path_to_taxdb_location, "-C", BLAST_DATABASE_DIR],
                     shell=False)
        returncode = proc.wait(timeout=600)
        if returncode != 0:
            raise SubprocessError
        logger.info("INFO:DONE")
        return returncode
    except TimeoutExpired as e:
        logger.warning("ERROR:TIMEOUT EXPIRED DURING DOWNLOAD OF TAXONOMY DATABASE: {}".format(e))
        logger.info(
            "INFO:IF YOU HAVE NO STABLE INTERNET CONNECTION TRY TO RESTART THE WEBSERVER ONCE YOU HAVE A BETTER CONNECTION")
        logger.info("INFO:YOU CAN MANUALLY LOAD THE TAXONOMY DATABASE INTO THE DATABASE FOLDER")

        pid = proc.pid
        parent = psutil.Process(pid)

        for child in parent.children(recursive=True):
            child.kill()
        parent.kill()
        raise TimeoutExpired(e)
    except SubprocessError as e:
        logger.warning("ERROR:WGET RESULTED IN AN ERROR: {}".format(e))
        logger.info("INFO:YOU CAN MANUALLY LOAD THE TAXONOMY DATABASE INTO THE DATABASE FOLDER")

        pid = proc.pid
        parent = psutil.Process(pid)

        for child in parent.children(recursive=True):
            child.kill()
        parent.kill()
        raise SubprocessError(e)
    except Exception as e:
        raise Exception(e)

''' get_species_taxids_into_file

spawns a child process via subprocess Popen and waits for its returncode 
it can be used to create a taxonomic node file in the project media folder

:param taxonomic_node
    :type int
:param taxids_filename
    :type str
:returns Popen.wait(timeout=200) returncode
'''


@shared_task
def write_species_taxids_into_file(taxonomic_node, taxids_filename):
    # full path in docker: blast/reciprocal_blast/celery_blast/media/taxonomic_node_files/
    filepath_species_taxids = os.getcwd() + '/' + TAXONOMIC_NODES + taxids_filename
    logger.info('invoking get_spexies_taxids.sh script (e-direct tool) with parameter -t and input node {},'
                'output is redirected into {}'.format(taxonomic_node, filepath_species_taxids))
    # invoke the get_species_taxids.sh script and redirect ouput into file
    try:
        returncode = 1
        # iteration over all possible taxonomic nodes
        with open(filepath_species_taxids, 'w') as taxfile:
            for node in taxonomic_node:
                e_direct_process = Popen(
                    ['get_species_taxids.sh', '-t', str(node)],
                    stdout=taxfile,
                    stderr=subSTDOUT
                )
                # communicate with subprocess : https://docs.python.org/3/library/subprocess.html#subprocess.Popen.communicate
                # wait for process to terminate and set returncode attribute
                try:
                    logger.info(
                        'INFO:waiting for popen instance {} to finish with timeout set to {}'.format(
                            e_direct_process.pid, 200))
                    returncode = e_direct_process.wait(timeout=200)
                    logger.info('returncode : {}'.format(returncode))

                except TimeoutExpired as e:
                    logger.info('ERROR:timeout expired - trying to kill process {}'.format(e_direct_process.pid))
                    e_direct_process.kill()
                    raise Exception(
                        'exception during waiting for popen instance : {} \n\t returncode of popen.wait : {}'.format(e,
                                                                                                                     returncode))
        return returncode
    except SubprocessError as e:
        logger.info('subprocess throwed exception: {}'.format(e))
        raise Exception(
            'exception occured during invocation of:\n\t get_species_taxids_into_file function : {}'.format(e))


# TODO documentation
@shared_task(bind=True)
def execute_reciprocal_blast_project(self, project_id):
    try:
        snakemake_working_dir = BLAST_PROJECT_DIR + str(project_id) + '/'
        snakemake_config_file = BLAST_PROJECT_DIR + str(project_id) + '/snakefile_config'
        snakefile_dir = 'static/snakefiles/reciprocal_blast/Snakefile'

        progress_recorder = ProgressRecorder(self)
        progress_recorder.set_progress(0, 100, "started process")

        try:
            update_blast_project_with_task_result_model(project_id, str(self.request.id))
        except Exception as e:
            logger.warning('couldnt update blastproject with exception : {}'.format(e))
            raise Exception('couldnt update blastproject with exception : {}'.format(e))

        try:
            logger.info('INFO:trying to start snakemake reciprocal BLAST workflow')
            progress_recorder.set_progress(25, 100, 'PROGRESS')

            '''
            #snakemake --snakefile '../../../static/snakefiles/reciprocal_blast/Snakefile' --cores 1 --configfile 'media/blast_project/1/snakefile_config --directory 'media/blast_project/1'
            cmd = 'snakemake --snakefile {} --wms-monitor {} --cores 1 --configfile {} --directory {} --keep-incomplete -q'.format(
                snakefile_dir, settings.PANOPTES_IP, snakemake_config_file, snakemake_working_dir
            )
            reciprocal_blast_snakemake = Popen(
                cmd, shell=True
            )       
            '''

            reciprocal_blast_snakemake = Popen(
                ['snakemake',
                 '--snakefile', snakefile_dir,
                 '--cores', '1',
                 '--configfile', snakemake_config_file,
                 '--directory', snakemake_working_dir,
                 ], shell=False)  # -q

            progress_recorder.set_progress(50, 100, "PROGRESS")

            logger.info('INFO:waiting for popen instance {} to finish with timeout set to {}'.format(
                reciprocal_blast_snakemake.pid, 604800))
            returncode = reciprocal_blast_snakemake.wait(timeout=settings.SUBPROCESS_TIME_LIMIT)  # 66 min 604800 = 7d
            logger.info('returncode : {}'.format(returncode))
            if (returncode != 0):
                logger.warning('subprocess Popen reciprocal BLAST resulted in an error!')
                progress_recorder.set_progress(0, 100, "FAILURE")
                raise Exception('Popen hasnt succeeded ...')

            progress_recorder.set_progress(100, 100, "SUCCESS")
            logger.info('INFO:creating external tools model')
            create_external_tools_after_snakemake_workflow_finishes(project_id, 'local')
            logger.info('INFO:update phylo and msa task with id of the reciprocal BLAST')
            external_tools = ExternalTools.objects.get_external_tools_based_on_project_id(project_id, 'local')
            external_tools.update_for_all_query_sequences_msa_task(str(self.request.id))
            external_tools.update_for_all_query_sequences_phylo_task(str(self.request.id))

            return returncode

        except TimeoutExpired as e:
            logger.info('ERROR:timeout expired - trying to kill process {}'.format(reciprocal_blast_snakemake.pid))
            if 'reciprocal_blast_snakemake' in locals():
                pid = reciprocal_blast_snakemake.pid
                parent = psutil.Process(pid)
                for child in parent.children(recursive=True):
                    child.kill()
                parent.kill()

            raise Exception(
                'ERROR:exception during waiting for popen reciprocal blast task: {} \n\t returncode of popen.wait : {}'.format(
                    e, returncode))

        except SubprocessError as e:
            logger.info('ERROR:subprocess throwed exception: {}'.format(e))
            if 'reciprocal_blast_snakemake' in locals():
                pid = reciprocal_blast_snakemake.pid
                parent = psutil.Process(pid)
                for child in parent.children(recursive=True):
                    child.kill()
                parent.kill()
            raise Exception(
                'ERROR:exception occured during invokation of:\n\t reciprocal blast snakefile : {}'.format(e))

    except SoftTimeLimitExceeded as e:
        logger.info("ERROR:Reciprocal BLAST reached Task Time Limit")
        if 'reciprocal_blast_snakemake' in locals():
            pid = reciprocal_blast_snakemake.pid
            parent = psutil.Process(pid)
            for child in parent.children(recursive=True):
                child.kill()
            parent.kill()
        raise Exception("ERROR Reciprocal BLAST reached Task Time Limit")

@shared_task(bind=True)
def execute_remote_reciprocal_blast_project(self, project_id):
    try:
        snakemake_working_dir = REMOTE_BLAST_PROJECT_DIR + str(project_id) + '/'
        snakemake_config_file = REMOTE_BLAST_PROJECT_DIR + str(project_id) + '/snakefile_config'
        snakefile_dir = 'static/snakefiles/reciprocal_blast/remote_searches/Snakefile'

        progress_recorder = ProgressRecorder(self)
        progress_recorder.set_progress(0, 100, "started process")

        try:
            update_remote_blast_project_with_task_result_model(project_id, str(self.request.id))
        except Exception as e:
            logger.warning('couldnt update blastproject with exception : {}'.format(e))
            raise Exception('couldnt update blastproject with exception : {}'.format(e))

        try:
            logger.info('INFO:trying to start snakemake reciprocal BLAST workflow')
            progress_recorder.set_progress(25, 100, 'PROGRESS')

            reciprocal_blast_snakemake = Popen(
                ['snakemake',
                 '--snakefile', snakefile_dir,
                 '--cores', '1',
                 '--configfile', snakemake_config_file,
                 '--directory', snakemake_working_dir,
                 ], shell=False)  # -q

            progress_recorder.set_progress(50, 100, "PROGRESS")

            logger.info('INFO:waiting for popen instance {} to finish with timeout set to {}'.format(
                reciprocal_blast_snakemake.pid, 604800))
            returncode = reciprocal_blast_snakemake.wait(timeout=settings.SUBPROCESS_TIME_LIMIT)  # 66 min 604800 = 7d
            logger.info('returncode : {}'.format(returncode))
            if (returncode != 0):
                logger.warning('subprocess Popen reciprocal BLAST resulted in an error!')
                progress_recorder.set_progress(0, 100, "FAILURE")
                raise Exception('Popen hasnt succeeded ...')

            progress_recorder.set_progress(100, 100, "SUCCESS")
            logger.info('INFO:creating external tools model')
            create_external_tools_after_snakemake_workflow_finishes(project_id,'remote')
            logger.info('INFO:update phylo and msa task with id of the reciprocal BLAST')
            external_tools = ExternalTools.objects.get_external_tools_based_on_project_id(project_id, 'remote')
            external_tools.update_for_all_query_sequences_msa_task(str(self.request.id))
            external_tools.update_for_all_query_sequences_phylo_task(str(self.request.id))

            return returncode

        except TimeoutExpired as e:
            logger.info('ERROR:timeout expired - trying to kill process {}'.format(reciprocal_blast_snakemake.pid))
            if 'reciprocal_blast_snakemake' in locals():
                pid = reciprocal_blast_snakemake.pid
                parent = psutil.Process(pid)
                for child in parent.children(recursive=True):
                    child.kill()
                parent.kill()

            raise Exception(
                'ERROR:exception during waiting for popen reciprocal blast task: {} \n\t returncode of popen.wait : {}'.format(
                    e, returncode))

        except SubprocessError as e:
            logger.info('ERROR:subprocess throwed exception: {}'.format(e))
            if 'reciprocal_blast_snakemake' in locals():
                pid = reciprocal_blast_snakemake.pid
                parent = psutil.Process(pid)
                for child in parent.children(recursive=True):
                    child.kill()
                parent.kill()
            raise Exception(
                'ERROR:exception occured during invokation of:\n\t reciprocal blast snakefile : {}'.format(e))

    except SoftTimeLimitExceeded as e:
        logger.info("ERROR:Reciprocal BLAST reached Task Time Limit")
        if 'reciprocal_blast_snakemake' in locals():
            pid = reciprocal_blast_snakemake.pid
            parent = psutil.Process(pid)
            for child in parent.children(recursive=True):
                child.kill()
            parent.kill()
        raise Exception("ERROR Reciprocal BLAST reached Task Time Limit")

# TODO documentation
@shared_task(bind=True)
def execute_makeblastdb_with_uploaded_genomes(self, database_id, path_to_database, taxmap_file=None,
                                              taxonomic_node=None):
    try:
        progress_recorder = ProgressRecorder(self)
        progress_recorder.set_progress(0, 100, "started process")
        logger.info("INFO:database_id : {}, path_to_database : {}".format(database_id, path_to_database))
        try:
            update_blast_database_with_task_result_model(database_id, self.request.id)
        except Exception as e:
            logger.warning('ERROR:couldnt update blastdatabase with exception : {}'.format(e))
            raise Exception('ERROR:couldnt update blastdatabase with exception : {}'.format(e))

        try:
            logger.info('INFO:trying to start makeblastdb for formatting blast database : {}'.format(path_to_database))
            progress_recorder.set_progress(25, 100, 'PROGRESS')
            if taxmap_file != None:
                logger.info('INFO:execution with provided taxmap_file')
                makeblastdb_popen = Popen(
                    ['makeblastdb',
                     '-in', path_to_database,
                     '-out', path_to_database,
                     '-taxid_map', BLAST_DATABASE_DIR + str(database_id) + '/acc_taxmap.table',
                     '-dbtype', 'prot',
                     '-input_type', 'fasta',
                     '-parse_seqids'
                     ], shell=False)

            elif taxonomic_node != None:
                logger.info('INFO:execution with provided taxonomic_node')
                # logger.info('INFO:makeblastdb -in {} -out {} -taxid {} -dbtype {} -input_type {} -parse_seqids'.format(
                #    path_to_database, path_to_database, str(taxonomic_node),'prot','fasta'
                # ))
                makeblastdb_popen = Popen(
                    ['makeblastdb',
                     '-in', path_to_database,
                     '-out', path_to_database,
                     '-taxid', str(taxonomic_node),
                     '-dbtype', 'prot',
                     '-input_type', 'fasta',
                     '-parse_seqids'
                     ], shell=False)
            else:
                logger.info("ERROR:there is no taxmap file or taxonomic node ...\n")
                raise SubprocessError(
                    'ERROR:no subprocess to execute : taxmap_file : {} , taxonomic_node : {}...'.format(taxmap_file,
                                                                                                        taxonomic_node))
            progress_recorder.set_progress(90, 100, 'PROGRESS')

            logger.info(
                'INFO:waiting for popen instance {} to finish with timeout set to {}'.format(makeblastdb_popen.pid,
                                                                                             settings.SUBPROCESS_TIME_LIMIT))
            returncode = makeblastdb_popen.wait(timeout=settings.SUBPROCESS_TIME_LIMIT)  # 66 min 604800 = 7d
            if returncode > 0:
                logger.info('ERROR:returncode : {}'.format(returncode))
                progress_recorder.set_progress(0, 100, "FAILURE")
                raise Exception('ERROR:error during makeblastdb attempt ...')
            else:
                logger.info('INFO:returncode : {}'.format(returncode))
                progress_recorder.set_progress(100, 100, "SUCCESS")
            return returncode

        except TimeoutExpired as e:
            logger.info('ERROR:timeout expired ... trying to kill process {}'.format(makeblastdb_popen.pid))
            makeblastdb_popen.kill()
            raise Exception(
                'exception during waiting for popen instance : {} \n\t returncode of popen.wait > 0'.format(e))

        except SubprocessError as e:
            logger.warning('ERROR:error during execution of makeblastdb with exception : {}'.format(e))
            makeblastdb_popen.kill()
            raise Exception(
                'ERROR:exception in Subprocess Popen call : {}'.format(e))

    except SoftTimeLimitExceeded as e:
        logger.info("ERROR:makeblastdb process reached Task Time Limit")
        makeblastdb_popen.kill()
        raise Exception("ERROR:makeblastdb process reached Task Time Limit")


'''calculate_database_statistics_task
    
    This function calculates the database statistics based on the reciprocal_results_with_taxonomy.csv
    and the database csv files. It produces several output files. 
    The {database_name}_with_taxonomic_information.csv file for the underlying forward database, 
    {taxonomic_unit}_database_statistics_normalized.csv, {taxonomic_unit}_database_statistics.csv, 
    
    :param self
        :type TaskResult object
    :param project_id
        :type int
    :param taxonomic_units
        :type list[str]
    
'''


@shared_task(bind=True)
def calculate_database_statistics_task(self, project_id: int, user_email: str, taxonomic_units: list):
    try:
        progress_recorder = ProgressRecorder(self)
        progress_recorder.set_progress(0, 100, "STARTED")
        logger.info("INFO:started database statistics for project: {}".format(project_id))
        try:
            update_blast_project_with_database_statistics_task_result_model(project_id, self.request.id)
            progress_recorder.set_progress(25, 100, "PROGRESS")


        except Exception as e:
            logger.warning('ERROR:couldnt update blast project with exception : {}'.format(e))
            raise Exception('ERROR:couldnt update blast project with exception : {}'.format(e))

        logfile = BLAST_PROJECT_DIR + str(project_id) + '/log/calculate_database_statistics.log'
        if os.path.isdir(BLAST_PROJECT_DIR + str(project_id) + '/log'):
            logger.info("INFO:Starting database statistics task")
            calculate_database_statistics(project_id, logfile=logfile, user_email=user_email,
                                          taxonomic_units=taxonomic_units)

        else:
            logger.warning("WARNING: cant write {}".format(logfile))

        progress_recorder.set_progress(100, 100, "SUCCESS")
        logger.info("DONE:calculating database statistics for project: {}".format(project_id))

    except SoftTimeLimitExceeded as e:
        logger.info("ERROR:database statistics calculation reached Task Time Limit")
        raise Exception("ERROR:database statistics calculation reached Task Time Limit")
    except Exception as e:
        raise Exception("ERROR: unknown exception occurred: {}".format(e))


'''download_and_decompress_cdd_database OBSOLETE

    This task can be used to refresh the CDD database which is per default loaded before
    web server startup.

'''


@shared_task(bind=True)
def download_and_decompress_cdd_database(self):
    logger.info("INFO:NO CDD DATBASE\n"
                "INFO:trying to download CDD database from: {}".format(CDD_DATABASE_URL))

    try:
        cdd_ftp_path = CDD_DATABASE_URL
        if os.path.isdir(BLAST_DATABASE_DIR + "CDD/") == False:
            os.mkdir(BLAST_DATABASE_DIR + "CDD/")
        path_to_cdd_location = BLAST_DATABASE_DIR
        path_to_cdd_location = path_to_cdd_location + 'Cdd_LE.tar.gz'
        logger.info("INFO:DOWNLOADING CONSVERED DOMAIN DATABASE INTO {}".format(BLAST_DATABASE_DIR + "CDD/"))
        proc = Popen(["wget", cdd_ftp_path, "-q", "-O", path_to_cdd_location], shell=False)
        returncode = proc.wait(timeout=3000)
        if returncode != 0:
            raise SubprocessError
        logger.info("INFO:EXTRACTING CONSERVED DOMAIN DATABASE")

        proc = Popen(
            ["tar", "-zxvf", path_to_cdd_location, "-C", "/blast/reciprocal_blast/media/databases/CDD/"],
            shell=False)
        returncode = proc.wait(timeout=800)
        if returncode != 0:
            raise SubprocessError
        os.remove(path_to_cdd_location)  # remove the zip of cdd database
        logger.info("INFO:DONE DOWNLOADING CONSERVED DOMAIN DATABASE")
    except TimeoutExpired as e:
        logger.warning("ERROR:TIMEOUT EXPIRED DURING DOWNLOAD OF CONSERVED DOMAIN DATABASE: {}".format(e))
        logger.info(
            "INFO:IF YOU HAVE NO STABLE INTERNET CONNECTION TRY TO RESTART THE WEBSERVER ONCE YOU HAVE A BETTER CONNECTION")
        logger.info("INFO:YOU CAN MANUALLY LOAD THE DOMAIN DATABASE DATABASE INTO THE DATABASE/CDD FOLDER")

        if 'proc' in locals():
            pid = proc.pid
            parent = psutil.Process(pid)

            for child in parent.children(recursive=True):
                child.kill()
            parent.kill()
        else:
            logger.warning("WARNING:CHECK FOR UNFINISHED PROCESSES OR RESTART THE WEB-SERVER")
        raise TimeoutExpired(e)
    except SubprocessError as e:
        logger.warning("ERROR:WGET RESULTED IN AN ERROR: {}".format(e))
        logger.info("INFO:YOU CAN MANUALLY LOAD THE CONSERVED DOMAIN DATABASE INTO THE DATABASE/CDD FOLDER")

        if 'proc' in locals():
            pid = proc.pid
            parent = psutil.Process(pid)

            for child in parent.children(recursive=True):
                child.kill()
            parent.kill()
        else:
            logger.warning("WARNING:CHECK FOR UNFINISHED PROCESSES OR RESTART THE WEB-SERVER")
        raise SubprocessError(e)
    except Exception as e:
        logger.warning("ERROR:UNSUSPECTED ERROR OCCURRED: {}".format(e))
        logger.warning("WARNING:YOU CANT USE THE CDD EXTENSION")
        logger.warning("INFO:YOU CAN MANUALLY LOAD THE CONSERVED DOMAIN DATABASE INTO THE DATABASE/CDD FOLDER")
        raise Exception(e)


