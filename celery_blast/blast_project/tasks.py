#The best practice is to create a common logger for all of your tasks at the top of your module:

import os
from subprocess import Popen, PIPE as subPIPE, STDOUT as subSTDOUT, SubprocessError, TimeoutExpired
from django.conf import settings
from celery import shared_task
from celery.utils.log import get_task_logger
from celery_progress.backend import ProgressRecorder

from .py_django_db_services import update_blast_project_with_task_result_model, update_blast_database_with_task_result_model, create_external_tools_after_snakemake_workflow_finishes
#logger for celery worker instances
logger = get_task_logger(__name__)

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
    #full path in docker: blast/reciprocal_blast/celery_blast/media/taxonomic_node_files/
    filepath_species_taxids = os.getcwd() + '/media/taxonomic_node_files/' + taxids_filename
    logger.info('invoking get_spexies_taxids.sh script (e-direct tool) with parameter -t and input node {},'
                'output is redirected into {}'.format(taxonomic_node, filepath_species_taxids))
    #invoke the get_species_taxids.sh script and redirect ouput into file
    try:
        e_direct_process = Popen(
            ['get_species_taxids.sh','-t', str(taxonomic_node)],
            stdout=open(filepath_species_taxids,'w'),
            stderr=subSTDOUT
        )

    except SubprocessError as e:
        logger.info('subprocess throwed exception: {}'.format(e))
        raise Exception('exception occured during invokation of:\n\t get_species_taxids_into_file function : {}'.format(e))

    #communicate with subprocess : https://docs.python.org/3/library/subprocess.html#subprocess.Popen.communicate
    #wait for process to terminate and set returncode attribute
    try:
        logger.info('waiting for popen instance {} to finish with timeout set to {}'.format(e_direct_process.pid,200))
        returncode = e_direct_process.wait(timeout=200)
        logger.info('returncode : {}'.format(returncode))
        return returncode

    except TimeoutExpired as e:
        logger.info('timeout expired ... trying to kill process {}'.format(e_direct_process.pid))
        e_direct_process.kill()

        raise Exception('exception during waiting for popen instance : {} \n\t returncode of popen.wait : {}'.format(e,returncode))

#TODO documentation
@shared_task(bind=True)
def execute_reciprocal_blast_project(self,project_id):
    snakemake_working_dir = 'media/blast_projects/' + str(project_id) + '/'
    snakemake_config_file = 'media/blast_projects/' + str(project_id) + '/snakefile_config'
    snakefile_dir = 'static/snakefiles/reciprocal_blast/Snakefile'

    progress_recorder = ProgressRecorder(self)
    progress_recorder.set_progress(0, 100,"started process")

    try:
        update_blast_project_with_task_result_model(project_id, str(self.request.id))
    except Exception as e:
        logger.warning('couldnt update blastproject with exception : {}'.format(e))
        raise Exception('couldnt update blastproject with exception : {}'.format(e))
    try:
        logger.info('trying to start snakemake reciprocal BLAST workflow')
        progress_recorder.set_progress(25,100,'PROGRESS')
        #snakemake --snakefile '../../../static/snakefiles/reciprocal_blast/Snakefile' --cores 1 --configfile 'media/blast_project/13/snakefile_config --directory 'media/blast_project/13'
        reciprocal_blast_snakemake = Popen(
            ['snakemake',
             '--snakefile',snakefile_dir,
             '--wms-monitor',settings.PANOPTES_IP,
             '--cores','1',
             '--configfile',snakemake_config_file,
             '--directory',snakemake_working_dir,
             '--keep-incomplete'], shell=False, stdout=subPIPE, stderr=subSTDOUT)
        progress_recorder.set_progress(50, 100, "PROGRESS")

    except SubprocessError as e:
        logger.info('subprocess throwed exception: {}'.format(e))
        raise Exception(
            'exception occured during invokation of:\n\t reciprocal blast snakefile : {}'.format(e))

    try:
        logger.info('waiting for popen instance {} to finish with timeout set to {}'.format(reciprocal_blast_snakemake.pid, 604800))
        returncode = reciprocal_blast_snakemake.wait(timeout=604800) #66 min 604800 = 7d
        logger.info('returncode : {}'.format(returncode))
        if (returncode != 0):
            logger.warning('subprocess Popen reciprocal BLAST resulted in an error!')
            progress_recorder.set_progress(0, 100, "FAILURE")
            raise Exception('Popen hasnt succeeded ...')

        progress_recorder.set_progress(100, 100, "SUCCESS")
        create_external_tools_after_snakemake_workflow_finishes(project_id)

        return returncode

    except TimeoutExpired as e:
        logger.info('timeout expired ... trying to kill process {}'.format(reciprocal_blast_snakemake.pid))
        reciprocal_blast_snakemake.kill()

        raise Exception(
            'exception during waiting for popen instance : {} \n\t returncode of popen.wait : {}'.format(e, returncode))

@shared_task(bind=True)
def execute_makeblastdb_with_uploaded_genomes(self,database_id,path_to_database,taxmap_file=None, taxonomic_node=None):
    progress_recorder = ProgressRecorder(self)
    progress_recorder.set_progress(0, 100,"started process")
    logger.info("database_id : {}, path_to_database : {}".format(database_id,path_to_database))
    try:
        update_blast_database_with_task_result_model(database_id,self.request.id)
    except Exception as e:
        logger.warning('couldnt update blastdatabase with exception : {}'.format(e))
        raise Exception('couldnt update blastdatabase with exception : {}'.format(e))

    try:
        logger.info('trying to start makeblastdb for formatting blast database : {}'.format(path_to_database))
        progress_recorder.set_progress(25,100,'PROGRESS')
        if taxmap_file != None:
            logger.info('execution with provided taxmap_file')
            makeblastdb_popen = Popen(
                ['makeblastdb',
                 '-in',path_to_database,
                 '-out',path_to_database,
                 '-taxid_map','media/databases/'+str(database_id)+'/acc_taxmap.table',
                 '-dbtype','prot',
                 '-input_type','fasta',
                 '-parse_seqids'
                ], shell=False, stdout=subPIPE, stderr=subSTDOUT)
        elif taxonomic_node != None:
            logger.info('execution with provided taxonomic_node')
            logger.info('makeblastdb -in {} -out {} -taxid {} -dbtype {} -input_type {} -parse_seqids'.format(
                path_to_database, path_to_database, str(taxonomic_node),'prot','fasta'
            ))
            makeblastdb_popen = Popen(
                ['makeblastdb',
                 '-in',path_to_database,
                 '-out',path_to_database,
                 '-taxid',str(taxonomic_node),
                 '-dbtype','prot',
                 '-input_type','fasta',
                 '-parse_seqids'
                ], shell=False, stdout=subPIPE, stderr=subSTDOUT)
        else:
            raise SubprocessError('no subprocess to execute : taxmap_file : {} , taxonomic_node : {}...'.format(taxmap_file,taxonomic_node))
        progress_recorder.set_progress(90,100,'PROGRESS')

    except SubprocessError as e:
        logger.warning('error during execution of makeblastdb with exception : {}'.format(e))

    try:
        logger.info(
            'waiting for popen instance {} to finish with timeout set to {}'.format(makeblastdb_popen.pid,
                                                                                    4000))
        returncode = makeblastdb_popen.wait(timeout=604800)  # 66 min 604800 = 7d
        if returncode > 0:
            logger.info('ERROR: returncode : {}'.format(returncode))
            progress_recorder.set_progress(0, 100, "FAILURE")
            raise Exception('error during makeblastdb attempt ...')
        else:
            logger.info('returncode : {}'.format(returncode))
            progress_recorder.set_progress(100, 100, "SUCCESS")
        return returncode

    except TimeoutExpired as e:
        logger.info('timeout expired ... trying to kill process {}'.format(makeblastdb_popen.pid))
        makeblastdb_popen.kill()

        raise Exception('exception during waiting for popen instance : {} \n\t returncode of popen.wait > 0'.format(e))