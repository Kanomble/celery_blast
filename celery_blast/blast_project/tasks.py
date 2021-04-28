#The best practice is to create a common logger for all of your tasks at the top of your module:

import os
from subprocess import Popen, PIPE as subPIPE, STDOUT as subSTDOUT, SubprocessError, TimeoutExpired

from celery import shared_task
from celery.utils.log import get_task_logger
from celery_progress.backend import ProgressRecorder

from .py_django_db_services import update_blast_project_with_task_result_model
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
        logger.info('trying to start snakemake recipcoal BLAST workflow')
        #snakemake --snakefile 'static/snakefiles/reciprocal_blast/Snakefile' --cores 1 --configfile 'media/blast_projects/4/snakefile_config' --directory 'media/blast_projects/4/'
        reciprocal_blast_snakemake = Popen(
            ['snakemake',
             '--snakefile',snakefile_dir,
             '--wms-monitor','http://172.23.0.5:5000',
             '--cores','1','--configfile',
             snakemake_config_file,'--directory',
             snakemake_working_dir], shell=False, stdout=subPIPE, stderr=subSTDOUT)
        progress_recorder.set_progress(50, 100, "executed snakemake")

    except SubprocessError as e:
        logger.info('subprocess throwed exception: {}'.format(e))
        raise Exception(
            'exception occured during invokation of:\n\t reciprocal blast snakefile : {}'.format(e))

    try:
        logger.info('waiting for popen instance {} to finish with timeout set to {}'.format(reciprocal_blast_snakemake.pid, 4000))
        returncode = reciprocal_blast_snakemake.wait(timeout=4000) #66 min
        logger.info('returncode : {}'.format(returncode))
        progress_recorder.set_progress(100, 100, "SUCCESS")
        return returncode

    except TimeoutExpired as e:
        logger.info('timeout expired ... trying to kill process {}'.format(reciprocal_blast_snakemake.pid))
        reciprocal_blast_snakemake.kill()

        raise Exception(
            'exception during waiting for popen instance : {} \n\t returncode of popen.wait : {}'.format(e, returncode))