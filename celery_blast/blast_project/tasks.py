#The best practice is to create a common logger for all of your tasks at the top of your module:

import os
import subprocess

from celery import shared_task
from celery.utils.log import get_task_logger

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
        e_direct_process = subprocess.Popen(
            ['get_species_taxids.sh','-t', str(taxonomic_node)],
            stdout=open(filepath_species_taxids,'w'),
            stderr=subprocess.STDOUT
        )

    except subprocess.SubprocessError as e:
        logger.info('subprocess throwed exception: {}'.format(e))
        raise Exception('exception occured during invokation of:\n\t get_species_taxids_into_file function : {}'.format(e))

    #communicate with subprocess : https://docs.python.org/3/library/subprocess.html#subprocess.Popen.communicate
    #wait for process to terminate and set returncode attribute
    try:
        logger.info('waiting for popen instance {} to finish with timeout set to {}'.format(e_direct_process.pid,200))
        returncode = e_direct_process.wait(timeout=200)
        logger.info('returncode : {}'.format(returncode))
        return returncode

    except subprocess.TimeoutExpired as e:
        logger.info('timeout expired ... trying to kill process {}'.format(e_direct_process.pid))
        e_direct_process.kill()

        raise Exception('exception during waiting for popen instance : {} \n\t returncode of popen.wait : {}'.format(e,returncode))
