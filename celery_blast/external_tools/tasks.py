import json

from celery import shared_task
from celery.utils.log import get_task_logger
import requests
from .models import ExternalTools
from celery_progress.backend import ProgressRecorder
from subprocess import Popen, PIPE as subPIPE, STDOUT as subSTDOUT, SubprocessError, TimeoutExpired
from .py_services import check_if_target_sequences_are_available, check_if_msa_file_is_available

#logger for celery worker instances
logger = get_task_logger(__name__)

@shared_task(bind=True)
def execute_multiple_sequence_alignment_new(self, project_id, query_sequence_id):
    try:
        logger.info("trying to start mafft multiple sequence alignment")
        progress_recorder = ProgressRecorder(self)
        progress_recorder.set_progress(0, 100, "PROGRESS")


        external_tools = ExternalTools.objects.get_external_tools_based_on_project_id(project_id)
        external_tools.update_query_sequences_msa_task(query_sequence_id, str(self.request.id))
        logger.info("updated query sequence model with taskresult instance : {}".format(str(self.request.id)))

        path_to_project = 'data/blast_projects/' + str(project_id) + '/'
        path_to_query_file = path_to_project + query_sequence_id + '/target_sequences.faa'

        target_sequence_status = check_if_target_sequences_are_available(path_to_query_file)
        if target_sequence_status == 0:
            #mafft invocation with default settings
            output = path_to_project + query_sequence_id + '/target_sequences.msa'
            msa_task = Popen(
                ['mafft',path_to_query_file,'>',output], shell=False, stdout=subPIPE, stderr=subSTDOUT)
            logger.info(
                'waiting for popen instance {} to finish with timeout set to {}'.format(msa_task.pid,
                                                                                        4000))
            returncode = msa_task.wait(4000)
            if returncode != 0:
                raise Exception("Popen hasnt succeeded, returncode != 0: {}".format(returncode))
            else:
                return 0
        elif target_sequence_status == 1:
            raise FileNotFoundError("query file with targets for msa does not exist!")
        elif target_sequence_status == 2:
            raise Exception("not enough target sequences")

    except Exception as e:
        raise Exception("[-] Couldnt perform multiple sequence alignment task with Exception: {}".format(e))

@shared_task(bind=True)
def execute_phylogenetic_tree_building_new(self,project_id,query_sequence_id):
    try:
        path_to_msa_file=""
        check_if_msa_file_is_available(path_to_msa_file)
    except Exception as e:
        raise Exception("[-] Couldnt perform phylogenetic tree task with Exception: {}".format(e))
#TODO documentation - refactoring with logger statements
@shared_task(bind=True)
def execute_multiple_sequence_alignment(self, project_id, query_sequence_id):
    try:
        logger.info("trying to execute mafft multiple sequence alignment per request to bioinformatic tools container")
        progress_recorder = ProgressRecorder(self)
        progress_recorder.set_progress(0, 100, "PROGRESS")


        external_tools = ExternalTools.objects.get_external_tools_based_on_project_id(project_id)
        external_tools.update_query_sequences_msa_task(query_sequence_id, str(self.request.id))
        logger.info("updated query sequence model with taskresult instance : {}".format(str(self.request.id)))
        url = "http://tools:5001/perform_simple_msa/" + str(project_id) +'/' + query_sequence_id
        response = requests.post(url, json={"query_sequence_id":query_sequence_id, "project_id":project_id})
        logger.info("received response from : {} with status_code : {}".format(url,response.status_code))
        if response.status_code == 200:
            progress_recorder.set_progress(100, 100, "SUCCESS")
            return 0
        elif response.status_code == 400:
            raise Exception("[-] status code for http request to bioinformatic tools container with url : {} is 400".format(url))
        elif response.status_code == 500:
            raise Exception(
                "[-] status code for http request to bioinformatic tools container with url : {} is 500".format(url))
        else:
            raise Exception("[-] error in response from bioinformatic tools container with url : {} status code is : {}"
                            .format(url,response.status_code))
    except Exception as e:
        raise Exception("[-] couldnt perform simple msa task with exception : {}".format(e))

@shared_task(bind=True)
def execute_phylogenetic_tree_building(self,project_id,query_sequence_id):
    try:
        logger.info("trying to execute fasttree phylogenetic tree construction per request to bioinformatic tools container")

        progress_recorder = ProgressRecorder(self)
        progress_recorder.set_progress(0,100,"PROGRESS")

        external_tools = ExternalTools.objects.get_external_tools_based_on_project_id(project_id)
        logger.info("cheking if msa task succeeded for query sequence : {}".format(query_sequence_id))
        if external_tools.check_if_msa_task_is_completed(query_sequence_id):
            external_tools.update_query_sequences_phylo_task(query_sequence_id,str(self.request.id))
            logger.info("updated query sequence model with taskresult instance : {}".format(str(self.request.id)))

            url = "http://tools:5001/perform_phylo_task/" + str(project_id) + '/' + query_sequence_id
            response = requests.post(url, json={"query_sequence_id": query_sequence_id, "project_id": project_id})
            logger.info("received response from : {} with status_code : {}".format(url, response.status_code))

            if response.status_code == 200:
                progress_recorder.set_progress(100, 100, "SUCCESS")
                return 0
            elif response.status_code == 400:
                raise Exception(
                    "[-] status code for http request to bioinformatic tools container with url : {} is 400".format(url))
            elif response.status_code == 500:
                raise Exception(
                    "[-] status code for http request to bioinformatic tools container with url : {} is 500".format(url))
            else:
                raise Exception("[-] error in response from bioinformatic tools container with url : {} status code is : {}"
                                .format(url, response.status_code))
        else:
            raise Exception("[-] error in performing fasttree task, there is no multiple sequence alignment file for your query sequence : {}".format(query_sequence_id))

    except Exception as e:
        raise Exception("[-] couldnt perform phylogenetic tree building with exception : {}".format(e))

@shared_task(bind=True)
def execute_multiple_sequence_alignment_for_all_query_sequences(self, project_id):
    try:
        logger.info("trying to execute mafft multiple sequence alignment for all query sequences per request to bioinformatic tools container")

        progress_recorder = ProgressRecorder(self)
        progress_recorder.set_progress(0, 100, "PROGRESS")

        external_tools = ExternalTools.objects.get_external_tools_based_on_project_id(project_id)
        external_tools.update_for_all_query_sequences_msa_task(str(self.request.id))
        logger.info("updated multiple query sequence models with taskresult instance : {}".format(str(self.request.id)))

        query_sequence_ids = [qseq.query_accession_id for qseq in external_tools.query_sequences.get_queryset()]
        json_query_sequence_ids = json.dumps(query_sequence_ids)

        url = "http://tools:5001/perform_simple_msa_with_all_qseqs/" + str(project_id)
        response = requests.post(url,json={'query_sequence_ids':json_query_sequence_ids,'project_id':project_id})
        logger.info("received response from : {} with status_code : {}".format(url,response.status_code))

        if response.status_code == 200:
            progress_recorder.set_progress(100, 100, "SUCCESS")
            return 0
        elif response.status_code == 400:
            raise Exception(
                "[-] status code for http request to bioinformatic tools container with url : {} is 400".format(url))
        elif response.status_code == 500:
            raise Exception(
                "[-] status code for http request to bioinformatic tools container with url : {} is 500".format(url))
        else:
            raise Exception("[-] error in response from bioinformatic tools container with url : {} status code is : {}"
                            .format(url, response.status_code))
    except Exception as e:
        raise Exception("[-] couldnt perform msa for all query sequences with exception : {}".format(e))

#TODO documentation - fasttree call just for sequences with SUCCESS of MSA task ...
@shared_task(bind=True)
def execute_fasttree_phylobuild_for_all_query_sequences(self, project_id):
    try:
        logger.info("trying to execute fasttree phylogenetic tree construction for all query sequences per request to bioinformatic tools container")

        progress_recorder = ProgressRecorder(self)
        progress_recorder.set_progress(0, 100, "PROGRESS")

        external_tools = ExternalTools.objects.get_external_tools_based_on_project_id(project_id)

        external_tools.update_for_all_query_sequences_phylo_task(str(self.request.id))
        logger.info("updated multiple query sequence models with taskresult instance : {}".format(str(self.request.id)))

        query_sequence_ids = []
        for qseq in external_tools.query_sequences.get_queryset():
            logger.info("check if msa task succeeded for query sequence : {}".format(qseq.query_accession_id))
            if external_tools.check_if_msa_task_is_completed(qseq.query_accession_id):
                query_sequence_ids.append(qseq.query_accession_id)
                logger.info("\tmsa task succeeded ... added query sequence to target list for phylogenetic tree construction")

        json_query_sequence_ids = json.dumps(query_sequence_ids)

        url = "http://tools:5001/perform_fasttree_phylobuild_with_all_qseqs/" + str(project_id)
        response = requests.post(url,json={'query_sequence_ids':json_query_sequence_ids,'project_id':project_id})
        logger.info("received response from : {} with status_code : {}".format(url,response.status_code))

        if response.status_code == 200:
            progress_recorder.set_progress(100, 100, "SUCCESS")
            return 0
        elif response.status_code == 400:
            raise Exception(
                "[-] status code for http request to bioinformatic tools container with url : {} is 400".format(url))
        elif response.status_code == 500:
            raise Exception(
                "[-] status code for http request to bioinformatic tools container with url : {} is 500".format(url))
        else:
            raise Exception("[-] error in response from bioinformatic tools container with url : {} status code is : {}"
                            .format(url, response.status_code))
    except Exception as e:
        raise Exception("[-] couldnt perform phylogenetic tree construction for all query sequences with exception : {}".format(e))

