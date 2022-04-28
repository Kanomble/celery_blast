import json


from celery import shared_task
from celery.utils.log import get_task_logger
import requests
from .models import ExternalTools
from celery_progress.backend import ProgressRecorder
from subprocess import Popen, PIPE as subPIPE, STDOUT as subSTDOUT, SubprocessError, TimeoutExpired
from .py_services import check_if_target_sequences_are_available, check_if_msa_file_is_available, create_html_output_for_newicktree
from .entrez_search_service import execute_entrez_search, create_random_filename, save_entrez_search_model

#logger for celery worker instances
logger = get_task_logger(__name__)

@shared_task(bind=True)
def entrez_search_task(self,database:str,entrez_query:str,user_id:int):
    try:
        logger.info("trying to start entrez search")
        progress_recorder = ProgressRecorder(self)
        progress_recorder.set_progress(0, 100, "PROGRESS")

        randomly_generated_filename = create_random_filename() + '.tsf'
        esearch_output_filepath = 'media/esearch_output/result_dataframe_' + randomly_generated_filename
        save_entrez_search_model(database=database,
                                 entrez_query=entrez_query,
                                 file_name=esearch_output_filepath,
                                 task_result_id=self.request.id,
                                 user_id=user_id)
        returncode = execute_entrez_search(database, entrez_query, esearch_output_filepath)

        if returncode != 0:
            raise Exception("Popen hasnt succeeded, returncode != 0: {}".format(returncode))
        else:
            progress_recorder.set_progress(100, 100, "SUCCESS")
            return 0
    except Exception as e:
        raise Exception("[-] Couldnt perform entrez search with exception: {}".format(e))


'''
execute_multiple_sequence_alignment 
'''
@shared_task(bind=True)
def execute_multiple_sequence_alignment(self, project_id, query_sequence_id):
    try:
        logger.info("trying to start mafft multiple sequence alignment")
        progress_recorder = ProgressRecorder(self)
        progress_recorder.set_progress(0, 100, "PROGRESS")
        external_tools = ExternalTools.objects.get_external_tools_based_on_project_id(project_id)
        external_tools.update_query_sequences_msa_task(query_sequence_id, str(self.request.id))
        logger.info("updated query sequence model with taskresult instance : {}".format(str(self.request.id)))
        path_to_project = 'media/blast_projects/' + str(project_id) + '/'
        path_to_query_file = path_to_project + query_sequence_id + '/target_sequences.faa'

        target_sequence_status = check_if_target_sequences_are_available(path_to_query_file)
        if target_sequence_status == 0:
            #mafft invocation with default settings
            path_to_mafft_output = path_to_project + query_sequence_id + '/target_sequences.msa'
            cmd = "mafft {} > {}".format(path_to_query_file, path_to_mafft_output)
            msa_task = Popen(cmd, shell=True)
            logger.info(
                'waiting for popen instance {} to finish with timeout set to {}'.format(msa_task.pid,
                                                                                        40000))
            progress_recorder.set_progress(20, 100, "PROGRESS")

            returncode = msa_task.wait(40000)
            if returncode != 0:
                raise Exception("Popen hasnt succeeded, returncode != 0: {}".format(returncode))
            else:
                progress_recorder.set_progress(100, 100, "SUCCESS")
                return 0
        elif target_sequence_status == 1:
            raise FileNotFoundError("query file with targets for msa does not exist!")
        elif target_sequence_status == 2:
            raise Exception("not enough target sequences")

    except Exception as e:
        raise Exception("[-] Couldnt perform multiple sequence alignment task with Exception: {}".format(e))

@shared_task(bind=True)
def execute_phylogenetic_tree_building(self,project_id,query_sequence_id):
    try:
        logger.info("trying to execute fasttree phylogenetic tree construction per request to bioinformatic tools container")

        progress_recorder = ProgressRecorder(self)
        progress_recorder.set_progress(0,100,"PROGRESS")

        external_tools = ExternalTools.objects.get_external_tools_based_on_project_id(project_id)
        logger.info("cheking if msa task succeeded for query sequence : {}".format(query_sequence_id))

        path_to_project = 'media/blast_projects/' + str(project_id) + '/'
        path_to_msa_file = path_to_project + query_sequence_id + '/target_sequences.msa'
        path_to_fasttree_output = path_to_project + query_sequence_id + '/target_sequences.nwk'
        msa_status = check_if_msa_file_is_available(path_to_msa_file)
        if msa_status == 0 and external_tools.check_if_msa_task_is_completed(query_sequence_id):

            external_tools.update_query_sequences_phylo_task(query_sequence_id,str(self.request.id))
            logger.info("updated query sequence model with taskresult instance : {}".format(str(self.request.id)))
            progress_recorder.set_progress(20, 100, "PROGRESS")
            cmd = "fasttree -lg {} > {}".format(path_to_msa_file, path_to_fasttree_output)
            phylo_task = Popen(cmd,shell=True)
            progress_recorder.set_progress(30, 100, "PROGRESS")
            logger.info(
                'waiting for popen instance {} to finish with timeout set to {}'.format(phylo_task.pid,
                                                                                        40000))
            returncode = phylo_task.wait(4000)
            if returncode != 0:
                raise Exception("Popen hasnt succeeded, returncode != 0: {}".format(returncode))

            returncode = create_html_output_for_newicktree(path_to_fasttree_output,project_id, query_sequence_id)
            if returncode != 0:
                raise Exception("HTML building hasnt succeeded, returncode != 0: {}".format(returncode))

            else:
                progress_recorder.set_progress(100, 100, "SUCCESS")
                return 0
        elif msa_status == 1:
            raise FileNotFoundError("msa file does not exist!")
    except Exception as e:
        raise Exception("[-] Couldnt perform phylogenetic tree task with Exception: {}".format(e))

@shared_task(bind=True)
def execute_multiple_sequence_alignment_for_all_query_sequences(self, project_id):
    try:
        logger.info("trying to start mafft multiple sequence alignment for multiple query sequence targets")
        progress_recorder = ProgressRecorder(self)
        progress_recorder.set_progress(0, 100, "PROGRESS")
        external_tools = ExternalTools.objects.get_external_tools_based_on_project_id(project_id)
        external_tools.update_for_all_query_sequences_msa_task(str(self.request.id))
        logger.info("updated multiple query sequence models with taskresult instance : {}".format(str(self.request.id)))
        progress_recorder.set_progress(19, 100, "PROGRESS")

        query_sequence_ids = [qseq.query_accession_id for qseq in external_tools.query_sequences.get_queryset()]

        path_to_project = 'media/blast_projects/' + str(project_id) + '/'

        progress = 80 / len(query_sequence_ids)
        counter = 1
        for qseqid in query_sequence_ids:
            path_to_query_file = path_to_project + qseqid + '/target_sequences.faa'
            path_to_mafft_output = path_to_project + qseqid + '/target_sequences.msa'

            target_sequence_status = check_if_target_sequences_are_available(path_to_query_file)

            if target_sequence_status == 0:
                cmd = "mafft {} > {}".format(path_to_query_file, path_to_mafft_output)
                msa_task = Popen(cmd, shell=True)
                logger.info(
                    'waiting for popen instance {} to finish with timeout set to {}'.format(msa_task.pid,
                                                                                            40000))
                returncode = msa_task.wait(40000)
                progress = int(progress * counter)
                if progress <= 80:
                    logger.info(
                        'progress of multiple sequence alignment task for all query sequences set to {}'.format(progress))
                    progress_recorder.set_progress(progress, 100, "PROGRESS")
                    counter += 1

                if returncode != 0:
                    raise Exception("Popen hasnt succeeded, returncode != 0: {}".format(returncode))
            elif target_sequence_status == 1:
                raise FileNotFoundError("query file with targets for msa does not exist!")
            elif target_sequence_status == 2:
                continue

        progress_recorder.set_progress(100, 100, "PROGRESS")
    except Exception as e:
        raise Exception("[-] Couldnt perform multiple sequence alignment task for multiple query sequences with Exception: {}".format(e))

@shared_task(bind=True)
def execute_fasttree_phylobuild_for_all_query_sequences(self, project_id):
    try:
        logger.info("trying to execute fasttree phylogenetic tree construction for all query sequences")
        progress_recorder = ProgressRecorder(self)
        progress_recorder.set_progress(0,100,"PROGRESS")
        external_tools = ExternalTools.objects.get_external_tools_based_on_project_id(project_id)

        external_tools.update_for_all_query_sequences_phylo_task(str(self.request.id))
        logger.info("updated multiple query sequence models with taskresult instance : {}".format(str(self.request.id)))
        path_to_project = 'media/blast_projects/' + str(project_id) + '/'

        query_sequence_ids = []
        for qseq in external_tools.query_sequences.get_queryset():
            logger.info("check if msa task succeeded for query sequence : {}".format(qseq.query_accession_id))
            if external_tools.check_if_msa_task_is_completed(qseq.query_accession_id):
                query_sequence_ids.append(qseq.query_accession_id)
                logger.info("\tmsa task succeeded ... added query sequence to target list for phylogenetic tree construction")
        progress = 80/len(query_sequence_ids)
        counter = 1
        for qseqid in query_sequence_ids:
            path_to_msa_file = path_to_project + qseqid + '/target_sequences.msa'
            path_to_fasttree_output = path_to_project + qseqid + '/target_sequences.nwk'
            msa_status = check_if_msa_file_is_available(path_to_msa_file)
            if msa_status == 0 and external_tools.check_if_msa_task_is_completed(qseqid):

                external_tools.update_query_sequences_phylo_task(qseqid, str(self.request.id))
                logger.info("updated query sequence model with taskresult instance : {}".format(str(self.request.id)))

                cmd = "fasttree -lg {} > {}".format(path_to_msa_file, path_to_fasttree_output)
                phylo_task = Popen(cmd, shell=True)
                logger.info(
                    'waiting for popen instance {} to finish with timeout set to {}'.format(phylo_task.pid,
                                                                                            40000))
                returncode = phylo_task.wait(4000)

                progress = int(progress * counter)
                if progress <= 80:
                    logger.info(
                        'progress of fasttree task for all query sequences set to {}'.format(
                            progress))
                    progress_recorder.set_progress(progress, 100, "PROGRESS")
                    counter += 1

                if returncode != 0:
                    raise Exception("Popen hasnt succeeded, returncode != 0: {}".format(returncode))
            elif msa_status == 1:
                raise FileNotFoundError("msa file does not exist!")
    except Exception as e:
        raise Exception("[-] couldnt perform fasttree task for all query sequences with exception: {}".format(e))
