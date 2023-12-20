import os
from subprocess import Popen, SubprocessError, TimeoutExpired, check_output

import pandas as pd

import psutil
from blast_project.py_django_db_services import update_external_tool_with_cdd_search, get_project_by_id, \
    update_blast_project_with_database_statistics_selection_task_result_model, update_domain_database_task_result_model
from celery import shared_task
from celery.exceptions import SoftTimeLimitExceeded
from celery.utils.log import get_task_logger
from celery_progress.backend import ProgressRecorder
from django.conf import settings

from .entrez_search_service import execute_entrez_search, create_random_filename, save_entrez_search_model, \
    download_esearch_protein_fasta_files, \
    update_entrezsearch_with_download_task_result, download_by_organism
from .models import ExternalTools, DomainDatabase
from .py_cdd_domain_search import produce_bokeh_pca_plot, write_domain_corrected_fasta_file
from .genbank_download_clinker_synteny import extract_assembly_ftp_paths_from_reciprocal_result_entries, \
    download_genbank_files, write_new_genbank_file, delete_sliced_genbank_files
from .py_services import check_if_target_sequences_are_available, check_if_msa_file_is_available, \
    slice_cdd_domain_corrected_fasta_file, slice_bw_query_fasta_file
from celery_blast.settings import BLAST_PROJECT_DIR, BLAST_DATABASE_DIR, CDD_DATABASE_URL
from os.path import isdir
from os import listdir, mkdir, remove
from refseq_transactions.tasks import download_refseq_assembly_summary
# logger for celery worker instances
logger = get_task_logger(__name__)


@shared_task(bind=True)
def calculate_phylogeny_based_on_database_statistics_selection(self, project_id:int, accession_data:dict):
    try:
        logger.info("starting phylogenetic inference based on selected RBHs from the database statistics plot")
        progress_recorder = ProgressRecorder(self)
        progress_recorder.set_progress(0, 100, "PROGRESS")
        logger.info("updating blast project with selection constrained cdd task")
        blast_project = get_project_by_id(project_id)
        update_blast_project_with_database_statistics_selection_task_result_model(project_id, str(self.request.id))
        logger.info("working with post data ...")
        accessions = []
        for ind in accession_data.keys():
            accessions.append(accession_data[ind])
            logger.info(accession_data[ind])

        path_to_query_subdir = blast_project.get_project_dir()
        path_to_fasta_file = path_to_query_subdir + '/bw_queries.faa'
        logger.info("trying to slice bw_queries.faa with accession list")
        path_to_sliced_fasta_file = slice_bw_query_fasta_file(path_to_query_subdir, path_to_fasta_file, accessions)
        logger.info("done slicing, created output: {}".format(path_to_sliced_fasta_file))
        progress_recorder.set_progress(15,100,"PROGRESS")
        logger.info("starting multiple sequence alignment")

        path_to_mafft_output = path_to_query_subdir + "/selection_sliced_mafft.msa"
        cmd = "mafft {} > {}".format(path_to_sliced_fasta_file, path_to_mafft_output)
        msa_task = Popen(cmd, shell=True)
        logger.info(
            'waiting for popen instance {} (mafft) to finish with timeout set to {}'.format(msa_task.pid,
                                                                                    40000))
        progress_recorder.set_progress(20, 100, "PROGRESS")

        returncode = msa_task.wait(40000)
        if returncode != 0:
            raise Exception("Popen hasnt succeeded, returncode != 0: {}".format(returncode))
        progress_recorder.set_progress(70, 100, "SUCCESS")
        logger.info("done calculating multiple sequence alignment with mafft")

        logger.info("starting phylogenetic inference with fasttree")
        path_to_fasttree_output = path_to_query_subdir + "/selection_sliced_phylogeny.nwk"
        cmd = "fasttree -lg {} > {}".format(path_to_mafft_output, path_to_fasttree_output)

        phylo_task = Popen(cmd, shell=True)
        progress_recorder.set_progress(71,100, "PROGRESS")
        logger.info(
            'waiting for popen instance {} (fasttree) to finish with timeout set to {}'.format(phylo_task.pid,
                                                                                    40000))
        returncode = phylo_task.wait(40000)
        if returncode != 0:
            raise Exception("Popen hasnt succeeded, returncode != 0: {}".format(returncode))
        logger.info("done with phylogenetic tree inference")
        logger.info("starting to generate html file with phylogeny ...")

        logger.info("generating table for shiptv metadata")

        reciprocal_result_table = pd.read_csv(path_to_query_subdir + '/reciprocal_results_with_taxonomy.csv' , index_col=0)
        reciprocal_result_table = reciprocal_result_table[reciprocal_result_table.sacc_transformed.isin(accessions)]
        path_to_rbh_table = path_to_query_subdir + '/selection_sliced_rbh_table.tsf'
        reciprocal_result_table.to_csv(path_to_rbh_table, sep="\t", index=None)

        path_to_html_tree = path_to_query_subdir + '/selection_sliced_phylogeny.html'
        shiptv_task = "shiptv --newick {} --metadata {} --output-html {}".format(path_to_fasttree_output,
                                                                                 path_to_rbh_table,
                                                                                 path_to_html_tree)
        shiptv_task = Popen(shiptv_task, shell=True)
        returncode = shiptv_task.wait(40000)
        if returncode != 0:
            raise Exception("Popen hasnt succeeded, returncode != 0: {}".format(returncode))
        logger.info("DONE")
        progress_recorder.set_progress(99,100,"PROGRESS")
    except Exception as e:
        raise Exception("[-] ERROR during inference of selection constrained phylogenetic inference on CDD bokeh plot"
                        " with exception: {}".format(e))



@shared_task(bind=True)
def calculate_phylogeny_based_on_selection(self, project_id:int, query_sequence:str, accession_data:dict):
    try:
        logger.info("starting phylogenetic inference based on selected RBHs from the CDD domain plot")
        progress_recorder = ProgressRecorder(self)
        progress_recorder.set_progress(0, 100, "PROGRESS")
        logger.info("updating external tools with selection constrained cdd task")
        external_tools = ExternalTools.objects.get_external_tools_based_on_project_id(project_id)
        external_tools.update_selection_constrained_CDD_phylogenetic_inference(query_sequence, str(self.request.id))
        blast_project = get_project_by_id(project_id)
        logger.info("working with post data ...")
        accessions = []
        for ind in accession_data.keys():
            accessions.append(accession_data[ind])
            logger.info(accession_data[ind])

        path_to_query_subdir = blast_project.get_project_dir() + '/' + query_sequence
        path_to_fasta_file = path_to_query_subdir + '/domain_corrected_target_sequences.faa'
        logger.info("trying to slice domain_corrected_target_sequences.faa with accession list")
        path_to_sliced_fasta_file = slice_cdd_domain_corrected_fasta_file(path_to_query_subdir, path_to_fasta_file, accessions)
        logger.info("done slicing, created output: {}".format(path_to_sliced_fasta_file))
        progress_recorder.set_progress(15,100,"PROGRESS")
        logger.info("starting multiple sequence alignment")

        path_to_mafft_output = path_to_query_subdir + "/selection_sliced_domain_mafft.msa"
        cmd = "mafft {} > {}".format(path_to_sliced_fasta_file, path_to_mafft_output)
        msa_task = Popen(cmd, shell=True)
        logger.info(
            'waiting for popen instance {} (mafft) to finish with timeout set to {}'.format(msa_task.pid,
                                                                                    40000))
        progress_recorder.set_progress(20, 100, "PROGRESS")

        returncode = msa_task.wait(40000)
        if returncode != 0:
            raise Exception("Popen hasnt succeeded, returncode != 0: {}".format(returncode))
        progress_recorder.set_progress(70, 100, "SUCCESS")
        logger.info("done calculating multiple sequence alignment with mafft")

        logger.info("starting phylogenetic inference with fasttree")
        path_to_fasttree_output = path_to_query_subdir + "/selection_sliced_domain_phylogeny.nwk"
        cmd = "fasttree -lg {} > {}".format(path_to_mafft_output, path_to_fasttree_output)

        phylo_task = Popen(cmd, shell=True)
        progress_recorder.set_progress(71,100, "PROGRESS")
        logger.info(
            'waiting for popen instance {} (fasttree) to finish with timeout set to {}'.format(phylo_task.pid,
                                                                                    40000))
        returncode = phylo_task.wait(40000)
        if returncode != 0:
            raise Exception("Popen hasnt succeeded, returncode != 0: {}".format(returncode))
        logger.info("done with phylogenetic tree inference")
        logger.info("starting to generate html file with phylogeny ...")
        path_to_rbh_table = path_to_query_subdir + '/rbh_table.tsf'
        path_to_html_tree = path_to_query_subdir + '/selection_sliced_domain_phylogeny.html'
        shiptv_task = "shiptv --newick {} --metadata {} --output-html {}".format(path_to_fasttree_output,
                                                                                 path_to_rbh_table, path_to_html_tree)
        shiptv_task = Popen(shiptv_task, shell=True)
        returncode = shiptv_task.wait(40000)
        if returncode != 0:
            raise Exception("Popen hasnt succeeded, returncode != 0: {}".format(returncode))
        logger.info("DONE")
        progress_recorder.set_progress(99,100,"PROGRESS")
    except Exception as e:
        raise Exception("[-] ERROR during inference of selection constrained phylogenetic inference on CDD bokeh plot"
                        " with exception: {}".format(e))

#TODO documentation + exception handling ...
'''synteny_calculation_task

'''
@shared_task(bind=True)
def synteny_calculation_task(self, project_id:int, query_sequence:str, rbh_dict:dict):
    try:
        logger.info("starting synteny calculation task for query sequence: {} in project {}".format(query_sequence, project_id))
        progress_recorder = ProgressRecorder(self)
        progress_recorder.set_progress(0, 100, "PROGRESS")

        returncode = delete_sliced_genbank_files(project_id, query_sequence)
        if returncode == 1:
            logger.warning("there is no directory for the given project_id: {} and query_sequence: {}".format(project_id,query_sequence))
            return 1
        external_tools = ExternalTools.objects.get_external_tools_based_on_project_id(project_id)
        external_tools.update_query_sequences_synteny_calculation_task(query_sequence, str(self.request.id))
        blast_project = get_project_by_id(project_id)

        result_data_path = blast_project.get_project_dir() + '/reciprocal_results.csv'
        target_rbhs = []
        for key in rbh_dict.keys():
            target_rbhs.append(rbh_dict[key])
        forward_database = blast_project.project_forward_database
        database_table_path = forward_database.path_to_database_file + '/'
        database_table_path += forward_database.get_pandas_table_name()

        database_genbank_file_path = forward_database.path_to_database_file + '/synteny_analysis'
        sequence_id_to_ftp_path_dict = extract_assembly_ftp_paths_from_reciprocal_result_entries(database_table_path,
                                                                                                 result_data_path,
                                                                                                 target_rbhs, project_id,
                                                                                                 query_sequence)
        progress_recorder.set_progress(20, 100, "PROGRESS")

        if isdir(database_genbank_file_path) == False:
            logger.info("creating synteny analysis folder for forward database: {} with filepath: {}".format(forward_database.database_name, forward_database.path_to_database_file))
            mkdir(database_genbank_file_path)

        logger.info("starting download of genbank files ...")
        genbank_filelist = download_genbank_files(sequence_id_to_ftp_path_dict,database_genbank_file_path, project_id)
        if len(genbank_filelist) == 0:
            logger.warning("there are genbank files left")
        else:
            progress_recorder.set_progress(75, 100, "PROGRESS")

            logger.info("done downloading starting to slice genbank file entries for synteny analysis")
            returncode = write_new_genbank_file(sequence_id_to_ftp_path_dict,
                                                database_genbank_file_path, 15, query_sequence, project_id)

            logger.info("trying to execute clinker")
            working_directory = BLAST_PROJECT_DIR + str(project_id) + '/' + query_sequence + '/'
            logger.info("executing clinker in working directory: {}".format(working_directory))
            cmd = "clinker "+working_directory+'*.gbk'+" -p "+working_directory+'clinker_result_plot.html'+' -i 0.25'
            proc = Popen(cmd, shell=True)
            returncode = proc.wait(timeout=2000)
            progress_recorder.set_progress(99, 100, "PROGRESS")
            if returncode != 0:
                logger.warning("error during donwload of the cdd database")
                raise SubprocessError

        logger.info("done")
        return 0
    except Exception as e:
        raise Exception("[-] ERROR during synteny calculation task with exception: {}".format(e))


'''check_database_integrity

    This function uses the blastdbcheck program to verify the database integrity. Before it starts
    verification it checks if the database directory exists. The database that has to be checked is specified
    by following code: path_to_database + database_name (e.g.: 'media/databases/1/' + 'high_quality_db.database').
    
    :param path_to_database
        :type str
    :param database_name
        :type str
    :returns valid
        :type boolean
'''


@shared_task(bind=True)
def check_database_integrity(self, path_to_database:str, database_name:str):
    try:
        progress_recorder = ProgressRecorder(self)
        progress_recorder.set_progress(0, 100, "PROGRESS")
        if isdir(path_to_database) == False:
            valid = False
        else:
            database_to_check = path_to_database + database_name
            verification = check_output(['blastdbcheck','-db', database_to_check])
            verification = str(verification)
            if "Result=SUCCESS" in verification and "No errors reported" in verification:
                valid = True
            else:
                valid = False
        progress_recorder.set_progress(99, 100, "PROGRESS")
        return valid
    except Exception as e:
        raise Exception("[-] ERROR during BLAST database integrity check with exception: {}".format(e))

'''setup_cathi_download_cdd_refseq_genbank_assembly_files

    Function for downloading the cdd database. The function checks if a directory exists and 
    removes previously downloaded files. It then downloads the CDD database from: CDD_DATABASE_URL
    After successfully downloading the CDD database, this function downloads the RefSeq and GenBank 
    assembly summary files.

'''


@shared_task(bind=True)
def setup_cathi_download_cdd_refseq_genbank_assembly_files(self):
    try:
        progress_recorder = ProgressRecorder(self)
        progress_recorder.set_progress(0, 100, "PROGRESS")

        domain_database_query_set = DomainDatabase.objects.all()

        if len(domain_database_query_set) != 1:
            DomainDatabase.objects.all().delete()
            domain_database_model = DomainDatabase(domain_database_loaded=False)
        else:
            domain_database_model = domain_database_query_set[0]
            domain_database_model.domain_database_loaded = False

        logger.info("trying to update domaindatabase model instance: {}".format(domain_database_model))
        update_domain_database_task_result_model(domain_database_model.id, self.request.id)


        cdd_path = BLAST_DATABASE_DIR + 'CDD/'
        if isdir(cdd_path) == False:
            logger.info("creating cdd directory in: {}".format(cdd_path))
            mkdir(cdd_path)
        else:
            logger.info("cdd directory exists")

        # cleaning up cdd directory
        for file in listdir(cdd_path):
            file_to_remove = cdd_path + file
            remove(file_to_remove)

        cdd_ftp_path = CDD_DATABASE_URL
        path_to_cdd_location = BLAST_DATABASE_DIR
        path_to_cdd_location = path_to_cdd_location + 'Cdd_LE.tar.gz'
        download_directory = BLAST_DATABASE_DIR + 'CDD'
        logger.info("downloading CDD into {}".format(download_directory))
        proc = Popen(["wget", cdd_ftp_path, "-q", "-O", path_to_cdd_location], shell=False)
        returncode = proc.wait(timeout=4000)
        progress_recorder.set_progress(70, 100, "PROGRESS")
        if returncode != 0:
            logger.warning("error during donwload of the cdd database")
            raise SubprocessError
        logger.info("successfully downloaded the cdd database")
        proc = Popen(
            ["tar", "-zxvf", path_to_cdd_location, "-C", download_directory + '/'],
            shell=False)
        returncode = proc.wait(timeout=2000)
        progress_recorder.set_progress(80, 100, "PROGRESS")
        if returncode != 0:
            logger.warning("error during extraction of the cdd database")
            raise SubprocessError
        logger.info("successfully extracted the cdd database")
        remove(path_to_cdd_location)

        domain_database_query_set = DomainDatabase.objects.all()
        domain_database_model = domain_database_query_set[0]
        domain_database_model.domain_database_loaded = True
        domain_database_model.save()

        try:
            download_refseq_assembly_summary("RefSeq")
            download_refseq_assembly_summary("GenBank")
        except Exception as e:
            logger.warning("error during download of the RefSeq and/or GenBank assembly summary files.")
            raise Exception("[-] ERROR during download of RefSeq and/or GenBank summary files with exception: {}.".format(e))
    except TimeoutExpired as e:
        logger.warning("error during downloading the cdd database, timeout expired with exception: {}".format(e))
        raise Exception("[-] ERROR timeout expired for downloading the cdd database with exception: {}".format(e))

    except Exception as e:
        raise Exception("[-] ERROR during downloading and decompressing the CDD database with exception: {}".format(e))

@shared_task(bind=True)
def download_entrez_search_associated_protein_sequences(self, search_id: int):
    try:
        progress_recorder = ProgressRecorder(self)

        logger.info("trying to download protein accessions with entrez")
        progress_recorder.set_progress(0, 100, "PROGRESS")

        returncode = update_entrezsearch_with_download_task_result(search_id, self.request.id)
        if returncode == 0:
            logger.info("starting download")
            returncode = download_esearch_protein_fasta_files(search_id)
            if returncode == 0:
                progress_recorder.set_progress(100, 100, "SUCCESS")
                return 0
            elif returncode == 1:
                logger.info("error in downloading protein sequences")
                return 0
        elif returncode == 1:
            logger.info("file already exists")
            return 0
        else:
            raise Exception(
                "Error during saving taskresult instance to download_task_result field of entrezsearch with id: {}".format(
                    search_id))

    except SoftTimeLimitExceeded:
        if 'returncode' in locals():
            if returncode != 0:
                logger.info("soft time limit exceeded for process with pid : {}".format(returncode))
        else:
            raise Exception("soft time limit exceeded for entrez downloading task ...")

    except Exception as e:
        raise Exception("[-] an error occurred during downloading fasta files with entrez: {}".format(e))


@shared_task(bind=True)
def entrez_search_task(self, database: str,number_records: str, entrez_query: str, user_id: int):
    try:
        logger.info("trying to start entrez search")
        progress_recorder = ProgressRecorder(self)
        progress_recorder.set_progress(0, 100, "PROGRESS")

        randomly_generated_filename = create_random_filename() + '.tsf'
        logger.info("filename for search output: {}".format(randomly_generated_filename))
        esearch_output_filepath = settings.ESEARCH_OUTPUT + 'result_dataframe_' + randomly_generated_filename

        entrez_search = save_entrez_search_model(database=database,
                                                 entrez_query=entrez_query,
                                                 file_name=esearch_output_filepath,
                                                 task_result_id=self.request.id,
                                                 user_id=user_id)

        returncode = execute_entrez_search(database, number_records, entrez_query, esearch_output_filepath, entrez_search)

        if returncode != 0:
            raise Exception("Popen hasnt succeeded, returncode != 0: {}".format(returncode))
        else:
            progress_recorder.set_progress(100, 100, "SUCCESS")
            entrez_search.update_paper_entries()
            return 0
    except SoftTimeLimitExceeded:
        if 'returncode' in locals():
            if returncode != 0:
                logger.info("soft time limit exceeded for process with pid : {}".format(returncode))
        else:
            raise Exception("soft time limit exceeded for entrez search task ...")
    except Exception as e:
        raise Exception("[-] Couldnt perform entrez search with exception: {}".format(e))


@shared_task(bind=True)
def download_organism_protein_sequences_task(self, search_id: int, organism_download: str, email: str):
    try:
        logger.info("trying to start organism sequence download")
        progress_recorder = ProgressRecorder(self)
        progress_recorder.set_progress(0, 100, "PROGRESS")

        returncode = download_by_organism(search_id, organism_download, email)

        if returncode != 0:
            logger.info("couldnt perform download task of organisms with exception : {}".format(returncode))
            return returncode
        else:
            progress_recorder.set_progress(100, 100, "SUCCESS")
            return 0


    except Exception as e:
        raise Exception("[-] couldnt perform download of organisms with exception: {}".format(e))


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
        path_to_project = settings.BLAST_PROJECT_DIR + str(project_id) + '/'
        path_to_query_file = path_to_project + query_sequence_id + '/target_sequences.faa'

        target_sequence_status = check_if_target_sequences_are_available(path_to_query_file)
        if target_sequence_status == 0:
            # mafft invocation with default settings
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
def execute_phylogenetic_tree_building(self, project_id, query_sequence_id):
    try:
        logger.info(
            "trying to execute fasttree phylogenetic tree construction per request to bioinformatic tools container")

        progress_recorder = ProgressRecorder(self)
        progress_recorder.set_progress(0, 100, "PROGRESS")

        external_tools = ExternalTools.objects.get_external_tools_based_on_project_id(project_id)
        logger.info("cheking if msa task succeeded for query sequence : {}".format(query_sequence_id))

        path_to_project = settings.BLAST_PROJECT_DIR + str(project_id) + '/'
        path_to_msa_file = path_to_project + query_sequence_id + '/target_sequences.msa'
        path_to_fasttree_output = path_to_project + query_sequence_id + '/target_sequences.nwk'
        msa_status = check_if_msa_file_is_available(path_to_msa_file)
        if msa_status == 0 and external_tools.check_if_msa_task_is_completed(query_sequence_id):

            external_tools.update_query_sequences_phylo_task(query_sequence_id, str(self.request.id))
            logger.info("updated query sequence model with taskresult instance : {}".format(str(self.request.id)))
            progress_recorder.set_progress(20, 100, "PROGRESS")
            cmd = "fasttree -lg {} > {}".format(path_to_msa_file, path_to_fasttree_output)
            phylo_task = Popen(cmd, shell=True)
            progress_recorder.set_progress(30, 100, "PROGRESS")
            logger.info(
                'waiting for popen instance {} to finish with timeout set to {}'.format(phylo_task.pid,
                                                                                        40000))
            returncode = phylo_task.wait(40000)
            if returncode != 0:
                raise Exception("Popen hasnt succeeded, returncode != 0: {}".format(returncode))

        elif msa_status == 1:
            raise FileNotFoundError("msa file does not exist!")

        progress_recorder.set_progress(100, 100, "SUCCESS")

        return 0
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

        path_to_project = settings.BLAST_PROJECT_DIR + str(project_id) + '/'

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
                        'progress of multiple sequence alignment task for all query sequences set to {}'.format(
                            progress))
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
        raise Exception(
            "[-] Couldnt perform multiple sequence alignment task for multiple query sequences with Exception: {}".format(
                e))


@shared_task(bind=True)
def execute_fasttree_phylobuild_for_all_query_sequences(self, project_id):
    try:
        logger.info("trying to execute fasttree phylogenetic tree construction for all query sequences")
        progress_recorder = ProgressRecorder(self)
        progress_recorder.set_progress(0, 100, "PROGRESS")
        external_tools = ExternalTools.objects.get_external_tools_based_on_project_id(project_id)

        external_tools.update_for_all_query_sequences_phylo_task(str(self.request.id))
        logger.info("updated multiple query sequence models with taskresult instance : {}".format(str(self.request.id)))
        path_to_project = settings.BLAST_PROJECT_DIR + str(project_id) + '/'

        query_sequence_ids = []
        for qseq in external_tools.query_sequences.get_queryset():
            logger.info("check if msa task succeeded for query sequence : {}".format(qseq.query_accession_id))
            if external_tools.check_if_msa_task_is_completed(qseq.query_accession_id):
                query_sequence_ids.append(qseq.query_accession_id)
                logger.info(
                    "\tmsa task succeeded ... added query sequence to target list for phylogenetic tree construction")
        progress = 80 / len(query_sequence_ids)
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
                returncode = phylo_task.wait(40000)

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


'''cdd_domain_search_with_rbhs_task
        
    This function executes a rpsblast for searching in the CDD database. Settings for the rpsblast are user specific 
    and derive from the underlying form instance.
    
    :param self
        :type TaskResult
    :param project_id
        :type int
    :param rps_blast_task_data -> form data 
        :type dict[str] = value
'''


@shared_task(bind=True)
def cdd_domain_search_with_rbhs_task(self, project_id: int, rps_blast_task_data: dict):
    try:

        progress_recorder = ProgressRecorder(self)
        progress_recorder.set_progress(0, 100, "STARTED")
        target_query = rps_blast_task_data['query_sequence']
        logger.info("INFO:started domain search in the cdd database for project: {} and sequence {}".format(project_id,
                                                                                                            target_query))
        try:

            # project_id: int, query_sequence: str, task_id: int
            update_external_tool_with_cdd_search(project_id, target_query, self.request.id)
            progress_recorder.set_progress(25, 100, "PROGRESS")

        except Exception as e:
            logger.warning("ERROR: couldnt update blast project with exception: {}".format(e))
            raise Exception("ERROR: couldnt update blast project with exception: {}".format(e))

        path_to_project = settings.BLAST_PROJECT_DIR + str(project_id) + '/'
        path_to_cdd_db = settings.CDD_DIR
        path_to_query_file = path_to_project + target_query + '/target_sequences.faa'
        path_to_cdd_domain_search_output = path_to_project + target_query + '/cdd_domains.tsf'

        logger.info("INFO:preparing POPEN cmd for cdd search")

        proc = Popen(['rpsblast', '-query', path_to_query_file,
                      '-db', path_to_cdd_db,
                      "-outfmt", "6 qseqid qlen sacc slen qstart qend sstart send bitscore evalue pident",
                      "-out", path_to_cdd_domain_search_output,
                      '-evalue', str(rps_blast_task_data['rps_e_value']),
                      '-num_threads', str(rps_blast_task_data['rps_num_threads']),
                      '-max_hsps', str(rps_blast_task_data['rps_max_hsps']),
                      '-num_alignments', str(rps_blast_task_data['rps_num_alignments'])], shell=False)
        progress_recorder.set_progress(30, 100, "PROGRESS")
        logger.info("INFO:waiting for POPEN process with id: {} to finish".format(proc.pid))

        returncode = proc.wait(timeout=settings.SUBPROCESS_TIME_LIMIT)
        progress_recorder.set_progress(50, 100, "PROGRESS")

        if returncode != 0:
            raise SubprocessError

        logger.info("INFO:performing PCA analysis and build interactive visualization ...")
        produce_bokeh_pca_plot(project_id, target_query, taxonomic_unit='phylum')
        progress_recorder.set_progress(60, 100, "PROGRESS")

        logger.info("INFO:starting to write domain corrected fasta file ...")
        write_domain_corrected_fasta_file(project_id=project_id, query_sequence=target_query)
        progress_recorder.set_progress(65, 100, "PROGRESS")

        logger.info("INFO:starting mafft multiple sequence alignment with inferred domain segments ...")
        execute_domain_multiple_sequence_alignment(project_id, target_query)
        progress_recorder.set_progress(80, 100, "PROGRESS")

        logger.info("INFO:starting fasttree phylogenetic tree reconstruction with inffered domain msa ...")
        execute_phylogenetic_tree_building_with_domains(project_id, target_query)
        progress_recorder.set_progress(99, 100, "PROGRESS")

        logger.info("INFO:DONE")
        return returncode

    except SoftTimeLimitExceeded as e:
        logger.info("ERROR:CDD search reached Task Time Limit")
        if 'proc' in locals():
            pid = proc.pid
            # TODO check if process with pid exists
            parent = psutil.Process(pid)
            for child in parent.children(recursive=True):
                child.kill()
            parent.kill()
        raise Exception("ERROR CDD search reached Task Time Limit")

    except TimeoutExpired as e:
        logger.warning("ERROR:TIMEOUT EXPIRED DURING CDD SEARCH: {}".format(e))
        if 'proc' in locals():
            pid = proc.pid
            parent = psutil.Process(pid)

            for child in parent.children(recursive=True):
                child.kill()
            parent.kill()

    except SubprocessError as e:
        logger.warning("ERROR: POPEN CALL for cdd domain search resulted in exception: {}".format(e))
        if 'proc' in locals():

            pid = proc.pid
            parent = psutil.Process(pid)

            for child in parent.children(recursive=True):
                child.kill()
            parent.kill()
    except Exception as e:
        raise Exception(
            "ERROR: unknown exception occurred during cdd_domain_search_with_rbhs_task with exception: {}".format(e))


'''execute_domain_multiple_sequence_alignment
    
    This function executes a MAFFT multiple sequence alignment with overlap corrected target sequences.
    The overlap correction is based on the write_domain_corrected_fasta_file function, which creates a fasta file
    based on domain sequence segments.
    
    :param project_id
        :type int
    :param query_sequence_id
        :type str 
'''


@shared_task()
def execute_domain_multiple_sequence_alignment(project_id: int, query_sequence_id: str):
    try:
        logger.info("INFO: starting mafft multiple sequence alignment with RBHs of {} and inferred domains".format(
            query_sequence_id))
        path_to_project = settings.BLAST_PROJECT_DIR + str(project_id) + '/'
        path_to_query_file = path_to_project + query_sequence_id + '/domain_corrected_target_sequences.faa'

        target_sequence_status = check_if_target_sequences_are_available(path_to_query_file)
        if target_sequence_status == 0:
            # mafft invocation with default settings
            path_to_mafft_output = path_to_project + query_sequence_id + '/domain_corrected_target_sequences.msa'
            cmd = "mafft {} > {}".format(path_to_query_file, path_to_mafft_output)
            msa_task = Popen(cmd, shell=True)
            logger.info(
                'INFO: waiting for popen mafft instance {} to finish with timeout set to {}'.format(msa_task.pid,
                                                                                                    40000))

            returncode = msa_task.wait(40000)
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


'''execute_phylogenetic_tree_building_with_domains
    
    Function for phylogenetic tree reconstruction with the msa file of the overlap corrected target sequences.
    This function can only get executed after the MAFFT msa completed.
    
    :param project_id
        :type int
    :param query_sequence_id
        :type str
    
'''


@shared_task()
def execute_phylogenetic_tree_building_with_domains(project_id: int, query_sequence_id: str):
    try:
        logger.info("INFO: starting phylogenetic tree reconstruction with fasttree")

        logger.info("cheking if msa task succeeded for query sequence : {}".format(query_sequence_id))

        path_to_project = settings.BLAST_PROJECT_DIR + str(project_id) + '/'
        path_to_msa_file = path_to_project + query_sequence_id + '/domain_corrected_target_sequences.msa'
        path_to_fasttree_output = path_to_project + query_sequence_id + '/domain_corrected_target_sequences.nwk'
        msa_status = check_if_msa_file_is_available(path_to_msa_file)
        if msa_status == 0:

            cmd = "fasttree -lg {} > {}".format(path_to_msa_file, path_to_fasttree_output)
            phylo_task = Popen(cmd, shell=True)
            logger.info(
                'INFO: waiting for fasttree popen instance {} to finish with timeout set to {}'.format(phylo_task.pid,
                                                                                                       40000))
            returncode = phylo_task.wait(40000)
            if returncode != 0:
                raise Exception("Popen hasnt succeeded, returncode != 0: {}".format(returncode))
        elif msa_status == 1:
            raise FileNotFoundError("msa file does not exist!")
        return 0
    except Exception as e:
        raise Exception("[-] Couldnt perform phylogenetic tree task with Exception: {}".format(e))
