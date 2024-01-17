import os
import random
import shutil
import string
from subprocess import Popen, TimeoutExpired, SubprocessError

import psutil
from Bio import Entrez
from celery_blast.settings import ESEARCH_OUTPUT
from django.conf import settings
from django.contrib.auth.models import User
from django.db import transaction, IntegrityError
from django_celery_results.models import TaskResult

from .models import EntrezSearch

'''execute_entrez_search

    This function executes the search requested by the user. It is utilized in the external_tools/tasks.py
    script within the entrez_search_task function.
    
    :param database - from form
        :type str
    
    :param number_records
        :type str
    
    :param entrez_query - from form
        :type str
        
    :param output_filepath - s. ESEARCH_OUTPUT
        :type str
        
    :param entrez_search
        :type EntrezSearch
        
    :returns returncode
        :type int
'''
def execute_entrez_search(database: str, number_records: str, entrez_query: str, output_filepath: str, entrez_search: EntrezSearch) -> int:
    # starts an Esearch which downloads the requested search with the selected database and saves it in a file
    # for adding  more databases, the columns need to be added here, in the models.py get_pandas_table function and in forms.py to the EntrezSearchForm class
    try:
        xtract_format = {}
        xtract_format['pubmed'] = 'Id PubDate Source Author Title ELocationID'
        xtract_format['protein'] = 'Id Caption Title Organism'
        xtract_format['assembly'] = 'Id AssemblyName AssemblyStatus Organism Taxid'
        xtract_format['cdd'] = "Id Title Subtitle Abstract"
        xtract_format['protfam'] = "Id DispMethod DispReviewLevel string"

        cmd = 'esearch -db {} -query "{}" | efetch -format docsum -start 1 -stop {} | xtract -pattern DocumentSummary -sep "\t" -sep ": "  -element {} > {}'.format(
            database, entrez_query,number_records, xtract_format[database], output_filepath)

        process = Popen(cmd, shell=True)
        try:
            returncode = process.wait(timeout=20000)

            if returncode != 0:
                with transaction.atomic():
                    if os.path.isfile(entrez_search.file_name):
                        os.remove(entrez_search.file_name)
                    entrez_search.delete()

            return returncode
        except TimeoutExpired as e:
            process.kill()
            returncode = 'searchtime expired {}'.format(e)
            return returncode
    except Exception:
        return 1

'''download_by_organism

    This function is executed within the entrez search dashboard details page.
    It downloads protein fasta files for the specified organism via biopython.
    
    :param search_id
        :type int
    
    :param organism_download
        :type str
    
    :param email
        :type str
        
    :returns returncode
        :type int

'''
def download_by_organism(search_id: int, organism_download: str, email: str) -> int:
    # uses biopython entrez tool to download fasta files of a selected organism in an entrezsearch paper and saves it in a file
    try:
        if os.path.isdir(ESEARCH_OUTPUT + str(search_id)) == False:
            os.mkdir(ESEARCH_OUTPUT + str(search_id))

        with open(ESEARCH_OUTPUT + str(search_id) + "/download_by_organism.log", 'w') as logfile:
            logfile.write("INFO:starting to download proteins for the selected organism: {}\n".format(organism_download))
            entrez_search = get_entrezsearch_object_with_entrezsearch_id(search_id)

            organism_target_df = entrez_search.get_pandas_table()

            organism = organism_download
            filtered_organism_target_df = organism_target_df[organism_target_df.Organism == organism]
            filtered_organism_target_id_list = filtered_organism_target_df["Id"].to_list()
            Entrez.email = email  # use the email of the user
            fasta_file = ESEARCH_OUTPUT + str(search_id) + "/" + str(organism) + ".faa"
            logfile.write("INFO:starting to write fasta file: {}\n".format(fasta_file))
            with open(fasta_file,'w') as output:  # could use the name of the file used as input

                end = len(filtered_organism_target_id_list)
                begin = 0
                logfile.write("INFO:total number of ids: {}\nINFO:splitting ids into chunks of 500\n".format(end))
                step = 500
                steps = 500
                while begin < end:
                    if step >= end:
                        step = end
                    splitted_ids = filtered_organism_target_id_list[begin:step]
                    # range(X) tries for biopython calls
                    for attempt in range(10):
                        try:
                            handle = Entrez.efetch(id=splitted_ids, db="protein", retmode="xml")
                            record = Entrez.read(handle)
                            handle.close()
                        except Exception as e:
                            if attempt == 9:
                                logfile.write(
                                    "ERROR:biopython inference of proteins failed with exception {}\n".format(e))
                        else:
                            for rec in record:
                                output.write('>' + rec['GBSeq_locus'] + ' ' + rec['GBSeq_definition'] + "\n")
                                output.write(rec['GBSeq_sequence'] + "\n")
                            output.close()
                            break
                    begin += steps
                    step += steps
        return 0
    except Exception:
        return 1


'''download_esearch_protein_fasta_files
    
    This function is part of the download_entrez_search_associated_protein_sequences task in external_tools/tasks.py.
    The function downloads the protein sequences listed in the entrez search output. It uses the efetch tool for 
    the downloading process.
    
    If the search was based on the pubmed database, associated protein sequences are downloaded via the elink and 
    efetch tool. 
    
    The output of this function is a sequence file and a table containing the listed identifier. Both files 
    are created with a random character string.
    
    :param search_id
        :type int
        
    :returns returncode
        :type int
'''
def download_esearch_protein_fasta_files(search_id: int) -> int:
    # downloads a fasta file of the associated protein sequences of an enztezsearch if the entrez database is protein or pubmed
    try:
        entrez_search = get_entrezsearch_object_with_entrezsearch_id(search_id)

        database = entrez_search.database
        entrez_query = entrez_search.entrez_query

        file_path = entrez_search.file_name
        file_random_number = file_path.split("_")[-1].split(".")[0]

        target_fasta_file_path = ESEARCH_OUTPUT + 'target_fasta_file_' + str(file_random_number) + '.faa'
        sequence_list_file_path = ESEARCH_OUTPUT + 'sequence_list_' + str(file_random_number) + '.table'

        if os.path.isfile(target_fasta_file_path):
            if entrez_search.download_task_result != None:
                if entrez_search.download_task_result.status == 'SUCCESS':
                    return 1

        elif entrez_search.download_task_result != None:
            if entrez_search.download_task_result.status == 'SUCCESS':
                return 1
        if database == "protein":

            pandas_table = entrez_search.get_pandas_table()
            sequence_ids = list(pandas_table['Caption'])
            with open(sequence_list_file_path, 'w') as seq_id_file:
                for seqid in sequence_ids:
                    seq_id_file.write("{}\n".format(seqid))

            cmd = 'efetch -db protein -format fasta -input {} > {}'.format(sequence_list_file_path,
                                                                           target_fasta_file_path)
            process = Popen(cmd, shell=True)
            returncode = process.wait(timeout=settings.SUBPROCESS_TIME_LIMIT)

            with transaction.atomic():
                entrez_search.fasta_file_name = target_fasta_file_path
                entrez_search.save()

            return returncode
        elif database == "pubmed":
            cmd = 'esearch -db pubmed -query "{}" | elink -target protein | efetch -format fasta -start 1 -stop 100 > {}'.format(
                entrez_query, target_fasta_file_path)
            process = Popen(cmd, shell=True)
            returncode = process.wait(timeout=settings.SUBPROCESS_TIME_LIMIT)

            with transaction.atomic():
                entrez_search.fasta_file_name = target_fasta_file_path
                entrez_search.save()
            return returncode
        else:
            raise Exception("There is no download option for the selected database: {}".format(database))
    # catch subprocess exceptions
    except TimeoutExpired:
        # delete all child processes
        if 'process' in locals():
            pid = process.pid
            parent = psutil.Process(pid)
            for child in parent.children(recursive=True):
                child.kill()
            parent.kill()  # this is not enough need to kill all child processes
            if os.path.isfile(target_fasta_file_path):
                os.remove(target_fasta_file_path)
            return pid
        else:
            return 1
    except SubprocessError as e:
        raise Exception(
            "[-] Couldnt download protein fasta files for esearch {} on protein database and exception: {}".format(
                search_id, e))


'''update_entrezsearch_with_download_task_result
    
    Updates the EntrezSearch database model object with the associated celery TaskResult database model.
    
    :param search_id - identifier for the EntrezSearch object
        :type int
    
    :param task_id - identifier for the associated TaskResult object
        :type int

'''
def update_entrezsearch_with_download_task_result(search_id: int, task_id: int) -> int:
    try:
        entrez_search = get_entrezsearch_object_with_entrezsearch_id(search_id)

        file_path = entrez_search.file_name
        file_random_number = file_path.split("_")[-1].split(".")[0]
        target_fasta_file_path = ESEARCH_OUTPUT + 'target_fasta_file_' + str(file_random_number) + '.faa'

        if os.path.isfile(target_fasta_file_path):
            if entrez_search.download_task_result != None:
                if entrez_search.download_task_result.status == 'SUCCESS' or entrez_search.download_task_result.status == 'PROGRESS':
                    return 1
        elif entrez_search.download_task_result != None:
            if entrez_search.download_task_result.status != 'FAILURE':
                return 1

        with transaction.atomic():
            task_result = TaskResult.objects.get(task_id=task_id)
            entrez_search = get_entrezsearch_object_with_entrezsearch_id(search_id)
            entrez_search.download_task_result = task_result
            entrez_search.save()
        return 0
    except Exception as e:
        raise Exception(
            "[-] Couldnt update entrezsearch with taskresult instance for downloading protein accessions with exception: {}".format(
                e))


# TODO not in use
'''execute_entrez_efetch_fasta_files

    This function is currently not in use.

'''
def execute_entrez_efetch_fasta_files(database: str, entrez_query: str, output_filepath: str) -> int:
    try:
        if database == 'pubmed':
            cmd = 'esearch -db pubmed -query "{}" | elink -target protein | efilter -source refseq | efetch -format fasta > {}'.format(
                entrez_query, output_filepath)
        elif database == 'protein':
            cmd = 'esearch -db protein -query "{}" | efilter -source refseq | efetch -format fasta > {}'.format(
                entrez_query, output_filepath
            )
        else:
            raise Exception
        process = Popen(cmd, shell=True)

        returncode = process.wait(timeout=settings.SUBPROCESS_TIME_LIMIT)
        return returncode

    except TimeoutExpired as e:
        if 'process' in locals():
            pid = process.pid
            parent = psutil.Process(pid)
            for child in parent.children(recursive=True):
                child.kill()
            parent.kill()  # this is not enough need to kill all child processes
            if os.path.isfile(output_filepath):
                os.remove(output_filepath)
            return pid
        else:
            return 1
    except Exception:
        return 1

'''create_random_filename

    This function is used to create filenames for the output of the entrez search.

'''
def create_random_filename() -> str:
    # creates a rondom file name and returns it
    randomly_generated_filename = ''.join(random.choices(string.ascii_uppercase + string.digits, k=10))
    esearch_output_filepath = 'media/esearch_output/result_dataframe_' + randomly_generated_filename
    while os.path.isfile(esearch_output_filepath) == True:
        randomly_generated_filename = ''.join(random.choices(string.ascii_uppercase + string.digits, k=10))
        esearch_output_filepath = 'media/esearch_output/result_dataframe_' + randomly_generated_filename

    return randomly_generated_filename


'''get_entrezsearch_object_with_entrezsearch_id
    
    Returns the associated EntrezSearch database model object.
    
    :param search_id
        :type int
        
    :returns entrez_search
        :type EntrezSearch

'''
def get_entrezsearch_object_with_entrezsearch_id(search_id: int) -> EntrezSearch:
    # gets an entrez search database row based on a search_id and returns it
    try:
        entrez_search = EntrezSearch.objects.get(id=search_id)
        return entrez_search
    except Exception as e:
        raise Exception("[-] There is no entrez_search object with id: {} exception: {}".format(search_id, e))


# TODO also remove sequence_list_RNDNumber
'''delete_esearch_by_id
    
    Deletes the associated EntrezSearch object and all relevant files.
    
    :param search_id
        :type int
    
    :returns returncode
        :type int

'''
def delete_esearch_by_id(search_id: int) -> int:
    # deletes entrezsearch associated files and taskresult entries based on a search_id
    # returns 0 if it worked or 1 if it did not
    try:
        with transaction.atomic():

            entrez_search = EntrezSearch.objects.get(id=search_id)
            task_db_id = EntrezSearch.objects.values('search_task_result_id').filter(id=search_id)
            task_db_id = task_db_id[0]['search_task_result_id']
            task_db = TaskResult.objects.get(id=task_db_id)
            organism_db_file_name = ESEARCH_OUTPUT + str(search_id)

            if os.path.isdir(organism_db_file_name):
                shutil.rmtree(organism_db_file_name)

            if entrez_search != None:
                if os.path.isfile(entrez_search.file_name):
                    os.remove(entrez_search.file_name)
                if entrez_search.fasta_file_name != None:
                    if os.path.isfile(entrez_search.fasta_file_name):
                        os.remove(entrez_search.fasta_file_name)
                entrez_search.delete()
                task_db.delete()

                return 0
            else:
                return 1
    except Exception as e:
        raise IntegrityError(
            "ERROR during deletion of entrez_search with id : {} with exception: {}".format(search_id, e))


'''save_entrez_search_model

    Wrapper function to create the EntrezSearch model isntance from the Django form data.
    
    :param database
        :type str
    :param entrez_query
        :type str
    :param file_name - with random ending
        :type str
    :param task_result_id
        :type int
    :param user_id
        :type int
    
    :returns esearch_object - created EntrezSearch model
        :type EntrezSearch
'''
def save_entrez_search_model(database: str, entrez_query: str, file_name: str, task_result_id: int,
                             user_id: int) -> EntrezSearch:
    # saves users entrezsearch input into a datbase with extra information
    # returns the database row
    try:
        with transaction.atomic():
            user = User.objects.get(id=user_id)
            esearch_object = EntrezSearch.edirect_objects.create_entrez_search(database=database,
                                                                               entrez_query=entrez_query,
                                                                               file_path=file_name,
                                                                               search_task_result=task_result_id,
                                                                               entrez_user=user)
            return esearch_object
    except Exception as e:
        raise IntegrityError(
            "[-] An error occcurred during saving the edirect search object into the database: {}".format(e))
