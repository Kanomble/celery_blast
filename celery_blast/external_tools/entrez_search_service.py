import os
import subprocess
import string
import random
from django.db import transaction, IntegrityError
from .models import EntrezSearch
from Bio import Entrez
from django.contrib.auth.models import User
from django_celery_results.models import TaskResult

def execute_entrez_search(database: str, entrez_query: str, output_filepath: str, entrez_search: EntrezSearch) -> int:
    try:
        xtract_format = {}
        xtract_format['pubmed'] = 'Id PubDate Source Author Title ELocationID'
        xtract_format['protein'] = 'Id Caption Title Organism'

        cmd = 'esearch -db {} -query "{}" | efetch -format docsum | xtract -pattern DocumentSummary -sep "\t" -element {} > {}'.format(
            database, entrez_query, xtract_format[database], output_filepath)

        process = subprocess.Popen(cmd, shell=True)
        try:
            returncode = process.wait(timeout=20000)

            if returncode != 0:
                with transaction.atomic():
                    if os.path.isfile(entrez_search.file_name):
                        os.remove(entrez_search.file_name)
                    entrez_search.delete()

            return returncode
        except subprocess.TimeoutExpired as e:
            process.kill()
            returncode = 'searchtime expired {}'.format(e)
            return returncode
    except Exception:
        return 1

#TODO implementation
def download_esearch_protein_fasta_files(search_id:int):
    try:
        entrez_search = get_entrezsearch_object_with_entrezsearch_id(search_id)

        file_path = entrez_search.file_name
        file_random_number = file_path.split("_")[-1].split(".")[0]

        target_fasta_file_path = 'media/esearch_output/target_fasta_file_' + str(file_random_number) + '.faa'
        if os.path.isfile(target_fasta_file_path):
            if entrez_search.download_task_result != None:
                if entrez_search.download_task_result.status == 'SUCCESS':
                    return 1

        elif entrez_search.download_task_result != None:
            if entrez_search.download_task_result.status == 'SUCCESS':
                return 1

        Entrez.email = entrez_search.entrez_user.email
        pandas_table = entrez_search.get_pandas_table()
        sequence_ids = list(pandas_table['Caption'])

        if len(sequence_ids) > 500:
            sequence_ids = sequence_ids[0:499]

        handle = Entrez.efetch(db="protein", id=sequence_ids, retmode="xml")
        record = Entrez.read(handle)
        handle.close()

        with open(target_fasta_file_path,'w') as target_fasta_file:
            for rec in record:
                target_fasta_file.write('>' + rec['GBSeq_locus'] + ' ' + rec['GBSeq_definition'] + "\n")
                target_fasta_file.write(rec['GBSeq_sequence'] + "\n")

        #update entrezsearch instance with fastafile path
        with transaction.atomic():
            entrez_search.fasta_file_name = target_fasta_file_path
            entrez_search.save()

        return 0
    except Exception as e:
        raise Exception("[-] Couldnt download protein fasta files for esearch {} on protein database and exception: {}".format(search_id,e))

#TODO documentation
def update_entrezsearch_with_download_task_result(search_id,task_id):
    try:
        entrez_search = get_entrezsearch_object_with_entrezsearch_id(search_id)

        file_path = entrez_search.file_name
        file_random_number = file_path.split("_")[-1].split(".")[0]
        target_fasta_file_path = 'media/esearch_output/target_fasta_file_' + str(file_random_number) + '.faa'

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
        raise Exception("[-] Couldnt update entrezsearch with taskresult instance for downloading protein accessions with exception: {}".format(e))

#TODO implementation documentation
def execute_entrez_efetch_fasta_files(database: str, entrez_query: str, output_filepath: str) -> int:
    try:
        if database == 'pubmed':
            cmd = 'esearch -db pubmed -query "{}" | elink -target protein | efilter -source refseq | efetch -format fasta > {}'.format(
                entrez_query, output_filepath)
        elif database == 'protein':
            cmd = 'esearch -db protein -query "{}" | efilter -source refseq | efetch -format fasta > {}'.format(
                entrez_query,output_filepath
            )
        else:
            raise Exception
        process = subprocess.Popen(cmd, shell=True)
        try:
            returncode = process.wait(timeout=20000)
            return returncode
        except subprocess.TimeoutExpired as e:
            process.kill()
            returncode = 'searchtime expired {}'.format(e)
            return returncode
    except Exception:
        return 1


def create_random_filename() -> str:
    randomly_generated_filename = ''.join(random.choices(string.ascii_uppercase + string.digits, k = 10))
    esearch_output_filepath = 'media/esearch_output/result_dataframe_' + randomly_generated_filename
    while os.path.isfile(esearch_output_filepath) == True:
        randomly_generated_filename = ''.join(random.choices(string.ascii_uppercase + string.digits, k=10))
        esearch_output_filepath = 'media/esearch_output/result_dataframe_' + randomly_generated_filename

    return randomly_generated_filename

def get_entrezsearch_object_with_entrezsearch_id(search_id):
    try:
        entrez_search = EntrezSearch.objects.get(id=search_id)
        return entrez_search
    except Exception as e:
        raise Exception("[-] There is no entrez_search object with id: {} exception: {}".format(search_id, e))

#TODO delete associated TaskResult object
def delete_esearch_by_id(search_id):
    try:
        with transaction.atomic():
            entrez_search = EntrezSearch.objects.get(id=search_id)
            if entrez_search != None:
                if os.path.isfile(entrez_search.file_name):
                    os.remove(entrez_search.file_name)
                if entrez_search.fasta_file_name != None:
                    if os.path.isfile(entrez_search.fasta_file_name):
                        os.remove(entrez_search.fasta_file_name)
                entrez_search.delete()
                return 0
            else:
                return 1
    except Exception as e:
        raise IntegrityError("ERROR during deletion of entrez_search with id : {} with exception: {}".format(search_id,e))

def save_entrez_search_model(database:str,entrez_query:str,file_name:str, task_result_id:int, user_id:int) -> EntrezSearch:
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
        raise IntegrityError("[-] An error occcurred during saving the edirect search object into the database: {}".format(e))

