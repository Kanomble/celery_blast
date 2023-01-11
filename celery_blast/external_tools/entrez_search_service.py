import os
import psutil
import subprocess
import string
import random
import shutil
from django.db import transaction, IntegrityError
from .models import EntrezSearch
from django.contrib.auth.models import User
from django_celery_results.models import TaskResult
from django.conf import settings
from Bio import Entrez

def execute_entrez_search(database: str, entrez_query: str, output_filepath: str, entrez_search: EntrezSearch) -> int:
    #starts an Esearch which downloads the requested search with the selected database and saves it in a file
    # for adding  more databases, the columns need to be added here, in the models.py get_pandas_table function and in forms.py to the EntrezSearchForm class
    try:
        xtract_format = {}
        xtract_format['pubmed'] = 'Id PubDate Source Author Title ELocationID'
        xtract_format['protein'] = 'Id Caption Title Organism'
        xtract_format['assembly'] = 'Id AssemblyName AssemblyStatus Organism Taxid'
        xtract_format['cdd'] = "Id Title Subtitle Abstract"
        xtract_format['protfam'] = "Id DispMethod DispReviewLevel string"


        cmd = 'esearch -db {} -query "{}" | efetch -format docsum | xtract -pattern DocumentSummary -sep "\t" -sep ": "  -element {} > {}'.format(
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

def download_by_organism(search_id: int, organism_download: str, email: str) -> int:
    #uses biopython entrez tool to download fasta files of a selected organism in an entrezsearch paper and saves it in a file
    try:
        entrez_search = get_entrezsearch_object_with_entrezsearch_id(search_id)

        organism_target_df = entrez_search.get_pandas_table()

        organism = organism_download
        filtered_organism_target_df = organism_target_df[organism_target_df.Organism == organism]
        filtered_organism_target_id_list = filtered_organism_target_df["Id"].to_list()
        Entrez.email = email  # use the email of the user
        if os.path.isdir("media/esearch_output/"+ str(search_id))== False:
            os.mkdir("media/esearch_output/"+ str(search_id))
        output = open("media/esearch_output/"+ str(search_id)+"/"+str(organism)+".faa", 'w')  # could use the name of the file used as input

        end = len(filtered_organism_target_id_list)
        begin = 0
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
                            "ERROR:inference of taxonomic informations failed for query sequence dataframe of {} with exception {}\n".format(
                                query, e))
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



#TODO implementation
def download_esearch_protein_fasta_files(search_id:int) -> int:
    #downloads a fasta file of the associated protein sequences of an enztezsearch if the entrez database is protein or pubmed
    try:
        entrez_search = get_entrezsearch_object_with_entrezsearch_id(search_id)

        database = entrez_search.database
        entrez_query = entrez_search.entrez_query

        file_path = entrez_search.file_name
        file_random_number = file_path.split("_")[-1].split(".")[0]

        target_fasta_file_path = 'media/esearch_output/target_fasta_file_' + str(file_random_number) + '.faa'
        sequence_list_file_path = 'media/esearch_output/sequence_list_' + str(file_random_number) + '.table'

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
            with open(sequence_list_file_path,'w') as seq_id_file:
                for seqid in sequence_ids:
                    seq_id_file.write("{}\n".format(seqid))


            cmd = 'efetch -db protein -format fasta -input {} > {}'.format(sequence_list_file_path,target_fasta_file_path)
            process = subprocess.Popen(cmd, shell=True)
            returncode = process.wait(timeout=settings.SUBPROCESS_TIME_LIMIT)

            with transaction.atomic():
                entrez_search.fasta_file_name = target_fasta_file_path
                entrez_search.save()

            return returncode
        if database == "pubmed":

            cmd = 'esearch -db pubmed -query "{}" | elink -target protein | efetch -format fasta > {}'.format(entrez_query,target_fasta_file_path)
            process = subprocess.Popen(cmd, shell=True)
            returncode = process.wait(timeout=settings.SUBPROCESS_TIME_LIMIT)

            with transaction.atomic():
                entrez_search.fasta_file_name = target_fasta_file_path
                entrez_search.save()
            return returncode

    #catch either subprocess
    except subprocess.TimeoutExpired:
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
    except subprocess.SubprocessError as e:
        raise Exception("[-] Couldnt download protein fasta files for esearch {} on protein database and exception: {}".format(search_id,e))

#TODO documentation
def update_entrezsearch_with_download_task_result(search_id: int,task_id: int) -> int:
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

        returncode = process.wait(timeout=settings.SUBPROCESS_TIME_LIMIT)
        return returncode

    except subprocess.TimeoutExpired as e:
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


def create_random_filename() -> str:
    #creates a rondom file name and returns it
    randomly_generated_filename = ''.join(random.choices(string.ascii_uppercase + string.digits, k = 10))
    esearch_output_filepath = 'media/esearch_output/result_dataframe_' + randomly_generated_filename
    while os.path.isfile(esearch_output_filepath) == True:
        randomly_generated_filename = ''.join(random.choices(string.ascii_uppercase + string.digits, k=10))
        esearch_output_filepath = 'media/esearch_output/result_dataframe_' + randomly_generated_filename

    return randomly_generated_filename

def get_entrezsearch_object_with_entrezsearch_id(search_id: int) -> int:
    #gets an entrez search database row based on a search_id and returns it
    try:
        entrez_search = EntrezSearch.objects.get(id=search_id)
        return entrez_search
    except Exception as e:
        raise Exception("[-] There is no entrez_search object with id: {} exception: {}".format(search_id, e))


def delete_esearch_by_id(search_id: int):
    #deletes en entrazsearch assoiciated files and taskresult entry based on a search_id
    #returns 0 if it worked or 1 if it did not
    try:
        with transaction.atomic():

            entrez_search = EntrezSearch.objects.get(id=search_id)
            task_db_id = EntrezSearch.objects.values('search_task_result_id').filter(id=search_id)
            task_db_id = task_db_id[0]['search_task_result_id']
            task_db = TaskResult.objects.get(id=task_db_id)
            organism_db_file_name = "media/esearch_output/"+str(search_id)

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
        raise IntegrityError("ERROR during deletion of entrez_search with id : {} with exception: {}".format(search_id,e))

def save_entrez_search_model(database:str,entrez_query:str,file_name:str, task_result_id:int, user_id:int) -> EntrezSearch:
    #saves users entrezsearch input into a datbase with extra information
    #returns the database row
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