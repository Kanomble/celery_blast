import os
import pandas as pd
from django.conf import settings
import json

'''slice_cdd_domain_corrected_fasta_file
    
    This function is part of the selection constrained CDD phylogenetic inference task pipeline.
    
    :param path_to_fasta_file
        :type str
    :param accessions
        :type list[str]
    
    :returns selection_sliced_domain_fasta_file
        :type str
'''
def slice_cdd_domain_corrected_fasta_file(path_to_query_subdir:str, path_to_fasta_file:str, accessions:list)->str:
    try:
        if os.path.isdir(path_to_query_subdir) == False:
            raise NotADirectoryError("{} is not a directory".format(path_to_query_subdir))
        elif os.path.isfile(path_to_fasta_file) == False:
            raise FileNotFoundError("{} is not a file".format(path_to_fasta_file))
        else:
            targets = 0
            output_file_path = path_to_query_subdir + '/selection_sliced_domain_fasta.faa'
            with open(path_to_fasta_file, 'r') as input_file:
                with open(output_file_path, 'w') as output_file:
                    for line in input_file.readlines():
                        if line.startswith(">"):
                            header = line.split(">")[1].rstrip()
                            if header in accessions:
                                output_file.write(">"+header+"\n")
                                switch = True
                                targets += 1
                            else:
                                switch = False
                        else:
                            if switch == True:
                                output_file.write(line.rstrip()+"\n")
            if targets < len(accessions):
                raise Exception("some accessions are missing in the target file ... ")
            return output_file_path
    except Exception as e:
        raise Exception("[-] ERROR couldnt slice domain corrected fasta file with Exception: {}".format(e))

'''read_query_sequence_tbh_table
    
    This function is executed after pressing the calculate synteny button within the synteny dashboard.
    It loads the rbh table for the specified query sequence id as a json dictionary.
    The json dictionary will be rendered by the DataTables library.
    
    :param project_id
        :type int
    :param qseqid 
        :type str
    
    :returns json
        :type json
'''
def read_query_sequence_rbh_table(project_id:int, qseqid:str):
    try:
        target_table_path = settings.BLAST_PROJECT_DIR + str(project_id) + '/' + qseqid + '/rbh_table.tsf'
        if os.path.isfile(target_table_path) == False:
            raise Exception("{} does not exist!".format(target_table_path))
        else:
            table = pd.read_csv(target_table_path, sep="\t", header=0)
            json_records = table.reset_index().to_json(orient='records')
            json_data = json.loads(json_records)
            return json_data
    except Exception as e:
        raise Exception("[-] ERROR during reading of query sequence rbh table:"
                        " {} - {} with exception: {}".format(project_id, qseqid, e))

def get_list_of_query_sequence_folder(project_id):
    path_to_project = settings.BLAST_PROJECT_DIR + str(project_id)
    try:
        qseqids = [pathname for pathname in os.listdir(path_to_project) if
                   os.path.isdir(os.path.join(path_to_project, pathname))]
        if '.snakemake' in qseqids:
            qseqids.remove('.snakemake')
        return qseqids
    except OSError as e:
        raise Exception("[-] couldnt perform query id folder listing with exception : {}".format(e))


def get_msa_files_from_folder_list(project_id, folders):
    try:
        qseqid_to_msa = {}
        path_to_folder = settings.BLAST_PROJECT_DIR + str(project_id)
        for folder in folders:
            qseqid_to_msa[folder] = False
            target_path = path_to_folder + '/' + folder
            files = os.listdir(target_path)
            if "target_sequences.msa" in files:
                qseqid_to_msa[folder] = True
        return qseqid_to_msa
    except Exception as e:
        raise Exception("[-] couldnt identify if qseq folder contains msa file or not, exception : {}".format(e))


'''
check if there are > 2 target sequences available (for msa task).
    
    :returns 0 if task can proceed, 1 if file was not found and 2 if there are not enough sequences
'''


def check_if_target_sequences_are_available(path_to_query_file: str) -> int:
    try:
        if os.path.isfile(path_to_query_file):
            count = 1
            with open(path_to_query_file, 'r') as query_file:
                for line in query_file.readlines():
                    if line.startswith(">"):
                        count += 1
            if count >= 2:
                return 0
            else:  # not enough target sequences
                return 2
        else:  # FileNotFound
            return 1
    except Exception as e:
        raise Exception(
            "[-] error during checking amount of sequence targets for msa task with exception: {}".format(e))


'''
check if multiple sequence file is available (for phylo task)
'''


def check_if_msa_file_is_available(path_to_msa_file: str) -> int:
    try:
        if os.path.isfile(path_to_msa_file):
            return 0
        else:
            return 1
    except Exception as e:
        raise Exception("[-] error during checking if msa file exists with exception: {}".format(e))


'''delete_cdd_search_output

    This function is executed in delete_cdd_domain_search_view of this app. It 
    parses through the potential output of the cdd search task an deletes all associated files.
    
    :param query_sequence
        :type str
    :param project_id
        :type int
        
    :returns returncode
        :type int
'''


def delete_cdd_search_output(query_sequence: str, project_id: int) -> int:
    try:
        path_to_project_dir = settings.BLAST_PROJECT_DIR + str(project_id) + '/' + query_sequence + '/'
        cdd_table = path_to_project_dir + 'cdd_domains.tsf'
        bokeh_plot = path_to_project_dir + 'pca_bokeh_domain_plot.html'
        domain_corrected_faa_file = path_to_project_dir + 'domain_corrected_target_sequences.faa'
        msa_file = path_to_project_dir + 'domain_corrected_target_sequences.msa'
        tree_file = path_to_project_dir + 'domain_corrected_domain_corrected_target_sequences.nwk'
        filelist = [cdd_table, bokeh_plot, domain_corrected_faa_file, msa_file, tree_file]
        for file in filelist:
            if os.path.isfile(file):
                os.remove(file)
        return 0
    except Exception as e:
        raise Exception("[-] error during deletion of cdd search output with exception: {}".format(e))


'''check_if_cdd_search_can_get_executed
    
    This function ensures that the calculation of principal components based on the CDD search results is useful.
    It checks how many domains are present in the query sequence and how many target sequences have been inferred.
    If there are not more than two domain, the PCA after the CDD search wouldnt yield any significant results.
    
    :param query_sequence
        :type str
    :param project_id
        :type int
    
    :returns returncode - 1 = dont execute PCA plotting, 0 = execute PCA plotting
        :type int
'''


def check_if_cdd_search_can_get_executed(query_sequence: str, project_id: int) -> int:
    try:
        path_to_project_dir = settings.BLAST_PROJECT_DIR + str(project_id) + '/'
        path_to_query_domains = path_to_project_dir + 'query_domains.tsf'
        sequence_id_file = path_to_project_dir + query_sequence + '/' + 'target_sequence_ids.txt'
        if os.path.isfile(sequence_id_file) == False:
            raise FileNotFoundError("[-] there is no file for path: {}".format(sequence_id_file))
        if os.path.isfile(path_to_query_domains) == False:
            raise FileNotFoundError("[-] there is no file for path: {}".format(path_to_query_domains))

        with open(sequence_id_file, 'r') as seqidfile:
            number_of_queries = len(seqidfile.readlines())

        if number_of_queries <= 2:
            return 1
        header = "qseqid qlen sacc slen qstart qend sstart send qseq sseq bitscore evalue pident stitle".split(" ")
        cdd_queries = pd.read_table(path_to_query_domains, header=None)
        cdd_queries.columns = header
        cdd_queries['transformed_qseqid'] = cdd_queries.qseqid.apply(lambda x: x.split(".")[0])
        number_of_domains = len(cdd_queries[cdd_queries['transformed_qseqid'] == query_sequence])
        if number_of_domains <= 2:
            return 1
        return 0
    except Exception as e:
        raise Exception("[-] error checking if a PCA for the CDD search result "
                        "dataframe for: {} is practicable, with exception: {}".format(query_sequence, e))


'''get_html_results
    This function is used in the phylogenetic dashboard, it returns
    the plain HTML as a string. 

    :param project_id
        :type int
    :param filename - html result filename
        :type str
    :param html_result_path
        :type str

    :returns data - string representation of a pandas html table
        :type list[str]
'''


def get_html_results(project_id: int, filename: str, html_result_path=settings.BLAST_PROJECT_DIR) -> list:
    try:
        with open(html_result_path + str(project_id) + "/" + filename) as res:
            data = res.readlines()
        return data
    except Exception as e:
        raise FileNotFoundError("[-] ERROR: Couldn't read file {} with Exception: {}".format(filename, e))