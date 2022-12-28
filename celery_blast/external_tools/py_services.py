import os
from django.conf import settings
import pandas as pd

def create_html_output_for_newicktree(path_to_fasttree_output, project_id, query_accession):
    try:
        if os.path.isfile(path_to_fasttree_output):
            path_to_html_phylogeny = path_to_fasttree_output.split(".nwk")[0] + '.html'
            if os.path.isdir('static/images/result_images/' + str(project_id) + '/'+query_accession):
                path_to_static_html_phylogeny = 'static/images/result_images/' + str(project_id) + '/'+query_accession+'/target_sequences.html'
            else:
                os.mkdir('static/images/result_images/' + str(project_id) + '/'+query_accession)
                path_to_static_html_phylogeny = 'static/images/result_images/' + str(project_id) + '/'+query_accession+'/target_sequences.html'

            return 0
        else:
            return 1
    except Exception as e:
        raise Exception("[-] couldnt create html table for the newick file: {}".format(path_to_fasttree_output))

def get_list_of_query_sequence_folder(project_id):
    path_to_project = settings.BLAST_PROJECT_DIR + str(project_id)
    try:
        qseqids = [pathname for pathname in os.listdir(path_to_project) if os.path.isdir(os.path.join(path_to_project,pathname))]
        if '.snakemake' in qseqids:
            qseqids.remove('.snakemake')
        return qseqids
    except OSError as e:
        raise Exception("[-] couldnt perform query id folder listing with exception : {}".format(e))

def get_msa_files_from_folder_list(project_id,folders):
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
            else: #not enough target sequences
                return 2
        else: #FileNotFound
            return 1
    except Exception as e:
        raise Exception("[-] error during checking amount of sequence targets for msa task with exception: {}".format(e))

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
def delete_cdd_search_output(query_sequence:str, project_id:int)->int:
    try:
        path_to_project_dir = settings.BLAST_PROJECT_DIR + str(project_id) + '/' + query_sequence + '/'
        cdd_table = path_to_project_dir + 'cdd_domains.tsf'
        bokeh_plot = path_to_project_dir + 'pca_bokeh_domain_plot.html'
        if os.path.isfile(cdd_table):
            os.remove(cdd_table)
        if os.path.isfile(bokeh_plot):
            os.remove(bokeh_plot)
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
def check_if_cdd_search_can_get_executed(query_sequence:str, project_id:int)->int:
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

