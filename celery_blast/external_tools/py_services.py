import os
from .phy2html import create_html_tree

def create_html_output_for_newicktree(path_to_fasttree_output, project_id, query_accession):
    try:
        if os.path.isfile(path_to_fasttree_output):
            path_to_html_phylogeny = path_to_fasttree_output.split(".nwk")[0] + '.html'
            if os.path.isdir('static/images/result_images/' + str(project_id) + '/'+query_accession):
                path_to_static_html_phylogeny = 'static/images/result_images/' + str(project_id) + '/'+query_accession+'/target_sequences.html'
            else:
                os.mkdir('static/images/result_images/' + str(project_id) + '/'+query_accession)
                path_to_static_html_phylogeny = 'static/images/result_images/' + str(project_id) + '/'+query_accession+'/target_sequences.html'

            html_table_list = create_html_tree(path_to_fasttree_output,path_to_html_phylogeny)
            html_table_list = create_html_tree(path_to_fasttree_output, path_to_static_html_phylogeny)
            return 0
        else:
            return 1
    except Exception as e:
        raise Exception("[-] couldnt create html table for the newick file: {}".format(path_to_fasttree_output))

def get_list_of_query_sequence_folder(project_id):
    path_to_project = 'media/blast_projects/' + str(project_id)
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
        path_to_folder = 'media/blast_projects/' + str(project_id)
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