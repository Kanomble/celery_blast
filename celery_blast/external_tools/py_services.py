import os

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