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
