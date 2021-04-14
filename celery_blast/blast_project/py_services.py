from os import listdir

'''list_taxonomic_files

utilization in create_taxonomic_file_view
returns a list of all files in the media/taxonomic_node_files folder that end with .taxids

'''
def list_taxonomic_files():
    try:
        files_in_taxonomic_node_files = listdir('media/taxonomic_node_files/')
        return [file for file in files_in_taxonomic_node_files if file.endswith('.taxids')]
    except Exception as e:
        raise Exception('exception ocurred in blast_project/py_services.list_taxonomic_files : {}'.format(e))