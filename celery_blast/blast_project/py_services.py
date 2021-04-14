from os import listdir

'''list_taxonomic_files
utilization in create_taxonomic_file_view, returns a list of all files in the media/taxonomic_node_files folder
'''
def list_taxonomic_files():
    try:
        return listdir('media/taxonomic_node_files/')
    except Exception as e:
        raise Exception('exception ocurred in blast_project/py_services.list_taxonomic_files : {}'.format(e))