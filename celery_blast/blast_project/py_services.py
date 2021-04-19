from .models import BlastDatabase
from os.path import isdir
from os import mkdir, listdir
from shutil import rmtree
from django.db import IntegrityError, transaction

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


''' delete_refseqgenome_and_associated_directories_by_id

    deletes the blastdatabase entry by id
    removes the blastdatabase directory
    
    raises an integrity error if exception occurs
    
    :param database_id
        :type int
'''
def delete_blastdb_and_associated_directories_by_id(database_id):
    try:
        with transaction.atomic():
            blastdatabase = BlastDatabase.objects.get(id=database_id)
            if isdir('media/databases/' + str(database_id)):
                rmtree('media/databases/' + str(database_id))
            blastdatabase.delete()
    except Exception as e:
        raise IntegrityError("couldnt delete refseqgenome entry : {}".format(e))

''' create_blastdatabase_directory
    
    creates an directory with the name database_id in the /media/databases folder 
    
    :param database_id
        :type int
'''
def create_blastdatabase_directory(database_id):
    try:
        mkdir('media/databases/' + str(database_id))
        return 'media/databases/' + str(database_id)
    except Exception as e:
        raise IntegrityError(
            'something went wrong during database directory creation: {}'.format(e))


def upload_file(project_file, destination):
    try:
        with open(destination, 'wb+') as dest:
            for chunk in project_file.chunks():
                dest.write(chunk)
    except Exception as e:
        raise IntegrityError(
            'exception during file upload of : {} : exception : {}'.format(project_file.name,e))

