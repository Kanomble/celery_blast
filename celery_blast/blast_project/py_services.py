from .models import BlastDatabase, BlastProject
from os.path import isdir, isfile
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
        raise IntegrityError("couldnt delete blast database entry : {}".format(e))

#TODO documentation
def delete_project_and_associated_directories_by_id(project_id):
    try:
        with transaction.atomic():
            project = BlastProject.objects.get(id=project_id)
            if isdir('media/blast_projects/' + str(project_id)):
                rmtree('media/blast_projects/' + str(project_id))
            if isdir('static/images/result_images/'+str(project_id)):
                rmtree('static/images/result_images/'+str(project_id))
            project.delete()
    except Exception as e:
        raise IntegrityError("couldnt delete blast project entry : {}".format(e))


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

#TODO documentation
def upload_file(project_file, destination):
    try:
        with open(destination, 'wb+') as dest:
            for chunk in project_file.chunks():
                dest.write(chunk)
    except Exception as e:
        raise IntegrityError(
            'exception during file upload of : {} : exception : {}'.format(project_file.name,e))

#TODO documentation
#loads the reciprocal results table that is written with one of the last rules in the snakefiles
def get_html_results(project_id,filename):
    try:
        with open("media/blast_projects/"+str(project_id)+"/"+filename) as res:
            data = res.readlines()
        return data
    except Exception as e:
        raise FileNotFoundError("Couldn't read file {} with Exception: {}".format(e))

#TODO documentation
def html_table_exists(project_id,filename):
    if(isfile("media/blast_projects/"+str(project_id)+"/"+filename)):
        return True
    else:
        return False