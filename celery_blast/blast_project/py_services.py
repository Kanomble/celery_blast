from .models import BlastProject
from refseq_transactions.models import BlastDatabase
from os.path import isdir, isfile
from os import mkdir, listdir
from shutil import rmtree
from django.db import IntegrityError, transaction
from blast_project import py_biopython as pyb

'''check_if_taxdb_exists
    
    This function is used before project creation to check if the taxonomy database exists in the database directory.
    
    :returns True or False
        :type boolean
'''
def check_if_taxdb_exists()->bool:
    if isfile('media/databases/taxdb.btd') and isfile('media/databases/taxdb.bti'):
        return True
    else:
        return False

def check_if_file_exists(file_path)->bool:
    return isfile(file_path)


'''list_taxonomic_files

utilization in create_taxonomic_file_view and refseqdatabaseform 
returns a list of all files and their corresponding total line length in the media/taxonomic_node_files folder that end with .taxids

'''
def list_taxonomic_files():
    try:
        files_in_taxonomic_node_files = listdir('media/taxonomic_node_files/')
        files = []
        length = []
        for file in files_in_taxonomic_node_files:
            lines = 0
            with open('media/taxonomic_node_files/'+file) as f:
                for line in f:
                    lines = lines + 1
            if file.endswith('.taxids'):
                files.append(file)
                length.append(lines)
        #[file for file in files_in_taxonomic_node_files if file.endswith('.taxids')]
        return files, length
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
def create_blastdatabase_directory(database_id,database_filepath='media/databases/'):
    try:
        mkdir(database_filepath + str(database_id))
        return database_filepath + str(database_id)
    except Exception as e:
        raise IntegrityError(
            'something went wrong during database directory creation: {}'.format(e))

'''upload_file
    simple function for uploading a file to a specific server side location, defined by destination
'''
def upload_file(project_file, destination:str):
    try:
        with open(destination, 'wb+') as dest:
            for chunk in project_file.chunks():
                dest.write(chunk)
    except Exception as e:
        raise IntegrityError(
            'exception during file upload of : {} : exception : {}'.format(project_file.name,e))


'''get_html_results
    This function is used in blast_project.views.load_reciprocal_result_html_table_view, it returns
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
#loads the reciprocal results table that is written with one of the last rules in the snakefiles
def get_html_results(project_id:int,filename:str,html_result_path="media/blast_projects/")->list:
    try:
        with open(html_result_path+str(project_id)+"/"+filename) as res:
            data = res.readlines()
        return data
    except Exception as e:
        raise FileNotFoundError("[-] ERROR: Couldn't read file {} with Exception: {}".format(filename,e))

'''html_table_exists
    Checks if the html table exists.
'''
def html_table_exists(project_id,filename,html_result_path="media/blast_projects/"):
    if(isfile(html_result_path+str(project_id)+"/"+filename)):
        return True
    else:
        return False

#TODO documentation
def concatenate_genome_fasta_files_in_db_dir(path_to_database,database_title,genome_files):
    try:
        database_name = database_title.replace(' ', '_').upper() + '.database'
        with open(path_to_database+database_name,'w') as dbfile:
            for file in genome_files:
                with open(path_to_database+file, 'r') as gfile:
                    for line in gfile.readlines():
                        dbfile.write(line)
    except Exception as e:
        raise IntegrityError('couldnt concatenate database files : {}'.format(e))
