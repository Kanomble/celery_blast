from .models import BlastDatabase, BlastProject
from os.path import isdir, isfile
from os import mkdir, listdir
from shutil import rmtree
from django.db import IntegrityError, transaction
from blast_project import py_biopython as pyb
import pandas as pd

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

#TODO documentation
def write_pandas_table_for_one_genome_file(blast_database,organism_name,assembly_level,taxonomic_node,assembly_accession):
    try:

        path_to_database = 'media/databases/' + str(blast_database.id) + '/'

        if assembly_accession == None:
            assembly_accession = 'not provided'
        if organism_name == None:
            organism_name = 'not provided'

        pandas_table_file = open(path_to_database + blast_database.get_pandas_table_name(), 'w')
        pandas_table_file.write(',assembly_accession,organism_name,taxid,species_taxid,assembly_level,ftp_path\n')
        pandas_table_file.write('0,{},{},{},{},{},{}\n'.format(
            assembly_accession,organism_name,taxonomic_node,taxonomic_node,assembly_level,'uploaded genome'
        ))
        pandas_table_file.close()
    except Exception as e:
        raise IntegrityError("Couldnt write pandas dataframe for your uploaded genome file, with exception : {}".format(e))


#TODO documentation
def write_pandas_table_for_multiple_uploaded_files(blast_database, genomes_to_organism_and_taxid_dict):
    try:
        path_to_database = 'media/databases/' + str(blast_database.id)+'/'
        with open(path_to_database + blast_database.get_pandas_table_name(), 'w') as  pandas_table_file:
            pandas_table_file.write(',assembly_accession,organism_name,taxid,species_taxid,assembly_level,ftp_path\n')
            for line_index, key in enumerate(list(genomes_to_organism_and_taxid_dict.keys())):
                pandas_table_file.write(str(line_index) + ',')
                pandas_table_file.write(key + ',')
                pandas_table_file.write(genomes_to_organism_and_taxid_dict[key][0] + ',')
                pandas_table_file.write(genomes_to_organism_and_taxid_dict[key][1] + ',' + genomes_to_organism_and_taxid_dict[key][1] + ',')
                pandas_table_file.write("Chromosome"+ ',uploaded genome\n')
        return 0
    except Exception as e:
        raise IntegrityError('couldnt write database table : {}'.format(e))


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



#TODO documentation
def write_pandas_table_for_uploaded_genomes(blast_database,
                                            assembly_accessions_file,
                                            assembly_levels_file,
                                            organisms_file,
                                            user_email):
    try:
        path_to_database = 'media/databases/' + str(blast_database.id)+'/'

        taxonomic_nodes,organisms = get_list_of_taxonomic_nodes_based_on_organisms_file(organisms_file,user_email)

        assembly_accessions = []
        assembly_levels = []

        if assembly_accessions_file != None:
            for line in assembly_accessions_file:
                line = line.decode().rstrip()
                if line != '' and line != '\r':
                    assembly_accessions.append(line)

            if assembly_levels_file != None:
                for line in assembly_levels_file:
                    line = line.decode().rstrip()
                    if line != '' and line != '\r':
                        assembly_levels.append(line)
            else:
                for line in range(len(taxonomic_nodes)):
                    assembly_levels.append('not provided')

        elif assembly_accessions_file == None:
            for line in range(len(taxonomic_nodes)):
                assembly_accessions.append('not provided')
            if assembly_levels_file != None:
                for line in assembly_levels_file:
                    line = line.decode().rstrip()
                    if line != '' and line != '\r':
                        assembly_levels.append(line)
            else:
                for line in range(len(taxonomic_nodes)):
                    assembly_levels.append('not provided')

        pandas_table_file = open(path_to_database+blast_database.get_pandas_table_name(),'w')
        pandas_table_file.write(',assembly_accession,organism_name,taxid,species_taxid,assembly_level,ftp_path\n')
        for line_index in range(len(taxonomic_nodes)):
            pandas_table_file.write(str(line_index)+',')
            pandas_table_file.write(assembly_accessions[line_index]+',')
            pandas_table_file.write(organisms[line_index]+',')
            pandas_table_file.write(taxonomic_nodes[line_index]+','+taxonomic_nodes[line_index]+',')
            pandas_table_file.write(assembly_levels[line_index]+',uploaded genome\n')
        pandas_table_file.close()

    except Exception as e:
        raise IntegrityError('couldnt write database table : {}'.format(e))

#TODO documentation
def get_list_of_taxonomic_nodes_based_on_organisms_file(organisms_file,user_email):
    try:
        taxids = []
        organisms = []
        for line in organisms_file:
            line = line.decode().rstrip()

            if line != '':
                taxid = pyb.get_species_taxid_by_name(user_email,line)
                taxids.append(taxid)
                organisms.append(line)
        return taxids, organisms
    except Exception as e:
        raise IntegrityError('couldnt translate organism names into taxonomic nodes with exception : {}'.format(e))
