import os

import pandas as pd

from .models import BlastProject, BlastDatabase, AssemblyLevels, BlastSettings
from external_tools.models import ExternalTools
from .py_services import create_blastdatabase_directory,concatenate_genome_fasta_files_in_db_dir, upload_file, write_pandas_table_for_uploaded_genomes, write_pandas_table_for_one_genome_file,write_pandas_table_for_multiple_uploaded_files, pyb
from django_celery_results.models import TaskResult
from django.db import IntegrityError, transaction
from pandas import read_csv, Series

'''py_django_db_services

This script provides functions that should serve as a layer between the database and other services.

'''

#TODO documentation
def create_external_tools_after_snakemake_workflow_finishes(project_id):
    try:
        with transaction.atomic():
            external_tool = ExternalTools.objects.create_external_tools(project_id)
        return external_tool
    except IntegrityError as e:
        raise IntegrityError("[-] couldnt create external tools model with exception : {}".format(e))

#following functions are utilized in the dashboard_view
''' get_users_blast_projects
    
    returns a query-set of blast_projects from the logged in user
    
    :param userid (request.user.id)
        :type int
    :returns django.db.models.query.QuerySet of BlastProjects
'''
def get_users_blast_projects(userid):
    return BlastProject.objects.get_blast_projects_by_userid(userid)

''' get_all_blast_databases

    :returns all blastdatabases as a query-set
'''
def get_all_blast_databases():
    return BlastDatabase.objects.all()


#TODO documentation
def get_project_by_id(project_id):
    return BlastProject.objects.get(id=project_id)

''' create_project_from_form

Interacts with the BlastProject model. Calls the BlastProjectManager by using the "objects" field of the BlastProject
model. Executes the "create_blast_project" function of the BlastProjectManager, which saves the BlastProject into
the postgresql database. The input for this function is maintained by the py_project_creation.py script. 

    :param valid_project_form
        :type django form object
    :param user
        :type django user object (model)
    :param fw_settings
        :type django model
    :param bw_settings
        :type django model
    :param query_sequence_filename
        :type string
    :return blast_project
        :type django model (BlastProject)
'''
def create_project_from_form(valid_project_form,user,fw_settings,bw_settings,query_sequence_filename):
    try:
        blast_project = BlastProject.objects.create_blast_project(
            project_title=valid_project_form.cleaned_data['project_title'],
            search_strategy='blastp',
            project_query_sequences=query_sequence_filename,
            project_user=user,
            project_forward_settings=fw_settings,
            project_backward_settings=bw_settings,
            project_forward_database=valid_project_form.cleaned_data['project_forward_database'],
            project_backward_database=valid_project_form.cleaned_data['project_backward_database'],
            species_name_for_backward_blast=valid_project_form.cleaned_data['species_name_for_backward_blast']
        )
        return blast_project
    except Exception as e:
        raise IntegrityError('couldnt create blast project with exception : {}'.format(e))
#TODO documentation
def create_blast_settings_from_form(fwOrBw,valid_settings_form):
    try:
        if(fwOrBw == 'fw'):
            blast_settings = BlastSettings.objects.create(e_value=valid_settings_form.cleaned_data['fw_e_value'],
                                            word_size=valid_settings_form.cleaned_data['fw_word_size'],
                                            num_alignments=valid_settings_form.cleaned_data['fw_num_alignments'],
                                            max_target_seqs=valid_settings_form.cleaned_data['fw_max_target_seqs'],
                                            num_threads=valid_settings_form.cleaned_data['fw_num_threads'],
                                            max_hsps=valid_settings_form.cleaned_data['fw_max_hsps'])
        elif(fwOrBw == 'bw'):
            blast_settings = BlastSettings.objects.create(e_value=valid_settings_form.cleaned_data['bw_e_value'],
                                            word_size=valid_settings_form.cleaned_data['bw_word_size'],
                                            num_alignments=valid_settings_form.cleaned_data['bw_num_alignments'],
                                            max_target_seqs=valid_settings_form.cleaned_data['bw_max_target_seqs'],
                                            num_threads=valid_settings_form.cleaned_data['bw_num_threads'],
                                            max_hsps=valid_settings_form.cleaned_data['bw_max_hsps'])
        else:
            raise IntegrityError('fwOrBw is wrong ...')

        return blast_settings
    except Exception as e:
        raise IntegrityError('something went wrong during creation of blastsettings with Exception : {}'.format(e))

''' update_blast_project_with_task_result_model

everytime a celery reciprocal blast task gets executed, particularly the function execute_reciprocal_blast, the TaskResult object gets 
connected to the BlastProject model. This function is responsible for this linkage, thus it loads the relevant BlastProject and TaskResult by their 
corresponding IDs and saves them back to the database.

    :param project_id
        :type int
    :param task_id
        :type int
'''
def update_blast_project_with_task_result_model(project_id,task_id):
    try:
        blast_project = BlastProject.objects.get(id=project_id)
        taskresult = TaskResult.objects.get(task_id=task_id)
        blast_project.project_execution_snakemake_task = taskresult
        blast_project.save()
    except Exception as e:
        raise IntegrityError('problem during updating of blastproject model with task result instance exception : {}'.format(e))

'''update_blast_database_with_task_result_model

see above description. Similar to update_blast_project_with_task_result_model just for the celery database creation task, particularly the functions
execute_makeblastdb_with_uploaded_genomes and download_blast_databases_based_on_summary_file.

    :param database_id
        :type int
    :param task_id
        :type int

'''
def update_blast_database_with_task_result_model(database_id,task_id):
    try:
        blastdb = BlastDatabase.objects.get(id=database_id)
        taskresult = TaskResult.objects.get(task_id=task_id)
        blastdb.database_download_and_format_task = taskresult
        blastdb.save()
    except Exception as e:
        raise IntegrityError('problem during updating of database model with task result instance exception : {}'.format(e))

'''get_database_by_id
    
wrapper for the model function objects.get, return a BlastDatabase object with the provided key.

    :param database_id
        :type int

'''
def get_database_by_id(database_id):
    try:
        blastdb=BlastDatabase.objects.get(id=database_id)
        return blastdb
    except Exception as e:
        raise IntegrityError('there is no database with this {} id : {}'.format(database_id,e))

'''get_all_succeeded_databases

This functions uses the BlastDatabase model manager and executes the custom function get_databases_with_succeeded_tasks.
As the name says, it returns all BlastDatabase models with the 'SUCCESS' entry in their corresponding TaskResult objects. 

'''
def get_all_succeeded_databases():
    return BlastDatabase.objects.get_databases_with_succeeded_tasks()

#TODO documentation
def create_and_save_refseq_database_model(database_name,database_description,assembly_levels,assembly_entries,attached_taxonomic_file=None):
    try:

        #create model refseq genome objects (s. models.py file)
        #path_to_database_file = 'media/' + 'databases/' + 'refseq_databases/' + database_description.replace(' ','_').upper() + '.database.faa'
        if attached_taxonomic_file != None:
            blast_database = BlastDatabase.objects.create(
                database_name=database_name,
                database_description=database_description,
                assembly_entries=assembly_entries,
                attached_taxonomic_node_file=attached_taxonomic_file)
        else:
            blast_database = BlastDatabase.objects.create(
                database_name=database_name,
                database_description=database_description,
                assembly_entries=assembly_entries)


        #get all associated assembly levels (max. 4)
        assembly_levels_models = AssemblyLevels.objects.filter(assembly_level__in=assembly_levels)

        for assembly_level in assembly_levels_models:
            blast_database.assembly_levels.add(assembly_level)

        blast_database.path_to_database_file = 'media/databases/' + str(blast_database.id)
        blast_database.save()
        return blast_database
    except Exception as e:
        raise IntegrityError('couldnt save refseq genome model into database with exception : {}'.format(e))

#TODO implementation documentation
def save_uploaded_genomes_into_database(database_title,database_description,genome_file,assembly_entries,
                                     assembly_level,taxonomic_node,user_email,assembly_accession=None,
                                     organism_name=None,taxmap_file=None,organism_file=None,
                                     assembly_accession_file=None,assembly_level_file=None):

    try:
        blast_database = BlastDatabase.objects.create(database_name=database_title,
                                                      database_description=database_description,
                                                      assembly_entries=assembly_entries,
                                                      uploaded_files=True)
        #blast_database.path_to_database_file = 'media/databases/' + str(blast_database.id)
        assembly_levels = AssemblyLevels.objects.filter(assembly_level__contains=assembly_level)
        for assembly_lvl in assembly_levels:
            blast_database.assembly_levels.add(assembly_lvl)

        if taxmap_file != None:
            create_database_directory_and_upload_files(blast_database,
                                                       genome_file,
                                                       taxmap_file=taxmap_file)

            write_pandas_table_for_uploaded_genomes(blast_database,
                                                    assembly_accession_file,
                                                    assembly_level_file,
                                                    organism_file,
                                                    user_email)
        elif taxmap_file == None:
            create_database_directory_and_upload_files(blast_database,
                                                       genome_file)

            write_pandas_table_for_one_genome_file(blast_database,
                                                   organism_name,
                                                   assembly_level,
                                                   taxonomic_node,
                                                   assembly_accession)

        blast_database.path_to_database_file = "media/databases/"+str(blast_database.id)
        blast_database.save()
        return blast_database
    except Exception as e:
        raise IntegrityError('couldnt save uploaded genome model into database with exception : {}'.format(e))

#TODO documentation
def save_uploaded_multiple_file_genomes_into_database(cleaned_data_multiple_files, amount_of_entries, user_email):
    #1st save db model
    #2nd upload files
    #3rd write accession table
    #4th concatenate files
    #5th makeblastdb cmd

    database_title = cleaned_data_multiple_files['database_title']
    database_description = cleaned_data_multiple_files['database_description']

    try:
        blast_database = BlastDatabase.objects.create(database_name=database_title,
                                                      database_description=database_description,
                                                      assembly_entries=amount_of_entries,
                                                      uploaded_files=True)
        #blast_database.path_to_database_file = 'media/databases/' + str(blast_database.id)
        assembly_levels = AssemblyLevels.objects.filter(assembly_level__contains="Chromosome")
        for assembly_lvl in assembly_levels:
            blast_database.assembly_levels.add(assembly_lvl)

        #creation of database directory:
        create_blastdatabase_directory(database_id=blast_database.id)
        path_to_database = 'media/databases/' + str(blast_database.id) + '/'

        genomes_to_organism_and_taxid_dict = {}
        with open(path_to_database+'acc_taxmap.table','w') as taxmap_file:
            for index in range(int(amount_of_entries)):
                file = 'genome_file_field_{}'.format(index)
                organism = 'organism_name_{}'.format(index)

                organism = cleaned_data_multiple_files[organism]
                taxid = pyb.get_species_taxid_by_name(user_email,organism)
                file = cleaned_data_multiple_files[file]

                upload_file(file,path_to_database + file.name)

                with open(path_to_database + file.name,'r') as current_genome_file:
                    for line in current_genome_file.readlines():
                        if line.startswith(">"):
                            acc = line.split(" ")[0].split(">")[1]
                            taxmap_file.write(acc+"\t"+str(taxid)+"\n")

                genomes_to_organism_and_taxid_dict[file.name] = [organism,taxid]

        write_pandas_table_for_multiple_uploaded_files(blast_database,genomes_to_organism_and_taxid_dict)

        blast_database.path_to_database_file = "media/databases/"+str(blast_database.id)

        genome_files = list(genomes_to_organism_and_taxid_dict.keys())
        concatenate_genome_fasta_files_in_db_dir(path_to_database, database_title, genome_files)
        blast_database.save()
        return blast_database
    except Exception as e:
        raise IntegrityError('couldn save multiple files genome model into database with exception : {}'.format(e))

#TODO documentation
def create_database_directory_and_upload_files(blast_database,genome_file,taxmap_file=None):
    try:

        path_to_database = 'media/databases/' + str(blast_database.id) + '/'

        create_blastdatabase_directory(database_id=blast_database.id)

        if taxmap_file != None:
            upload_file(taxmap_file, path_to_database+'acc_taxmap.table')
        upload_file(genome_file, path_to_database + blast_database.database_name.replace(' ','_').upper()+'.database')
    except Exception as e:
        raise IntegrityError('couldnt upload genome or taxmap file into database directory with exception : {}'.format(e))

#TODO documentation
def check_if_taxid_is_in_database(database_id, taxonomic_node):
    path_to_database = 'media/databases/' + str(database_id) + '/'
    database = get_database_by_id(database_id)
    pandas_table_file = path_to_database + database.get_pandas_table_name()
    df = read_csv(pandas_table_file, header=0, index_col=0)
    boolean = int(taxonomic_node) in list(df['taxid'])
    return boolean

'''check_if_sequences_are_in_databases
This function is used in forms.py for the reciprocal BLAST project creation. 
It validates 
'''
def check_if_sequences_are_in_database(database_id, sequences):
    path_to_database = 'media/databases/' + str(database_id) + '/'

    taxmap_files = os.listdir(path_to_database)
    taxmap_files = [file for file in taxmap_files if file.endswith('table')]

    pandas_table_file = path_to_database + taxmap_files[0]

    #TODO
    #what to do if database is too big?
    df = read_csv(pandas_table_file, header=None, sep="\t")
    df.columns = ['AccessionId', 'TaxId']
    df = df['AccessionId'].map(lambda acc: acc.split(".")[0])

    to_compare = Series(sequences)
    to_compare = to_compare[~to_compare.isin(df)]
    print(df)
    if len(to_compare) != 0:
        return list(to_compare)
    else:
        return True
