import os
from .models import BlastProject, BlastSettings
from refseq_transactions.models import BlastDatabase, AssemblyLevels
from external_tools.models import ExternalTools, QuerySequences
from blast_project import py_biopython as pyb
from .py_services import create_blastdatabase_directory,concatenate_genome_fasta_files_in_db_dir, upload_file
from django_celery_results.models import TaskResult
from django.db import IntegrityError, transaction
from pandas import read_csv, Series
from shutil import rmtree
from celery_blast.settings import BLAST_DATABASE_DIR

'''update_external_tool_with_cdd_search

    This script uses model based functions of the ExternalTools model to update the QuerySequence model with 
    the associated celery cdd search TaskResult model. This can then be used to track the status of the cdd search.
    
    :param project_id
        :type int
    :param query_sequence
        :type str
    :param task_id
        :type int

'''
def update_external_tool_with_cdd_search(project_id:int,query_sequence:str,task_id:int):
    try:
        external_tool = ExternalTools.objects.get_external_tools_based_on_project_id(project_id)
        external_tool.update_query_sequences_cdd_search_task(query_sequence,task_id)
    except Exception as e:
        raise Exception("[-] ERROR couldnt update query sequence model with cdd search task with exception: {}".format(e))


'''get_query_sequence_of_external_tools
    
    This function returns the query_sequence_model of the specified query_sequence. 
    The query_sequence is the identifier of a protein used in the blast projects.
    
    :param project_id
        :type int
    :param query_sequence
        :type str
    
    :return query_sequence_model
        :type django model
'''
def get_query_sequence_of_external_tools(project_id:int,query_sequence:str):
    try:
        external_tool = ExternalTools.objects.get_external_tools_based_on_project_id(project_id)
        query_sequence_model = QuerySequences.objects.filter(external_tool_for_query_sequence=external_tool,
                                      query_accession_id=query_sequence)
        return query_sequence_model
    except Exception as e:
        raise Exception("[-] ERROR couldnt get query sequence model instance from external tools model with exception: {}".format(e))

'''get_reciprocal_result_target_fasta_files

    This function returns a list of filepaths to the fasta files of RBHs.
    
    :param project_id
        :type int
    
    :return target_fasta_files
        :type list[str]
    :return target_queries
        :type list[str]
'''
def get_reciprocal_result_target_fasta_files_and_queries(project_id: int):
    try:

        blast_project = get_project_by_id(project_id)
        project_dir = blast_project.get_project_dir()
        project_dir += '/'

        queries = blast_project.get_list_of_query_sequences()

        target_fasta_files = []
        target_queries = []
        if os.path.isdir(project_dir):
            for query in queries:
                if os.path.isdir(project_dir + query):
                    if os.path.isfile(project_dir + query + '/target_sequences.faa'):
                        target_fasta_files.append(project_dir + query + '/target_sequences.faa')
                        target_queries.append(query)
                else:
                    raise NotADirectoryError("{} is not a directory".format(project_dir + query))
        else:
            raise NotADirectoryError("{} is not a directory".format(project_dir))

        return target_fasta_files, target_queries
    except Exception as e:
        raise Exception("[-] ERROR during get_reciprocal_result_target_fasta_files with exception: {}".format(e))

'''py_django_db_services

This script provides functions that should serve as a layer between the database and other services.

'''

#TODO documentation
def create_external_tools_after_snakemake_workflow_finishes(project_id:int)->ExternalTools:
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
def get_users_blast_projects(userid:int):
    return BlastProject.objects.get_blast_projects_by_userid(userid)

''' get_all_blast_databases

    :returns all blastdatabases as a query-set
'''
def get_all_blast_databases():
    return BlastDatabase.objects.all()

'''delete_failed_or_unknown_databases
    
    This function parses the directory of all active databases by using the function get_all_databases().
    If a directory without a corresponding database.id exists, it will get removed.
    
    :returns 0 by success
        :type int
    :returns 1 by failure
        :type int
    
'''
def delete_failed_or_unknown_databases():
    try:
        # check if there are database directories that do not reside in the postgres database
        databases = get_all_blast_databases()
        ids = [int(database.id) for database in databases]
        for database_id in os.listdir(BLAST_DATABASE_DIR):
            try:
                identifier = int(database_id)
                if identifier not in ids:
                    if os.path.isdir(BLAST_DATABASE_DIR + str(identifier)):
                        rmtree(BLAST_DATABASE_DIR + str(identifier) + '/')
            except:
                continue
        return 0
    except:
        return 1

'''get_project_by_id

    Returns the associated BlastProject model instance.

    :param project_id
        :type int
        
    :returns BlastProject
        :type blast_project.models.BlastProject
'''
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
def create_project_from_form(valid_project_form,user,fw_settings,bw_settings,query_sequence_filename,filepath='media/blast_projects/'):
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
            species_name_for_backward_blast=valid_project_form.cleaned_data['species_name_for_backward_blast'],
            filepath=filepath
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

'''update_blast_project_with_database_statistics_task_result_model
    
    see above description. Similar to update_blast_project_with_task_result_model function but for the database statistics task.
    
    :param project_id
        :type int
    :param task_id
        :type int
'''
def update_blast_project_with_database_statistics_task_result_model(project_id:int,task_id:int):
    try:
        blast_project = BlastProject.objects.get(id=project_id)
        taskresult = TaskResult.objects.get(task_id=task_id)
        blast_project.project_database_statistics_task = taskresult
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

        blast_database.path_to_database_file = BLAST_DATABASE_DIR + str(blast_database.id)
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

        blast_database.path_to_database_file = BLAST_DATABASE_DIR+str(blast_database.id)
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
    #5th save database into postgresql

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
        create_blastdatabase_directory(database_id=blast_database.id, database_filepath=BLAST_DATABASE_DIR)
        path_to_database = BLAST_DATABASE_DIR + str(blast_database.id) + '/'

        genomes_to_organism_and_taxid_dict = {}
        with open(path_to_database+'acc_taxmap.table','w') as taxmap_file:
            for index in range(int(amount_of_entries)):
                file = 'genome_file_field_{}'.format(index)
                organism = 'organism_name_{}'.format(index)

                organism = cleaned_data_multiple_files[organism]
                #take the first taxid
                taxid = pyb.get_species_taxid_by_name(user_email,organism)[0]
                file = cleaned_data_multiple_files[file]

                upload_file(file,path_to_database + file.name)

                with open(path_to_database + file.name,'r') as current_genome_file:
                    for line in current_genome_file.readlines():
                        if line.startswith(">"):
                            acc = line.split(" ")[0].split(">")[1]
                            taxmap_file.write(acc+"\t"+str(taxid)+"\n")

                genomes_to_organism_and_taxid_dict[file.name] = [organism,taxid]

        write_pandas_table_for_multiple_uploaded_files(blast_database,genomes_to_organism_and_taxid_dict)

        blast_database.path_to_database_file = BLAST_DATABASE_DIR+str(blast_database.id)

        genome_files = list(genomes_to_organism_and_taxid_dict.keys())
        concatenate_genome_fasta_files_in_db_dir(path_to_database, database_title, genome_files)
        blast_database.save()
        return blast_database
    except Exception as e:
        raise IntegrityError('couldn save multiple files genome model into database with exception : {}'.format(e))

#TODO documentation
def create_database_directory_and_upload_files(blast_database,genome_file,taxmap_file=None):
    try:

        path_to_database = BLAST_DATABASE_DIR + str(blast_database.id) + '/'
        create_blastdatabase_directory(database_id=blast_database.id,database_filepath=BLAST_DATABASE_DIR)

        if taxmap_file != None:
            upload_file(taxmap_file, path_to_database+'acc_taxmap.table')
        upload_file(genome_file, path_to_database + blast_database.database_name.replace(' ','_').upper()+'.database')
    except Exception as e:
        raise IntegrityError('couldnt upload genome or taxmap file into database directory with exception : {}'.format(e))

'''check_if_taxid_is_in_databaase
    
    Is used during ProjectCreationForm validation. Checks if the selected database
    contains entries for the specified taxonomic node. Usually there is just one
    taxonomic node but some species do possess more than one node (e.g. Bacillus subtilis). Therefore
    the function returns true if just one taxonomic node is in the backward database. 
    
    :param database_id
        :type int
    :param taxonomic_nodes - derived by get_species_taxid_by_name from .py_biopython
        :type list[int]
    
    :returns True OR not_in_database
        :type boolean OR list[int]
'''
def check_if_taxid_is_in_database(database_id:int, taxonomic_nodes:list):
    #path_to_database = 'media/databases/' + str(database_id) + '/'
    database = get_database_by_id(database_id)
    #path now available for testing
    path_to_database = database.path_to_database_file + '/'
    pandas_table_file = path_to_database + database.get_pandas_table_name()
    df = read_csv(pandas_table_file, header=0, index_col=0)
    not_in_database=[]
    #sometimes there is more than one taxonomic node for an organism
    #the counter resembles the number of corresponding organisms, therefor
    #at the end of validation the counter must be one int smaller than the length
    #of taxonomic_nodes
    counter = len(taxonomic_nodes)
    for node in taxonomic_nodes:
        #and int(node) not in list(df['species_taxid']) is not allowed due to database formatting
        if int(node) not in list(df['taxid']):
            not_in_database.append(node)
        else:
            counter -= 1
    if counter == len(taxonomic_nodes):
        return False
    else:
        return True

'''check_if_sequences_are_in_databases
    This function is used in forms.py for the reciprocal BLAST project creation. 
    The function checks if the query sequence resides in the backward BLAST database.
    
    :params database_id
        :type int
    
    :params sequences - Accession IDs of the query sequences
        :type list[str]
        
    :returns True if all sequences reside in the database OR list[str] with sequences that do not reside in the database
        :type boolean OR list
'''
def check_if_sequences_are_in_database(database_id:int, sequences:list):
    #path_to_database = 'media/databases/' + str(database_id) + '/'
    #correct path for testing
    database = get_database_by_id(database_id)
    path_to_database = database.path_to_database_file + '/'

    taxmap_files = os.listdir(path_to_database)
    taxmap_files = [file for file in taxmap_files if file.endswith('.table')]

    to_compare = Series(sequences)
    for index,taxfile in enumerate(taxmap_files):
        pandas_taxmap_table_file = path_to_database + taxmap_files[index]

        #what to do if database is too big? -> acc_taxmap files are always restricted to contain the accession ids of maximal 500 genome entries
        df = read_csv(pandas_taxmap_table_file, header=None, sep="\t")
        df.columns = ['AccessionId', 'TaxId']
        df = df['AccessionId'].map(lambda acc: acc.split(".")[0])
        to_compare = to_compare[~to_compare.isin(df)]
        #no need for parsing additional taxmap files if the sequences reside in the current one
        if len(to_compare) == 0:
            return True

    if len(to_compare) != 0:
        return list(to_compare)

    elif len(to_compare) == 0:
        return True


'''write_pandas_table_for_one_genome_file

    This function is used in blast_project.py_django_db_services.save_uploaded_genomes_into_database.
'''
def write_pandas_table_for_one_genome_file(blast_database, organism_name, assembly_level, taxonomic_node,
                                           assembly_accession):
    try:
        # 'media/databases/' + str(blast_database.id) +
        path_to_database = str(blast_database.path_to_database_file) + '/'

        if assembly_accession == None:
            assembly_accession = 'not provided'
        if organism_name == None:
            organism_name = 'not provided'

        with open(path_to_database + blast_database.get_pandas_table_name(), 'w') as pandas_table_file:
            pandas_table_file.write(',assembly_accession,organism_name,taxid,species_taxid,assembly_level,ftp_path\n')
            pandas_table_file.write('0,{},{},{},{},{},{}\n'.format(
                assembly_accession, organism_name, taxonomic_node, taxonomic_node, assembly_level, 'uploaded genome'
            ))

    except Exception as e:
        raise IntegrityError(
            "[-] ERROR: Couldnt write pandas dataframe for your uploaded genome file, with exception : {}".format(e))

#TODO documentation
def write_pandas_table_for_multiple_uploaded_files(blast_database, genomes_to_organism_and_taxid_dict):
    try:
        path_to_database = BLAST_DATABASE_DIR + str(blast_database.id)+'/'
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

'''get_list_of_taxonomic_nodes_based_on_organism_file
    
    This function is executed in the write_pandas_table_for_uploaded_genomes. It converts the 
    provided organism names to taxonomic ids by taking the first taxid occurence of the biopython call.
    If no taxid is found a integrity error is raised. This should not happend due to Form validation.
    
    :param organisms_file
        :type django.core.files.uploadedfile.TemporaryUploadedFile
    :param user_email
        :type str
    
    :returns taxids, organisms
        :type list[int], list[str]
        
'''
def get_list_of_taxonomic_nodes_based_on_organisms_file(organisms_file,user_email:str):
    try:
        print(type(organisms_file))
        taxids = []
        organisms = []
        for line in organisms_file:
            line = line.decode().rstrip()

            if line != '':
                taxid = pyb.get_species_taxid_by_name(user_email,line)
                taxids.append(taxid[0])
                organisms.append(line)
        return taxids, organisms
    except Exception as e:
        raise IntegrityError('couldnt translate organism names into taxonomic nodes with exception : {}'.format(e))

'''write_pandas_table_for_uploaded_genomes
    
    This function writes the pandas table for the BlastDatabase instance.
    The table is written into media/databases/BlastDatabase.id. Provided organism names
    are translated to taxonomic identifier. Assembly levels and assembly accessions are optional
    and if not provided the table entries are filled with the note: "not provided".
    The columns of the resulting pandas table are: index,assembly_accession,organism_name,taxid,species_taxid,assembly_level,ftp_path.
    The ftp_path column is filled with the string: "uploaded genome".
    If an exception occurs the function raises an integrity error.
    
    :param blast_database
        :type refseq_transactions.models.BlastDatabase
    :param assembly_accessions_file
        :type django.core.files.uploadedfile.TemporaryUploadedFile
    :param assembly_levels_file
        :type django.core.files.uploadedfile.TemporaryUploadedFile
    :param organisms_file
        :type django.core.files.uploadedfile.TemporaryUploadedFile
    :param user_email
        :type str
'''
def write_pandas_table_for_uploaded_genomes(blast_database:BlastDatabase,
                                            assembly_accessions_file,
                                            assembly_levels_file,
                                            organisms_file,
                                            user_email:str):
    try:
        path_to_database = BLAST_DATABASE_DIR + str(blast_database.id)+'/'

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