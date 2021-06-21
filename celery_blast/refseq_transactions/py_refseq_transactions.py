from blast_project.models import BlastDatabase
from blast_project.py_services import create_blastdatabase_directory, upload_file
from blast_project.py_django_db_services import create_and_save_refseq_database_model, get_database_by_id
from .py_services import write_pandas_table_to_project_dir, transform_data_table_to_json_dict, filter_duplicates_by_ftp_path
from os.path import isfile
import pandas as pd
from django.db import IntegrityError, transaction

''' 
transactions with models (manager)
'''
def get_databases_without_tasks():
    return BlastDatabase.objects.get_databases_without_executed_tasks()

def get_databases_with_tasks():
    return BlastDatabase.objects.get_databases_with_executed_tasks()

def get_downloaded_databases():
    return BlastDatabase.objects.get_databases_with_succeeded_tasks()

def get_failed_tasks():
    return BlastDatabase.objects.get_databases_with_failed_tasks()

def get_databases_in_progress():
    return BlastDatabase.objects.get_databases_with_task_on_progress()

#TODO this function needs to get refactored maybe put some functionality in the form and model for BlastDatabase
''' create_blastdatabase_table_and_directory
    
    Processes the form data of the create_blast_database_model_and_directory view POST request.
    It checks wether a database is limited by taxonomy. 
    Therefore the valid form can inherit a taxonomy file or can use a present taxonomy file (media/taxonomic_node_files).
    If a filehandle is provided it gets uploaded into the media/taxonomic_node_files folder.
    Then the assembly summary file with the filepath media/refseq_summary_file/assembly_summary_refseq.txt gets transformed into 
    a pandas table with the read_current_assembly_summary_with_pandas(assembly_levels) function. Assembly levels come from the form.
    The table gets either limited by taxonomy or not. A BlastDatabase instance is created and saved into the database with the 
    create_and_save_refseq_database_model(...) function. With the create_blastdatabase_directory(database_id) and 
    write_pandas_table_to_project_dir(...) function, a database directory with the associated summary file is created.
    
    :param valid_blastdatabase_form
        :django.form.is_valid() == True
    :returns None
'''
def create_blastdatabase_table_and_directory(valid_blastdatabase_form):
    try:
        with transaction.atomic():
            database_name = valid_blastdatabase_form.cleaned_data['database_name']
            database_description = valid_blastdatabase_form.cleaned_data['database_description']
            assembly_levels = valid_blastdatabase_form.cleaned_data['assembly_levels']

            if valid_blastdatabase_form.cleaned_data.get('taxid_file', False):
                # upload taxonomic information file
                taxid_file = valid_blastdatabase_form.cleaned_data['taxid_file']
                taxid_file_path = 'media/' + 'taxonomic_node_files/' + taxid_file.name
                upload_file(taxid_file, taxid_file_path)

                # get data from refseq_summary_file and limit by assembly_levels
                refseq_table = read_current_assembly_summary_with_pandas(assembly_levels)
                # limit file by taxonomy
                taxonomy_table = read_taxonomy_table(taxid_file.name)
                filtered_table = filter_table_by_taxonomy(refseq_table, taxonomy_table)

                filtered_table = filter_duplicates_by_ftp_path(filtered_table)


                if len(filtered_table) == 0:
                    raise IntegrityError("the database doesnt contain any entries, pls apply an other filter method!")

                # create new refseq database model and save it into the database
                new_blastdb = create_and_save_refseq_database_model(
                    database_name=database_name,
                    database_description=database_description,
                    assembly_levels=assembly_levels,
                    assembly_entries=len(filtered_table),
                    attached_taxonomic_file=taxid_file_path)

                # create directory in media/databases/
                refseq_database_table_path = create_blastdatabase_directory(new_blastdb.id)

                write_pandas_table_to_project_dir(refseq_database_table_path,
                                                  filtered_table,
                                                  database_name)

            elif valid_blastdatabase_form.cleaned_data.get('taxid_uploaded_file', False):
                taxid_file = valid_blastdatabase_form.cleaned_data['taxid_uploaded_file']
                taxid_file_path = 'media/' + 'taxonomic_node_files/' + taxid_file

                refseq_table = read_current_assembly_summary_with_pandas(assembly_levels)
                taxonomy_table = read_taxonomy_table(taxid_file)
                filtered_table = filter_table_by_taxonomy(refseq_table, taxonomy_table)
                filtered_table = filter_duplicates_by_ftp_path(filtered_table)

                if len(filtered_table) == 0:
                    raise IntegrityError("the database doesnt contain any entries, pls apply an other filter method!")

                    # create new refseq database model and save it into the database
                new_blastdb = create_and_save_refseq_database_model(
                    database_name=database_name,
                    database_description=database_description,
                    assembly_levels=assembly_levels,
                    assembly_entries=len(filtered_table),
                    attached_taxonomic_file=taxid_file_path)
                # create directory in media/databases
                refseq_database_table_path = create_blastdatabase_directory(new_blastdb.id)

                write_pandas_table_to_project_dir(refseq_database_table_path,
                                                  filtered_table,
                                                  database_name)
            else:
                filtered_table = read_current_assembly_summary_with_pandas(assembly_levels)
                filtered_table = filter_duplicates_by_ftp_path(filtered_table)
                new_blastdb = create_and_save_refseq_database_model(
                    database_name=database_name,
                    database_description=database_description,
                    assembly_levels=assembly_levels,
                    assembly_entries=len(filtered_table))

                refseq_database_table_path = create_blastdatabase_directory(new_blastdb.id)

                write_pandas_table_to_project_dir(refseq_database_table_path,
                                                  filtered_table,
                                                  database_name)
    except Exception as e:
        raise IntegrityError('couldnt create blastdatabase : {}'.format(e))


''' read_current_assembly_summary_with_pandas
    
    Transforms the blast database refseq summary file into a pandas dataframe.
    It then filters that dataframe against the assembly_levels that have been chosen by the user and
    relevant fields. Thereby it transforms the ftp_path field into a real ftp path for the assembly
    protein sequence files. The returned dataframe has six columns: 
    ['assembly_accession', 'organism_name', 'taxid', 'species_taxid','assembly_level', 'ftp_path'].
    
    :param assembly_levels
        :type list with str objects
    :returns pandas dataframe
'''
def read_current_assembly_summary_with_pandas(assembly_levels):
    summary_file_path = 'media/databases/refseq_summary_file/assembly_summary_refseq.txt'
    if(isfile(summary_file_path) == False):
        raise ValueError('assembly summary file does not exist!')

    #function for changing the ftp_header in the pandas table
    def set_protein_assembly_file(ftp_path):
        protein_genome = ftp_path.split('/')[-1:][0]
        protein_genome = ftp_path + '/' + str(protein_genome) + '_protein.faa.gz'
        return protein_genome

    #init parsing refseq table with pandas
    try:
        refseq_table = pd.read_table(summary_file_path, skiprows=[0, 1], header=None, usecols=range(22))
        header = ["assembly_accession", "bioproject", "biosample", "wgs_master", "refseq_category", "taxid",
                  "species_taxid", "organism_name", "infraspecific_name", "isolate", "version_status", "assembly_level",
                  "release_type", "genome_rep", "seq_rel_date", "asm_name", "submitter", "gbrs_paired_asm",
                  "paired_asm_comp", "ftp_path", "excluded_from_refseq", "relation_to_type_material"]
        refseq_table.columns = header
    except Exception as e:
        raise ValueError("exception during pandas parsing of assembly_summary_refseq.txt file ...\n\tException: {}".format(e))

    #extract necessary data fields: assembly number, names, taxids and the correct ftp_filepath for downloading with gzip
    try:
        refseq_table = refseq_table[['assembly_accession', 'organism_name', 'taxid', 'species_taxid','assembly_level', 'ftp_path']]
        #python lambda function applied to each row in the dataframe
        refseq_table['ftp_path'] = refseq_table['ftp_path'].apply(lambda row: set_protein_assembly_file(row))

        pandas_genome_level_dataframes = []
        for genome_level in assembly_levels:
            pandas_genome_level_dataframes.append(refseq_table[refseq_table['assembly_level'] == genome_level])

        desired_refseq_genomes_dataframe = pd.concat(pandas_genome_level_dataframes)
        return desired_refseq_genomes_dataframe

    except Exception as e:
        raise ValueError("exception during extraction of smaller dataframe from refseq_table dataframe ...\n\tException: {}".format(e))


''' read_taxonomy_table

    Transforms a taxonomic node file into a pandas dataframe (with one column).
    
    :param taxfilename
        :type str (filename)
    :returns pandas dataframe
'''
def read_taxonomy_table(taxfilename):
    filepath = 'media/taxonomic_node_files/' + taxfilename
    if (isfile(filepath) == False):
        raise ValueError("there is no taxonomy file called: {}".format(filepath))
    taxonomy_file = pd.read_table(filepath, header=None)
    # species_taxid and taxid should normally be interchangeable, the species_taxid may inherit more informations
    # to current strain (have a look at the README description of the refseq summary file)
    taxonomy_file.columns = ['taxid']
    return taxonomy_file

''' filter_table_by_taxonomy
    
    uses the pandas.Dataframe.merge() method to create a dataframe filtered on the species_taxid column.
    
    :param refseq_table
        :type pandas dataframe
    :param taxonomy_table
        :type pandas dataframe
    :returns pandas dataframe
'''
def filter_table_by_taxonomy(refseq_table, taxonomy_table):
    #on field can be changed to on=['taxid']
    return refseq_table.merge(taxonomy_table, how='inner', on=['taxid'])

''' read_database_table_by_database_id_and_return_json
    
    This function is called by the ajax call in the datatable_blast_database_details.html template.
    
    :param database_id
        :type int
    :returns jsonarray
'''
def read_database_table_by_database_id_and_return_json(database_id):
    blastdb = get_database_by_id(database_id)
    tablefile_name = blastdb.database_name.replace(' ', '_').upper()
    print("[++++++]",tablefile_name)
    table = pd.read_csv(blastdb.path_to_database_file + '/' + tablefile_name,header=0,index_col=0)
    print(table)
    json = transform_data_table_to_json_dict(table)
    return json

