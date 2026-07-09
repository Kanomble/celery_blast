# functions for the creation of BLAST databases which are triggered via POST and GET requests in the views.py file of this package
from os.path import isfile
from random import choices
from string import digits, ascii_uppercase

import pandas as pd
from blast_project.py_django_db_services import create_and_save_refseq_database_model, get_database_by_id
from blast_project.py_services import create_blastdatabase_directory, upload_file
from blast_project.tasks import write_species_taxids_into_file
from django.db import IntegrityError, transaction

from .models import BlastDatabase
from .py_services import write_pandas_table_to_project_dir, transform_data_table_to_json_dict, \
    filter_duplicates_by_ftp_path
from celery_blast.settings import REFSEQ_ASSEMBLY_FILE, TAXONOMIC_NODES, BLAST_DATABASE_DIR

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

def get_all_databases():
    return BlastDatabase.objects.all()

def get_databases_in_progress():
    return BlastDatabase.objects.get_databases_with_task_on_progress()


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
        
    :returns returncode - 0
        :type int - 0
'''
def create_blastdatabase_table_and_directory(valid_blastdatabase_form):
    try:
        with transaction.atomic():
            database_name = valid_blastdatabase_form.cleaned_data['database_name']
            database_description = valid_blastdatabase_form.cleaned_data['database_description']
            assembly_levels = valid_blastdatabase_form.cleaned_data['assembly_levels']
            assembly_summary_file = valid_blastdatabase_form.cleaned_data['database_summary_file']



            if assembly_summary_file == "GenBank":
                assembly_summary_file = REFSEQ_ASSEMBLY_FILE + 'assembly_summary_genbank.txt'
            elif assembly_summary_file == "RefSeq":
                assembly_summary_file = REFSEQ_ASSEMBLY_FILE + 'assembly_summary_refseq.txt'

            if isfile(assembly_summary_file) == False:
                raise Exception("[-] ERROR there is no assembly_summary_file for {} present".format(assembly_summary_file))

            if len(assembly_levels) == 0:
                assembly_levels = ['Chromosome', 'Scaffold', 'Complete Genome', 'Contig']

            if valid_blastdatabase_form.cleaned_data.get('taxid_file', False):
                # upload taxonomic information file
                taxid_file = valid_blastdatabase_form.cleaned_data['taxid_file']
                taxid_file_path = TAXONOMIC_NODES + taxid_file.name
                upload_file(taxid_file, taxid_file_path)

                # get data from refseq_summary_file and limit by assembly_levels
                refseq_table = read_current_assembly_summary_with_pandas(assembly_levels, assembly_summary_file)
                # limit file by taxonomy
                taxonomy_table = read_taxonomy_table(taxid_file.name)
                filtered_table = filter_table_by_taxonomy(refseq_table, taxonomy_table)

                filtered_table = filter_duplicates_by_ftp_path(filtered_table)

                if len(filtered_table) == 0:
                    raise IntegrityError("the database doesnt contain any entries, pls apply another filter method!")

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
                taxid_file_path = TAXONOMIC_NODES + taxid_file

                refseq_table = read_current_assembly_summary_with_pandas(assembly_levels, assembly_summary_file)
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
            elif valid_blastdatabase_form.cleaned_data.get('taxid_text_field', False):
                taxid_list = valid_blastdatabase_form.cleaned_data['taxid_text_field']
                randomly_generated_filename = ''.join(choices(ascii_uppercase + digits, k=10))
                taxid_filename = taxid_list[0].lower() + '_' + randomly_generated_filename + '_node_file.taxids'
                taxid_file_path = TAXONOMIC_NODES + taxid_filename
                # biopython does not translate higher taxonomic nodes into genus/species level nodes
                try:
                    write_species_taxids_into_file(taxid_list, taxid_filename)
                except Exception as e:
                    raise Exception(
                        "[-] ERROR creating taxonomic_node file for pandas dataframe parsing with exception: {}".format(
                            e))

                refseq_table = read_current_assembly_summary_with_pandas(assembly_levels, assembly_summary_file)
                taxonomy_table = read_taxonomy_table(taxid_filename)
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
                # create directory in BLAST_DATABASE_DIR
                refseq_database_table_path = create_blastdatabase_directory(new_blastdb.id)

                write_pandas_table_to_project_dir(refseq_database_table_path,
                                                  filtered_table,
                                                  database_name)
            else:
                filtered_table = read_current_assembly_summary_with_pandas(assembly_levels, assembly_summary_file)
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

        return 0
    except Exception as e:
        raise IntegrityError('couldnt create blastdatabase : {}'.format(e))


''' read_current_assembly_summary_with_pandas
    
    Transforms the blast database refseq summary file into a pandas dataframe.
    It then filters that dataframe against the assembly_levels that have been chosen by the user and
    relevant fields. Thereby it transforms the ftp_path field into a real ftp path for the assembly
    protein sequence files. The returned dataframe has six columns: 
    ['assembly_accession', 'organism_name', 'taxid', 'species_taxid','assembly_level', 'ftp_path'].
    
    :param assembly_levels
        :type list[str]
    :returns desired_refseq_genomes_dataframe
        :type pd.DataFrame
'''
def read_current_assembly_summary_with_pandas(
    assembly_levels: list,
    summary_file_path: str
) -> pd.DataFrame:

    if not isfile(summary_file_path):
        raise ValueError("Assembly summary file does not exist: {}".format(summary_file_path))

    def set_protein_assembly_file(ftp_path):
        try:
            if isinstance(ftp_path, str) and ftp_path:
                protein_genome = ftp_path.rstrip("/").split("/")[-1]
                return ftp_path.rstrip("/") + "/" + protein_genome + "_protein.faa.gz"
            return ftp_path
        except Exception as e:
            raise ValueError(
                "Problem during parsing the ftp_path column in the assembly summary file"
            ) from e

    try:
        with open(summary_file_path, "r", encoding="utf-8") as rfile:
            rfile.readline()
            header = rfile.readline().lstrip("#").rstrip("\n").split("\t")

        refseq_table = pd.read_table(
            summary_file_path,
            sep="\t",
            skiprows=2,
            header=None,
            names=header,
            dtype=str,
            keep_default_na=False,
        )

        if "assembly_level" in refseq_table.columns and assembly_levels:
            refseq_table = refseq_table[
                refseq_table["assembly_level"].isin(assembly_levels)
            ]

        if "ftp_path" in refseq_table.columns:
            refseq_table["ftp_path"] = refseq_table["ftp_path"].apply(
                set_protein_assembly_file
            )

        return refseq_table[["assembly_accession","organism_name","taxid","species_taxid","assembly_level","ftp_path"]]

    except pd.errors.ParserError as e:
        raise ValueError(
            "Could not parse assembly summary file as a valid tab-separated table. "
            "The file likely contains a malformed line or was corrupted during download. "
            "Original pandas error: {}".format(e)
        ) from e

    except Exception as e:
        raise ValueError(
            "Exception during pandas parsing of assembly summary file: {}".format(e)
        ) from e

''' read_taxonomy_table

    Transforms a taxonomic node file into a pandas dataframe (with one column).
    
    :param taxfilename
        :type str (filename)
    :returns pandas dataframe
'''
def read_taxonomy_table(taxfilename: str) -> pd.DataFrame:
    filepath = TAXONOMIC_NODES + taxfilename
    if (isfile(filepath) == False):
        raise ValueError("there is no taxonomy file called: {}".format(filepath))
    taxonomy_file = pd.read_table(filepath, header=None)
    # species_taxid and taxid should normally be interchangeable, the species_taxid may inherit more information
    # to current strain (have a look at the README description of the refseq summary file)
    taxonomy_file.columns = ['taxid']
    taxonomy_file = taxonomy_file.astype({"taxid": str})
    return taxonomy_file


'''read_taxonomy_list
    
    Transforms a list of taxonomic nodes into a pandas dataframe (with one column).
    
    :param taxid_list
        :type list[int]
    
    :returns taxid_df
        :type pd.DataFrame
        
'''
def read_taxonomy_list(taxid_list: list) -> pd.DataFrame:
    try:
        taxid_df = pd.DataFrame({'taxid': taxid_list})
        taxid_df = taxid_df.astype({"taxid": str})
        return taxid_df
    except Exception as e:
        raise Exception("[-] error creating taxid dataframe with exception: {}".format(e))


''' filter_table_by_taxonomy
    
    uses the pandas.Dataframe.merge() method to create a dataframe filtered on the species_taxid column.
    
    :param refseq_table
        :type pandas dataframe
    :param taxonomy_table
        :type pandas dataframe
    :returns pandas dataframe
'''
def filter_table_by_taxonomy(refseq_table, taxonomy_table):
    # on field can be changed to on=['taxid']
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
    table = pd.read_csv(blastdb.path_to_database_file + '/' + tablefile_name, header=0, index_col=0)
    json = transform_data_table_to_json_dict(table)
    return json


'''create_blastdb_dir_and_table_based_on_user_selection
    
    This function triggers processing of user selected proteomes.
    Selection can be done within the datatable_blast_database_detail.html file.
    
    :param form_data_dict
        :type dict
        
    :returns database_id
        :type int
'''
def create_blastdb_dir_and_table_based_on_user_selection(form_data_dict:dict)->int:
    try:

        try:
            keys = list(form_data_dict.keys())
            values = [form_data_dict[key].split(',') for key in keys]
        except Exception as e:
            raise Exception("[-] ERROR creating pandas dataframe from dictionary with exception: {}".format(e))

        selected_proteome_table = pd.DataFrame(values)
        selected_proteome_table.columns = ["assembly_accession","organism_name","taxid","species_taxid","assembly_level","ftp_path"]
        assembly_levels = list(selected_proteome_table.assembly_level.unique())

        database_name = "selected_proteomes_database"
        counter = 1
        databases = get_all_databases()
        for db in databases:
            if database_name in db.database_name:
                counter += 1
        database_name = database_name + "_" + str(counter)

        database_description = "selected proteomes database"
        new_blastdb = create_and_save_refseq_database_model(
                    database_name=database_name,
                    database_description=database_description,
                    assembly_levels=assembly_levels,
                    assembly_entries=len(selected_proteome_table))
        database_table_path = create_blastdatabase_directory(new_blastdb.id)

        write_pandas_table_to_project_dir(database_table_path,
                                          selected_proteome_table,
                                          database_name)

        return new_blastdb.id
    except Exception as e:
        raise Exception("[-] ERROR processing selected proteomes with exception: {}".format(e))


# TODO not in use
'''check_for_db_updates

    This function searches for BLAST database updates if the previous download and formatting task 
    finished successfully.
    
    :param blast_database_id
        :type int
'''


def check_for_db_updates(blast_database_id: int):
    try:
        # retrieve database model
        blast_db = get_database_by_id(blast_database_id)
        # check if previous task finished
        if blast_db.database_download_and_format_task.status == 'SUCCESS' and blast_db.uploaded_files == False:
            # retrieve taxonomic node file
            taxid_filepath = blast_db.attached_taxonomic_node_file
            taxonomic_nodes = read_taxonomy_table(taxid_filepath)
            # retrieving assembly level information
            ass_levels = []
            for level in blast_db.assembly_levels.all():
                ass_levels.append(level.assembly_level)
            # retrieve taxonomic information
            assembly_summary = read_current_assembly_summary_with_pandas(assembly_levels=ass_levels)
            # filter assembly summary with taxonomic information
            filtered_table = filter_table_by_taxonomy(assembly_summary, taxonomic_nodes)
            # reading current blast_db table
            blast_db_table_path = blast_db.path_to_database_file + '/' + blast_db.get_pandas_table_name()

        # download and format task: progress or failure -> return 1
        else:
            return 1
        return 0
    except Exception as e:
        raise Exception("[-] ERROR updating blast database with exception: {}".format(e))
