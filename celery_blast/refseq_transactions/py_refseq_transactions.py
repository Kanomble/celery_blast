from blast_project.models import BlastDatabase
from blast_project.py_services import create_blastdatabase_directory, upload_file
from blast_project.py_django_db_services import create_and_save_refseq_database_model
from .py_services import write_pandas_table_to_project_dir
from os.path import isfile
import pandas as pd
from django.db import IntegrityError, transaction

def get_databases_without_tasks():
    return BlastDatabase.objects.get_databases_without_executed_tasks()

def get_databases_with_tasks():
    return BlastDatabase.objects.get_databases_with_executed_tasks()

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

                if len(filtered_table) == 0:
                    raise IntegrityError("the database doesnt contain any entries, pls apply an other filter method!")

                # create new refseq database model and save it into the database
                new_blastdatabase = create_and_save_refseq_database_model(
                    database_name=database_name,
                    database_description=database_description,
                    assembly_levels=assembly_levels,
                    assembly_entries=len(filtered_table),
                    attached_taxonomic_file=taxid_file_path)

                # create directory in media/databases/
                refseq_database_table_path = create_blastdatabase_directory(new_blastdatabase.id)

                write_pandas_table_to_project_dir(refseq_database_table_path,
                                                  filtered_table,
                                                  database_name)

            elif valid_blastdatabase_form.cleaned_data.get('taxid_uploaded_file', False):
                taxid_file = valid_blastdatabase_form.cleaned_data['taxid_uploaded_file']
                taxid_file_path = 'media/' + 'taxonomic_node_files/' + taxid_file

                refseq_table = read_current_assembly_summary_with_pandas(assembly_levels)
                taxonomy_table = read_taxonomy_table(taxid_file)
                filtered_table = filter_table_by_taxonomy(refseq_table, taxonomy_table)

                if len(filtered_table) == 0:
                    raise IntegrityError("the database doesnt contain any entries, pls apply an other filter method!")

                    # create new refseq database model and save it into the database
                new_blastdatabase = create_and_save_refseq_database_model(
                    database_name=database_name,
                    database_description=database_description,
                    assembly_levels=assembly_levels,
                    assembly_entries=len(filtered_table),
                    attached_taxonomic_file=taxid_file_path)
                # create directory in media/databases/refseq_databases
                refseq_database_table_path = create_blastdatabase_directory(new_blastdatabase.id)

                write_pandas_table_to_project_dir(refseq_database_table_path,
                                                  filtered_table,
                                                  database_name)
            else:
                filtered_table = read_current_assembly_summary_with_pandas(assembly_levels)
                new_blastdatabase = create_and_save_refseq_database_model(
                    database_name=database_name,
                    database_description=database_description,
                    assembly_levels=assembly_levels,
                    assembly_entries=len(filtered_table))

                refseq_database_table_path = create_blastdatabase_directory(new_blastdatabase.id)

                write_pandas_table_to_project_dir(refseq_database_table_path,
                                                  filtered_table,
                                                  database_name)
    except Exception as e:
        raise IntegrityError('couldnt create blastdatabase : {}'.format(e))


def read_current_assembly_summary_with_pandas(refseq_level_checklist):
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
        refseq_table = pd.read_table(summary_file_path, skiprows=[0, 1], header=None)
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
        refseq_table['ftp_path'] = refseq_table['ftp_path'].apply(lambda row: set_protein_assembly_file(row))

        pandas_genome_level_dataframes = []
        for genome_level in refseq_level_checklist:
            pandas_genome_level_dataframes.append(refseq_table[refseq_table['assembly_level'] == genome_level])

        desired_refseq_genomes_dataframe = pd.concat(pandas_genome_level_dataframes)

        #tuple list for dropdown menu, not implemented yet
        #html_input_list = tuple(zip(refseq_table['assembly_accession'], refseq_table['organism_name']))
    except Exception as e:
        raise ValueError("exception during extraction of smaller dataframe from refseq_table dataframe ...\n\tException: {}".format(e))
    return desired_refseq_genomes_dataframe

def read_taxonomy_table(taxfilename):
    filepath = 'media/taxonomic_node_files/' + taxfilename
    if (isfile(filepath) == False):
        raise ValueError("there is no taxonomy file called: {}".format(filepath))
    taxonomy_file = pd.read_table(filepath, header=None)
    # species_taxid and taxid should normally be interchangeable, the species_taxid may inherit more informations
    # to current strain (have a look at the README description of the refseq summary file)
    taxonomy_file.columns = ['species_taxid']
    return taxonomy_file

def filter_table_by_taxonomy(refseq_table, taxonomy_table):
    # species_taxid
    return refseq_table.merge(taxonomy_table, how='inner', on=['species_taxid'])