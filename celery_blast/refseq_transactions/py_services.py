from os.path import isfile
from django.db import IntegrityError
from blast_project.py_django_db_services import get_database_by_id
import pandas as pd
import json
''' refseq_file_exists
    
    checks if the refseq_summary_file exists in the desired media directory
    
    :returns
        :type boolean
'''
def refseq_file_exists():
    return isfile('media/databases/refseq_summary_file/assembly_summary_refseq.txt')


def filter_duplicates_by_ftp_path(pandas_table):
    try:
        pandas_table = pandas_table[pandas_table['ftp_path'].duplicated() == False]
        if(len(pandas_table) == 0):
            raise Exception('there are no entries in the pandas table')
        return pandas_table
    except Exception as e:
        raise IntegrityError('couldnt filter pandas table by duplicates Exception : {}'.format(e))

#TODO documentation
def write_pandas_table_to_project_dir(blastdatabase_path, pandas_table, database_name):
    try:
        pandas_table_filepath = blastdatabase_path + '/' + database_name.replace(' ', '_').upper()
        pandas_table.to_csv(pandas_table_filepath)
    except Exception as e:
        raise IntegrityError('couldnt write pandas table to refseq genome directory: {}'.format(e))

#TODO documentation
def transform_data_table_to_json_dict(df):
    json_records = df.reset_index().to_json(orient='records')
    data = []
    data = json.loads(json_records)
    return data

#TODO documentation
def write_blastdatabase_snakemake_configfile(database_id,task_id):
    try:
        bdb_summary_table_name = get_database_by_id(database_id).get_pandas_table_name()
        configfile_path = 'media/databases/' + str(database_id) + '/snakefile_config'
        with open(configfile_path,'w') as cf_file:
            cf_file.write('db_summary: '+"\""+bdb_summary_table_name+"\""+"\n")
            cf_file.write('task_id: '+"\""+task_id+"\""+"\n")
    except Exception as e:
        raise IntegrityError('exception during writing snakemake configfile with exception : {}'.format(e))