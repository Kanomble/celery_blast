from os.path import isfile
from django.db import IntegrityError
import pandas as pd
import json
''' refseq_file_exists
    
    checks if the refseq_summary_file exists in the desired media directory
    
    :returns
        :type boolean
'''
def refseq_file_exists():
    return isfile('media/databases/refseq_summary_file/assembly_summary_refseq.txt')

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
