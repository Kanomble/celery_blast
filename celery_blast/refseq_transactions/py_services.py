from os.path import isfile, isdir
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

#TODO documentation
def get_database_download_and_formatting_task_result(database_id):
    try:
        blastdb = get_database_by_id(database_id)
        results = blastdb.database_download_and_format_task.result
        try:
            progress = json.loads(results)['percent']
            return progress
        except Exception as e:
            return 'DONE'
    except Exception as e:
        raise Exception("couldnt read progress with exception : {}".format(e))


#TODO documentation
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


#TODO documentation (blastdatabase model - without taks_id ..)
def write_blastdatabase_snakemake_configfile(database_id,task_id):
    try:
        bdb_summary_table_name = get_database_by_id(database_id).get_pandas_table_name()
        configfile_path = 'media/databases/' + str(database_id) + '/snakefile_config'
        with open(configfile_path,'w') as cf_file:
            cf_file.write('db_summary: '+"\""+bdb_summary_table_name+"\""+"\n")
            cf_file.write('task_id: '+"\""+task_id+"\""+"\n")
    except Exception as e:
        raise IntegrityError('exception during writing snakemake configfile with exception : {}'.format(e))

#TODO documentation
def read_database_download_and_format_logfile(database_id):
    try:
        database_file_path = 'media/databases/'+str(database_id)
        snakemake_logfile_path = database_file_path + '/snakemake_progress.log'
        if(isdir(database_file_path)):
            if(isfile(snakemake_logfile_path)):
                fd = open(snakemake_logfile_path, 'r')
                progress = round(float(fd.readlines()[-1]), 3)
                fd.close()
                return progress
            else:
                return 0
        else:
            raise Exception
    except Exception as e:
        raise Exception("unable to track progress out of logfile ...")

#TODO documentation
def get_bdb_summary_table_name(database_id):
    try:
        return get_database_by_id(database_id).get_pandas_table_name()
    except Exception as e:
        raise Exception('couldnt read blast database pandas table name with exception : {} and id : {}'.format(e,database_id))

#TODO documentation
#SETUP UTILITY FUNCTIONS
def get_ftp_paths_and_taxids_from_summary_file(database_id):
    try:
        bdb_summary_table_name = get_database_by_id(database_id).get_pandas_table_name()
        filepath = 'media/databases/' + str(database_id) + '/' + bdb_summary_table_name
        dataframe = pd.read_table(filepath,header=0,index_col=0,delimiter=",")
        #there shouldnt be any duplicates
        dataframe = dataframe[dataframe['ftp_path'].duplicated() == False]
        return dict(zip(dataframe['ftp_path'],dataframe['taxid']))
    except Exception as e:
        raise Exception('couldnt read database summary table with exception : {}'.format(e))

#TODO documentation
def write_progress_database_transactions(database_id, progress):
    try:
        filepath = 'media/databases/' + str(database_id) + '/task_progress.log'
        progress_log = open(filepath,'a')
        progress_log.write(str(progress)+'\n')
        progress_log.close()
        return progress
    except Exception as e:
        raise Exception("error writing logfile exception : {}".format(e))

#TODO documentation
def get_current_progress_database_transactions(database_id):
    try:
        filepath = 'media/databases/' + str(database_id) + '/task_progress.log'
        fd = open(filepath,'r')
        progress = round(float(fd.readlines()[-1]),3)
        fd.close()
        return progress
    except Exception as e:
        raise Exception("error reading logfile exception : {}".format(e))