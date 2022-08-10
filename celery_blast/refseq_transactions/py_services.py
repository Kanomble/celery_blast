#service function for the BLAST database formatting and download procedures as well as for the view functions of this package
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

'''get_database_download_and_formatting_task_result_progress
    
    This function is used within an ajax_call that sends a json response to the client.
    The json object contains the "progress" key with the celery task progress variable, that is defined by
    the ProgressRecorder class.
    
    :param database_id
        :type int
    
    :returns progress
        :type str
'''
def get_database_download_and_formatting_task_result_progress(database_id:int)->str:
    try:
        blastdb = get_database_by_id(database_id)
        results = blastdb.database_download_and_format_task.result
        progress = json.loads(results)['percent']
        return progress
    except Exception as e:
        raise Exception("couldnt read progress with exception : {}".format(e))


'''filter_duplicates_by_ftp_path
    This function deletes duplicate ftp_path entries in the BLAST database pandas dataframe.
    The function is executed during database table formatting in the create_blastdatabase_table_and_directory function. 
    
    :params database_table
        :type pd.DataFrame
        
    :returns database_table
        :type pd.DataFrame
'''
def filter_duplicates_by_ftp_path(database_table:pd.DataFrame)->pd.DataFrame:
    try:
        database_table = database_table[database_table['ftp_path'].duplicated() == False]
        if(len(database_table) == 0):
            raise Exception('there are no entries in the pandas table')
        return database_table
    except Exception as e:
        raise IntegrityError('couldnt filter pandas table by duplicates Exception : {}'.format(e))

'''write_pandas_table_to_project_dir
    This function converts the BLAST database pandas dataframe to csv file, which is saved in the 
    corresponding database directory. The csv file is named after the database_name parameter, in uppercase characters,
    and with whitespace converted to an underscore character.
    The function is executed during database table formatting in the create_blastdatabase_table_and_directory function.
    
    :params blastdatabase_path
        :type str
        
    :params database_table
        :type pd.DataFrame
    
    :params database_name
        :type str
'''
def write_pandas_table_to_project_dir(blastdatabase_path:str, database_table:pd.DataFrame, database_name:str):
    try:
        database_table_filepath = blastdatabase_path + '/' + database_name.replace(' ', '_').upper()
        database_table.to_csv(database_table_filepath)
    except Exception as e:
        raise IntegrityError('couldnt write pandas table to refseq genome directory: {}'.format(e))

'''transform_data_table_to_json_dict

    This function is executed in the read_database_table_by_database_id_and_return_json function, 
    it is used for pandas dataframe conversion to an json object. 
    The json object is used in a view ajax function for the HTML table, that is in turn controlled by the
    DataTables package.
    
    :params df - BLAST database dataframe
        :type pd.DataFrame
        
    :returns data
        :type json object
'''
def transform_data_table_to_json_dict(df:pd.DataFrame):
    json_records = df.reset_index().to_json(orient='records')
    data = json.loads(json_records)
    return data

'''get_bdb_summary_table_name

    This function is executed in the celery task download_blast_databases_based_on_summary_file.
    The function returns the BLAST database pandas dataframe name, which is saved within the postgresql database. 
    
    :params database_id
        :type int
        
    :returns BLAST database pandas table name 
        :type str
'''
def get_bdb_summary_table_name(database_id):
    try:
        return get_database_by_id(database_id).get_pandas_table_name()
    except Exception as e:
        raise Exception('couldnt read blast database pandas table name with exception : {} and id : {}'.format(e,database_id))

'''get_ftp_paths_and_taxids_from_summary_file
    
    This function is executed in the celery task download_blast_databases_based_on_summary_file.
    The function returns a dictionary composed of genome ftp_paths (to the NCBI refseq ftp-server) and
    the corresponding taxonomic identifier. 
    
    :params database_id
        :type int
        
    :returns dictionary_ftp_paths_taxids
        :type dict - dict[ftp_path,taxid]

'''
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

#TODO this function is currently not used, in the near future it can be used to fill a custom logfile
def write_progress_database_transactions(database_id:int, progress:str)->str:
    try:
        filepath = 'media/databases/' + str(database_id) + '/task_progress.log'
        progress_log = open(filepath,'a')
        progress_log.write(str(progress)+'\n')
        progress_log.close()
        return progress
    except Exception as e:
        raise Exception("error writing logfile exception : {}".format(e))

#TODO this function is currently not used, in the near future it can be used to fill a custom logfile
def get_current_progress_database_transactions(database_id):
    try:
        filepath = 'media/databases/' + str(database_id) + '/task_progress.log'
        fd = open(filepath,'r')
        progress = round(float(fd.readlines()[-1]),3)
        fd.close()
        return progress
    except Exception as e:
        raise Exception("error reading logfile exception : {}".format(e))