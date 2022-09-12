from blast_project import py_django_db_services as py_db_service
import pandas as pd
import matplotlib.pyplot as plt
from os.path import isfile, isdir
from os import remove, listdir
from django.db import IntegrityError

#TODO extend this function to produce not only class based statistics
'''calculate_database_statistics

    This function calculates database statistics based on the reciprocal results of a finished pipeline project.
    This function depends on the database csv file, therefore it is not included within the Snakemake pipeline.

    :param project_id
        :type int

    :returns df
        :type pd.core.frame.DataFrame
'''
def calculate_database_statistics(project_id: int)->pd.core.frame.DataFrame:
    try:
        project = py_db_service.get_project_by_id(project_id)
        forward_db = project.project_forward_database
        path_to_db_csv = forward_db.path_to_database_file
        db_name = forward_db.database_name.replace(' ','_').upper()
        path_to_db_csv = path_to_db_csv + '/' + db_name
        path_to_reciprocal_results = 'media/blast_projects/' + str(project_id) + '/reciprocal_results_with_taxonomy.csv'
        if isfile(path_to_reciprocal_results) is False or isfile(path_to_db_csv) is False:
            raise FileNotFoundError

        result_data = pd.read_csv(path_to_reciprocal_results,index_col=0)
        db_df=pd.read_csv(path_to_db_csv,index_col=0)

        #database is based on refseq genome files --> filtering based on assembly_accessions and scomnames
        if forward_db.uploaded_files is False:
            acc = lambda ids: ids.split('_')[2] + '_' + ids.split('_')[3]
            result_data['assembly_accession'] = result_data['sseqid']

            result_data['assembly_accession'] = result_data['assembly_accession'].map(acc)
            class_counts = []
            for query in result_data['qseqid'].unique():
                df = result_data[result_data['qseqid'] == query]
                if len(df) > 0:
                    df = df.drop_duplicates(subset=['assembly_accession'], keep='first')
                    m_df = db_df.merge(df, on='assembly_accession')
                    m_df = m_df.drop_duplicates(subset=['scomnames'], keep='first')
                    values = m_df['class'].value_counts()
                    values = values.rename(query)
                    class_counts.append(values)

        #database is based on custom uploaded files --> assembly accessions are not accurate, therfore filtering is based on taxids and scomnames
        else:
            class_counts = []
            result_data.rename(columns={'staxids': 'species_taxid'}, inplace=True)
            for query in result_data['qseqid'].unique():
                df = result_data[result_data['qseqid'] == query]
                if len(df) > 0:
                    db_df = db_df.drop_duplicates(subset=["assembly_accession"], keep='first')
                    db_df = db_df.drop_duplicates(subset=["species_taxid"], keep='first')

                    m_df = db_df[['species_taxid', 'organism_name']].merge(df, on='species_taxid')
                    m_df = m_df.drop_duplicates(subset=['scomnames'], keep='first')
                    values = m_df['class'].value_counts()
                    values = values.rename(query)
                    class_counts.append(values)

        if len(class_counts) > 0:
            df = pd.DataFrame(class_counts)
            df = df.fillna(0)
            df = df.transpose()
            df_path = 'media/blast_projects/' + str(project_id)+ '/database_statistics.csv'
            df.to_csv(df_path)
            for column in df.columns:
                df[column].plot.bar(figsize=(12, 8), title=column, grid=True)#TODO directly saving?
                plot_path = 'media/blast_projects/' + str(project_id)+ '/' + str(column) + '/' + 'database_result_statistics.png'
                plt.tight_layout()
                plt.savefig(plot_path)
        else:
            df=pd.DataFrame()
        return df

    except Exception as e:
        raise Exception("ERROR during calculation of database statistics with exception: {}".format(e))

'''get_database_statistics_task_status
    
    This function returns the status of the database statistic task.
    Status not executed: NOTEXEC
    Status ongoing: PROGRESS
    Status finished: SUCCESS
    Status failed: FAILURE
    
    :param project_id
        :type int
        
    :return status
        :type str
'''
def get_database_statistics_task_status(project_id:int)->str:
    try:
        project = py_db_service.get_project_by_id(project_id)
        if(project.project_database_statistics_task):
            task_status = str(project.project_database_statistics_task.status)
        else:
            task_status = "NOTEXEC"
        return task_status
    except Exception as e:
        raise Exception("ERROR during database statistics task status query with exception: {}".format(e))

'''delete_database_statistics_task_and_output
    
    This function removes the output generated by the calculate_database_statistics_task.
    Returns if function was successfully executed, else exception is thrown.
    
    :param project_id
        :type int
    
    :return returncode
        :type int
'''
def delete_database_statistics_task_and_output(project_id:int)->int:
    try:

        project = py_db_service.get_project_by_id(project_id)
        if(project.project_database_statistics_task):
            project.project_database_statistics_task.delete()
        else:
            raise Exception("ERROR there is no database statistics task for this project ...")

        path_to_project='media/blast_projects/'+str(project_id)
        if isdir(path_to_project):
            path_to_database_statistics='media/blast_projects/' + str(project_id)+ '/database_statistics.csv'
            if isfile(path_to_database_statistics):
                remove(path_to_database_statistics)
            for directory in listdir(path_to_project):
                directory=path_to_project+'/'+directory
                if isdir(directory):
                    path_to_database_statistics_image=directory+'/'+'database_result_statistics.png'
                    if(isfile(path_to_database_statistics_image)):
                        remove(path_to_database_statistics_image)
        return 0
    except IntegrityError as e:
        raise IntegrityError("ERROR during deletion of database statistics task result object with exception: {}".format(e))
    except Exception as e:
        raise Exception("ERROR during removing database statistics output with exception: {}".format(e))