from blast_project import py_django_db_services as py_db_service
import pandas as pd
import matplotlib.pyplot as plt
from os.path import isfile

'''calculate_database_statistics

    this function calculates database statistics based on the reciprocal results of a finished pipeline object

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

        result_data = pd.read_csv(path_to_reciprocal_results,
                                  index_col=0)
        db_df=pd.read_csv(path_to_db_csv,index_col=0)
        db_df.columns=['Assembly Accession','Organism Name','Taxid','Species Taxid', 'Assembly Level', 'FTP Path']
        if forward_db.uploaded_files is False:
            print("[+] LETS GO")
            acc = lambda ids: ids.split('_')[2] + '_' + ids.split('_')[3]
            result_data['Assembly Accession'] = result_data['sseqid']

            result_data['Assembly Accession'] = result_data['Assembly Accession'].map(acc)
            class_counts = []
            for query in result_data['qseqid'].unique():
                df = result_data[result_data['qseqid'] == query]
                if len(df) > 0:
                    df = df.drop_duplicates(subset=['Assembly Accession'], keep='first')
                    m_df = db_df.merge(df, on='Assembly Accession')
                    m_df = m_df.drop_duplicates(subset=['scomnames'], keep='first')
                    values = m_df['class'].value_counts()
                    values = values.rename(query)
                    class_counts.append(values)

            if len(class_counts) > 0:
                df = pd.DataFrame(class_counts)
                df = df.fillna(0)
                df = df.transpose()
                for column in df.columns:
                    fig = df[column].plot.bar(figsize=(12, 8), title=column, grid=True)
                    plot_path = 'media/blast_projects/' + str(project_id)+ '/' + str(column) + '/' + 'database_result_statistics.png'
                    plt.tight_layout()
                    plt.savefig(plot_path)

            print("[+] DONE")
        else:
            print("[+] DO SOMETHING WITH CUSTOM UPLOADED DATABASES")

            return 0

    except Exception as e:
        raise Exception("ERROR during calculation of database statistics with exception: {}".format(e))