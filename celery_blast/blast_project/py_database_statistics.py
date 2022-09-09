from blast_project import py_django_db_services as py_db_service
import pandas as pd
from os.path import isfile

'''calculate_database_statistics

    this function calculates database statistics based on the reciprocal results of a finished pipeline object

    :param project_id
        :type int

    :returns df
        :type pd.core.frame.DataFrame
'''
def calculate_database_statistics(project_id: int) -> pd.core.frame.DataFrame:
    try:
        project = py_db_service.get_project_by_id(project_id)
        forward_db = project.project_forward_database
        path_to_db_csv = forward_db.path_to_database_file
        db_name = forward_db.database_name.replace(' ','_').upper()
        path_to_db_csv = path_to_db_csv + '/' + db_name
        if isfile(path_to_db_csv):
            return 0
        else:
            return 1

    except Exception as e:
        raise Exception("ERROR during calculation of database statistics with exception: {}".format(e))