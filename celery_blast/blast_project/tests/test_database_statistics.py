from django.test import TestCase, tag
import pandas as pd
from blast_project.py_database_statistics import add_taxonomic_information_to_db, \
                                                extract_taxonomic_information
from os.path import isfile

class DatabaseStatisticsTest(TestCase):

    path_to_test_files = 'testfiles/database_statistics/'
    path_to_log_files = 'testfiles/database_statistics/log/'

    def test_add_taxonomic_information_to_db(self):
        logfile_path = self.path_to_log_files + 'add_taxonomic_information_to_db.log'
        database_dataframe = pd.read_csv(self.path_to_test_files + 'HIGH_QUALITY_CYANOBACTERIA_DATABASE', index_col=0)
        taxonomy_df = add_taxonomic_information_to_db('lukas.becker@hhu.de',
                                                      logfile_path,
                                                      list(database_dataframe['taxid'].unique()))

        self.assertTrue(isfile(logfile_path))
        self.assertTrue(len(taxonomy_df) == len(list(database_dataframe['taxid'].unique())))

    def test_extract_taxonomic_information(self):
        logfile_path = self.path_to_log_files + 'extract_taxonomic_information.log'
        database_dataframe = pd.read_csv(self.path_to_test_files + 'HIGH_QUALITY_CYANOBACTERIA_DATABASE', index_col=0)
        reciprocal_results_dataframe = pd.read_csv(self.path_to_test_files + 'reciprocal_results_with_taxonomy.csv', index_col=0)
        tax_counts = extract_taxonomic_information(logfile_path, False,
                                                   reciprocal_results_dataframe,database_dataframe,'order')

        for query_list in tax_counts:
            if query_list.name == "WP_010872350":
                synechococcus = query_list['Synechococcales']
                nostocales = query_list['Nostocales']
                chroococcales = query_list['Chroococcales']

        self.assertTrue(len(tax_counts) == 5)
        self.assertTrue(type(tax_counts) == list)
        self.assertTrue(type(tax_counts[0]) == pd.core.series.Series)
        self.assertTrue(synechococcus == 97)
        self.assertTrue(nostocales == 69)
        self.assertTrue(chroococcales == 24)