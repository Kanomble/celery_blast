from blast_project import py_django_db_services as py_db_service
import pandas as pd
import matplotlib.pyplot as plt
from os.path import isfile, isdir
from os import remove, listdir
from django.db import IntegrityError
from Bio import Entrez
import json

'''add_taxonomic_information_to_db

    This function adds taxonomic information to the database dataframe.
    It uses the biopython entrez interface.

    :param user_email
        :type str
    :param logfile
        :type str
    :param taxids
        :type list

    :returns db_df
        :type pandas.core.frame.DataFrame
'''
def add_taxonomic_information_to_db(user_email: str,logfile:str, taxids:list) -> pd.core.frame.DataFrame:
    try:
        with open(logfile,'w') as log:
            Entrez.email = user_email

            taxid = []
            taxonomy = []
            genus = []
            superfamily = []
            family = []
            order = []
            classt = []
            phylum = []

            # looping over taxids in db df
            # looping steps for every 500 taxonomic identifier
            end = len(taxids)
            begin = 0
            step = 500
            steps = 500

            log.write("INFO:Starting Analysis from 0 to {}\n".format(end))
            while begin < end:
                if step >= end:
                    step = end
                splitted_ids = taxids[begin:step]
                for attempt in range(10):
                    try:
                        handle = Entrez.efetch(id=splitted_ids, db="taxonomy", retmode="xml")
                        record = Entrez.read(handle)
                        handle.close()
                    except Exception as e:
                        if attempt == 9:
                            raise Exception

                    else:
                        for i in range(len(record)):
                            taxonomy.append(record[i]['ScientificName'])
                            if 'AkaTaxIds' in record[i].keys():
                                for akaid in record[i]['AkaTaxIds']:
                                    if int(akaid) in splitted_ids:
                                        taxid.append(akaid)
                                        log.write("\tINFO: AkaTaxIds detected: {}\n".format(akaid))
                                        break
                                else:
                                    taxid.append(record[i]['TaxId'])
                            else:
                                taxid.append(record[i]['TaxId'])
                            for j in record[i]['LineageEx']:
                                if j['Rank'] == 'genus':
                                    genus.append(j['ScientificName'])
                                if j['Rank'] == 'superfamily':
                                    superfamily.append(j['ScientificName'])
                                if j['Rank'] == 'family':
                                    family.append(j['ScientificName'])
                                if j['Rank'] == 'order':
                                    order.append(j['ScientificName'])
                                if j['Rank'] == 'class':
                                    classt.append(j['ScientificName'])
                                if j['Rank'] == 'phylum':
                                    phylum.append(j['ScientificName'])

                            if (len(taxonomy) != len(genus)):
                                genus.append('unknown')
                            if (len(taxonomy) != len(superfamily)):
                                superfamily.append('unknown')
                            if (len(taxonomy) != len(family)):
                                family.append('unknown')
                            if (len(taxonomy) != len(order)):
                                order.append('unknown')
                            if (len(taxonomy) != len(classt)):
                                classt.append('unknown')
                            if (len(taxonomy) != len(phylum)):
                                phylum.append('unknown')

                            if len(record) != len(splitted_ids):
                                missing_ids = [m_taxid for m_taxid in splitted_ids if m_taxid not in taxid]
                                for m_taxid in missing_ids:
                                    log.write("WARNING: problem during fetching of taxonomic information for: {}\n".format(m_taxid))
                                    taxid.append(m_taxid)
                                    taxonomy.append('unknown')
                                    genus.append('unknown')
                                    superfamily.append('unknown')
                                    family.append('unknown')
                                    order.append('unknown')
                                    classt.append('unknown')
                                    phylum.append('unknown')
                        break

                log.write("INFO: Done with chunk: {} - {}\n".format(begin, step))
                begin += steps
                step += steps

            columns = ["taxid", 'organism_name_taxdb', 'genus', 'family', 'superfamily', 'order', 'class', 'phylum']
            tax_db = pd.DataFrame([taxid, taxonomy, genus, family, superfamily, order, classt, phylum])
            tax_db = tax_db.transpose()
            tax_db.columns = columns
            return tax_db
    except Exception as e:
        log.write("ERROR:problem during fetching and appending taxonomic information with exception: {}\n".format(e))
        raise Exception("ERROR during database table extension with taxonomic information with excpetion: {}".format(e))


'''extract_taxonomic_information

    This function loops over the queries in the reciprocal result dataframe. 
    Each query specific result dataframe is getting merged with the dataframe of the database that inherits the taxonomic informations.
    Therefore it

    :param logfile
        :type str
    :param uploaded
        :type bool
    :param rbh_df
        :type pd.core.frame.DataFrame
    :param db_df
        :type pd.core.frame.DataFrame
    :param taxonomic_unit
        :type str

    :returns tax_counts
        :type list

'''
def extract_taxonomic_information(logfile: str, uploaded: bool, rbh_df: pd.core.frame.DataFrame,
                                  db_df: pd.core.frame.DataFrame, taxonomic_unit: str) -> list:
    with open(logfile, 'w') as log:
        try:
            tax_counts = []
            db_df = db_df.drop_duplicates(subset=["assembly_accession"], keep='first')

            if uploaded is False:
                log.write("INFO:working on downloaded refseq database entries\n")
                acc = lambda ids: ids.split('_')[2] + '_' + ids.split('_')[3]
                rbh_df['assembly_accession'] = rbh_df['sseqid']
                rbh_df['assembly_accession'] = rbh_df['assembly_accession'].map(acc)
                for query in rbh_df['qseqid'].unique():
                    log.write("\tINFO:extracting taxonomic information for {}\n".format(query))
                    df = rbh_df[rbh_df['qseqid'] == query]
                    if len(df) > 0:
                        df = df.drop_duplicates(subset=['assembly_accession'], keep='first')
                        m_df = db_df.merge(df, on='assembly_accession')
                        m_df = m_df.drop_duplicates(subset=['scomnames'], keep='first')
                        if taxonomic_unit + '_x' in m_df.columns:
                            values = m_df[taxonomic_unit + '_x'].value_counts()
                        else:
                            values = m_df[taxonomic_unit].value_counts()
                        values = values.rename(query)
                        tax_counts.append(values)

            # database is based on custom uploaded files --> assembly accessions are not accurate, therfore filtering is based on taxids and scomnames
            else:
                log.write("INFO:working on uploaded database entries\n")
                rbh_df.rename(columns={'staxids': 'species_taxid'}, inplace=True)
                db_df = db_df.drop_duplicates(subset=["species_taxid"], keep='first')

                for query in rbh_df['qseqid'].unique():
                    log.write("\tINFO:extracting taxonomic information for {}\n".format(query))
                    df = rbh_df[rbh_df['qseqid'] == query]
                    if len(df) > 0:
                        m_df = db_df.merge(df, on='species_taxid')
                        m_df = m_df.drop_duplicates(subset=['scomnames'], keep='first')
                        if taxonomic_unit + '_x' in m_df.columns:
                            values = m_df[taxonomic_unit + '_x'].value_counts()
                        else:
                            values = m_df[taxonomic_unit].value_counts()
                        values = values.rename(query)
                        tax_counts.append(values)
            log.write("DONE\n")
            return tax_counts
        except Exception as e:
            log.write("ERROR:{}\n".format(e))
            raise Exception("[-] ERROR in extract_taxonomic_information with exception: {}".format(e))

'''tax_counts_to_db_statistic_tables

    This function produces pandas dataframes with statistical information of the reciprocal best hits 
    and the specified taxonomic unit. It writes two CSV files in the relevant project directory.
    The first file inherits all hits for the specified taxonomic unit and the query sequences, the second
    dataframe inherits normalized values of the first table. Normalization is based on the number of 
    entries with the specified taxonomic unit within the database table.
    
    :param logfile
        :type str
    :param project_id
        :type int
    :param db_df
        :type pd.core.frame.DataFrame
    :param tax_counts
        :type list
    :param taxonomic_unit
        :type str
    
    :returns (df, normalized_df)
        :type (pd.core.frame.DataFrame, pd.core.frame.DataFrame)

'''
def tax_counts_to_db_statistic_tables(logfile: str, project_id: int, db_df: pd.core.frame.DataFrame,
                                      tax_counts: list, taxonomic_unit: str) -> tuple:
    with open(logfile, 'w') as log:
        try:
            log.write("INFO:starting statistical inference for {}\n".format(taxonomic_unit))
            if len(tax_counts) > 0:

                percentage = db_df['class'].value_counts()
                percentage = percentage.rename('Database')
                tax_counts.append(percentage)

                log.write("INFO:trying to write table with {} values for each query\n".format(taxonomic_unit))
                df = pd.DataFrame(tax_counts)
                df = df.fillna(0)
                df = df.transpose()
                df_path = 'media/blast_projects/' + str(project_id) + '/database_statistics.csv'
                df.to_csv(df_path)
                log.write("INFO:produced {} table\n".format(df_path))

                log.write("INFO:trying to normalize values\n")
                f = lambda x, y: (100 / y) * x
                normalizing_results = tax_counts
                normalizing_results.append(percentage)
                normalized_df = pd.DataFrame(normalizing_results)
                normalized_df = normalized_df.transpose().fillna(0)
                normalized_table = []

                for query in normalized_df.columns:
                    if query != 'Database':
                        temp_res = f(normalized_df[query], percentage)
                        temp_res.name = query
                        normalized_table.append(temp_res)

                log.write(
                    "INFO:trying to write normalization table with percentage of {} values for each query\n".format(
                        taxonomic_unit))
                normalized_df = pd.DataFrame(normalized_table)
                normalized_df = normalized_df.fillna(0)
                df_path = 'media/blast_projects/' + str(project_id) + '/database_statistics_normalized.csv'
                normalized_df.to_csv(df_path)
                log.write("INFO:produced {} table\n".format(df_path))
            else:
                log.write("INFO:Dataframes are empty\n")
                df = pd.DataFrame()
                normalized_df = pd.DataFrame()
            return df, normalized_df
        except Exception as e:
            log.write("ERROR:{}".format(e))
            raise Exception("[-] ERROR in tax_counts_to_db_statistics_tables function with exception {}".format(e))

#TODO extend this function to produce not only class based statistics
'''calculate_database_statistics

    This function calculates database statistics based on the reciprocal results of a finished pipeline project.
    This function depends on the database csv file, therefore it is not included within the Snakemake pipeline.

    :param project_id
        :type int
    :param logfile
        :type str
    :param user_email
        :type str
        
    :returns df
        :type pd.core.frame.DataFrame
'''
def calculate_database_statistics(project_id: int,logfile:str,user_email:str)->pd.core.frame.DataFrame:
    with open(logfile, 'w') as log:
        try:

            log.write("INFO:Starting to calculate database statistics\n")

            path_to_project = 'media/blast_projects/' + str(project_id)
            project = py_db_service.get_project_by_id(project_id)
            forward_db = project.project_forward_database
            path_to_db_csv = forward_db.path_to_database_file
            db_name = forward_db.database_name.replace(' ','_').upper()
            path_to_db_csv = path_to_db_csv + '/' + db_name
            new_database_name = path_to_db_csv + '_with_taxonomic_information.csv'

            path_to_reciprocal_results = path_to_project + '/reciprocal_results_with_taxonomy.csv'
            if isfile(path_to_reciprocal_results) is False or isfile(path_to_db_csv) is False:
                log.write("ERROR: {} or {} file not found\n".format(path_to_db_csv, path_to_reciprocal_results))
                raise FileNotFoundError

            result_data = pd.read_csv(path_to_reciprocal_results,index_col=0)

            if(isfile(new_database_name)):
                log.write("INFO:database: {} exists\n".format(new_database_name))
                db_df = pd.read_csv(new_database_name, index_col=0)
                log.write("INFO:Done loading result and database dataframe\n")
            else:
                db_df = pd.read_csv(path_to_db_csv,index_col=0)
                log.write("INFO:Done loading result and database dataframe\n")
                log.write("INFO:Trying to add taxonomic information to dataframe ...\n")
                add_tax_logfile=path_to_project+'/log/add_taxonomic_information_to_db.log'
                tax_df = add_taxonomic_information_to_db(user_email, add_tax_logfile, list(db_df['taxid'].unique()))
                tax_df['taxid'] = tax_df['taxid'].astype('int64')
                db_df = db_df.merge(tax_df, on='taxid')
                log.write("INFO:DONE with adding taxonomic information\n")
                log.write("INFO:trying to write database csv file into project directory\n")
                db_df.to_csv(new_database_name)
                log.write("INFO:DONE writing database csv file\n")

            log.write("INFO:starting function extract_taxonomic_information ...")
            logfile_tax_count_function=path_to_project+'log/extract_taxonomic_information.log'
            tax_counts=extract_taxonomic_information(logfile_tax_count_function, forward_db.uploaded_files, result_data, db_df, 'class')
            log.write("INFO:Done extracting list of taxonomic informations\n")

            log.write("INFO:Starting to produce output tables and graphs\n")
            logfile_tax_to_db_stat_function = path_to_project + 'log/tax_counts_to_db_statistics_tables.log'
            df, normalized_df = tax_counts_to_db_statistic_tables(logfile_tax_to_db_stat_function, project_id, db_df, tax_counts, 'class')
            log.write("DONE\n")
            return df

        except Exception as e:
            log.write("ERROR: problems during calculation of database statistics {}\n".format(e))
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
    :param logfile
        :type str
    :return returncode
        :type int
'''
def delete_database_statistics_task_and_output(project_id:int,logfile:str)->int:
    try:
        with open(logfile,'w') as log:
            log.write("INFO:trying to delete database statistics output and task\n")
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
                    log.write("INFO:removed: {}\n".format(path_to_database_statistics))
                for directory in listdir(path_to_project):
                    directory=path_to_project+'/'+directory
                    if isdir(directory):
                        path_to_database_statistics_image=directory+'/'+'database_result_statistics.png'
                        if(isfile(path_to_database_statistics_image)):
                            remove(path_to_database_statistics_image)
                            log.write("INFO:removed: {}\n".format(path_to_database_statistics_image))
            log.write("DONE\n")
            return 0
    except IntegrityError as e:
        log.write("ERROR: couldnt delete database statistics task with exception: {}\n".format(e))
        raise IntegrityError("ERROR during deletion of database statistics task result object with exception: {}".format(e))
    except Exception as e:
        log.write("ERROR: couldnt delete database statistics task with exception: {}\n".format(e))
        raise Exception("ERROR during removing database statistics output with exception: {}".format(e))




'''transform_normalized_database_table_to_json

'''
def transform_normalized_database_table_to_json(project_id:int):
    try:

        df_path = 'media/blast_projects/' + str(project_id) + '/database_statistics_normalized.csv'
        if isfile(df_path):
            table = pd.read_csv(df_path,header=0,index_col=0)
            json_records = table.reset_index().to_json(orient='records')
            data = json.loads(json_records)
            return data
        else:
            raise FileNotFoundError
    except FileNotFoundError:
        raise Exception("[-] File {} does not exist!".format('media/blast_projects/' + str(df_path)))
    except Exception as e:
        raise Exception("[-] Couldnt transform normalized database statistics to json with exception: {}".format(e))

