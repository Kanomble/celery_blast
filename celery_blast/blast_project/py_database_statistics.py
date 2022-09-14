from blast_project import py_django_db_services as py_db_service
import pandas as pd
import matplotlib.pyplot as plt
from os.path import isfile, isdir
from os import remove, listdir
from django.db import IntegrityError
from Bio import Entrez

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
                                        log.write("\tINFO: AkaTaxIds detected: {}".format(akaid))
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
            return tax_db, taxid
    except Exception as e:
        log.write("ERROR:problem during fetching and appending taxonomic information with exception: {}\n".format(e))
        raise Exception("ERROR during database table extension with taxonomic information with excpetion: {}".format(e))


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
    try:
        with open(logfile,'w') as log:
            log.write("INFO:Starting to calculate database statistics\n")

            path_to_project='media/blast_projects/' + str(project_id)

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
                db_df = pd.read_csv(new_database_name, index_col=0)
                log.write("INFO:Done loading result and database dataframe\n")
            else:
                db_df=pd.read_csv(path_to_db_csv,index_col=0)
                log.write("INFO:Done loading result and database dataframe\n")
                log.write("INFO:Trying to add taxonomic information to dataframe ...\n")
                add_tax_logfile=path_to_project+'/log/add_taxonomic_information_to_db.log'
                tax_df=add_taxonomic_information_to_db(user_email, add_tax_logfile, list(db_df['taxid'].unique()))
                tax_df['taxid'] = tax_df['taxid'].astype('int64')
                db_df = db_df.merge(tax_df, on='taxid')
                log.write("INFO:DONE with adding taxonomic information\n")

                log.write("INFO:trying to write database csv file into project directory\n")

                db_df.to_csv(new_database_name)
                log.write("INFO:DONE writing database csv file\n")

            #database is based on refseq genome files --> filtering based on assembly_accessions and scomnames
            if forward_db.uploaded_files is False:

                log.write("INFO:working on downloaded refseq database entries\n")

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
                log.write("INFO:working on uploaded database entries\n")
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

            log.write("INFO:Done merging database and result dataframes\n")
            log.write("INFO:Starting to produce output tables and graphs\n")
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
                log.write("INFO:Dataframe is empty\n")
                df=pd.DataFrame()
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

