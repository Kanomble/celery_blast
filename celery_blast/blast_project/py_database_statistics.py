import json
from os import remove, listdir
from os.path import isfile, isdir

import altair as alt
# output_file-to save the layout in file, show-display the layout , output_notebook-to configure the default output state  to generate the output in jupytor notebook.
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from Bio import Entrez
from blast_project.models import BlastProject
from bokeh.core.enums import MarkerType
# output_file-to save the layout in file, show-display the layout , output_notebook-to configure the default output state  to generate the output in jupytor notebook.
from bokeh.io import output_file, save
# This is for creating layout
from bokeh.layouts import column, gridplot, row
# ColumnDataSource makes selection of the column easier and Select is used to create drop down
from bokeh.models import ColumnDataSource, Select, Spinner, MultiSelect, ColorPicker, RangeSlider, DataTable, \
    TableColumn
# To create intractive plot we need this to add callback method.
from bokeh.models import CustomJS, Legend, Button
# Figure objects have many glyph methods that can be used to draw vectorized graphical glyphs. example of glyphs-circle, line, scattter etc.
from bokeh.plotting import figure
from django.db import IntegrityError

plt.rcParams['legend.fontsize'] = 10
from bokeh.palettes import Spectral #inferno, viridis, magma,
# from matplotlib.ticker import MaxNLocator
from random import shuffle
from celery_blast.settings import BLAST_PROJECT_DIR

'''get_project_by_id

    Returns the associated BlastProject model instance.
    This function resides here due to circular import problems.
    
    :param project_id
        :type int
    
    :returns BlastProject
        :type blast_project.models.BlastProject
'''
def get_project_by_id(project_id):
    return BlastProject.objects.get(id=project_id)

'''add_taxonomic_information_to_db

    This function adds taxonomic information to the database dataframe.
    It uses the biopython entrez interface. First it eliminates all duplicate taxonomic ids, then 
    it parses the remaining taxonomic ids in chunks of 500 unique entries. Those batches of 500 unique entries
    are used within the Biopython Entrez.efetch function to query the taxonomy database. 
    After collecting all necessary information within the LineageEx element in the retrieved record, it checks the data
    for completeness. If the data is not complete, it adds an "unknown" field to the list which is used to build
    the taxonomy pandas dataframe that is returned.

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
                        log.write("\tWARNING:exception occured: {}\n".format(e))
                        if attempt == 9:
                            raise Exception

                    else:
                        for i in range(len(record)):
                            taxonomy.append(record[i]['ScientificName'])
                            #if taxonomic identifiers are renamed there is no field TaxId but AkaTaxIds
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
                                log.write("\tINFO:appending unknown for object: {} of record: {}\n".format(i,round(steps/step,1)))
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
            log.write("DONE\n")
            return tax_db
    except Exception as e:
        log.write("ERROR:problem during fetching and appending taxonomic information with exception: {}\n".format(e))
        raise Exception("ERROR during database table extension with taxonomic information with excpetion: {}".format(e))


'''extract_taxonomic_information

    This function loops over the queries in the reciprocal result dataframe. 
    Each query specific result dataframe is getting merged with the dataframe of the database that inherits the taxonomic information.
    The function filters the dataframe for unique entries within the specified taxonomic_unit and it returns a list with pandas series 
    as entries. The pandas series contains the number of unique occurrences of the specified taxonomic node and the respective query sequences.

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
        :type list[pd.core.series.Series]

'''
def extract_taxonomic_information(logfile: str, uploaded: bool, rbh_df: pd.core.frame.DataFrame,
                                  db_df: pd.core.frame.DataFrame, taxonomic_unit: str) -> list:
    with open(logfile, 'w') as log:
        try:
            tax_counts = []
            db_df = db_df.drop_duplicates(subset=["assembly_accession"], keep='first')
            log.write("INFO:starting to extract taxonomic information from reciprocal result dataframe ...\n")

            # database is based on downloaded refseq files therfore the assembly_accession column can be used to filter for duplicates
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
                    else:
                        raise Exception("After filtering, there are no database entries left.")

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
                    else:
                        raise Exception("After filtering, there are no database entries left.")
            log.write("DONE\n")
            return tax_counts
        except Exception as e:
            log.write("ERROR:{}\n".format(e))
            raise Exception("[-] ERROR in extract_taxonomic_information with exception: {}".format(e))

'''tax_counts_to_db_statistic_tables

    This function produces pandas dataframes with statistical information of the reciprocal best hits 
    and the specified taxonomic unit in comparison with the underlying database.
    It writes two CSV files in the relevant project directory.
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

                #counting the number of {taxonomic_unit} database entries
                percentage = db_df[taxonomic_unit].value_counts()
                percentage = percentage.rename('Database')
                tax_counts.append(percentage)

                log.write("INFO:trying to write table with {} values for each query\n".format(taxonomic_unit))
                df = pd.DataFrame(tax_counts)
                df = df.fillna(0)
                df = df.transpose()
                df_path = BLAST_PROJECT_DIR + str(project_id) + '/' + taxonomic_unit + '_database_statistics.csv'
                df.to_csv(df_path)
                log.write("INFO:produced {} table\n".format(df_path))

                #normalization based on number of entries in the database
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
                df_path = BLAST_PROJECT_DIR + str(project_id) + '/' + taxonomic_unit + '_database_statistics_normalized.csv'
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

'''calculate_database_statistics

    This function calculates database statistics based on the reciprocal results of a finished pipeline project.
    This function depends on the database csv file, therefore it is not included within the Snakemake pipeline.

    :param project_id
        :type int
    :param logfile
        :type str
    :param user_email
        :type str
    :param taxonomic_units
        :type list[str]
        
    :returns 0
        :type int
'''
def calculate_database_statistics(project_id: int,logfile:str,user_email:str, taxonomic_units:list)->int:
    with open(logfile, 'w') as log:
        try:
            log.write("INFO:defining taxonomic units")
            log.write("INFO:Starting to calculate database statistics\n")

            path_to_project = BLAST_PROJECT_DIR + str(project_id)
            project = get_project_by_id(project_id)
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

                #addition of taxonomic information via the biopython entrez module in add_taxonomic_information_to_db
                log.write("INFO:Trying to add taxonomic information to dataframe ...\n")
                add_tax_logfile=path_to_project+'/log/add_taxonomic_information_to_db.log'
                tax_df = add_taxonomic_information_to_db(user_email, add_tax_logfile, list(db_df['taxid'].unique()))
                tax_df['taxid'] = tax_df['taxid'].astype('int64')
                db_df = db_df.merge(tax_df, on='taxid')
                log.write("INFO:DONE with adding taxonomic information\n")

                log.write("INFO:trying to write database csv file into project directory\n")
                db_df.to_csv(new_database_name)
                log.write("INFO:DONE writing database csv file\n")

            #extraction of taxonomic information
            for taxonomic_unit in taxonomic_units:
                #result_data[taxonomic_unit].astype = str

                #path_to_project = BLAST_PROJECT_DIR + str(project_id)
                normalized_df_filepath = path_to_project + '/' + taxonomic_unit + '_database_statistics_normalized.csv'
                df_filepath = path_to_project + '/' + taxonomic_unit + '_database_statistics.csv'
                if isfile(normalized_df_filepath) is False or isfile(df_filepath) is False:
                    log.write("INFO:starting function extract_taxonomic_information ...\n")
                    logfile_tax_count_function = path_to_project + '/log/' + taxonomic_unit + '_extract_taxonomic_information.log'
                    tax_counts=extract_taxonomic_information(logfile_tax_count_function, forward_db.uploaded_files, result_data, db_df, taxonomic_unit)
                    log.write("INFO:Done extracting list of taxonomic informations\n")

                    log.write("INFO:Starting to produce output tables and graphs\n")
                    logfile_tax_to_db_stat_function = path_to_project + '/log/' + taxonomic_unit + '_tax_counts_to_db_statistics_tables.log'
                    df, normalized_df = tax_counts_to_db_statistic_tables(logfile_tax_to_db_stat_function, project_id, db_df, tax_counts, taxonomic_unit)
                    log.write("INFO:DONE extraction of {} database statistics\n".format(taxonomic_unit))
                else:
                    log.write("INFO:the database statistics dataframes do exist, skipping creation procedure\n")
                    #df = pd.read_csv(df_filepath,index_col=0,header=0)
                    normalized_df = pd.read_csv(normalized_df_filepath,index_col=0,header=0)

                #database_statistics_to_altair_plots(project_id: int, taxonomic_unit: str, full_df: pd.DataFrame, normalized_df: pd.DataFrame, logfile: str)
                log.write("INFO:starting to produce altair plots for database statistics dataframes\n")
                logfile_altair_plots = path_to_project + '/log/' + taxonomic_unit + '_database_statistics_to_altair_plots.log'
                database_statistics_to_altair_plots(project_id,taxonomic_unit,result_data,normalized_df,logfile_altair_plots)

            log.write("INFO:starting to produce interactive bokeh result plot\n")
            create_bokeh_plots(result_df=result_data, database=db_df,
                               taxonomic_unit='phylum',project_id=project_id)
            log.write("DONE\n")
            return 0

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
        project = get_project_by_id(project_id)
        if(project.project_database_statistics_task):
            task_status = str(project.project_database_statistics_task.status)
        else:
            task_status = "NOTEXEC"
        return task_status
    except Exception as e:
        raise Exception("ERROR during database statistics task status query with exception: {}".format(e))

#TODO remove also static plots
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
            project = get_project_by_id(project_id)
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

            #path_to_static_dir = "static/images/result_images/" + str(project_id) + "/"
            #path_to_altair_plot = path_to_static_dir + taxonomic_unit + "_altair_plot_normalized.html"
            log.write("DONE\n")
            return 0
    except IntegrityError as e:
        log.write("ERROR: couldnt delete database statistics task with exception: {}\n".format(e))
        raise IntegrityError("ERROR during deletion of database statistics task result object with exception: {}".format(e))
    except Exception as e:
        log.write("ERROR: couldnt delete database statistics task with exception: {}\n".format(e))
        raise Exception("ERROR during removing database statistics output with exception: {}".format(e))


'''transform_normalized_database_table_to_json
    
    This function reads the taxonomic_unit specific normalized database statistics dataframe from the project directory
    and returns a json object, that is used for feeding the datatables ajax calls in the database statistics dashboard. 
    
    :param project_id
        :type int
    :param taxonomic_unit
        :type str
        
    :returns data
        :type list[json_dicts]
'''
def transform_normalized_database_table_to_json(project_id:int,taxonomic_unit:str)->list:
    try:
        df_path = BLAST_PROJECT_DIR + str(project_id) + '/' + taxonomic_unit + '_database_statistics_normalized.csv'
        if isfile(df_path):
            table = pd.read_csv(df_path,header=0,index_col=0)
            number_queries = len(table.index)
            table = table.round(2)
            table = table.transpose()[(table == 0.0).sum() != number_queries]
            if(table.shape[1] > table.shape[0]):
                table = table.transpose()
            json_records = table.reset_index().to_json(orient='records')
            data = json.loads(json_records)
            return data
        else:
            raise FileNotFoundError
    except FileNotFoundError:
        raise Exception("[-] File {} does not exist!".format('media/blast_projects/' + str(df_path)))
    except Exception as e:
        raise Exception("[-] Couldnt transform normalized database statistics to json with exception: {}".format(e))

'''database_statistics_to_altair_plots
#OBSOLETE ??
'''
def database_statistics_to_altair_plots(project_id:int,taxonomic_unit:str,result_data:pd.DataFrame,normalized_df:pd.DataFrame,logfile:str):
    try:
        with open(logfile,'w') as log:
            path_to_static_dir = "static/images/result_images/" + str(project_id) + "/"
            log.write("INFO:checking if static dir: {} exists\n".format(path_to_static_dir))
            path_to_altair_plot = path_to_static_dir + taxonomic_unit + "_altair_plot_normalized.html"
            if isdir(path_to_static_dir):
                log.write("INFO:static directory exists, producing altair plots\n")
                number_queries = len(normalized_df.index)
                normalized_df = normalized_df.transpose()[(normalized_df == 0.0).sum() != number_queries]
                normalized_df = normalized_df.transpose()

                if len(normalized_df.columns) >= 2 and len(normalized_df.columns) < 10:

                    brush = alt.selection(type='interval')
                    bitscore = alt.Chart(result_data).mark_point().encode(
                        x='bitscore',
                        y='pident',
                        tooltip=['qseqid', 'sacc_transformed', 'pident', 'bitscore'],
                        color=alt.condition(brush, taxonomic_unit, alt.value('lightgray'))
                    ).add_selection(
                        brush
                    )

                    bar = alt.Chart(result_data).mark_bar().encode(
                        y=taxonomic_unit,
                        color=taxonomic_unit,
                        x='count('+taxonomic_unit+')'
                    ).transform_filter(
                        brush
                    )
                    chart = bitscore & bar
                    chart.save(path_to_altair_plot)
                else:
                    number_queries = len(normalized_df.index)
                    normalized_df = normalized_df.transpose()[(normalized_df == 0.0).sum() != number_queries]
                    normalized_df = normalized_df.transpose()
                    #TODO this plot is wrong -> qseqid is not correctly integrated
                    if len(normalized_df.columns) <= 15 and len(normalized_df.columns) >= 2:
                        altair_df = pd.melt(normalized_df)
                        altair_df.columns = [taxonomic_unit, "Relative # of RBHs"]
                        chart = alt.Chart(altair_df).mark_bar().encode(
                            x=taxonomic_unit,
                            y="Relative # of RBHs",
                            color=taxonomic_unit
                        )
                        chart.save(path_to_altair_plot)
                    else:
                        input_dropdown = alt.binding_select(options=list(normalized_df.columns), name=taxonomic_unit)
                        selection = alt.selection_single(fields=[taxonomic_unit], bind=input_dropdown)
                        color = alt.condition(selection,
                                              alt.Color(taxonomic_unit+':N', legend=None),
                                              alt.value('lightgray'))

                        chart = alt.Chart(result_data).mark_point().encode(
                            x='bitscore',
                            y='pident',
                            tooltip=['qseqid', 'sacc_transformed', 'pident', 'bitscore',taxonomic_unit],
                            color=color
                        ).add_selection(
                            selection
                        ).interactive(
                        )
                        chart.save(path_to_altair_plot)

                log.write("INFO:done saving database statistic altair plot for taxonomic unit: {}\n".format(taxonomic_unit))
            else:
                log.write("ERROR:static directory does not exist\n")
                raise IsADirectoryError(path_to_static_dir)
    except Exception as e:
        raise Exception("[-] ERROR in producing altair plots for database statistics with exception: {}".format(e))

################################## BOKEH INTERACTIVE RESULTS ###########################################################

'''create_bokeh_plots
    
    This function serves as hub for building the interactive bokeh result plot. It triggers the create_linked_bokeh_plot
    function that produces the interactive bokeh result plot. A logfile of this procedure is written into the corresponding
    log directory.
    
    :param result_df - the reciprocal_results_with_taxonomy.csv dataframe
        :type pd.DataFrame
    :param database - the database dataframe with taxonomic information
        :type pd.DataFrame
    :param taxonomic_unit - the taxonomic unit that is used for creating the first selection
        :type str
    :param project_id
        :type int

'''
def create_bokeh_plots(result_df:pd.DataFrame,database:pd.DataFrame,taxonomic_unit: str, project_id: int):
        try:
            path_to_project = BLAST_PROJECT_DIR + str(project_id)
            logfile_bokeh_plots = path_to_project + '/log/' + 'create_bokeh_plots.log'
            create_linked_bokeh_plot(logfile_bokeh_plots, result_df, database, taxonomic_unit, project_id)
        except Exception as e:
            raise Exception("[-] ERROR during creation of bokeh plots with exception : {}".format(e))

'''create_color_palette_selection_callback
    
    This function is used for creating a custom javascript callback for the color palette menu.
    The user can choose from a MultiSelect widget, which color palette should be used to render the circles. 
    It overwrites the "color" column of the currently selected ColumnDataSource items.
    
    :param curr- the current selection
        :type bokeh.models.ColumnDataSource
    :param color_menu - the current color_menu based on taxonomy e.g. Phylum ...
        :type bokeh.models.Select
    :param taxonomy_table_callback_dict - unique taxonomic entries within the defined taxonomic units (e.g. Phylum ...)
        :type dict['str'] = list[str, ...]
        
    :returns c_palette_callback
        :type bokeh.models.CustomJS
    
'''
def create_color_palette_selection_callback(curr: ColumnDataSource, color_menu: Select,
                                            taxonomy_table_callback_dict: dict) -> CustomJS:
    try:
        palettes = Spectral
        c_palette_callback = CustomJS(
            args=dict(sc=curr, color_menu=color_menu, tax_menu=taxonomy_table_callback_dict, pals=palettes), code="""
                                    // the callback value is a number 3,4,5,6,7,8,9,10,11,12
                                    var call_back_object = cb_obj.value

                                    var unique_organisms = []


                                    for(var i = 0; i<sc.get_length();i++){
                                        if(unique_organisms.includes(sc.data[color_menu.value][i]) == false){
                                            unique_organisms.push(sc.data[color_menu.value][i])
                                        }
                                    }

                                    if(unique_organisms.length <= pals[call_back_object].length){
                                        var color_dict = {}
                                        for(var i = 0; i < unique_organisms.length; i++){
                                            if(i == pals[call_back_object].length){
                                                i = 0
                                            }
                                            color_dict[unique_organisms[i]] = pals[call_back_object][i]
                                        }
                                        sc.data['color'] = []
                                        for(var i = 0; i<sc.data[color_menu.value].length;i++){
                                            sc.data['color'].push(color_dict[sc.data[color_menu.value][i]])
                                        }
                                     sc.change.emit();
                                     }

        """)
        return c_palette_callback
    except Exception as e:
        raise Exception("[-] ERROR couldnt create color palette selection callback with exception: {}".format(e))


'''create_color_palette_selection
    
    This function is used to pass the number of colored items to the custom javascript callback function for the
    colorpalette. The options are hardcoded based on the "Spectral" color palette. Spectral3 will color the circles in
    three different colors. The acutal coloring is done via changing the color column values in the curr ColumnDataSource via
    the create_color_palette_selection_callback function.
    
    :returns color_palette_menu
        :type bokeh.models.Select
'''
def create_color_palette_selection()->Select:
    try:

        options = [(str(val), "Spectral" + str(val)) for val in range(3, 12)]
        color_palette_menu = Select(options=options,
                                    value=str(3),
                                    title="Select a color palette")
        return color_palette_menu
    except Exception as e:
        raise Exception("[-] ERROR couldnt create color palette selection with exception: {}".format(e))

'''build_taxonomy_menu
    
    This function is used for creating the taxonomy MulitSelect menu for the different taxonomic units:
    Phylum, Class, Order, Family and Genus. For each of those items a menu is created. In order to pass selected items
    to callback functions, those menus are used as keys in a taxonomy menu dictionary. This allows accurate plot
    updates on certain changes. The selection is based on unique taxonomic entries within the bokeh_dataframe variable,
    which is an abstraction of the reciprocal_result_with_taxonomy.csv dataframe.    
    
    :param bokeh_dataframe
        :type pd.DataFrame
    :param taxonomic_unit
        :type str
    
    :returns tax_menu
        :type bokeh.models.MultiSelect
'''
def build_taxonomy_menu(bokeh_dataframe: pd.DataFrame, taxonomic_unit: str)->MultiSelect:
    try:
        unique_tax = list(bokeh_dataframe[taxonomic_unit].unique())
        if len(unique_tax) > 1:
            tax_menu = MultiSelect(options=unique_tax, value=[unique_tax[0], unique_tax[1]],
                                   title='Select: ' + taxonomic_unit.capitalize())
        else:
            tax_menu = MultiSelect(options=unique_tax, value=[unique_tax[0]],
                                   title='Select: ' + taxonomic_unit.capitalize())

        return tax_menu
    except Exception as e:
        raise Exception("[-] ERROR creating taxonomy menu for bokeh plot with exception: {}".format(e))

'''build_json_callback_for_taxonomy

    This function produces the custom javascript callback function for the taxonomy. It is similar to the 
    build_json_callback_for_taxonomy function in external_tools.py_cdd_domain_search.py but slightly more complex
    as there are more connected data sources.
    Based on this callback function the parameters column_dat and table_dat are getting changed.
    
    :param column_dat - represents the selected data
        :type bokeh.models.ColumnDataSource 
    :param static_dat - represents the full data
        :type bokeh.models.ColumnDataSource
    :param table_dat - represents the data used within the Bokeh DataTable
        :type bokeh.models.ColumnDataSource
    :param domains - list of the CDD accession in the bokeh DataTable (cdd_dataframe)
        :type list[str,...]
    :param taxonomic_unit - taxonomic unit where this callback is attached to
        :type str
    :param tax_selection - current taxonomic selection of the user: 
                           e.g. tax_selection_dict = {'class':class_menu,'order':order_menu,'family':family_menu,
                           'genus':genus_menu}
        :type dict[str] = bokeh.models.MultiSelect
    :param menu_qseqid = current query sequence selection
        :type bokeh.models.MultiSelect
    :param xaxis_menu
        :type bokeh.models.Select
    :param yaxis_menu
        :type bokeh.models.Select
    :param color_menu
        :type bokeh.models.Select
    :param color_dict - this variable will hold the color value for the taxonomy key (e.g. cyanobacteria = blue)
        :type dict
    
    :returns tax_menu_callback
        :type bokeh.models.CustomJS
'''
def build_json_callback_for_taxonomy(column_dat: ColumnDataSource, static_dat: ColumnDataSource,
                                     table_dat: ColumnDataSource, taxonomic_unit: str, tax_selection: dict,
                                     menu_qseqid: MultiSelect, xaxis_menu: Select, yaxis_menu: Select,
                                     color_menu: Select, color_dict: dict,
                                     taxonomy_table_callback_dict: dict) -> CustomJS:
    tax_menu_callback = CustomJS(args=dict(sc=column_dat,
                                           source=static_dat,
                                           tax_unit=taxonomic_unit,
                                           selected_taxonomy=tax_selection,
                                           menu_qseqids=menu_qseqid,
                                           table_data=table_dat,
                                           xaxis_menu=xaxis_menu, yaxis_menu=yaxis_menu,
                                           color_menu=color_menu, color_dict=color_dict,
                                           tax_dict=taxonomy_table_callback_dict), code="""

                        var call_back_object = cb_obj.value 
                        var tab_dict = {};
                        var tab_dict_static = {};
                        var tab_dict_org_counter = {};

                        for(var i = 0;i<table_data.get_length();i++){
                            tab_dict[table_data.data['# TaxName'][i]] = 0
                            tab_dict_org_counter[table_data.data['# TaxName'][i]] = 0
                            tab_dict_static[table_data.data['# TaxName'][i]] = table_data.data['# Different Organisms In DB'][i]                                                  
                        }


                        var unique_class = []
                        var unique_order = []
                        var unique_family = []
                        var unique_genus = []

                        let keys = Object.keys(sc.data)
                        for(var i = 0; i < keys.length; i++){
                            sc.data[keys[i]] = []
                        }

                        var taxid_arr = []
                        for(var i = 0; i < source.get_length(); i++){
                            for(var j = 0; j < call_back_object.length; j++){  
                                for(var k = 0; k < menu_qseqids.value.length; k++){
                                    if(source.data['qseqid'][i] == menu_qseqids.value[k]){
                                        if(source.data[tax_unit][i] == call_back_object[j]){

                                            if(unique_order.includes(source.data['order'][i]) == false){
                                                    unique_order.push(source.data['order'][i])
                                                }

                                            if(unique_class.includes(source.data['class'][i]) == false){
                                                unique_class.push(source.data['class'][i])
                                            }


                                            if(unique_family.includes(source.data['family'][i]) == false){
                                                unique_family.push(source.data['family'][i])
                                            }

                                            if(unique_genus.includes(source.data['genus'][i]) == false){
                                                unique_genus.push(source.data['genus'][i])
                                            }

                                            for(var x = 0; x < keys.length; x++){
                                                if((keys[x] != 'x') && (keys[x] != 'y') && (keys[x] != 'color')){
                                                    sc.data[keys[x]].push(source.data[keys[x]][i])

                                                }
                                            }

                                            sc.data['color'].push(color_dict[source.data[color_menu.value][i]])
                                            sc.data['x'].push(source.data[xaxis_menu.value][i])
                                            sc.data['y'].push(source.data[yaxis_menu.value][i])
                                            //tab_dict[call_back_object[j]]+=1
                                            for(var l = 0;l<tax_dict[color_menu.value].length;l++){
                                                if(source.data[color_menu.value][i] == tax_dict[color_menu.value][l]){
                                                    tab_dict[tax_dict[color_menu.value][l]]+=1
                                                    //hier muss ich wissen ob color_menu.value zu source.data[tax_unit][i] gehört
                                                    if(taxid_arr.includes(source.data['staxids'][i]) == false){
                                                        taxid_arr.push(source.data['staxids'][i])
                                                        tab_dict_org_counter[tax_dict[color_menu.value][l]]+=1
                                                    }                                           
                                                }                                           
                                            }


                                        }            
                                    }
                                }
                            }
                        }

                        for(var key in selected_taxonomy) {
                            if(key == 'order'){
                                selected_taxonomy[key].options = unique_order
                                selected_taxonomy[key].value = unique_order                                                        
                            }

                            if(key == 'class'){
                                selected_taxonomy[key].options = unique_class
                                selected_taxonomy[key].value = unique_class                                                        
                            }

                            if(key == 'family'){
                                selected_taxonomy[key].options = unique_family
                                selected_taxonomy[key].value = unique_family                                                        
                            }

                            if(key == 'genus'){
                                selected_taxonomy[key].options = unique_genus
                                selected_taxonomy[key].value = unique_genus                                                        
                            }
                        }

                        table_data.data['#RBHs'] = []
                        table_data.data['# Different Organisms In DB'] = []
                        table_data.data['# TaxName'] = []
                        table_data.data['# Different Organisms In Selection'] = []
                        table_data.data['index'] = []
                        var counter = 1
                        for(let key in tab_dict){
                            for(var l = 0;l<tax_dict[color_menu.value].length;l++){
                                if(key == tax_dict[color_menu.value][l]){
                                    table_data.data['#RBHs'].push(tab_dict[key])
                                    table_data.data['# Different Organisms In DB'].push(tab_dict_static[key])
                                    table_data.data['# TaxName'].push(key)
                                    table_data.data['# Different Organisms In Selection'].push(tab_dict_org_counter[key])
                                    table_data.data['index'].push(counter)
                                    counter += 1                               
                                }
                            }
                        }

                        table_data.change.emit();
                        sc.change.emit();
                        """)
    return tax_menu_callback


'''create_y_axis_menu
    
    
'''
def create_y_axis_menu(circle, axis, data_column):
    y_axis_menu = Select(options=['bitscore', 'pident', 'evalue', 'slen'],  # ,'slen'
                         value='pident',
                         title="Select Y axis elements")

    y_axis_menu_callback = CustomJS(args=dict(gl=circle, plot=axis, data=data_column), code='''
           var call_back_object = cb_obj.value;
           data.data['y'] = data.data[call_back_object]
           plot[1].axis_label = call_back_object;
           data.change.emit();

    ''')

    y_axis_menu.js_on_change('value', y_axis_menu_callback)
    return y_axis_menu

'''create_x_axis_menu

'''
def create_x_axis_menu(circle, axis, data_column):
    x_axis_menu = Select(options=['bitscore', 'pident', 'evalue', 'slen'],  # ,'slen'
                         value='bitscore',
                         title="Select X axis elements")

    x_axis_menu_callback = CustomJS(args=dict(gl=circle, plot=axis, data=data_column), code='''
           var call_back_object = cb_obj.value;
           data.data['x'] = data.data[call_back_object]
           plot[0].axis_label = call_back_object;
           data.change.emit();

    ''')

    x_axis_menu.js_on_change('value', x_axis_menu_callback)
    return x_axis_menu

'''create_color_and_marker_dictionaries_for_bokeh_dataframe

'''
def create_color_and_marker_dictionaries_for_bokeh_dataframe(result_data: pd.DataFrame) -> tuple:
    try:
        # prepare distinct colors for the specified taxonomic unit
        color_dict = {}
        for tax_unit in ['phylum', 'order', 'class', 'family', 'genus']:
            num_colors = len(result_data[tax_unit].unique())

            if num_colors > 256:
                clrs = sns.color_palette('pastel', n_colors=num_colors)
                clrs = clrs.as_hex()
                color_dict.update(dict(zip(result_data[tax_unit].unique(), clrs)))

            else:
                clrs = sns.color_palette('pastel', n_colors=num_colors)
                clrs = clrs.as_hex()
                color_dict.update(dict(zip(result_data[tax_unit].unique(), clrs)))  # magma(n)

        # prepare custom marker for each query sequence
        marker_dict = {}
        marker = list(MarkerType)
        # just use colorable marker types
        for m in ["x", "y", "dot", "dash", "cross", "asterisk"]:
            marker.remove(m)
        shuffle(marker)

        for i, query in enumerate(result_data['qseqid'].unique()):
            marker_dict[query] = marker[i % len(marker)]

        return color_dict, marker_dict
    except Exception as e:
        raise Exception("[-] ERROR creating marker and color data for RBH result plot with exception: {}".format(e))


'''create_qseqid_menu_callback

'''
def create_qseqid_menu_callback(Overall: ColumnDataSource, Curr: ColumnDataSource, DbData: ColumnDataSource,
                                taxonomic_unit: str, tax_menu: MultiSelect, xaxis_menu: MultiSelect,
                                yaxis_menu: MultiSelect,
                                color_menu: Select, color_dict: dict, taxonomy_menus: list) -> CustomJS:
    try:
        menu_qseqid_callback = CustomJS(args=dict(source=Overall, sc=Curr, table_data=DbData,
                                                  tax_unit=taxonomic_unit,
                                                  xaxis_menu=xaxis_menu, yaxis_menu=yaxis_menu,
                                                  color_menu=color_menu, color_dict=color_dict,
                                                  taxonomy_menus=taxonomy_menus), code="""
        var tab_dict = {};
        var tab_dict_static = {};
        var tab_dict_org_count = {};
        for(var i = 0;i<table_data.get_length();i++){
            tab_dict[table_data.data['# TaxName'][i]] = 0
            tab_dict_org_count[table_data.data['# TaxName'][i]] = 0
            tab_dict_static[table_data.data['# TaxName'][i]] = table_data.data['# Different Organisms In DB'][i]
        }
        var call_back_object = cb_obj.value


        let keys = Object.keys(sc.data)
        for(var i = 0; i < keys.length; i++){
            sc.data[keys[i]] = []
        }

        var taxid_arr = []
        for(var i = 0; i < source.get_length(); i++){
            for(var j = 0; j < call_back_object.length; j++){


                if(source.data['qseqid'][i] == call_back_object[j]){

                     if(taxonomy_menus[0].value.includes(source.data['phylum'][i]) == true){
                        if(taxonomy_menus[1].value.includes(source.data['class'][i]) == true){
                            if(taxonomy_menus[2].value.includes(source.data['order'][i]) == true){
                                if(taxonomy_menus[3].value.includes(source.data['family'][i]) == true){
                                     if(taxonomy_menus[4].value.includes(source.data['genus'][i]) == true){
                                       for(var x = 0; x < keys.length; x++){
                                            if((keys[x] != 'x') && (keys[x] != 'y') && (keys[x] != 'color')){
                                                sc.data[keys[x]].push(source.data[keys[x]][i])
                                            }
                                        }

                                        sc.data['color'].push(color_dict[source.data[color_menu.value][i]])
                                        sc.data['x'].push(source.data[xaxis_menu.value][i])
                                        sc.data['y'].push(source.data[yaxis_menu.value][i])
                                        tab_dict[source.data[tax_unit][i]]+=1
                                        if(taxid_arr.includes(source.data['staxids'][i]) == false){
                                            taxid_arr.push(source.data['staxids'][i])
                                            tab_dict_org_count[source.data[tax_unit][i]]+=1
                                        }
                                     }                               
                                }
                            }
                        }
                     }





                }
            }
        }
        table_data.data['#RBHs'] = []
        table_data.data['# Different Organisms In DB'] = []
        table_data.data['# TaxName'] = []
        table_data.data['# Different Organisms In Selection'] = []
        table_data.data['index'] = []
        var counter = 1
        for(let key in tab_dict){
            console.log(key)
            if(key != '# TaxName'){
                table_data.data['#RBHs'].push(tab_dict[key])
                table_data.data['# Different Organisms In DB'].push(tab_dict_static[key])
                table_data.data['# TaxName'].push(key)
                table_data.data['# Different Organisms In Selection'].push(tab_dict_org_count[key])
                table_data.data['index'].push(counter)
                counter += 1            
            }
        }

        console.log(table_data.data)
        table_data.change.emit();
        sc.change.emit();
        """)
        return menu_qseqid_callback
    except Exception as e:
        raise Exception("[-] ERROR creating the custom js callback for the qseqid menu with exception: {}".format(e))


'''create_initial_bokeh_database_data
    
'''
def create_initial_bokeh_database_data(database: pd.DataFrame, data_selection: pd.DataFrame, taxcount_df: pd.DataFrame,
                                       taxonomic_unit: str)->pd.DataFrame:
    try:
        # unique database entries
        db_df = pd.DataFrame(database[taxonomic_unit].value_counts())

        selection = pd.DataFrame(data_selection[taxonomic_unit].value_counts())

        db_df.columns = ['value']
        selection.columns = ['value']

        db_df[taxonomic_unit] = db_df.index
        selection[taxonomic_unit] = selection.index

        db_df.index = range(len(db_df))
        selection.index = range(len(selection))

        db_df = selection.merge(db_df, on=taxonomic_unit, how='outer')
        db_df = taxcount_df.merge(db_df, on=taxonomic_unit, how='outer')
        db_df = db_df[['value_x', 'value_y', 'value', taxonomic_unit]]
        db_df.columns = ['#RBHs', '# Different Organisms In DB', '# Different Organisms In Selection',
                         '# TaxName']

        db_df = db_df.fillna(0)

        return db_df
    except Exception as e:
        raise Exception("[-] ERROR creating database dataframe for bokeh RBH result plot with exception: {}".format(e))


'''create_initial_bokeh_result_data

'''
def create_initial_bokeh_result_data(result_data: pd.DataFrame, taxonomic_unit: str) -> tuple:
    try:
        # RBH result dataframe
        result_data = result_data.loc[:,
                      ['order', 'class', 'phylum', 'genus', 'family', 'bitscore', 'pident', 'stitle', 'scomnames',
                       'staxids', 'qseqid',
                       'sacc_transformed', 'slen']]  # ,'slen'
        result_data = result_data.sort_values(by=taxonomic_unit)

        color_dict, marker_dict = create_color_and_marker_dictionaries_for_bokeh_dataframe(result_data)

        # lambda functions for adding color and marker columns
        create_color_scheme = lambda value: color_dict[value]
        create_marker_scheme = lambda value: marker_dict[value]

        result_data['x'] = result_data['bitscore']
        result_data['y'] = result_data['pident']
        result_data['color'] = result_data[taxonomic_unit].apply(create_color_scheme)
        result_data['marker'] = result_data['qseqid'].apply(create_marker_scheme)

        return result_data, color_dict

    except Exception as e:
        raise Exception("[-] ERROR creating result dataframe for bokeh RBH result plot with exception: {}".format(e))


'''create_initial_bokeh_data_selection

'''
def create_initial_bokeh_data_selection(result_data: pd.DataFrame, taxonomic_unit: str):
    try:
        unique_tax = list(result_data[taxonomic_unit].unique())
        unique_qseqids = list(result_data['qseqid'].unique())

        if len(unique_tax) > 1:
            data_selection = result_data[
                (result_data[taxonomic_unit] == unique_tax[0]) | (result_data[taxonomic_unit] == unique_tax[1])
                ]
        else:
            data_selection = result_data[result_data[taxonomic_unit] == unique_tax[0]]  # prepare table dataframe

        if len(unique_qseqids) > 1:
            data_selection = data_selection[
                (data_selection['qseqid'] == unique_qseqids[0]) | (data_selection['qseqid'] == unique_qseqids[1])
                ]
        else:
            data_selection = data_selection[data_selection['qseqid'] == unique_qseqids[0]]

        taxcount_df = pd.DataFrame(data_selection.staxids.value_counts())
        taxcount_df['value'] = taxcount_df['staxids']
        taxcount_df['staxids'] = taxcount_df.index
        taxcount_df.index = pd.Index(range(len(taxcount_df)))
        taxid_to_taxonomic_unit = lambda taxid: \
        data_selection[data_selection.staxids == taxid][taxonomic_unit].unique()[0]
        taxcount_df[taxonomic_unit] = taxcount_df.staxids.apply(taxid_to_taxonomic_unit)
        taxcount_df = pd.DataFrame(taxcount_df[taxonomic_unit].value_counts())

        taxcount_df.columns = ['value']
        taxcount_df[taxonomic_unit] = taxcount_df.index
        taxcount_df.index = range(len(taxcount_df))

        return data_selection, taxcount_df
    except Exception as e:
        raise Exception(
            "[-] ERROR creating initial result dataframe selection for bokeh RBH result plot with exception: {}".format(
                e))

'''create_linked_bokeh_plot

    Main function for the interactive bokeh reciprocal result plots.

'''
def create_linked_bokeh_plot(logfile: str, result_data: pd.DataFrame, database: pd.DataFrame, taxonomic_unit: str,
                               project_id: int) -> int:
    try:
        with open(logfile, 'w') as log:
            #initializing output directories
            path_to_static_dir = "static/images/result_images/" + str(project_id) + "/"
            log.write("INFO:checking if static dir: {} exists\n".format(path_to_static_dir))
            #no static saving until now
            path_to_static_bokeh_plot = path_to_static_dir + "interactive_bokeh_plot.html"
            path_to_project_dir = BLAST_PROJECT_DIR + str(project_id) + '/' + "interactive_bokeh_plot.html"

            # create bokeh dataframes for plots and tables
            data_all, color_dict = create_initial_bokeh_result_data(result_data, taxonomic_unit)
            # selection subset for initial plot data
            data_selection, taxcount_df = create_initial_bokeh_data_selection(data_all, taxonomic_unit)
            db_df = create_initial_bokeh_database_data(database, data_selection, taxcount_df, taxonomic_unit)

            data_selection_phylum, taxcount_df_phylum = create_initial_bokeh_data_selection(data_all, 'phylum')
            db_df_phylum = create_initial_bokeh_database_data(database, data_selection, taxcount_df_phylum, 'phylum')
            data_selection_class, taxcount_df_class = create_initial_bokeh_data_selection(data_all, 'class')
            db_df_class = create_initial_bokeh_database_data(database, data_selection, taxcount_df_class, 'class')
            data_selection_order, taxcount_df_order = create_initial_bokeh_data_selection(data_all, 'order')
            db_df_order = create_initial_bokeh_database_data(database, data_selection, taxcount_df_order, 'order')
            data_selection_family, taxcount_df_family = create_initial_bokeh_data_selection(data_all, 'family')
            db_df_family = create_initial_bokeh_database_data(database, data_selection, taxcount_df_family, 'family')
            data_selection_genus, taxcount_df_genus = create_initial_bokeh_data_selection(data_all, 'genus')
            db_df_genus = create_initial_bokeh_database_data(database, data_selection, taxcount_df_genus, 'genus')

            # setup bokeh classes
            Overall = ColumnDataSource(data=data_all)
            Curr = ColumnDataSource(data=data_selection)
            DbData = ColumnDataSource(data=db_df)

            table_data_dict = {
                'phylum': ColumnDataSource(data=db_df_phylum),
                'class': ColumnDataSource(data=db_df_class),
                'order': ColumnDataSource(data=db_df_order),
                'family': ColumnDataSource(data=db_df_family),
                'genus': ColumnDataSource(data=db_df_genus)
            }

            # unique_tax = list(data_all[taxonomic_unit].unique())
            unique_qseqids = list(data_all['qseqid'].unique())

            # selection subset for initial plot data


            table = DataTable(source=DbData, width=390, height=275,
                              sizing_mode="scale_both", reorderable=True, sortable=True, fit_columns=True,
                              columns=[
                                  TableColumn(field='#RBHs', title='#RBHs'),
                                  TableColumn(field='# Different Organisms In DB', title='# Different Organisms In DB'),
                                  TableColumn(field='# Different Organisms In Selection',
                                              title='# Different Organisms In Selection'),
                                  TableColumn(field='# TaxName', title='# TaxName'),
                              ])

            qseq_values = []
            if len(unique_qseqids) > 1:
                qseq_values.append(unique_qseqids[0])
                qseq_values.append(unique_qseqids[1])
            else:
                qseq_values.append(unique_qseqids[0])

            menu_qseqids = MultiSelect(options=unique_qseqids, value=qseq_values,
                                       title='Select target query sequence')  # drop down menu

            # menu_callback = create_taxonomy_menu_callback(Overall,
            #                                     Curr,
            #                                     DbData,
            #                                     taxonomic_unit
            #                                     ,menu_qseqids)

            TOOLTIPS = [
                ("stitle", "@stitle"),
                ("bitscore,pident", "@bitscore, @pident"),
                ("sacc RBH to qseqid", "@sacc_transformed RBH to @qseqid "),
                ("scomname", "@scomnames"),
            ]

            # x_range=(0, result_data['bitscore'].max() + result_data['bitscore'].min())
            p = figure(x_axis_label='bitscore', y_axis_label='pident',
                       plot_height=700, plot_width=900,
                       tooltips=TOOLTIPS,
                       tools="lasso_select, reset,save, box_zoom,undo,redo,wheel_zoom, pan",
                       title="Reciprocal Best Hits Result Data",
                       )  # ,tools="box_select, reset" creating figure object

            p.add_layout(Legend(), 'left')

            circle = p.scatter(x='x', y='y', color='color', marker='marker', size=10, line_width=1,
                               line_color='black',
                               source=Curr, legend_field=taxonomic_unit)  # plotting the data using glyph circle

            p.legend.glyph_width = 40
            p.legend.glyph_height = 40

            color_menu = Select(options=['phylum', 'class', 'order', 'family', 'genus'],
                                value=taxonomic_unit, title="Select Legend Color")
            color_callback = CustomJS(
                args=dict(legend=p.legend.items[0], sc=Curr, source=Overall, color_dict=color_dict,
                          table_data=DbData, table_data_dict=table_data_dict, menu_qseqids=menu_qseqids), code='''

                                var tax_unit = cb_obj.value

                                var tab_dict = {};
                                var tab_dict_org_count = {}
                                var tab_dict_static = {}; //numbers dont change
                                for(var i = 0;i<table_data_dict[tax_unit].get_length();i++){
                                    tab_dict[table_data_dict[tax_unit].data['# TaxName'][i]] = 0
                                    tab_dict_org_count[table_data_dict[tax_unit].data['# TaxName'][i]] = 0
                                    tab_dict_static[table_data_dict[tax_unit].data['# TaxName'][i]] = table_data_dict[tax_unit].data['# Different Organisms In DB'][i]
                                }

                                legend.label = {'field':tax_unit}
                                var length = sc.get_length()
                                sc.data['color']=[]
                                for(var i = 0; i < length; i++){
                                    sc.data['color'].push(color_dict[sc.data[tax_unit][i]])
                                }

                                var taxid_arr = []
                                for(var i = 0; i < source.get_length(); i++){
                                    for(var k = 0; k < menu_qseqids.value.length; k++){
                                        if(source.data['qseqid'][i] == menu_qseqids.value[k]){
                                            if(taxid_arr.includes(source.data['staxids'][i]) == false){
                                                taxid_arr.push(source.data['staxids'][i])
                                                tab_dict_org_count[source.data[tax_unit][i]]+=1
                                            }
                                            tab_dict[source.data[tax_unit][i]]+=1
                                        }
                                    }
                                }

                                table_data.data['#RBHs'] = []
                                table_data.data['# Different Organisms In DB'] = []
                                table_data.data['# TaxName'] = []
                                table_data.data['# Different Organisms In Selection'] = []
                                table_data.data['index'] = []
                                var counter = 1
                                for(let key in tab_dict){
                                    table_data.data['#RBHs'].push(tab_dict[key])
                                    table_data.data['# Different Organisms In DB'].push(tab_dict_static[key])
                                    table_data.data['# TaxName'].push(key)
                                    table_data.data['# Different Organisms In Selection'].push(tab_dict_org_count[key])
                                    table_data.data['index'].push(counter)
                                    counter += 1
                                }
                                table_data.change.emit();
                                sc.change.emit();
            ''')
            color_menu.js_on_change('value', color_callback)

            phylum_menu = build_taxonomy_menu(data_all, 'phylum')
            class_menu = build_taxonomy_menu(data_all, 'class')
            order_menu = build_taxonomy_menu(data_all, 'order')
            family_menu = build_taxonomy_menu(data_all, 'family')
            genus_menu = build_taxonomy_menu(data_all, 'genus')

            unique_phylum = list(result_data['phylum'].unique())
            unique_class = list(result_data['class'].unique())
            unique_order = list(result_data['order'].unique())
            unique_family = list(result_data['family'].unique())
            unique_genus = list(result_data['genus'].unique())

            taxonomy_table_callback_dict = {
                'phylum': unique_phylum,
                'class': unique_class,
                'order': unique_order,
                'family': unique_family,
                'genus': unique_genus
            }

            tax_menus = [phylum_menu, class_menu, order_menu, family_menu, genus_menu]

            x_axis_menu = create_x_axis_menu(circle, p.axis, Curr)
            y_axis_menu = create_y_axis_menu(circle, p.axis, Curr)

            tax_selection_dict = {'class': class_menu, 'order': order_menu, 'family': family_menu, 'genus': genus_menu}
            phylum_menu_callback = build_json_callback_for_taxonomy(Curr, Overall, DbData, 'phylum',
                                                                    tax_selection_dict, menu_qseqids,
                                                                    x_axis_menu, y_axis_menu, color_menu,
                                                                    color_dict, taxonomy_table_callback_dict)
            phylum_menu.js_on_change('value', phylum_menu_callback)

            tax_selection_dict = {'order': order_menu, 'family': family_menu, 'genus': genus_menu}
            class_menu_callback = build_json_callback_for_taxonomy(Curr, Overall, DbData, 'class', tax_selection_dict,
                                                                   menu_qseqids, x_axis_menu, y_axis_menu, color_menu,
                                                                   color_dict, taxonomy_table_callback_dict)
            class_menu.js_on_change('value', class_menu_callback)

            tax_selection_dict = {'family': family_menu, 'genus': genus_menu}
            order_menu_callback = build_json_callback_for_taxonomy(Curr, Overall, DbData, 'order', tax_selection_dict,
                                                                   menu_qseqids, x_axis_menu, y_axis_menu,
                                                                   color_menu, color_dict, taxonomy_table_callback_dict)
            order_menu.js_on_change('value', order_menu_callback)

            tax_selection_dict = {'genus': genus_menu}
            family_menu_callback = build_json_callback_for_taxonomy(Curr, Overall, DbData, 'family', tax_selection_dict,
                                                                    menu_qseqids, x_axis_menu, y_axis_menu,
                                                                    color_menu, color_dict,
                                                                    taxonomy_table_callback_dict)
            family_menu.js_on_change('value', family_menu_callback)

            genus_menu_callback = build_json_callback_for_taxonomy(Curr, Overall, DbData, 'genus', {},
                                                                   menu_qseqids, x_axis_menu, y_axis_menu,
                                                                   color_menu, color_dict, taxonomy_table_callback_dict)
            genus_menu.js_on_change('value', genus_menu_callback)

            menu_qseqid_callback = create_qseqid_menu_callback(Overall,
                                                               Curr,
                                                               DbData,
                                                               taxonomic_unit,
                                                               phylum_menu,
                                                               x_axis_menu, y_axis_menu,
                                                               color_menu, color_dict, tax_menus)

            menu_qseqids.js_on_change('value', menu_qseqid_callback)

            range_slider = RangeSlider(start=0, end=result_data['bitscore'].max() + result_data['bitscore'].min(),
                                       value=(result_data['bitscore'].min(), result_data['bitscore'].max()), step=1,
                                       title="Bitscore Range Slider")

            circle_size_spinner = Spinner(title="Circle size",
                                          low=0, high=60, step=5,
                                          value=circle.glyph.size,
                                          width=200
                                          )

            line_size_spinner = Spinner(title="Circle line size",
                                        low=0, high=20, step=1,
                                        value=circle.glyph.line_width,
                                        width=200
                                        )
            line_color_picker = ColorPicker(color='black', title="Line Color")

            range_slider.js_link("value", p.x_range, "start", attr_selector=0)
            range_slider.js_link("value", p.x_range, "end", attr_selector=1)

            line_size_spinner.js_link("value", circle.glyph, "line_width")
            circle_size_spinner.js_link("value", circle.glyph, "size")
            line_color_picker.js_link('color', circle.glyph, 'line_color')

            selection_callback = CustomJS(
                args=dict(sc=Curr, source=Overall, table_data=DbData,
                          qseqids=menu_qseqids, color_menu=color_menu), code="""
                var call_back_object = cb_obj.indices
                var tax_unit = color_menu.value
                var tab_dict = {};
                var tab_dict_static = {};
                var tab_dict_org_counter = {};
                for(var i = 0;i<table_data.get_length();i++){
                    tab_dict[table_data.data['# TaxName'][i]] = 0
                    tab_dict_org_counter[table_data.data['# TaxName'][i]] = 0
                    tab_dict_static[table_data.data['# TaxName'][i]] = table_data.data['# Different Organisms In DB'][i]
                }
                var taxid_arr = []
                for(var i = 0; i < call_back_object.length; i++){
                    tab_dict[sc.data[tax_unit][call_back_object[i]]]+=1
                    if(taxid_arr.includes(sc.data['staxids'][[call_back_object[i]]]) == false){
                        taxid_arr.push(sc.data['staxids'][[call_back_object[i]]])
                        tab_dict_org_counter[sc.data[tax_unit][[call_back_object[i]]]]+=1
                    }  
                }
                table_data.data['#RBHs'] = []
                table_data.data['# Different Organisms In DB'] = []
                table_data.data['# TaxName'] = []
                table_data.data['# Different Organisms In Selection'] = []
                for(let key in tab_dict){
                    table_data.data['#RBHs'].push(tab_dict[key])
                    table_data.data['# Different Organisms In DB'].push(tab_dict_static[key])
                    table_data.data['# TaxName'].push(key)
                    table_data.data['# Different Organisms In Selection'].push(tab_dict_org_counter[key])
                }
                table_data.change.emit();
                """)

            Curr.selected.js_on_change('indices', selection_callback)

            download_selection_callback = CustomJS(args=dict(sc=Curr, tax_unit=taxonomic_unit), code="""
                var temp = []
                var csvFileData = []
                for(var i = 0; i < sc.selected.indices.length; i++){
                    temp = [sc.data['qseqid'][sc.selected.indices[i]],
                            sc.data['sacc_transformed'][sc.selected.indices[i]],
                            sc.data['staxids'][sc.selected.indices[i]]]
                    csvFileData.push(temp)
                }
                //define the heading for each row of the data  
                var csv = `qseqid,sacc,staxids\n`;  
                //merge the data with CSV  
                csvFileData.forEach(function(row) {  
                        csv += row.join(',');  
                        csv += `\n`;  
                });  
                var json = JSON.stringify(csv);
                var blob = new Blob([csv], {type: "octet/stream"});
                var url  = window.URL.createObjectURL(blob);
                window.location.assign(url);
            """)

            download_selection_button = Button(label="Download Selection")
            download_selection_button.js_on_click(download_selection_callback)

            color_palette = create_color_palette_selection()
            color_palette_callback = create_color_palette_selection_callback(Curr, color_menu,
                                                                             taxonomy_table_callback_dict)
            color_palette.js_on_change('value', color_palette_callback)

            grid = gridplot([[column(p),
                              column(menu_qseqids, row(circle_size_spinner, line_size_spinner),
                                     range_slider, download_selection_button, table, x_axis_menu, y_axis_menu),
                              column(phylum_menu, class_menu, order_menu, family_menu, genus_menu, color_menu,
                                     color_palette)]],
                            toolbar_location='right')

            output_file(filename=path_to_project_dir,
                      title="Reciprocal BLAST Interactive Result Plot".format(
                           taxonomic_unit))
            save(grid)

        return 0
    except Exception as e:
        raise Exception("ERROR in producing bokeh plots for database statistics with exception: {}".format(e))



