from blast_project import py_django_db_services as py_db_service
import pandas as pd
import altair as alt
from os.path import isfile, isdir
from os import remove, listdir
from django.db import IntegrityError
from Bio import Entrez
import json
import seaborn as sns
#output_file-to save the layout in file, show-display the layout , output_notebook-to configure the default output state  to generate the output in jupytor notebook.
from bokeh.io import output_file, save
#ColumnDataSource makes selection of the column easier and Select is used to create drop down
from bokeh.models import ColumnDataSource, Spinner, MultiSelect, ColorPicker, RangeSlider, DataTable, TableColumn
#Figure objects have many glyph methods that can be used to draw vectorized graphical glyphs. example of glyphs-circle, line, scattter etc.
from bokeh.plotting import figure
#To create intractive plot we need this to add callback method.
from bokeh.models import CustomJS, Button
#This is for creating layout
from bokeh.layouts import column, gridplot, row
from bokeh.core.enums import MarkerType
from random import shuffle

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

                #counting the number of {taxonomic_unit} database entries
                percentage = db_df[taxonomic_unit].value_counts()
                percentage = percentage.rename('Database')
                tax_counts.append(percentage)

                log.write("INFO:trying to write table with {} values for each query\n".format(taxonomic_unit))
                df = pd.DataFrame(tax_counts)
                df = df.fillna(0)
                df = df.transpose()
                df_path = 'media/blast_projects/' + str(project_id) + '/' + taxonomic_unit + '_database_statistics.csv'
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
                df_path = 'media/blast_projects/' + str(project_id) + '/' + taxonomic_unit + '_database_statistics_normalized.csv'
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

                path_to_project = 'media/blast_projects/' + str(project_id)
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
                    df = pd.read_csv(df_filepath,index_col=0,header=0)
                    normalized_df = pd.read_csv(normalized_df_filepath,index_col=0,header=0)

                #database_statistics_to_altair_plots(project_id: int, taxonomic_unit: str, full_df: pd.DataFrame, normalized_df: pd.DataFrame, logfile: str)
                log.write("INFO:starting to produce altair plots for database statistics dataframes\n")
                logfile_altair_plots = path_to_project + '/log/' + taxonomic_unit + '_database_statistics_to_altair_plots.log'
                database_statistics_to_altair_plots(project_id,taxonomic_unit,result_data,normalized_df,logfile_altair_plots)

                create_bokeh_plots(result_df=result_data,database=db_df,normalized_df=normalized_df,
                                   taxonomic_unit=taxonomic_unit,project_id=project_id)
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
        project = py_db_service.get_project_by_id(project_id)
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
        df_path = 'media/blast_projects/' + str(project_id) + '/' + taxonomic_unit + '_database_statistics_normalized.csv'
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

#TODO documentation
def create_bokeh_plots(result_df:pd.DataFrame,database:pd.DataFrame,normalized_df:pd.DataFrame,taxonomic_unit: str, project_id: int):
        try:
            number_queries = len(normalized_df.index)
            normalized_df = normalized_df.transpose()[(normalized_df == 0.0).sum() != number_queries]
            normalized_df = normalized_df.transpose()

            path_to_project = 'media/blast_projects/' + str(project_id)
            length_database = len(normalized_df.columns)
            logfile_bokeh_plots = path_to_project + '/log/' + taxonomic_unit + '_database_statistics_to_bokeh_plots.log'
            if length_database <= 15:
                create_linked_bokeh_plot(logfile_bokeh_plots,result_df,database,taxonomic_unit,project_id)
            else:
                create_unlinked_bokeh_plot(logfile_bokeh_plots,result_df, database,taxonomic_unit,project_id)
        except Exception as e:
            raise Exception("[-] ERROR during creation of bokeh plots with exception : {}".format(e))

#TODO documentation
def create_linked_bokeh_plot(logfile: str, result_data: pd.DataFrame, database: pd.DataFrame, taxonomic_unit: str,
                             project_id: int):
    try:
        with open(logfile, 'w') as log:

            path_to_static_dir = "static/images/result_images/" + str(project_id) + "/"
            log.write("INFO:checking if static dir: {} exists\n".format(path_to_static_dir))
            path_to_bokeh_plot = path_to_static_dir + taxonomic_unit + "_bokeh_plot.html"

            data_all = result_data.loc[:,
                       [taxonomic_unit, 'bitscore', 'pident', 'stitle', 'scomnames', 'staxids', 'qseqid',
                        'sacc_transformed']]
            data_all = data_all.sort_values(by=taxonomic_unit)
            num_colors = len(data_all[taxonomic_unit].unique())
            clrs = sns.color_palette('pastel', n_colors=num_colors)
            clrs = clrs.as_hex()
            color_dict = dict(zip(data_all[taxonomic_unit].unique(), clrs))
            create_color_scheme = lambda value: color_dict[value]

            data_all['color'] = data_all[taxonomic_unit].apply(create_color_scheme)
            data_selection = data_all.copy()

            Overall = ColumnDataSource(data=data_all)
            Curr = ColumnDataSource(data=data_selection)

            df = pd.DataFrame(data_all.value_counts(taxonomic_unit))
            df[taxonomic_unit] = df.index
            df.index = [x for x in range(0, len(df.values))]
            df.columns = ["value", taxonomic_unit]
            df['color'] = df[taxonomic_unit].apply(create_color_scheme)

            dbStatistics = pd.DataFrame(
                database[database[taxonomic_unit].isin(data_all[taxonomic_unit].unique())][
                    taxonomic_unit].value_counts())
            dbStatistics.columns = ["value"]
            dbStatistics[taxonomic_unit] = dbStatistics.index
            dbStatistics.index = [x for x in range(0, len(dbStatistics.values))]
            dbStatistics['color'] = dbStatistics[taxonomic_unit].apply(create_color_scheme)

            dbData = ColumnDataSource(data=dbStatistics)
            bCurr = ColumnDataSource(data=df)

            unique_qseqids = list(data_all['qseqid'].unique())
            menu_qseqids = MultiSelect(options=unique_qseqids, value=unique_qseqids,
                                       title='Select target query sequence')  # drop down menu
            unique_orders = list(data_all[taxonomic_unit].unique())
            menu = MultiSelect(options=unique_orders, value=unique_orders,
                               title='Select: ' + taxonomic_unit.capitalize())  # drop down menu

            callback = CustomJS(args=dict(
                source=Overall,
                sc=Curr,
                sbcDb=dbData,
                sbc=bCurr,
                color_code=color_dict,
                tax_unit=taxonomic_unit,
                menu_qseqids=menu_qseqids), code="""
                            var call_back_object = cb_obj.value
                            sc.data['bitscore']=[]
                            sc.data['pident']=[]
                            sc.data['color']=[]
                            sc.data[tax_unit]=[]
                            sc.data['index']=[]
                            sc.data['stitle']=[]
                            sc.data['scomnames']=[]
                            sc.data['staxids']=[]
                            sc.data['qseqid']=[]
                            sc.data['sacc_transformed']=[]
                            var dict = new Object()
                            var dict = {}
                            var indexDict = {}
                            var taxid_count = {}

                            for(var j = 0; j < call_back_object.length; j++){
                                dict[call_back_object[j]]=0
                                taxid_count[call_back_object[j]]=[]
                                indexDict[call_back_object[j]]=j+1
                            }


                            for(var i = 0; i < source.get_length(); i++){
                                for(var j = 0; j < call_back_object.length; j++){
                                    for(var k = 0; k < menu_qseqids.value.length; k++){
                                        if(source.data['qseqid'][i] == menu_qseqids.value[k]){
                                            if(source.data[tax_unit][i] == call_back_object[j]){
                                                sc.data['bitscore'].push(source.data['bitscore'][i])
                                                sc.data['pident'].push(source.data['pident'][i])
                                                sc.data['color'].push(source.data['color'][i])
                                                sc.data[tax_unit].push(source.data[tax_unit][i])
                                                sc.data['index'].push(source.data['index'][i])
                                                sc.data['stitle'].push(source.data['stitle'][i])
                                                sc.data['scomnames'].push(source.data['scomnames'][i])
                                                sc.data['staxids'].push(source.data['staxids'][i])
                                                sc.data['qseqid'].push(source.data['qseqid'][i])
                                                sc.data['sacc_transformed'].push(source.data['sacc_transformed'][i])
                                                dict[call_back_object[j]]+=1
                                                if(taxid_count[call_back_object[j]].includes(source.data['staxids'][i]) == false){
                                                        taxid_count[call_back_object[j]].push(source.data['staxids'][i])
                                                    }
                                            }                                        
                                        }
                                    }

                                }
                            }

                            sbc.data['value']=[]
                            sbc.data[tax_unit]=[]
                            sbc.data['color']=[]
                            sbc.data['index']=[]
                            for(let key in dict) {
                                sbc.data['value'].push(dict[key])
                                sbc.data[tax_unit].push(key)
                                sbc.data['color'].push(color_code[key])
                                sbc.data['index'].push(indexDict[key])
                            }

                            sbcDb.data['value']=[]
                            sbcDb.data[tax_unit]=[]
                            sbcDb.data['color']=[]
                            sbcDb.data['index']=[]
                            for(let key in taxid_count) {
                                sbcDb.data['value'].push(taxid_count[key].length)
                                sbcDb.data[tax_unit].push(key)
                                sbcDb.data['color'].push(color_code[key])
                                sbcDb.data['index'].push(indexDict[key])
                            }


                            sbcDb.change.emit();
                            sbc.change.emit();
                            sc.change.emit();
                            """)

            qseqid_callback = CustomJS(args=dict(
                source=Overall,
                sc=Curr,
                sbcDb=dbData,
                sbc=bCurr,
                color_code=color_dict,
                tax_unit=taxonomic_unit,
                menu_tax=menu), code="""
                            var call_back_object = cb_obj.value
                            sc.data['bitscore']=[]
                            sc.data['pident']=[]
                            sc.data['color']=[]
                            sc.data[tax_unit]=[]
                            sc.data['index']=[]
                            sc.data['stitle']=[]
                            sc.data['scomnames']=[]
                            sc.data['staxids']=[]
                            sc.data['qseqid']=[]
                            sc.data['sacc_transformed']=[]
                            var dict = new Object()
                            var dict = {}
                            var indexDict = {}
                            var taxid_count = {}

                            for(var j = 0; j < menu_tax.value.length; j++){
                                dict[menu_tax.value[j]]=0
                                taxid_count[menu_tax.value[j]]=[]
                                indexDict[menu_tax.value[j]]=j+1
                            }


                            for(var i = 0; i < source.get_length(); i++){
                                for(var j = 0; j < call_back_object.length; j++){
                                    for(var k = 0; k < menu_tax.value.length; k++){
                                        if(source.data[tax_unit][i] == menu_tax.value[k]){
                                            if(source.data['qseqid'][i] == call_back_object[j]){
                                                sc.data['bitscore'].push(source.data['bitscore'][i])
                                                sc.data['pident'].push(source.data['pident'][i])
                                                sc.data['color'].push(source.data['color'][i])
                                                sc.data[tax_unit].push(source.data[tax_unit][i])
                                                sc.data['index'].push(source.data['index'][i])
                                                sc.data['stitle'].push(source.data['stitle'][i])
                                                sc.data['scomnames'].push(source.data['scomnames'][i])
                                                sc.data['staxids'].push(source.data['staxids'][i])
                                                sc.data['qseqid'].push(source.data['qseqid'][i])
                                                sc.data['sacc_transformed'].push(source.data['sacc_transformed'][i])
                                                dict[menu_tax.value[k]]+=1
                                                if(taxid_count[menu_tax.value[k]].includes(source.data['staxids'][i]) == false){
                                                        taxid_count[menu_tax.value[k]].push(source.data['staxids'][i])
                                                    }
                                            }                                        
                                        }
                                    }

                                }
                            }

                            sbc.data['value']=[]
                            sbc.data[tax_unit]=[]
                            sbc.data['color']=[]
                            sbc.data['index']=[]
                            for(let key in dict) {
                                sbc.data['value'].push(dict[key])
                                sbc.data[tax_unit].push(key)
                                sbc.data['color'].push(color_code[key])
                                sbc.data['index'].push(indexDict[key])
                            }

                            sbcDb.data['value']=[]
                            sbcDb.data[tax_unit]=[]
                            sbcDb.data['color']=[]
                            sbcDb.data['index']=[]
                            for(let key in taxid_count) {
                                sbcDb.data['value'].push(taxid_count[key].length)
                                sbcDb.data[tax_unit].push(key)
                                sbcDb.data['color'].push(color_code[key])
                                sbcDb.data['index'].push(indexDict[key])
                            }


                            sbcDb.change.emit();
                            sbc.change.emit();
                            sc.change.emit();
                            """)
            TOOLTIPS = [
                ("stitle", "@stitle"),
                ("bitscore,pident", "@bitscore, @pident"),
                ("sacc RBH to qseqid", "@sacc_transformed RBH to @qseqid "),
                ("scomname", "@scomnames"),
            ]


            p = figure(x_axis_label='bitscore', y_axis_label='pident',
                       plot_height=850, plot_width=700,
                       tooltips=TOOLTIPS,
                       tools="box_select, reset, box_zoom, pan",
                       title="Reciprocal Best Hit - percent identity vs bitscore")  # ,tools="box_select, reset" creating figure object
            circle = p.circle(x='bitscore', y='pident', color='color', size=5, line_width=1, line_color='black',
                              source=Curr)
            menu.js_on_change('value', callback)
            menu_qseqids.js_on_change('value', qseqid_callback)
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

            b = figure(x_range=unique_orders, height=250, width=700,
                       title="Number of selected RBHs in the clade: {}".format(taxonomic_unit))
            b.vbar(x=taxonomic_unit, top='value', color='color', width=0.9, source=bCurr)

            b2 = figure(x_range=unique_orders, height=250, width=700,
                        title="Number of different organisms within the {} clade".format(taxonomic_unit))
            b2.vbar(x=taxonomic_unit, top='value', color='color', width=0.9, source=dbData)
            b.xaxis.major_label_orientation = "vertical"
            b2.xaxis.major_label_orientation = "vertical"

            bar_callback = CustomJS(args=dict(sc=Curr,
                                              sbc=bCurr,
                                              sbcDb=dbData,
                                              color_code=color_dict,
                                              tax_unit=taxonomic_unit,
                                              menu_qseqids=menu_qseqids), code="""
                                                var call_back_object = cb_obj.indices
                                                var dict = {}
                                                var taxid_count = {}
                                                var arr = []
                                                var indexDict = {}
                                                for(var j = 0; j < call_back_object.length; j++){
                                                    if(arr.includes(sc.data[tax_unit][call_back_object[j]]) == false){
                                                        arr.push(sc.data[tax_unit][call_back_object[j]])
                                                    }
                                                }

                                                for(var j = 0; j < arr.length; j++){
                                                    dict[arr[j]]=0
                                                    taxid_count[arr[j]]=[]
                                                    indexDict[arr[j]]=j+1
                                                }

                                                for(var i = 0; i < sc.get_length(); i++){
                                                    for(var j = 0; j < arr.length; j++){
                                                        for(var k = 0; k < menu_qseqids.value.length; k++){
                                                            if(sc.data['qseqid'][call_back_object[i]] == menu_qseqids.value[k])
                                                                if(arr[j] == sc.data[tax_unit][call_back_object[i]]){
                                                                    dict[arr[j]]+=1
                                                                    if(taxid_count[arr[j]].includes(sc.data['staxids'][call_back_object[i]]) == false){
                                                                        taxid_count[arr[j]].push(sc.data['staxids'][call_back_object[i]])
                                                                    }

                                                                }                                                        
                                                        }
                                                    }
                                                }

                                                sbc.data['value']=[]
                                                sbc.data[tax_unit]=[]
                                                sbc.data['color']=[]
                                                sbc.data['index']=[]
                                                for(let key in dict) {
                                                    sbc.data['value'].push(dict[key])
                                                    sbc.data[tax_unit].push(key)
                                                    sbc.data['color'].push(color_code[key])
                                                    sbc.data['index'].push(indexDict[key])
                                                }


                                                sbcDb.data['value']=[]
                                                sbcDb.data[tax_unit]=[]
                                                sbcDb.data['color']=[]
                                                sbcDb.data['index']=[]
                                                for(let key in taxid_count) {
                                                    sbcDb.data['value'].push(taxid_count[key].length)
                                                    sbcDb.data[tax_unit].push(key)
                                                    sbcDb.data['color'].push(color_code[key])
                                                    sbcDb.data['index'].push(indexDict[key])
                                                }


                                                sbcDb.change.emit();
                                                sbc.change.emit();
                                            """)
            Curr.selected.js_on_change('indices', bar_callback)

            line_size_spinner.js_link("value", circle.glyph, "line_width")
            circle_size_spinner.js_link("value", circle.glyph, "size")
            line_color_picker.js_link('color', circle.glyph, 'line_color')

            download_selection_callback = CustomJS(args=dict(sc=Curr), code="""
                var downloadable_items = []
                for(var i = 0; i < sc.selected.indices.length; i++){
                    downloadable_items.push(sc.data['sacc_transformed'][sc.selected.indices[i]])
                }

                var json = JSON.stringify(downloadable_items);
                var blob = new Blob([json],{type: "octet/stream"});
                var url = window.URL.createObjectURL(blob);
                window.location.assign(url);
            """)

            download_selection_button = Button(label="Download Selection")
            download_selection_button.js_on_click(download_selection_callback)

            grid = gridplot([[column(p), column(b, b2, menu, menu_qseqids, circle_size_spinner, download_selection_button)]],
                            toolbar_location='right', sizing_mode="stretch_both", merge_tools=True)
            output_file(filename=path_to_bokeh_plot, title="Interactive Graph Percent Identity vs. Bitscore linked to {} database entries".format(taxonomic_unit))
            save(grid)
        return 0
    except Exception as e:
        raise Exception("ERROR in producing bokeh plots for database statistics with exception: {}".format(e))


def create_unlinked_bokeh_plot(logfile:str,result_data: pd.DataFrame,database:pd.DataFrame, taxonomic_unit: str, project_id:int)->int:
    try:
        with open(logfile,'w') as log:
            path_to_static_dir = "static/images/result_images/" + str(project_id) + "/"
            # log.write("INFO:checking if static dir: {} exists\n".format(path_to_static_dir))
            path_to_bokeh_plot = path_to_static_dir + taxonomic_unit + "_bokeh_plot.html"

            # read database entries with taxonomy
            db_df = pd.DataFrame(database[taxonomic_unit].value_counts())
            db_df.columns = ['value']
            db_df[taxonomic_unit] = db_df.index
            db_df.index = range(len(db_df))

            data_all = result_data.loc[:,
                       [taxonomic_unit, 'bitscore', 'pident', 'stitle', 'scomnames', 'staxids', 'qseqid',
                        'sacc_transformed']]
            data_all = data_all.sort_values(by=taxonomic_unit)

            # prepare distinct colors for the specified taxonomic unit
            num_colors = len(data_all[taxonomic_unit].unique())
            clrs = sns.color_palette('pastel', n_colors=num_colors)
            clrs = clrs.as_hex()
            color_dict = dict(zip(data_all[taxonomic_unit].unique(), clrs))

            # prepare custom marker for each query sequence
            marker_dict = {}
            marker = list(MarkerType)
            for m in ["x", "y", "dot", "dash", "cross", "asterisk"]:
                marker.remove(m)
            shuffle(marker)

            for i, query in enumerate(data_all['qseqid'].unique()):
                marker_dict[query] = marker[i % len(marker)]

            create_color_scheme = lambda value: color_dict[value]
            create_marker_scheme = lambda value: marker_dict[value]

            # plot and the menu is linked with each other by this callback function
            unique_tax = list(data_all[taxonomic_unit].unique())
            unique_qseqids = list(data_all['qseqid'].unique())

            data_all['color'] = data_all[taxonomic_unit].apply(create_color_scheme)
            data_all['marker'] = data_all['qseqid'].apply(create_marker_scheme)

            # selection subset for initial plot data
            if len(unique_tax) > 1:

                data_selection = data_all[
                    (data_all[taxonomic_unit] == unique_tax[0]) | (data_all[taxonomic_unit] == unique_tax[1])]
                menu = MultiSelect(options=unique_tax, value=[unique_tax[0], unique_tax[1]],
                                   title='Select: ' + taxonomic_unit.capitalize())  # drop down menu

            else:
                data_selection = data_all[data_all[taxonomic_unit] == unique_tax[0]]
                menu = MultiSelect(options=unique_tax, value=[unique_tax[0]],
                                   title='Select: ' + taxonomic_unit.capitalize())  # drop down menu

            # prepare table dataframe
            taxcount_df = pd.DataFrame(data_selection.staxids.value_counts())
            taxcount_df['value'] = taxcount_df['staxids']
            taxcount_df['staxids'] = taxcount_df.index
            taxcount_df.index = pd.Index(range(len(taxcount_df)))
            # data_selection.merge(taxcount_df,on='staxids', how='outer')
            taxid_to_taxonomic_unit = lambda taxid: data_selection[data_selection.staxids == taxid][taxonomic_unit].unique()[0]
            taxcount_df[taxonomic_unit] = taxcount_df.staxids.apply(taxid_to_taxonomic_unit)
            taxcount_df = pd.DataFrame(taxcount_df[taxonomic_unit].value_counts())

            db_df = pd.DataFrame(database[taxonomic_unit].value_counts())
            selection = pd.DataFrame(data_selection[taxonomic_unit].value_counts())

            db_df.columns = ['value']
            selection.columns = ['value']
            taxcount_df.columns = ['value']

            db_df[taxonomic_unit] = db_df.index
            selection[taxonomic_unit] = selection.index
            taxcount_df[taxonomic_unit] = taxcount_df.index

            db_df.index = range(len(db_df))
            selection.index = range(len(selection))
            taxcount_df.index = range(len(taxcount_df))

            db_df = selection.merge(db_df, on=taxonomic_unit, how='outer')
            db_df = taxcount_df.merge(db_df, on=taxonomic_unit, how='outer')
            db_df = db_df[['value_x', 'value_y', 'value', taxonomic_unit]]
            db_df.columns = ['#RBHs', '# Different Organisms In DB', '# Different Organisms In Selection',
                             taxonomic_unit]
            db_df = db_df.fillna(0)
            # setup bokeh classes
            Overall = ColumnDataSource(data=data_all)
            Curr = ColumnDataSource(data=data_selection)
            DbData = ColumnDataSource(data=db_df)

            table = DataTable(source=DbData, width=390, height=275,
                              sizing_mode="scale_both", reorderable=True, sortable=True, fit_columns=True,
                              columns=[
                                  TableColumn(field='#RBHs', title='#RBHs'),
                                  TableColumn(field='# Different Organisms In DB', title='# Different Organisms In DB'),
                                  TableColumn(field='# Different Organisms In Selection',
                                              title='# Different Organisms In Selection'),
                                  TableColumn(field=taxonomic_unit, title=taxonomic_unit.capitalize()),
                              ])

            menu_qseqids = MultiSelect(options=unique_qseqids, value=unique_qseqids,
                                       title='Select target query sequence')  # drop down menu

            menu_callback = CustomJS(args=dict(source=Overall, sc=Curr, table_data=DbData,
                                               tax_unit=taxonomic_unit, menu_qseqids=menu_qseqids), code="""

            var tab_dict = {};
            var tab_dict_org_count = {}
            var tab_dict_static = {}; //numbers dont change
            for(var i = 0;i<table_data.get_length();i++){
                tab_dict[table_data.data[tax_unit][i]] = 0
                tab_dict_org_count[table_data.data[tax_unit][i]] = 0
                tab_dict_static[table_data.data[tax_unit][i]] = table_data.data['# Different Organisms In DB'][i]
            }


            var call_back_object = cb_obj.value
            console.log("menu_qseqids: ", menu_qseqids.value)
            sc.data['bitscore']=[]
            sc.data['pident']=[]
            sc.data['color']=[]
            sc.data['marker']=[]
            sc.data[tax_unit]=[]
            sc.data['index']=[]
            sc.data['stitle']=[]
            sc.data['scomnames']=[]
            sc.data['staxids']=[]
            sc.data['qseqid']=[]
            sc.data['sacc_transformed']=[]

            var taxid_arr = []
            for(var i = 0; i < source.get_length(); i++){
                for(var j = 0; j < call_back_object.length; j++){  
                    for(var k = 0; k < menu_qseqids.value.length; k++){
                        if(source.data['qseqid'][i] == menu_qseqids.value[k]){
                            if(source.data[tax_unit][i] == call_back_object[j]){
                                sc.data['bitscore'].push(source.data['bitscore'][i])
                                sc.data['pident'].push(source.data['pident'][i])
                                sc.data['color'].push(source.data['color'][i])
                                sc.data['marker'].push(source.data['marker'][i])
                                sc.data[tax_unit].push(source.data[tax_unit][i])
                                sc.data['index'].push(source.data['index'][i])
                                sc.data['stitle'].push(source.data['stitle'][i])
                                sc.data['scomnames'].push(source.data['scomnames'][i])
                                sc.data['staxids'].push(source.data['staxids'][i])
                                sc.data['qseqid'].push(source.data['qseqid'][i])
                                sc.data['sacc_transformed'].push(source.data['sacc_transformed'][i])
                                tab_dict[call_back_object[j]]+=1
                                if(taxid_arr.includes(source.data['staxids'][i]) == false){
                                    taxid_arr.push(source.data['staxids'][i])
                                    tab_dict_org_count[source.data[tax_unit][i]]+=1
                                }

                            }            
                        }
                    }
                }
            }

            table_data.data['#RBHs'] = []
            table_data.data['# Different Organisms In DB'] = []
            table_data.data[tax_unit] = []
            table_data.data['# Different Organisms In Selection'] = []
            for(let key in tab_dict){
                table_data.data['#RBHs'].push(tab_dict[key])
                table_data.data['# Different Organisms In DB'].push(tab_dict_static[key])
                table_data.data[tax_unit].push(key)
                table_data.data['# Different Organisms In Selection'].push(tab_dict_org_count[key])
            }


            table_data.change.emit();
            sc.change.emit();
            """)

            TOOLTIPS = [
                ("stitle", "@stitle"),
                ("bitscore,pident", "@bitscore, @pident"),
                ("sacc RBH to qseqid", "@sacc_transformed RBH to @qseqid "),
                ("scomname", "@scomnames"),
            ]


            p = figure(x_axis_label='bitscore', y_axis_label='pident',
                       plot_height=700, plot_width=700,
                       tooltips=TOOLTIPS,
                       tools="box_select, reset, box_zoom, pan", title="Number of RBHs - pident vs bitscore",
                       x_range=(0, result_data['bitscore'].max() + result_data['bitscore'].min())
                       )  # ,tools="box_select, reset" creating figure object
            circle = p.scatter(x='bitscore', y='pident', color='color', marker='marker', size=5, line_width=1,
                               line_color='black',
                               source=Curr, legend_field=taxonomic_unit)  # plotting the data using glyph circle

            menu.js_on_change('value', menu_callback)  # calling the function on change of selection

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

            menu_qseqid_callback = CustomJS(args=dict(source=Overall, sc=Curr, table_data=DbData,
                                                      tax_unit=taxonomic_unit, tax_selection=menu), code="""


            var tab_dict = {};
            var tab_dict_static = {};
            var tab_dict_org_count = {};
            for(var i = 0;i<table_data.get_length();i++){
                tab_dict[table_data.data[tax_unit][i]] = 0
                tab_dict_org_count[table_data.data[tax_unit][i]] = 0
                tab_dict_static[table_data.data[tax_unit][i]] = table_data.data['# Different Organisms In DB'][i]
            }

            var call_back_object = cb_obj.value
            sc.data['bitscore']=[]
            sc.data['pident']=[]
            sc.data['color']=[]
            sc.data['marker']=[]
            sc.data[tax_unit]=[]
            sc.data['index']=[]
            sc.data['stitle']=[]
            sc.data['scomnames']=[]
            sc.data['staxids']=[]
            sc.data['qseqid']=[]
            sc.data['sacc_transformed']=[]

            var taxid_arr = []

            for(var i = 0; i < source.get_length(); i++){
                for(var j = 0; j < call_back_object.length; j++){
                    for(var k = 0; k < tax_selection.value.length; k++){
                        if(source.data[tax_unit][i] == tax_selection.value[k]){
                            if(source.data['qseqid'][i] == call_back_object[j]){
                                sc.data['bitscore'].push(source.data['bitscore'][i])
                                sc.data['pident'].push(source.data['pident'][i])
                                sc.data['color'].push(source.data['color'][i])
                                sc.data['marker'].push(source.data['marker'][i])
                                sc.data[tax_unit].push(source.data[tax_unit][i])
                                sc.data['index'].push(source.data['index'][i])
                                sc.data['stitle'].push(source.data['stitle'][i])
                                sc.data['scomnames'].push(source.data['scomnames'][i])
                                sc.data['staxids'].push(source.data['staxids'][i])
                                sc.data['qseqid'].push(source.data['qseqid'][i])
                                sc.data['sacc_transformed'].push(source.data['sacc_transformed'][i])
                                tab_dict[tax_selection.value[k]]+=1
                                if(taxid_arr.includes(source.data['staxids'][i]) == false){
                                    taxid_arr.push(source.data['staxids'][i])
                                    tab_dict_org_count[source.data[tax_unit][i]]+=1
                                }

                            }
                        }
                    }

                }
            }


            table_data.data['#RBHs'] = []
            table_data.data['# Different Organisms In DB'] = []
            table_data.data[tax_unit] = []
            table_data.data['# Different Organisms In Selection'] = []
            for(let key in tab_dict){
                table_data.data['#RBHs'].push(tab_dict[key])
                table_data.data['# Different Organisms In DB'].push(tab_dict_static[key])
                table_data.data[tax_unit].push(key)
                table_data.data['# Different Organisms In Selection'].push(tab_dict_org_count[key])
            }

            table_data.change.emit();
            sc.change.emit();
            """)

            selection_callback = CustomJS(
                args=dict(sc=Curr, source=Overall, tax_unit=taxonomic_unit, table_data=DbData, menu=menu,
                          qseqids=menu_qseqids), code="""
                var call_back_object = cb_obj.indices

                var tab_dict = {};
                var tab_dict_static = {};
                var tab_dict_org_counter = {};
                for(var i = 0;i<table_data.get_length();i++){
                    tab_dict[table_data.data[tax_unit][i]] = 0
                    tab_dict_org_counter[table_data.data[tax_unit][i]] = 0
                    tab_dict_static[table_data.data[tax_unit][i]] = table_data.data['# Different Organisms In DB'][i]
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
                table_data.data[tax_unit] = []
                table_data.data['# Different Organisms In Selection'] = []
                for(let key in tab_dict){
                    table_data.data['#RBHs'].push(tab_dict[key])
                    table_data.data['# Different Organisms In DB'].push(tab_dict_static[key])
                    table_data.data[tax_unit].push(key)
                    table_data.data['# Different Organisms In Selection'].push(tab_dict_org_counter[key])
                }

                table_data.change.emit();
                """)

            Curr.selected.js_on_change('indices', selection_callback)

            menu_qseqids.js_on_change('value', menu_qseqid_callback)

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
                var csv = `qseqid,sacc,staxids,tax_unit\n`;  

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
            grid = gridplot([[column(p),
                              column(menu, menu_qseqids, row(circle_size_spinner, line_size_spinner),
                                     range_slider,download_selection_button, table)]],
                            toolbar_location='right')

            output_file(filename=path_to_bokeh_plot,
                        title="Interactive Graph Percent Identity vs. Bitscore linked to {} database entries".format(
                            taxonomic_unit))
            save(grid)

        return 0
    except Exception as e:
        raise Exception("ERROR in producing bokeh plots for database statistics with exception: {}".format(e))