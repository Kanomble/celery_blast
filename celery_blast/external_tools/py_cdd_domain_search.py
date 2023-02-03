import bokeh.models
import pandas as pd
from django.conf import settings
from sklearn.decomposition import PCA

#module to parse fasta files
from Bio import SeqIO
#output_file-to save the layout in file, show-display the layout , output_notebook-to configure the default output state  to generate the output in jupytor notebook.
from bokeh.io import output_file, save
#ColumnDataSource makes selection of the column easier and Select is used to create drop down
from bokeh.models import ColumnDataSource, MultiSelect,  DataTable, TableColumn, HTMLTemplateFormatter
#Figure objects have many glyph methods that can be used to draw vectorized graphical glyphs. example of glyphs-circle, line, scattter etc.
from bokeh.plotting import figure
#To create intractive plot we need this to add callback method.
from bokeh.models import CustomJS, Legend
#This is for creating layout
from bokeh.layouts import column, gridplot
import random

'''load_domain_query_data
    
    This function reads the domain dataframe for the query sequences and returns a dictionary with query sequence ids
    as keys and CDD accession hits as values.
    
    :param path_to_query
        :type str
    
    :return domain_dict
        :type dict[str] = list[str,...]
'''
def load_domain_query_data(path_to_query_domains:str)->dict:
    try:
        #sseq qseq stitle
        header = "qseqid qlen sacc slen qstart qend sstart send qseq sseq bitscore evalue pident stitle".split(" ")
        cdd_queries = pd.read_table(path_to_query_domains,header=None)
        cdd_queries.columns=header
        domain_dict = {}
        for qseq in cdd_queries.qseqid.unique():
            qseq_key = qseq.split(".")[0]
            domain_dict[qseq_key]=cdd_queries[cdd_queries.qseqid == qseq].sacc.unique()
        return domain_dict
    except Exception as e:
        raise Exception("[-] ERROR loading domain query dictionary with exception: {}".format(e))

'''load_domain_data
    
    This function reads the CDD domain dataframe for one particular query sequence, specified by qseqid_key.
    It creates a pandas dataframe with percent identity values for the RBHs based on all query sequence domains. 
    The resulting dataframe is a table with float values ranging from 0 to 100 representing the percent identity of the 
    query sequence domains of the RBHs in respect to the underlying CDD database.
    
    :param path_to_domains - e.g. 'media/blast_project/2/WP_087654321/cdd_domains.tsf'
        :type str
    :param qseqid_key - e.g. 'WP_087654321'
        :type str
    :param domain_dict - dictionary from function load_domain_query_data
        :type dict[str]=list[str, ...]
        
    :returns cdd_dataframe
        :type pd.DataFrame
    
'''
def load_domain_data(path_to_domains: str, qseqid_key: str, domain_dict: dict) -> pd.DataFrame:
    try:
        header = "qseqid qlen sacc slen qstart qend sstart send bitscore evalue pident".split(" ")
        cdd_output = pd.read_table(path_to_domains, header=None)
        cdd_output.columns = header
        transform_protids = lambda protid: protid.split(".")[0]
        cdd_output['transformed_qseqid'] = cdd_output.qseqid.apply(transform_protids)
        cdd_output['pident'] = cdd_output.pident.apply(lambda x: round(x / 100, 2))

        query_domain_df = pd.DataFrame(domain_dict[qseqid_key], columns=['sacc'])
        query_domain_df = query_domain_df.append({'sacc': 'additional_domains'}, ignore_index=True)

        header = list(query_domain_df['sacc'])
        cdd_dataframe = pd.DataFrame({qseqid_key: [1 for i in range(len(header))]}).transpose()
        cdd_dataframe.columns = header
        cdd_dataframe['additional_domains'] = 0

        apply_matrix_transform = lambda x: 1 if x != 0 else 0
        for query_seq in cdd_output.qseqid.unique():

            domain_df = cdd_output[cdd_output['qseqid'] == query_seq][['sacc', 'qseqid', 'pident']]
            domain_df = domain_df.drop_duplicates(subset='sacc').reset_index()
            domain_df = domain_df.drop('index', axis=1)

            #counting domains not present in the query sequence
            counter = 0
            for domain in domain_df.sacc:
                if domain not in header:
                    counter += 1

            merged_df = domain_df.merge(query_domain_df, on=['sacc'], how='right')
            merged_df = merged_df.fillna(0)
            merged_df['qseqid'] = merged_df['qseqid'].apply(apply_matrix_transform)
            index = merged_df[merged_df['sacc'] == 'additional_domains'].index[0]
            merged_df.loc[index, 'pident'] = counter

            merged_df = pd.DataFrame(merged_df['pident'])
            merged_df.columns = [query_seq]
            merged_df = merged_df.transpose()
            merged_df.columns = header
            cdd_dataframe = pd.concat([cdd_dataframe, merged_df])
        return cdd_dataframe
    except Exception as e:
        raise Exception("[-] ERROR loading domain dataframe with exception: {}".format(e))


'''add_color_column_to_dataframe

    This function adds a color column to the input dataframe based on 
    a tax column of this dataframe. The tax column is defined by the 
    provided taxonomic unit.

    :param result_df --> RBH result dataframe
        :type pd.DataFrame
    :param taxonomic_unit
        :type str

    :return result_df
        :type pd.DataFrame
    :return color_dict
        :type dict
'''
def add_color_column_to_dataframe(result_df: pd.DataFrame, taxonomic_unit) -> tuple:
    try:
        import matplotlib._color_data as mcd
        overlap = [name for name in mcd.CSS4_COLORS]
        overlap.remove("lightgrey")
        num_colors = len(result_df[taxonomic_unit].unique())
        # clrs = sns.color_palette('pastel', n_colors=num_colors)
        # clrs = clrs.as_hex()
        colors = []
        for i in range(num_colors):
            c = overlap[random.randint(0, len(overlap) - 1)]
            if c not in colors:
                colors.append(c)
            # color c is in colors --> one color for at least two different tax
            else:
                # select other random color
                c = overlap[random.randint(0, len(overlap) - 1)]
                if c not in colors:
                    colors.append(c)
                else:
                    colors.append('black')

        color_dict = dict(zip(result_df[taxonomic_unit].unique(), colors))  # colors
        create_color_scheme = lambda value: color_dict[value]
        result_df['color'] = result_df[taxonomic_unit].apply(create_color_scheme)

        return result_df, color_dict
    except Exception as e:
        raise Exception("[-] ERROR creating color column for result dataframe with exception: {}".format(e))


'''build_dataframe_for_bokeh

    This function transforms the result dataframe from the CDD domain search, the RBH
    result dataframe and the PCA result dataframe for bokeh ColumnData input.
    Additionally returns a list of columns for the bokeh DataTable.

    :param cdd_dataframe
        :type pandas.DataFrame
    :param pca_df 
        :type pandas.DataFrame
    :param selection
        :type pandas.DataFrame

    :return pca_df -> transformed dataframe
        :type pandas.DataFrame
    :return header -> CDD domains
        :type list[str]
'''
def build_dataframe_for_bokeh(cdd_dataframe: pd.DataFrame, pca_df: pd.DataFrame, selection: pd.DataFrame) -> tuple:
    try:
        cdd_dataframe = cdd_dataframe.reset_index()
        cdd_dataframe['sacc'] = list(cdd_dataframe['index'])
        cdd_dataframe = cdd_dataframe.drop("index", axis=1)
        cdd_dataframe['transformed_sacc'] = cdd_dataframe.sacc.apply(lambda x: x.split(".")[0])

        # get all PC0 -> PCX columns of the pca_df dataframe
        cols = []
        for col in pca_df.columns:
            try:
                int(col)
                col = "PC" + str(col)
                cols.append(col)
            except:
                cols.append(col)
                continue

        pca_df.columns = cols
        pca_df['sacc'] = list(selection.index)
        pca_df['bitscore'] = list(selection.bitscore)
        pca_df['pident'] = list(selection.pident)
        pca_df['order'] = list(selection.order)
        pca_df['evalue'] = list(selection.evalue)
        pca_df['family'] = list(selection.family)
        pca_df['genus'] = list(selection.genus)
        pca_df['phylum'] = list(selection.phylum)
        pca_df['class'] = list(selection['class'])
        pca_df['stitle'] = list(selection.stitle)

        pca_df['transformed_sacc'] = pca_df['sacc'].apply(lambda x: x.split(".")[0])
        pca_df = pca_df.merge(cdd_dataframe.loc[1:, :], on='sacc')

        # header for bokeh data table
        header = list(cdd_dataframe.columns)
        header.append('transformed_sacc_x')
        header.remove('transformed_sacc')

        return pca_df, header
    except Exception as e:
        raise Exception("[-] ERROR during creation of bokeh dataframe with exception: {}".format(e))


'''build_taxonomy_menu

    This function constructs a taxonomy menu for the bokeh plot.

    :param bokeh_dataframe
        :type pd.DataFrame
    :param taxonomic_unit
        :type str

    :return tax_menu
        :type MulitSelect -> bokeh
'''
def build_taxonomy_menu(bokeh_dataframe: pd.DataFrame, taxonomic_unit: str):
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


'''build_json_callback_for_selection

    This function renders the bokeh DataTable based on users lasso_tool, ... manual selection.
    The cb_obj.indices are the indices of the selected column_dat circles/points of the interactive plot.
     
    :param column_dat
        :type bokeh.models.ColumnDataSource
    :param table_dat
        :type bokeh.models.ColumnDataSource
    :param table_header
        :type list[str]
        
    :returns selection_callback
        :type bokeh.models.CustomJS
'''
def build_json_callback_for_selection(column_dat: ColumnDataSource, table_dat: ColumnDataSource,
                                      table_header: list) -> CustomJS:
    selection_callback = CustomJS(args=dict(sc=column_dat, table_data=table_dat, columns=table_header), code="""
            var call_back_object = cb_obj.indices


            table_data.data['transformed_sacc_x'] = []
            table_data.data['index'] = []
            table_data.data['sacc'] = []
            //producing empty table
            for(var i = 0; i < columns.length; i++){
                table_data.data[columns[i]] = []
            }

            for(var i = 0; i < call_back_object.length; i++){
                    for(var j = 0; j < columns.length; j++){
                        table_data.data[columns[j]].push(sc.data[columns[j]][call_back_object[i]])
                    }
                    table_data.data['transformed_sacc_x'].push(sc.data['transformed_sacc_x'][call_back_object[i]])
                    table_data.data['sacc'].push(sc.data['sacc'][call_back_object[i]])
                    table_data.data['index'].push(call_back_object[i])

            }

            table_data.change.emit();
            """)
    return selection_callback

'''build_json_callback_for_taxonomy
    
    This function produces the custom javascript callback function for the taxonomy. It is similar to the 
    build_json_callback_for_taxonomy function in blast_project.py_database_statistics.py.
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
    :param tax_selection - e.g. tax_selection_dict = {'class':class_menu,'order':order_menu,'family':family_menu,'genus':genus_menu}
        :type dict[str] = bokeh.models.MultiSelect
    
    :returns tax_menu_callback
        :type bokeh.models.CustomJS
'''
def build_json_callback_for_taxonomy(column_dat: ColumnDataSource, static_dat: ColumnDataSource,
                                     table_dat: ColumnDataSource, domains: list, taxonomic_unit: str,
                                     tax_selection: dict) -> CustomJS:
    tax_menu_callback = CustomJS(args=dict(sc=column_dat,
                                           source=static_dat,
                                           tax_unit=taxonomic_unit,
                                           columns=list(column_dat.data.keys()),
                                           domains=domains,
                                           selected_taxonomy=tax_selection,
                                           table_data=table_dat), code="""
                        var call_back_object = cb_obj.value                        
                        for(var i = 0;i < columns.length;i++){
                            sc.data[columns[i]]=[]
                        }

                        table_data.data['transformed_sacc_x'] = []
                        table_data.data['index'] = []
                        table_data.data['sacc'] = []

                        for(var i = 0; i < domains.length; i++){
                            table_data.data[domains[i]] = []
                        }

                        var unique_class = []
                        var unique_order = []
                        var unique_family = []
                        var unique_genus = []

                        for(var i = 0; i < source.get_length(); i++){
                            for(var j = 0; j < call_back_object.length; j++){
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


                                    for(var k = 0; k < columns.length;k++){
                                            sc.data[columns[k]].push(source.data[columns[k]][i])
                                        }

                                    for(var y = 0; y < domains.length; y++){
                                            table_data.data[domains[y]].push(source.data[domains[y]][i])
                                        }

                                    table_data.data['transformed_sacc_x'].push(source.data['transformed_sacc_x'][i])
                                    table_data.data['sacc'].push(source.data['sacc'][i])
                                    table_data.data['index'].push(i)
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


                        table_data.change.emit();
                        sc.change.emit();
                        """)
    return tax_menu_callback

'''build_table_html_formatter
    
    This function recolors the bokeh datatable, based on the pident values of the underlying CDD dataframe.
    It is used within the build_table_columns function, to render the datatable.
    
    :param value - pident
        :type str
        
    :returns formatter
        :type bokeh.models.HTMLTemplateFormatter
'''
def build_table_html_formatter(value:str)->bokeh.models.HTMLTemplateFormatter:
    condition_string=""" 
    if({val}<0.25)
        {{return("red")}}
    else if({val}>=0.25 && {val}<0.5)
        {{return("orange")}}
    else if({val}>=0.5 && {val}<0.75)
        {{return("lightgreen")}}
    else if({val}>=0.75)
        {{return("green")}}
    """.format(val=value)
    template="""
    <p style="color:<%=(function colorfromint(){{{cond}}}()) %>;"><%=value%></p>
            """.format(cond=condition_string)
    formatter = HTMLTemplateFormatter(template=template)
    return formatter


'''build_table_columns

    This function takes the CDD domain search dataframe as input and
    transforms the domain columns to TableColumn classes for the bokeh table.

    :param cdd_dataframe
        :type pandas.DataFrame


    :return table_columns
        :type TableColumn
    :return cdd_header
        :type list[str]
'''
def build_table_columns(cdd_dataframe: pd.DataFrame)->tuple:
    try:
        table_columns = []
        cdd_headers = []

        for col in cdd_dataframe.columns:
            col = col.replace(":", "_")
            if "CDD" in col:
                table_columns.append(
                    TableColumn(field=col, title=col, formatter=build_table_html_formatter(col))
                )
                cdd_headers.append(col)
            elif col == "transformed_sacc_x":
                table_columns.append(
                    TableColumn(field=col, title="Accession ID")
                )
            elif col == 'additional_domains':
                table_columns.append(TableColumn(field=col, title="Add. Domains"))
                cdd_headers.append(col)
        return table_columns, cdd_headers
    except Exception as e:
        raise Exception("[-] ERROR creating table column for bokeh table with exception: {}".format(e))


'''build_bokeh_plot

    This function produces an interactive bokeh plot for the visualization of 
    the CDD result dataframe.

    :param bokeh_dataframe -> PCA results and additional information
        :type pd.DataFrame
    :param domains -> list[str] for bokeh DataTable
        :type list
    :param taxonomic_unit
        :type str
    :param variances -> PCA variances
        :type list[str]
    :param query_sequence
        :type str
'''
def build_bokeh_plot(bokeh_dataframe: pd.DataFrame, domains: list, taxonomic_unit: str, variances: list,
                     query_sequence: str):
    try:
        # bokeh data preparation
        bokeh_dataframe = bokeh_dataframe.rename(
            columns=dict(zip([val for val in domains], [val.replace(":", "_") for val in domains])))
        columns = ['PC0', 'PC1', 'sacc', 'color', 'pident', 'bitscore', 'evalue', 'genus', 'family', 'order', 'phylum',
                   'class', 'stitle', 'additional_domains']
        domains = [domain.replace(":", "_") for domain in domains]
        columns.extend(domains)
        tab_dat = domains.copy()
        tab_dat.append('additional_domains')
        column_dat = ColumnDataSource(bokeh_dataframe[columns])
        static_data = ColumnDataSource(bokeh_dataframe[columns])
        table_dat = ColumnDataSource(bokeh_dataframe[tab_dat])
        table_columns, table_header = build_table_columns(bokeh_dataframe)

        TOOLTIPS = [
            ("family: ", "@family"),
            ("order: ", "@order"),
            ("genus: ", "@genus"),
            ("bitscore,pident: ", "@bitscore, @pident"),
            ("sacc: ", "@sacc"),
            ("title: ", "@stitle")
        ]

        # main figure properties
        p = figure(x_axis_label='PC1 with {}% captured variance'.format(variances[0]),
                   y_axis_label='PC2 with {}% captured variance'.format(variances[1]),
                   plot_height=700, plot_width=900,
                   tools="lasso_select, reset, save, box_zoom, undo, redo, wheel_zoom, pan",
                   tooltips=TOOLTIPS,
                   title="PCA on pidents of domains of {} RBHs".format(
                       query_sequence),
                   )

        p.add_layout(Legend(), 'left')
        p.legend.glyph_width = 40
        p.legend.glyph_height = 40

        # scatter plot
        circle = p.circle(x='PC0', y='PC1',
                          color='color', size=20, line_width=1, line_color='black',
                          source=column_dat,
                          legend_field=taxonomic_unit)

        # bokeh table
        table = DataTable(source=table_dat, width=400, height=275,
                          sizing_mode="stretch_both", reorderable=True, sortable=True, fit_columns=True,
                          columns=table_columns)

        # defining select funcionality
        selection_callback = build_json_callback_for_selection(column_dat, table_dat, table_header)  # circle
        column_dat.selected.js_on_change('indices', selection_callback)

        # taxonomy menus in bokeh plot
        phylum_menu = build_taxonomy_menu(bokeh_dataframe, 'phylum')
        class_menu = build_taxonomy_menu(bokeh_dataframe, 'class')
        order_menu = build_taxonomy_menu(bokeh_dataframe, 'order')
        family_menu = build_taxonomy_menu(bokeh_dataframe, 'family')
        genus_menu = build_taxonomy_menu(bokeh_dataframe, 'genus')

        # defining callback functions for the taxonomy menus
        tax_selection_dict = {'class': class_menu, 'order': order_menu, 'family': family_menu, 'genus': genus_menu}
        phylum_menu_callback = build_json_callback_for_taxonomy(column_dat, static_data, table_dat, table_header,
                                                                'phylum', tax_selection_dict)
        phylum_menu.js_on_change('value', phylum_menu_callback)

        tax_selection_dict = {'order': order_menu, 'family': family_menu, 'genus': genus_menu}
        class_menu_callback = build_json_callback_for_taxonomy(column_dat, static_data, table_dat, table_header,
                                                               'class', tax_selection_dict)
        class_menu.js_on_change('value', class_menu_callback)

        tax_selection_dict = {'family': family_menu, 'genus': genus_menu}
        order_menu_callback = build_json_callback_for_taxonomy(column_dat, static_data, table_dat, table_header,
                                                               'order', tax_selection_dict)
        order_menu.js_on_change('value', order_menu_callback)

        tax_selection_dict = {'genus': genus_menu}
        family_menu_callback = build_json_callback_for_taxonomy(column_dat, static_data, table_dat, table_header,
                                                                'family', tax_selection_dict)
        family_menu.js_on_change('value', family_menu_callback)

        genus_menu_callback = build_json_callback_for_taxonomy(column_dat, static_data, table_dat, table_header,
                                                               'genus', {})
        genus_menu.js_on_change('value', genus_menu_callback)

        return gridplot([[column(p), column(table), column(phylum_menu, class_menu, order_menu, family_menu, genus_menu)]], toolbar_location='right'), circle
    except Exception as e:
        raise Exception("[-] ERROR during creation of bokeh plots with exception: {}".format(e))

'''produce_bokeh_pca_plot
    
    This function is the controlling function for plotting the results of a principal component analysis based on the 
    CDD domain pident dataframe. It produces an interactive Bokeh plot.
    
    :param project_id
        :type int
    :param qseqid
        :type str
    :param taxonomic_unit - keyword
        :type str
    
    :returns return_value - 0 success, 1 failure
        :type int
'''
def produce_bokeh_pca_plot(project_id:int, qseqid: str,
                        taxonomic_unit='class')->int:
    try:
        path_to_output = settings.BLAST_PROJECT_DIR + str(project_id) + '/' + qseqid + '/pca_bokeh_domain_plot.html'
        path_to_query_domains = settings.BLAST_PROJECT_DIR + str(project_id) + '/query_domains.tsf'
        path_to_domains = settings.BLAST_PROJECT_DIR + str(project_id) + '/' + qseqid + '/cdd_domains.tsf'
        result_df_path = settings.BLAST_PROJECT_DIR + str(project_id) + '/reciprocal_results_with_taxonomy.csv'
        result_df = pd.read_csv(result_df_path, index_col=0)

        # add color column to result_df
        result_df, color_dict = add_color_column_to_dataframe(result_df, taxonomic_unit)

        selection = result_df[result_df['qseqid'] == qseqid][['sacc',
                                                              'color',
                                                              'pident',
                                                              'bitscore',
                                                              'evalue',
                                                              'genus',
                                                              'family',
                                                              'order',
                                                              'phylum',
                                                              'class',
                                                              'stitle']]
        selection = selection.reset_index()
        selection.index = selection['sacc']
        selection = selection.drop("index", axis=1)

        # load domain data
        domain_dict = load_domain_query_data(path_to_query_domains)
        cdd_dataframe = load_domain_data(path_to_domains, qseqid, domain_dict)

        # perform principal component analysis
        if len(cdd_dataframe.columns) > 1:
            selection = pd.merge(selection, cdd_dataframe, left_index=True, right_index=True)

            if len(selection.index) <= len(cdd_dataframe.columns):
                pca_selection = PCA(n_components=len(selection.index), svd_solver='full')
            else:
                pca_selection = PCA(n_components=len(cdd_dataframe.columns) - 1, svd_solver='full')
            #pca_selection = None
            cols = list(cdd_dataframe.columns)
            cols.remove('additional_domains')

            principal_components_selection = pca_selection.fit_transform(selection[cols])  # [final_df.columns]
            pca_df = pd.DataFrame(data=principal_components_selection)

            pca_df['color'] = list(selection['color'])

            # bokeh plotting
            variances = [round(pca_selection.explained_variance_ratio_[0] * 100, 3),
                         round(pca_selection.explained_variance_ratio_[1] * 100, 3)]
            bk_df, header = build_dataframe_for_bokeh(cdd_dataframe, pca_df, selection)
            grid, p = build_bokeh_plot(bk_df, header, taxonomic_unit, variances, qseqid)
            # save bokeh plot
            output_file(filename=path_to_output,
                      title="Principal Component Analysis of CDD-pidents".format(
                           taxonomic_unit))
            save(grid)
            return 0

        else:
            return 1

    except Exception as e:
        raise Exception("[-] ERROR with exception : {}".format(e))

########################################################################################################################

# returns false if both values are not in range
def check_range(x: int, y: int, boundary_x: int, boundary_y: int) -> bool:
    return boundary_x <= x <= boundary_y and boundary_x <= y <= boundary_y


# check current and potential new values for the the current boundaries and return updated values
def check_range_and_return_new_values(x: int, y: int, boundary_x: int, boundary_y: int) -> tuple:
    if boundary_x <= x <= boundary_y:
        new_x = boundary_x
        bool_x = False
    else:
        new_x = x
        bool_x = True

    if boundary_x <= y <= boundary_y:
        new_y = boundary_y
        bool_y = False
    else:
        bool_y = True
        new_y = y
    return (new_x, new_y), (bool_x, bool_y)


'''get_sequence_segments

    This function takes a pandas cdd dataframe for one protein sequence as input. It parses the qstart and qend
    columns and extracts overlapping parts within the sequence. It saves those overlapping parts in a list of tuples
    by the aid of check_range_and_return_new_values.

    :param qseqid_df
        :type pd.DataFrame

'''
def get_sequence_segments(qseqid_df: pd.DataFrame) -> list:
    try:
        segments = []
        # iterate over a sorted qstart qend set of the query sequence
        for i in sorted(set(zip(qseqid_df.qstart, qseqid_df.qend))):
            if len(segments) == 0:
                segments.append(list(i))
            else:
                curr_boundaries = segments[len(segments) - 1]
                # if one part of the segment is not within the current boundaries, update values or append a completly new segment
                if check_range(i[0], i[1], curr_boundaries[0], curr_boundaries[1]) == False:
                    values_new, bool_new = check_range_and_return_new_values(i[0], i[1], curr_boundaries[0],
                                                                             curr_boundaries[1])
                    if bool_new[0] == False:
                        segments.remove(curr_boundaries)
                        segments.append(list(values_new))
                    elif bool_new[1] == False:
                        segments.remove(curr_boundaries)
                        segments.append(list(values_new))

                    else:
                        segments.append(list(values_new))

        return segments
    except Exception as e:
        raise Exception(
            "[-] ERROR creating segment list for query sequence {} with exception {}".format(qseqid_df.qseqid[0], e))


'''write_domain_corrected_fasta_file

    This function executes the get_sequence_segments function for each unique query sequence within the cdd dataframe.
    It then parses the target_sequences.faa file within the query_sequence directory and writes a new fasta file based
    on the overlapping sequence segments. Writes a logfile to log/query_sequence/domain_corrected_fasta_file.log.

    :param project_id
        :type int
    :param query_sequence
        :type str

    :returns returncode
        :type int

'''
def write_domain_corrected_fasta_file(project_id: int, query_sequence: str) -> int:
    try:
        logfile_path = settings.BLAST_PROJECT_DIR + str(project_id) + '/log/' + query_sequence + '/domain_corrected_fasta_file.log'
        with open(logfile_path, "w") as logfile:
            logfile.write("INFO: Starting to produce multiple fasta file with overlapping domains.\n")
            cdd_filepath = settings.BLAST_PROJECT_DIR + str(project_id) + '/' + query_sequence + '/cdd_domains.tsf'
            fasta_filepath = settings.BLAST_PROJECT_DIR + str(project_id) + '/' + query_sequence + '/target_sequences.faa'
            domain_fasta_filepath = settings.BLAST_PROJECT_DIR + str(
                project_id) + '/' + query_sequence + '/domain_corrected_target_sequences.faa'
            # read cdd dataframe
            logfile.write("INFO: Loading dataframe with CDD hits.\n")
            df = pd.read_table(cdd_filepath, header=None)
            df.columns = "qseqid qlen sacc slen qstart qend sstart send bitscore evalue pident".split(" ")
            query_domains_range = {}
            logfile.write("INFO: Creating dictionary with unique query sequences (from CDD dataframe) as keys "
                          "and segments (list of integer values representing the overlap range) as values. \n")
            for qseqid in df.qseqid.unique():
                query_domains_range[qseqid] = get_sequence_segments(df[df['qseqid'] == qseqid])

            logfile.write("INFO: parsing target fasta file with biopython to slice sequences based on the segment list.\n")
            records = list(SeqIO.parse(fasta_filepath, "fasta"))
            # write new fasta file with only overlapping domain entries
            with open(domain_fasta_filepath, 'w') as domain_fasta:
                for rec in records:
                    identifier = rec.id
                    seq = ""
                    if identifier in list(query_domains_range.keys()):
                        for rng in query_domains_range[identifier]:
                            seq += str(rec.seq[rng[0]:rng[1]])

                        domain_fasta.write(">" + str(identifier) + "\n")
                        domain_fasta.write(seq + "\n")
                    else:
                        logfile.write("\tWARNING: for the sequence: {} no CDD have been found.\n".format(identifier))
            logfile.write("DONE")
        return 0
    except Exception as e:
        raise Exception("[-] ERROR producing domain corrected fasta file, with exception: {}".format(e))