import pandas as pd
import matplotlib.pyplot as plt
from bokeh.models import Button
import seaborn as sns
#output_file-to save the layout in file, show-display the layout , output_notebook-to configure the default output state  to generate the output in jupytor notebook.
from bokeh.io import output_file, save
#ColumnDataSource makes selection of the column easier and Select is used to create drop down
from bokeh.models import ColumnDataSource, Select, Spinner, MultiSelect, ColorPicker, RangeSlider, DataTable, TableColumn, HTMLTemplateFormatter
#Figure objects have many glyph methods that can be used to draw vectorized graphical glyphs. example of glyphs-circle, line, scattter etc.
from bokeh.plotting import figure
#To create intractive plot we need this to add callback method.
from bokeh.models import CustomJS, Legend
#This is for creating layout
from bokeh.layouts import column, gridplot, row
from bokeh.core.enums import MarkerType
plt.rcParams['legend.fontsize'] = 10
from bokeh.palettes import inferno, viridis, magma, Spectral
from random import shuffle


def bokeh_django_task_button(current_selection: ColumnDataSource, remote_or_local:str) -> CustomJS:
    task_selection_callback = CustomJS(args=dict(sc=current_selection, remote_or_local=remote_or_local), code="""
        var viewData = { 
            accessions : [],
            url : window.location.href
        };
        var jsonData = {};

        for(var i = 0; i < sc.selected.indices.length; i++){
            jsonData[i] = sc.data['sacc'][sc.selected.indices[i]]
        }

        viewData.accessions.push(jsonData);

        var json = JSON.stringify(viewData);
        
        var base_url = window.location.href
        base_url = base_url.split("/")[2]
        base_url = "http://" + base_url + "/external_tools/bokeh_database_task/" + remote_or_local
        $.ajax({
            type: "POST",
            url: base_url,
            data: json,
            success: function(result){
                setTimeout(function(){
                   window.location.reload();
                }, 2000);
            },
            dataType: "json"
        });


    """)
    return task_selection_callback

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


def create_initial_bokeh_result_data(result_data: pd.DataFrame, taxonomic_unit: str) -> tuple:
    try:
        # RBH result dataframe
        result_data = result_data.loc[:,
                      ['order', 'class', 'phylum', 'genus', 'family', 'bitscore', 'pident', 'stitle', 'scomnames',
                       'staxids', 'qseqid', 'sacc', 'slen']]  # ,'slen'
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


'''color_callback

    This function creates a js callback for changing the color based on the selection
    of a taxonomic_unit. The actual colors are saved within the color_dict, those colors
    are used to replace the color entry in the selection table and data data source.

    :param p_legend
        :type
    :param current_selection
        :type ColumnDataSource
    :param data
        :type ColumnDataSource
    :param color_dict
        :type dict
    :param menu_qseqids
        :type MultiSelect

    :returns callback
        :type CustomJS
'''


def create_color_callback(p_legend, current_selection: ColumnDataSource, data: ColumnDataSource, color_dict: dict,
                          menu_qseqids: MultiSelect) -> MultiSelect:
    callback = CustomJS(
        args=dict(legend=p_legend, sc=current_selection, source=data, color_dict=color_dict, menu_qseqids=menu_qseqids),
        code='''  
                    var tax_unit = cb_obj.value

                    legend.label = {'field':tax_unit}
                    var length = sc.get_length()
                    sc.data['color']=[]
                    for(var i = 0; i < length; i++){
                        sc.data['color'].push(color_dict[sc.data[tax_unit][i]])
                    }
                    sc.change.emit();
                ''')
    return callback


'''build_taxonomy_menu

    This function produces a MultiSelect widget that serves as taxonomy menu.

    :param bokeh_dataframe
        :type pd.DataFrame
    :param taxonomic_unit
        :type str

'''


def build_taxonomy_menu(bokeh_dataframe: pd.DataFrame, taxonomic_unit: str) -> MultiSelect:
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


def build_json_callback_for_taxonomy(column_dat: ColumnDataSource, static_dat: ColumnDataSource,
                                     taxonomic_unit: str, tax_selection: dict,
                                     menu_qseqid: MultiSelect, xaxis_menu: Select, yaxis_menu: Select,
                                     color_menu: Select, color_dict: dict,
                                     taxonomy_table_callback_dict: dict) -> CustomJS:
    tax_menu_callback = CustomJS(args=dict(sc=column_dat,
                                           source=static_dat,
                                           tax_unit=taxonomic_unit,
                                           selected_taxonomy=tax_selection,
                                           menu_qseqids=menu_qseqid,
                                           xaxis_menu=xaxis_menu, yaxis_menu=yaxis_menu,
                                           color_menu=color_menu, color_dict=color_dict,
                                           tax_dict=taxonomy_table_callback_dict), code="""

                        var call_back_object = cb_obj.value 


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
                                            for(var l = 0;l<tax_dict[color_menu.value].length;l++){
                                                if(source.data[color_menu.value][i] == tax_dict[color_menu.value][l]){
                                                    if(taxid_arr.includes(source.data['staxids'][i]) == false){
                                                        taxid_arr.push(source.data['staxids'][i])
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

                        sc.change.emit();
                        """)
    return tax_menu_callback


def create_color_palette_selection():
    try:
        # palettes = Spectral
        options = [(str(val), "Spectral" + str(val)) for val in range(3, 12)]
        color_palette_menu = Select(options=options,
                                    value=str(3),
                                    title="Select a color palette")
        return color_palette_menu
    except Exception as e:
        raise Exception("[-] ERROR couldnt create color palette selection")


def create_color_palette_selection_callback(curr: ColumnDataSource, color_menu: Select,
                                            taxonomy_table_callback_dict: dict) -> CustomJS:
    try:
        palettes = Spectral
        c_palette_callback = CustomJS(
            args=dict(sc=curr, color_menu=color_menu, tax_menu=taxonomy_table_callback_dict, pals=palettes), code="""
        // the callback value is a number 3,4,5,6,7,8,9,10,11,12
        var call_back_object = cb_obj.value

        var unique_organisms = []


        for(var i = 0; i<sc.get_length(); i++){
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
            for(var i = 0; i<sc.data[color_menu.value].length; i++){
                sc.data['color'].push(color_dict[sc.data[color_menu.value][i]])
            }
         sc.change.emit();
         }
        """)
        return c_palette_callback
    except Exception as e:
        raise Exception("[-] ERROR couldnt create color palette selection callback with exception: {}".format(e))


def create_qseqid_menu_callback(Overall: ColumnDataSource, Curr: ColumnDataSource,
                                xaxis_menu: Select, yaxis_menu: Select,
                                color_menu: Select, color_dict: dict, taxonomy_menus: list) -> CustomJS:
    try:
        menu_qseqid_callback = CustomJS(args=dict(source=Overall, sc=Curr,
                                                  xaxis_menu=xaxis_menu, yaxis_menu=yaxis_menu,
                                                  color_menu=color_menu, color_dict=color_dict,
                                                  taxonomy_menus=taxonomy_menus), code="""

        var tax_unit = color_menu.value;
        var call_back_object = cb_obj.value;

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
                                        if(taxid_arr.includes(source.data['staxids'][i]) == false){
                                            taxid_arr.push(source.data['staxids'][i])
                                        }
                                     }                               
                                }
                            }
                        }
                     }





                }
            }
        }

        sc.change.emit();
        """)
        return menu_qseqid_callback
    except Exception as e:
        raise Exception("[-] ERROR creating the custom js callback for the qseqid menu with exception: {}".format(e))


'''create_unlinked_bokeh_plot

    For documentation view the function definition for create_linked_bokeh_plot in the blast_project module
    py_database_statistics.py function.

'''


def create_unlinked_bokeh_plot(result_data: pd.DataFrame, taxonomic_unit: str) -> int:
    try:
        # change unknown entries to their corresponding higher taxonomic nodes
        result_data = result_data.copy()
        result_data.loc[result_data[result_data['class'] == "unknown"].index, 'class'] = \
            result_data[result_data['class'] == "unknown"]['phylum']
        result_data.loc[result_data[result_data['order'] == "unknown"].index, 'order'] = \
            result_data[result_data['order'] == "unknown"]['class']
        result_data.loc[result_data[result_data['family'] == "unknown"].index, 'family'] = \
            result_data[result_data['family'] == "unknown"]['order']
        result_data.loc[result_data[result_data['genus'] == "unknown"].index, 'genus'] = \
            result_data[result_data['genus'] == "unknown"]['family']

        # create bokeh dataframes for plot
        data_all, color_dict = create_initial_bokeh_result_data(result_data, taxonomic_unit)
        # selection subset for initial plot data
        data_selection, taxcount_df = create_initial_bokeh_data_selection(data_all, taxonomic_unit)

        ########

        # setup bokeh classes
        Overall = ColumnDataSource(data=data_all)
        Curr = ColumnDataSource(data=data_selection)

        # plot and the menu is linked with each other by this callback function
        unique_tax = list(data_all[taxonomic_unit].unique())
        unique_qseqids = list(data_all['qseqid'].unique())

        qseq_values = []
        if len(unique_qseqids) > 1:
            qseq_values.append(unique_qseqids[0])
            qseq_values.append(unique_qseqids[1])
        else:
            qseq_values.append(unique_qseqids[0])

        menu_qseqids = MultiSelect(options=unique_qseqids, value=qseq_values,
                                   title='Select target query sequence')  # drop down menu

        TOOLTIPS = [
            ("stitle", "@stitle"),
            ("bitscore,pident", "@bitscore, @pident"),
            ("sacc RBH to qseqid", "@sacc RBH to @qseqid "),
            ("scomname", "@scomnames"),
        ]

        # x_range=(0, result_data['bitscore'].max() + result_data['bitscore'].min())
        p = figure(x_axis_label='bitscore', y_axis_label='pident',
                   plot_height=700, plot_width=900,
                   tooltips=TOOLTIPS,
                   tools="lasso_select, reset,save, box_zoom,undo,redo,wheel_zoom, pan",
                   title="Number of RBHs - pident vs bitscore",
                   )  # ,tools="box_select, reset" creating figure object

        p.add_layout(Legend(), 'left')

        circle = p.scatter(x='x', y='y', color='color', marker='marker', size=10, line_width=1,
                           line_color='black',
                           source=Curr, legend_field=taxonomic_unit)  # plotting the data using glyph circle

        p.legend.glyph_width = 40
        p.legend.glyph_height = 40

        color_menu = Select(options=['phylum', 'class', 'order', 'family', 'genus'],
                            value=taxonomic_unit, title="Select Legend Color")
        color_callback = create_color_callback(p.legend.items[0], Curr, Overall, color_dict, menu_qseqids)
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
        phylum_menu_callback = build_json_callback_for_taxonomy(Curr, Overall, 'phylum',
                                                                tax_selection_dict, menu_qseqids,
                                                                x_axis_menu, y_axis_menu, color_menu,
                                                                color_dict, taxonomy_table_callback_dict)
        phylum_menu.js_on_change('value', phylum_menu_callback)

        tax_selection_dict = {'order': order_menu, 'family': family_menu, 'genus': genus_menu}
        class_menu_callback = build_json_callback_for_taxonomy(Curr, Overall, 'class', tax_selection_dict,
                                                               menu_qseqids, x_axis_menu, y_axis_menu, color_menu,
                                                               color_dict, taxonomy_table_callback_dict)
        class_menu.js_on_change('value', class_menu_callback)

        tax_selection_dict = {'family': family_menu, 'genus': genus_menu}
        order_menu_callback = build_json_callback_for_taxonomy(Curr, Overall, 'order', tax_selection_dict,
                                                               menu_qseqids, x_axis_menu, y_axis_menu,
                                                               color_menu, color_dict, taxonomy_table_callback_dict)
        order_menu.js_on_change('value', order_menu_callback)

        tax_selection_dict = {'genus': genus_menu}
        family_menu_callback = build_json_callback_for_taxonomy(Curr, Overall, 'family', tax_selection_dict,
                                                                menu_qseqids, x_axis_menu, y_axis_menu,
                                                                color_menu, color_dict,
                                                                taxonomy_table_callback_dict)
        family_menu.js_on_change('value', family_menu_callback)

        genus_menu_callback = build_json_callback_for_taxonomy(Curr, Overall, 'genus', {},
                                                               menu_qseqids, x_axis_menu, y_axis_menu,
                                                               color_menu, color_dict, taxonomy_table_callback_dict)
        genus_menu.js_on_change('value', genus_menu_callback)

        menu_qseqid_callback = create_qseqid_menu_callback(Overall,
                                                           Curr,
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

        download_selection_callback = CustomJS(args=dict(sc=Curr, tax_unit=taxonomic_unit), code="""
            var temp = []
            var csvFileData = []
            for(var i = 0; i < sc.selected.indices.length; i++){
                temp = [sc.data['qseqid'][sc.selected.indices[i]],
                        sc.data['sacc'][sc.selected.indices[i]],
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
            var file = new File([csv], "selection.csv" ,{type: "octet/stream"});
            var url = URL.createObjectURL(file);
            window.location.assign(url);
            URL.revokeObjectUrl(url);
        """)

        download_selection_button = Button(label="Download Selection")
        download_selection_button.js_on_click(download_selection_callback)

        bokeh_task_button_callback = bokeh_django_task_button(Curr, 'remote')
        bokeh_task_button = Button(label="Selection Constrained Phylogeny")
        bokeh_task_button.js_on_click(bokeh_task_button_callback)

        color_palette = create_color_palette_selection()

        color_palette_callback = create_color_palette_selection_callback(Curr, color_menu,
                                                                         taxonomy_table_callback_dict)
        color_palette.js_on_change('value', color_palette_callback)

        grid = gridplot([[column(p),
                          column(menu_qseqids, row(circle_size_spinner, line_size_spinner),
                                 range_slider, download_selection_button, bokeh_task_button, x_axis_menu, y_axis_menu),
                          column(phylum_menu, class_menu, order_menu, family_menu, genus_menu, color_menu,
                                 color_palette)]],
                        toolbar_location='right')

        output_file(filename=snakemake.output['bokeh_plot'],
                    title="BLAST Results".format(
                        taxonomic_unit))
        save(grid)

        return 0
    except Exception as e:
        raise Exception("ERROR in producing bokeh plots with exception: {}".format(e))

RETURNCODE=5
with open(snakemake.log['log'], 'w') as logfile:
    try:


        logfile.write("INFO:trying to load blast dataframe ...\n")
        result_df = pd.read_csv(snakemake.input['blast_results'], delimiter=",", index_col=0)

        logfile.write("INFO:loaded BLAST dataframe ...\n")


        if result_df['staxids'].dtype != 'int':
            logfile.write("INFO:changing staxids column datatype to int64 ... \n")
            result_df['staxids'] = result_df['staxids'].astype('int64')


        logfile.write("INFO:start producing interactive bokeh plots for target species families\n")

        # bokeh code
        try:
            bokeh_function = create_unlinked_bokeh_plot(result_df, 'family')
        except Exception as e:
            logfile.write("ERROR:exception was thrown during bokeh plot creation: {}\n".format(e))
            raise Exception("[-] couldnt create bokeh plot with exception: {}".format(e))

        logfile.write("INFO:finished bokeh procedure\n")

        logfile.write("DONE\n")
    except Exception as e:
        logfile.write("ERROR:{}\n".format(e))
        sys.exit(RETURNCODE)