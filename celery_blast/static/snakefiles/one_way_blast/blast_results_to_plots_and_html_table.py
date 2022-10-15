'''plot results and update taxonomic information

input: blast_results.table, config['query_sequence']
output: genus_bars.html, blast_results.html, blast_results_with_tax.table

Expands blast results table with tax information.
Creates altair plots based on tax information.
Produces blast results in html file.
'''

from Bio import Entrez
import pandas as pd
import altair as alt
import sys
import seaborn as sns
#output_file-to save the layout in file, show-display the layout , output_notebook-to configure the default output state  to generate the output in jupytor notebook.
from bokeh.io import output_file, save
#ColumnDataSource makes selection of the column easier and Select is used to create drop down
from bokeh.models import ColumnDataSource, Spinner, MultiSelect, ColorPicker, RangeSlider
#Figure objects have many glyph methods that can be used to draw vectorized graphical glyphs. example of glyphs-circle, line, scattter etc.
from bokeh.plotting import figure
#To create intractive plot we need this to add callback method.
from bokeh.models import CustomJS, Button
#This is for creating layout
from bokeh.layouts import column, gridplot

def create_unlinked_bokeh_plot(logfile: str,path_to_bokeh_plot:str,path_to_static:str, result_data: pd.DataFrame, taxonomic_unit: str) -> int:
    try:
        # log.write("INFO:checking if static dir: {} exists\n".format(path_to_static_dir))

        data_all = result_data.loc[:,
                   [taxonomic_unit, 'bitscore', 'pident', 'stitle', 'scomnames', 'staxids', 'qseqid',
                    'sacc']]
        data_all = data_all.sort_values(by=taxonomic_unit)
        num_colors = len(data_all[taxonomic_unit].unique())
        clrs = sns.color_palette('pastel', n_colors=num_colors)
        clrs = clrs.as_hex()
        color_dict = dict(zip(data_all[taxonomic_unit].unique(), clrs))
        create_color_scheme = lambda value: color_dict[value]

        # plot and the menu is linked with each other by this callback function
        unique_tax = list(data_all[taxonomic_unit].unique())
        unique_qseqids = list(data_all['qseqid'].unique())

        data_all['color'] = data_all[taxonomic_unit].apply(create_color_scheme)
        data_selection = data_all[
            (data_all[taxonomic_unit] == unique_tax[0]) | (data_all[taxonomic_unit] == unique_tax[1])]

        Overall = ColumnDataSource(data=data_all)
        Curr = ColumnDataSource(data=data_selection)
        # dbData=ColumnDataSource(data=dbData)

        TOOLTIPS = [
            ("stitle", "@stitle"),
            ("bitscore,pident", "@bitscore, @pident"),
            ("sacc RBH to qseqid", "@sacc RBH to @qseqid "),
            ("scomname", "@scomnames"),
        ]

        tax_menu = MultiSelect(options=unique_tax, value=[unique_tax[0], unique_tax[1]],
                               title='Select: ' + taxonomic_unit.capitalize())
        qseqid_menu = MultiSelect(options=unique_qseqids, value=unique_qseqids,
                                  title="Select target query sequences")

        tax_menu_callback = CustomJS(args=dict(source=Overall,
                                               sc=Curr,
                                               color_code=color_dict,
                                               tax_unit=taxonomic_unit,
                                               menu_qseqids=qseqid_menu), code="""
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
                    sc.data['sacc']=[]

                    for(var i = 0; i < source.get_length(); i++){
                        for(var j = 0; j < call_back_object.length; j++){
                            for(var k = 0; k < menu_qseqids.value.length; k++){
                                if(source.data[tax_unit][i] == call_back_object[j]){
                                    if(source.data['qseqid'][i] == menu_qseqids.value[k]){
                                        sc.data['bitscore'].push(source.data['bitscore'][i])
                                        sc.data['pident'].push(source.data['pident'][i])
                                        sc.data['color'].push(source.data['color'][i])
                                        sc.data[tax_unit].push(source.data[tax_unit][i])
                                        sc.data['index'].push(source.data['index'][i])
                                        sc.data['stitle'].push(source.data['stitle'][i])
                                        sc.data['scomnames'].push(source.data['scomnames'][i])
                                        sc.data['staxids'].push(source.data['staxids'][i])
                                        sc.data['qseqid'].push(source.data['qseqid'][i])
                                        sc.data['sacc'].push(source.data['sacc'][i])                                        
                                    }
                                }                                
                            }        
                        }
                    }


                    sc.change.emit();
                    """)
        qseqid_menu_callback = CustomJS(args=dict(source=Overall,
                                                  sc=Curr,
                                                  color_code=color_dict,
                                                  tax_unit=taxonomic_unit,
                                                  tax_selection=tax_menu), code="""
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
                    sc.data['sacc']=[]

                    for(var i = 0; i < source.get_length(); i++){
                        for(var j = 0; j < call_back_object.length; j++){
                            for(var k = 0; k < tax_selection.value.length; k++){
                                if(source.data[tax_unit][i] == tax_selection.value[k]){
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
                                        sc.data['sacc'].push(source.data['sacc'][i])                                        
                                    }
                                }                                
                            }        
                        }
                    }


                    sc.change.emit();
                    """)

        tax_menu.js_on_change('value', tax_menu_callback)  # calling the function on change of selection
        qseqid_menu.js_on_change('value', qseqid_menu_callback)

        p = figure(x_axis_label='bitscore', y_axis_label='pident',
                   plot_height=700, plot_width=700,
                   tooltips=TOOLTIPS,
                   tools="box_select, reset, box_zoom, pan", title="Number of RBHs - pident vs bitscore",
                   x_range=(0, result_data['bitscore'].max() + result_data['bitscore'].min())
                   )  # ,tools="box_select, reset" creating figure object
        circle = p.circle(x='bitscore', y='pident', color='color', size=5, line_width=1, line_color='black',
                          source=Curr,
                          legend_field=taxonomic_unit)  # plotting the data using glyph circle

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

        download_selection_callback = CustomJS(args=dict(sc=Curr), code="""
            var downloadable_items = []
            for(var i = 0; i < sc.selected.indices.length; i++){
                downloadable_items.push(sc.data['sacc'][sc.selected.indices[i]])
            }

            var json = JSON.stringify(downloadable_items);
            var blob = new Blob([json],{type: "octet/stream"});
            var url = window.URL.createObjectURL(blob);
            window.location.assign(url);
        """)

        download_selection_button = Button(label="Download Selection")
        download_selection_button.js_on_click(download_selection_callback)

        grid = gridplot(
            [[column(p),
              column(tax_menu, qseqid_menu, circle_size_spinner, line_size_spinner, line_color_picker, range_slider,
                     download_selection_button)]],
            toolbar_location='right', sizing_mode="stretch_both", merge_tools=True)

        output_file(filename=path_to_bokeh_plot,
                    title="Interactive Graph Percent Identity vs. Bitscore linked to {} database entries".format(
                        taxonomic_unit))
        save(grid)

        output_file(filename=path_to_static,
                    title="Interactive Graph Percent Identity vs. Bitscore linked to {} database entries".format(
                        taxonomic_unit))

        save(grid)

        return 0
    except Exception as e:
        raise Exception("ERROR in producing bokeh plots for database statistics with exception: {}".format(e))

def add_taxonomic_information_to_db(user_email:str,log:str,taxids:list)->pd.DataFrame:
    try:
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
                                log.write(
                                    "WARNING: problem during fetching of taxonomic information for: {}\n".format(m_taxid))
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

        columns = ["staxids", 'organism_name_taxdb', 'genus', 'family', 'superfamily', 'order', 'class', 'phylum']
        tax_db = pd.DataFrame([taxid, taxonomy, genus, family, superfamily, order, classt, phylum])
        tax_db = tax_db.transpose()
        tax_db.columns = columns
        return tax_db

    except Exception as e:
        log.write("ERROR:problem during fetching and appending taxonomic information with exception: {}\n".format(e))
        raise Exception("ERROR during database table extension with taxonomic information with excpetion: {}".format(e))


############################ MAIN SCRIPT ##############################

RETURNCODE=4
try:
    with open(snakemake.log['log'],'w') as logfile:
        logfile.write("INFO:starting to fetch taxonomic information...\n")
        #Entrez.email = snakemake.params['user_email']

        queries = {}
        queryfile = open(snakemake.input['query_file'], "r")
        for line in queryfile.readlines():
            if ">" in line:
                prot_id = line.split(">")[1].split(' ')[0]
                line = ' '.join(line.split(">")[1].split(' ')[1:]).rstrip()
                queries[prot_id] = line
        queryfile.close()

        df = pd.read_table(snakemake.input['blast_results'], delimiter="\t", header=None)
        df.columns = ["qseqid", "sseqid", "pident", "evalue", "bitscore", "qgi", "sgi", "sacc", "staxids", "sscinames", "scomnames",
                      "stitle"]
        unique_taxids = list(df["staxids"].unique())

        tax_df = add_taxonomic_information_to_db(snakemake.params['user_email'], logfile, unique_taxids)
        tax_df['staxids'] = tax_df['staxids'].astype('int64')
        result_df = df.merge(tax_df, on='staxids')


        dataframes = []

        logfile.write("INFO:parsing taxonomic information for queries\n")
        result_df.to_csv(snakemake.output['taxonomic_table'], sep='\t')
        logfile.write("INFO:start producing interactive bokeh plots for target species families\n")

        #bokeh code
        bokeh_function = create_unlinked_bokeh_plot(logfile,
                                                    snakemake.output['genus_bars'],
                                                    snakemake.params['genus_bars_static'],
                                                    result_df,'family')

        pd.set_option('colheader_justify', 'left')
        html_string = '''
        <html>
          <head>
            <title>BLAST Result Table</title>
            <!-- DataTables stylesheets-->
            <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.10.24/css/jquery.dataTables.css" crossorigin="anonymous">
            <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/select/1.3.2/css/select.dataTables.min.css" crossorigin="anonymous">
            <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/buttons/1.7.0/css/buttons.dataTables.min.css" crossorigin="anonymous">
          </head>
        
          <body>
            <div id="blast_results_table" style="display:none">
                {table}
            </div>
          </body>
        
            <script src="https://code.jquery.com/jquery-3.6.0.js" integrity="sha256-H+K7U5CnXl1h5ywQfKtSj8PCmoN9aaq30gDh27Xc0jk=" crossorigin="anonymous"></script>
            <!-- input scripts for DataTables: https://datatables.net/ -->
            <script type="text/javascript" charset="utf8" src="https://cdn.datatables.net/1.10.24/js/jquery.dataTables.js"></script>
            <script type="text/javascript" charset="utf8" src="https://cdn.datatables.net/select/1.3.2/js/dataTables.select.min.js"></script>
            <script type="text/javascript" charset="utf8" src="https://cdn.datatables.net/buttons/1.7.0/js/dataTables.buttons.min.js"></script>
            <script type="text/javascript" charset="utf8" src="https://cdn.datatables.net/buttons/1.7.0/js/buttons.html5.min.js"></script>
            <script>
            $(document).ready(function(){{
                var table = document.getElementsByTagName('table');
                table[0].id='myTable'
                $('#myTable').DataTable(
                    {{
                        dom: 'Bfrtip',
                        "lengthMenu": [ 100 ],
                        buttons: [
                            'copy',
                            'csv',
                            'selectAll',
                            'selectNone',
                            'selectRows'
                        ],
                        select: true
                    }}
                );
                var result_table = document.getElementById('blast_results_table');
                result_table.style.display = "block";
            }});
            </script>
        
        </html>
        '''

        # OUTPUT AN HTML FILE
        with open(snakemake.output['html_table'], 'w') as f:
            f.write(html_string.format(table=result_df.to_html(classes='mystyle')))
        logfile.write("DONE\n")
except Exception as e:
    logfile.write("ERROR:exception: {}\n".format(e))
    sys.exit(RETURNCODE)