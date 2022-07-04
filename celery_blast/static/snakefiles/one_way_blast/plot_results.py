'''plot results and update taxonomic information

input: blast_results.table, config['query_sequence']
output: genus_bars.html, blast_results.html, blast_results_with_tax.table

Expands blast results table with tax information.
Creates altair plots based on tax information.
Produces blast results in html file.
'''
#TODO fix print statements, add snakemake log

from Bio import Entrez
import pandas as pd
import altair as alt
from sys import exit

RETURNCODE=2
with open(snakemake.log['log'],'w') as logfile:
    try:
        Entrez.email = snakemake.params['user_email']

        #reading query file
        logfile.write("INFO:parsing query file for identifier and header lines\n")
        queries = {}
        queryfile = open(snakemake.input['query_file'], "r")
        for line in queryfile.readlines():
            if ">" in line:
                prot_id = line.split(">")[1].split(' ')[0]
                line = ' '.join(line.split(">")[1].split(' ')[1:]).rstrip()
                queries[prot_id] = line
        queryfile.close()

        logfile.write("INFO:reading BLAST result table\n")
        df = pd.read_table(snakemake.input['blast_results'], delimiter="\t", header=None)
        df.columns = ["qseqid", "sseqid", "pident", "evalue", "bitscore", "qgi", "sgi", "sacc", "staxids", "sscinames", "scomnames",
                      "stitle"]
        unique_queries = list(df["qseqid"].unique())
        dataframes = []

        #parsing result dataframe
        logfile.write("INFO:parsing BLAST table based on query identifier\n")
        for query in unique_queries:
            logfile.write("INFO:working on query: {}".format(query))
            dataframe = df.loc[df['qseqid'] == query].copy()
            dataframe['sacc'] = dataframe['sacc'].map(lambda protid: protid.split(".")[0])
            dataframe = dataframe.drop_duplicates(subset=['sacc'], keep="first")

            if len(dataframe) >= 4999:
                logfile.write("\tINFO:restrict dataframe size for altair - from {} to 5000 entries\n".format(
                    len(dataframe)))
                dataframe = dataframe[0:4999]


            logfile.write("\tINFO:reading taxonomic nodes of homologous sequences\n")
            staxids=[]
            for ids in list(dataframe['staxids']):
                if type(ids) == str:
                    staxids.append(ids.split(";")[0])
                elif type(ids) == int:
                    staxids.append(ids)

            result_record = []
            logfile.write("")
            end = len(dataframe[dataframe['qseqid'] == query])
            begin = 0
            step = 500
            steps = 500
            logfile.write("\tINFO:fetching taxonomic information with biopython\n")
            while begin < end:
                if step >= end:
                    step = end
                splitted_ids = staxids[begin:step]
                for attempt in range(10):
                    try:
                        handle = Entrez.efetch(id=splitted_ids, db="taxonomy", retmode="xml")
                        record = Entrez.read(handle)
                        handle.close()
                    except Exception as e:
                        if attempt == 9:
                            logfile.write("\tWARNING: ninth try for fetching taxonomic information with biopython, abort biopython call ...\n")
                            logfile.write("ERROR\n")
                            raise Exception

                    else:
                        for rec in record:
                            result_record.append(rec)


                        break
                begin += steps
                step += steps
            logfile.write("\tINFO:DONE")

            query_info = []
            taxonomy = []
            genus = []
            superfamily = []
            family = []
            logfile.write("INFO:attaching informations to result dataframe\n")
            for i in range(len(result_record)):
                query_info.append(queries[query])
                taxonomy.append(result_record[i]['ScientificName'])

                for j in result_record[i]['LineageEx']:
                    if j['Rank'] == 'genus':
                        genus.append(j['ScientificName'])
                    if j['Rank'] == 'superfamily':
                        superfamily.append(j['ScientificName'])
                    if j['Rank'] == 'family':
                        family.append(j['ScientificName'])

                if (len(taxonomy) != len(genus)):
                    genus.append('unknown')
                if (len(taxonomy) != len(superfamily)):
                    superfamily.append('unknown')
                if (len(taxonomy) != len(family)):
                    family.append('unknown')

            #print(len(genus) == len(taxonomy))
            # dataframe['taxonomic_name'] = taxonomy
            if (len(genus) == len(dataframe) and len(family) == len(dataframe) and len(superfamily) == len(dataframe) and len(
                    query_info) == len(dataframe)):
                # print("Yes!")
                dataframe['genus'] = genus
                dataframe['superfamily'] = superfamily
                dataframe['family'] = family
                dataframe['query_info'] = query_info
            else:
                # print("Nope!")
                break
            dataframes.append(dataframe)

        result_df = pd.concat(dataframes)
        logfile.write("INFO:writing taxonomic table\n")
        result_df.to_csv(snakemake.output['taxonomic_table'], sep='\t')

        #altair code
        logfile.write("INFO:using altair for taxonomy graphs\n")
        rf = result_df.drop_duplicates(subset='family', keep="first")
        bars = []
        selection = alt.selection_multi(fields=['family'])
        if len(dataframes) > 10:
            logfile.write("WARNING:restriction of altair plots to the first 10 query sequence entries\n")
        for df in dataframes[0:10]:
            bar = alt.Chart(df).mark_bar(tooltip=True).encode(
                alt.X("count()"),
                alt.Y("family"),
                color=alt.Color('genus', legend=None),
                tooltip=['count()', 'mean(bitscore)', 'mean(evalue)','genus'],
            ).transform_filter(selection).interactive().facet(facet='query_info')

            graph = bar
            bars.append(graph)

        make_selector = alt.Chart(rf).mark_rect().encode(y='family', color='family').add_selection(selection)
        graphics = make_selector | bars[0]
        if len(bars) > 1:
            for bar in bars[1:]:
                graphics |= bar

        graphics.save(snakemake.output['genus_bars'])
        graphics.save(snakemake.params['genus_bars_static'])

        #Requirements for altair
        charts_template = """
        <!DOCTYPE html>
        <html>
        <head>
          <style>
            #visCustom {{
              overflow-y: scroll;
              overflow-x: scroll;
              height: 50%;
              }}
          </style>
          <script src="https://cdn.jsdelivr.net/npm/vega@{vega_version}"></script>
          <script src="https://cdn.jsdelivr.net/npm/vega-lite@{vegalite_version}"></script>
          <script src="https://cdn.jsdelivr.net/npm/vega-embed@{vegaembed_version}"></script>
        </head>
        <body>
        
        <div id="visCustom"> 
            <div id="vis1"></div>
        </div>
        
        <script type="text/javascript">
          vegaEmbed('#vis1', {spec}).catch(console.error);
        </script>
        </body>
        </html>
        """
        with open(snakemake.output['genus_bars'], 'w') as f:
            f.write(charts_template.format(
                vega_version=alt.VEGA_VERSION,
                vegalite_version=alt.VEGALITE_VERSION,
                vegaembed_version=alt.VEGAEMBED_VERSION,
                spec=graphics.to_json(indent=None),
            ))
        with open(snakemake.params['genus_bars_static'], 'w') as f:
            f.write(charts_template.format(
                vega_version=alt.VEGA_VERSION,
                vegalite_version=alt.VEGALITE_VERSION,
                vegaembed_version=alt.VEGAEMBED_VERSION,
                spec=graphics.to_json(indent=None),
            ))


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
    except Exception as e:
        logfile.write("ERROR:{}\n".format(e))
        exit(RETURNCODE)