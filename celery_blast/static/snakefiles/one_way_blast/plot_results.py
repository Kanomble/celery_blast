from Bio import Entrez
import pandas as pd
import math
import altair as alt

Entrez.email = "lukas.becker@hhu.de"

queries = {}
queryfile = open(snakemake.input['query_file'], "r")
for line in queryfile.readlines():
    if ">" in line:
        prot_id = line.split(">")[1].split(' ')[0]
        line = ' '.join(line.split(">")[1].split(' ')[1:]).rstrip()
        queries[prot_id] = line
queryfile.close()

df = pd.read_csv(snakemake.input['blast_results'], delimiter="\t", header=None)
df.columns = ["qseqid", "sseqid", "evalue", "bitscore", "qgi", "sgi", "sacc", "staxids", "sscinames", "scomnames",
              "stitle"]
unique_queries = list(df["qseqid"].unique())
dataframes = []

for query in unique_queries:
    # print("processing : {}".format(query))
    dataframe = df.loc[df['qseqid'] == query].copy()
    dataframe['sacc'] = dataframe['sacc'].map(lambda protid: protid.split(".")[0])
    dataframe = dataframe.drop_duplicates(subset=['sacc'], keep="first")
    # print("Current length: {}".format(len(dataframe)))
    if len(dataframe) >= 4999:
        dataframe = dataframe[0:4999]
        print("New Length : {}".format(len(dataframe)))

    staxids=[]
    for ids in list(dataframe['staxids']):
        if type(ids) == str:
            staxids.append(ids.split(";")[0])
        elif type(ids) == int:
            staxids.append(ids)

    result_record = []

    end = len(dataframe[dataframe['qseqid'] == query])
    begin = 0
    step = 500
    steps = 500
    while begin < end:
        if step >= end:
            step = end
        # print("\t {} to {}".format(begin,step))
        splitted_ids = staxids[begin:step]
        for attempt in range(10):
            try:
                # print("length ids : {}".format(len(splitted_ids)))
                handle = Entrez.efetch(id=splitted_ids, db="taxonomy", retmode="xml")
                record = Entrez.read(handle)
                handle.close()
            except Exception as e:
                # print("attempt : {} queries : {}".format(attempt,query))
                if attempt == 9:
                    raise Exception

            else:
                for rec in record:
                    result_record.append(rec)
                # print("result record length : {}".format(len(result_record)))
                break
        begin += steps
        step += steps

    query_info = []
    taxonomy = []
    genus = []
    superfamily = []
    family = []
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

#cols = math.ceil(math.sqrt(len(result_df['qseqid'].unique())))
rf = result_df.drop_duplicates(subset='genus', keep="first")
bars = []
selection = alt.selection_multi(fields=['genus'])
for df in dataframes:
    # make_selector = alt.Chart(df).mark_rect().encode(y='genus', color='genus').add_selection(selection)
    bar = alt.Chart(df).mark_bar(tooltip=True).encode(
        alt.X("count()"),
        alt.Y("genus"),
        color=alt.Color('genus', legend=None),
        tooltip=['count()', 'mean(bitscore)', 'mean(evalue)'],
    ).transform_filter(selection).interactive().facet(facet='query_info')

    graph = bar
    bars.append(graph)

make_selector = alt.Chart(rf).mark_rect().encode(y='genus', color='genus').add_selection(selection)
graphics = make_selector | bars[0]
if len(bars) > 1:
    for bar in bars[1:]:
        graphics |= bar

graphics.save(snakemake.output['genus_bars'])
graphics.save(snakemake.params['genus_bars_static'])

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
# with open('fw_results.html', 'w') as f:
#    f.write(html_string.format(table=fw_res.to_html(classes='mystyle')))

# with open('bw_results.html', 'w') as f:
#    f.write(html_string.format(table=bw_res.to_html(classes='mystyle')))

with open(snakemake.output['html_table'], 'w') as f:
    f.write(html_string.format(table=result_df.to_html(classes='mystyle')))