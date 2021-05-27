#this script writes the RBHs identified by the reciprocal BLAST pipeline into an html table
import pandas as pd

rec_prot=pd.read_table(snakemake.input['rec_res'])
fw_res=pd.read_table(snakemake.input['fw_res'],header=None)
fw_res.columns=["qseqid", "sseqid", "evalue", "bitscore", "qgi", "sgi", "sacc", "staxids", "sscinames", "scomnames",
                  "stitle"]

fw_res['qseqid'] = fw_res['qseqid'].map(lambda line: line.split('.')[0])
fw_res['sacc'] = fw_res['sacc'].map(lambda line: line.split('.')[0])
rec_prot = rec_prot.rename(columns={"forward_genome_id": "sacc"})
rec_prot = rec_prot.rename(columns={"backward_genome_id": "qseqid"})
result_data = rec_prot.merge(fw_res,how='inner', on=['sacc','qseqid'])
#the backward blast is currently limited to output only the best match, but the best match can contain several hsps,
#thus it is possible that there are multiple lines of one qseqid present, which gets loaded by reading the dictionary for
#filtering reciprocal best hits
result_data = result_data.drop_duplicates('sacc', keep='first')
result_data = result_data.reset_index(drop=True)

'''
for i in range(0, len(result_data), 1):
    taxids = result_data.iat[i, 7]
    scientific_names = result_data.iat[i, 8]
    common_names = result_data.iat[i, 9]
'''

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

with open(snakemake.output['rec_html'], 'w') as f:
    f.write(html_string.format(table=result_data.to_html(classes='mystyle')))