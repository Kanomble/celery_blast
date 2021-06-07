#this script writes the RBHs identified by the reciprocal BLAST pipeline into an html table
import pandas as pd


result_data=pd.read_table(snakemake.input['res'],header=None)
result_data.columns=["qseqid", "sseqid", "evalue", "bitscore", "qgi", "sgi", "sacc", "staxids", "sscinames", "scomnames",
                  "stitle"]

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
    f.write(html_string.format(table=result_data.to_html(classes='mystyle')))