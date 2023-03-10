import pandas as pd
from sys import exit

ERRORCODE=17
with open(snakemake.log['log'] ,'w') as logfile:
    try:
        logfile.write("INFO:transforming query domain file into html table\n")
        result_data = pd.read_table(snakemake.input['query_domains'] ,header=None)
        columns = 'qseqid qlen sacc slen qstart qend sstart send qseq sseq bitscore evalue pident stitle'.split(" ")
        result_data.columns = columns
        result_data = result_data["qseqid qlen sacc slen qstart qend sstart send bitscore evalue pident stitle".split(" ")]

        logfile.write("INFO:preparing html string with DataTable CNNs ...\n")

        pd.set_option('colheader_justify', 'left')
        html_string = '''
        <html>
          <head>
            <title>Query Domain Table</title>
            <!-- DataTables stylesheets-->
            <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.10.24/css/jquery.dataTables.css" crossorigin="anonymous">
            <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/select/1.3.2/css/select.dataTables.min.css" crossorigin="anonymous">
            <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/buttons/1.7.0/css/buttons.dataTables.min.css" crossorigin="anonymous">
          </head>

          <body>
            <div id="query_domain_table" style="display:none">
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
                table[0].id='domainTable'
                $('#domainTable').DataTable(
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
                var result_table = document.getElementById('query_domain_table');
                result_table.style.display = "block";
            }});
            </script>

        </html>
        '''
        with open(snakemake.output['domain_html_table'], 'w') as f:
            f.write(html_string.format(table=result_data.to_html(classes='mystyle')))
        logfile.write("DONE\n")
    except Exception as e:
        logfile.write("ERROR:{}\n".format(e))
        exit(ERRORCODE)