#this script writes the RBHs identified by the reciprocal BLAST pipeline into an html table
import pandas as pd
import sys

ERRORCODE=7
with open(snakemake.log['log'],'w') as logfile:
    try:
        logfile.write("INFO:transforming reciprocal results csv files into html tables\n")
        result_data = pd.read_csv(snakemake.input['rec_res'],header=0,index_col=0)
        '''
        for i in range(0, len(result_data), 1):
            taxids = result_data.iat[i, 7]
            scientific_names = result_data.iat[i, 8]
            common_names = result_data.iat[i, 9]
        '''
        logfile.write("INFO:preparing html string ...\n")

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
        with open(snakemake.output['rec_html'], 'w') as f:
            f.write(html_string.format(table=result_data.to_html(classes='mystyle')))
        logfile.write("DONE\n")
    except Exception as e:
        logfile.write("ERROR:{}\n".format(e))
        sys.exit(ERRORCODE)
