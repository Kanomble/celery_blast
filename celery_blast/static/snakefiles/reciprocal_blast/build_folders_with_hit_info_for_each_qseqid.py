import pandas as pd
import matplotlib.pyplot as plt
from sys import exit
from os.path import isdir
from os import mkdir

ERRORCODE=8

with open(snakemake.log['log'],'w') as logfile:
    try:
        logfile.write("INFO:starting to write RBH summary file and directories for query sequences\n")
        logfile.write("INFO:loading reciprocal result csv file as pandas dataframe\n")
        result_df = pd.read_csv(snakemake.input['result_csv'],header=0, index_col=0)


        queries = []
        logfile.write("INFO:extracting query sequence identifier based on query file\n")
        with open(snakemake.input['query_file'],'r') as fhandle:
            for line in fhandle.readlines():
                if line[0] == '>':
                    query = line.split('>')[1].split(' ')[0]
                    if '.' in query:
                        query = query.split('.')[0]
                    logfile.write("\tINFO:extracted query: {}\n".format(query))
                    queries.append(query)

        logfile.write("INFO:looping over query sequences to identify amount of RBHs and to build sub-directories\n")
        hit_info_file = open(snakemake.output['hit_information'], 'w')
        hit_info_file.write("qseqid\thits\n")
        for query in queries:
            logfile.write("\tINFO:working with:{}\n".format(query))
            target_df = result_df[result_df['qseqid'] == query]
            hit_info_file.write("{}\t{}\n".format(query, len(target_df)))

            logfile.write("\tINFO:slicing target RBH dataframe by defined maximum sequences for MSA and phylo task: {}\n".format(
                snakemake.config['max_rbhs_for_phylo']
            ))
            target_df = target_df.sort_values(by='bitscore', ascending=False).sort_index().copy()
            logfile.write("\t\tINFO:constructing tab separated file ... \n")
            cols = ['sacc', 'qseqid', 'bitscore', 'pident', 'evalue', 'slen', 'query_info', 'phylum', 'class',
                    'order', 'family', 'genus']
            target_df_to_tab = target_df.loc[:,cols]
            # omit sequences not present within the multiple alignment file ...
            target_df_to_tab = target_df_to_tab[target_df_to_tab.query_info != "artificially added RBH - not in FW database"]

            logfile.write("\tINFO:number of sequences in RBH table: {}\n".format(len(target_df_to_tab)))
            target_df_to_tab.to_csv(query+'/rbh_table.tsf', index=False, sep="\t")

            logfile.write("\t\tINFO:constructing DataTable CDN html string for pandas html table ...\n")
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
            logfile.write("\t\tINFO:writing result_rbhs.html file for {}\n".format(query))

            result_rbh_html_filepath = query + '/' + "results_rbhs.html"
            with open(result_rbh_html_filepath, 'w') as f:
                f.write(html_string.format(table=target_df[['sacc_transformed','scomnames','staxids','pident','bitscore','evalue',"slen",'stitle']].to_html(classes='mystyle')))


            logfile.write("\t\tINFO:producing statistic plot for query sequence results\n")
            fig, ax = plt.subplots(2, 2)
            ax[0, 0].hist(target_df['pident'], edgecolor="black", color='yellow')
            ax[0, 0].set_title("percent identity")
            ax[0, 1].hist(target_df['bitscore'], edgecolor="black", color='orange')
            ax[0, 1].set_title("bitscore")
            ax[1, 0].scatter(y=target_df['evalue'], x=range(len(target_df['evalue'])), edgecolor="black", color='red')
            ax[1, 0].set_title("evalue")
            ax[1, 1].bar(height=target_df.groupby('scomnames').size().sort_values(ascending=False)[0:12],
                         x=range(1, len(target_df.groupby('scomnames').size().index[0:12]) + 1))
            ax[1, 1].set_xticks(range(1, len(target_df.groupby('scomnames').size().index[0:12]) + 1))
            ax[1, 1].set_xticklabels(target_df.groupby('scomnames').size().index[0:12], rotation=90)

            try:
                plt.tight_layout()
            except Exception as e:
                logfile.write("\t\tWARNING:can't use  the plt.tight_layout() function with exception: {}\n".format(e))

            logfile.write("\t\tINFO:saving plots to project directory\n")
            result_statistics=str(query)+ "/basic_statistics.png"
            plt.savefig(result_statistics, dpi=400)

            logfile.write("\t\tINFO:writing target ids into target_sequence_ids.txt\n")
            output_file_path = query + '/' + 'target_sequence_ids.txt'
            with open(output_file_path,'w') as output:
                sacc_list = list(target_df['sacc'].unique())
                for sacc in sacc_list:
                    if sacc.split(".")[0] not in queries:
                        output.write(sacc+'\n')

            plt.close('all')
        hit_info_file.close()
        logfile.write("DONE\n")
    except Exception as e:
        logfile.write("ERROR:{}\n".format(e))
        exit(ERRORCODE)