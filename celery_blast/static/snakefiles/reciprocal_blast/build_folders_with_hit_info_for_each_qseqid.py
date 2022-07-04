import pandas as pd
from sys import exit
ERRORCODE=8

with open(snakemake.log['log'],'w') as logfile:
    try:
        logfile.write("INFO:starting to write RBH summary file and directories for query sequences\n")
        result_df = pd.read_csv(snakemake.input['result_csv'],header=0, index_col=0)

        hit_info_file = open(snakemake.output['hit_information'],'w')
        queries = []
        logfile.write("INFO:extracting query sequence identifier based on query file\n")
        with open(snakemake.input['query_file'],'r') as fhandle:
            for line in fhandle.readlines():
                if line[0] == '>':
                    query = line.split('>')[1].split(' ')[0]
                    if '.' in query:
                        query = query.split('.')[0]
                    logfile.write("\textracted query: {}\n".format(query))
                    queries.append(query)

        logfile.write("INFO:looping over query sequences to identify amount of RBHs and to build directories\n")
        for query in queries:
            target_df = result_df[result_df['qseqid'] == query]
            sacc_list = list(target_df['sacc'].unique())

            output_file_path = query + '/' + 'target_sequence_ids.txt'

            output = open(output_file_path,'w')
            for sacc in sacc_list:
                output.write(sacc+'\n')
            output.close()

            hit_info_file.write("qseqid: {} hits: {}\n".format(query, len(target_df)))
            logfile.write("\tINFO:working with:{}\n".format(query))

        hit_info_file.close()
        logfile.write("DONE\n")
    except Exception as e:
        logfile.write("ERROR:{}\n".format(e))
        exit(ERRORCODE)