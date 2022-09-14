import pandas as pd

result_df = pd.read_csv(snakemake.input['results_csv'], header=0, index_col=0)
try:
    sacc_list = list(result_df['sacc'].unique())

    output = open(snakemake.output['ids'], 'w')
    for sacc in sacc_list:
        output.write(sacc + '\n')
    output.close()

except Exception as e:
    raise Exception("[-] Error during writing hit info as txt with Exception: {}".format(e))
'''
try:
    queries = []
    with open(snakemake.input['query_file'],'r') as fhandle:
        for line in fhandle.readlines():
            if line[0] == '>':
                query = line.split('>')[1].split(' ')[0]
                if '.' in query:
                    query = query.split('.')[0]
                print(query)
                queries.append(query)
    for query in queries:
        target_df = result_df[result_df['qseqid'] == query]
        #target_df = target_df.nlargest(max_hits, 'bitscore')
        sacc_list = list(target_df['sacc'].unique())
        output_file_path = query + '/' + 'target_sequence_ids.txt'
        output = open(output_file_path,'w')
        for sacc in sacc_list:
            output.write(sacc+'\n')
        output.close()
except Exception as e:
    raise Exception("[-] Error during building of folders with hit info with Exception: {}".format(e))

'''