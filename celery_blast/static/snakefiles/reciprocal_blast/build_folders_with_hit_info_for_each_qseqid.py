import pandas as pd

result_df = pd.read_csv(snakemake.input['result_csv'],header=0, index_col=0)

hit_info_file = open(snakemake.output['hit_information'],'w')
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
    sacc_list = list(target_df['sacc'].unique())

    output_file_path = query + '/' + 'target_sequence_ids.txt'

    output = open(output_file_path,'w')
    for sacc in sacc_list:
        output.write(sacc+'\n')
    output.close()

    hit_info_file.write("qseqid: {} hits: {}\n".format(query, len(target_df)))
hit_info_file.close()