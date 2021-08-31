import pandas as pd

result_df = pd.read_csv(snakemake.input['result_csv'],header=0, index_col=0)
hit_info_file = open(snakemake.output['hit_information'],'w')
for query in result_df['qseqid'].unique():

    target_df = result_df[result_df['qseqid'] == query]
    sacc_list = list(target_df['sacc'].unique())

    output_file_path = query + '/' + 'target_sequence_ids.txt'

    output = open(output_file_path,'w')
    for sacc in sacc_list:
        output.write(sacc+'\n')
    output.close()

    hit_info_file.write("qseqid: {} hits: {}\n".format(query, len(target_df)))

hit_info_file.close()