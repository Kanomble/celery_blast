import pandas as pd
import os

result_df = pd.read_csv(snakemake.input['result_csv'],header=0, index_col=0)

hit_info_file = open(snakemake.output['hit_information'],'w')
for qseqid in result_df['qseqid'].unique():
    try:
        os.mkdir(qseqid)
    except OSError as e:
        raise Exception("[-] Couldn't create folder for query sequence with exception : {}".format(e))

    target_df = result_df[result_df['qseqid'] == qseqid]
    sacc_list = list(target_df['sacc'].unique())

    output_file_path = qseqid + '/' + 'target_sequence_ids.txt'

    output = open(output_file_path,'w')
    for sacc in sacc_list:
        output.write(sacc+'\n')
    output.close()

    hit_info_file.write("qseqid: {} hits: {}\n".format(qseqid,len(target_df)))

hit_info_file.close()