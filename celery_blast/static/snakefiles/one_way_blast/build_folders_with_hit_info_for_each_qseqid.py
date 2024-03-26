import pandas as pd
x=3
#X is the amount of best hits we will allow
result_df = pd.read_csv(snakemake.input['results_csv'],header=0, index_col=0)
try:
    queries = []
    with open(snakemake.input['query_file'],'r') as fhandle:
        for line in fhandle.readlines():
            if line[0] == '>':
                prot_id = line.split(" ")[0].split(">")[1].strip()

                if "|" in prot_id:
                    prot_id = prot_id.split("|")[1]

                prot_id = prot_id.split(' ')[0].split(".")[0]
                queries.append(prot_id)

    for query in queries:
        target_df = result_df[result_df['qseqid'] == query]
        #target_df = target_df.nlargest(x, 'bitscore')
        sacc_list = list(target_df['sacc'].unique())

        output_file_path = query + '/' + 'target_sequence_ids.txt'

        output = open(output_file_path,'w')
        for sacc in sacc_list:
            output.write(sacc+'\n')
        output.close()
except Exception as e:
    raise Exception("[-] Error during building of folders with hit info with Exception: {}".format(e))