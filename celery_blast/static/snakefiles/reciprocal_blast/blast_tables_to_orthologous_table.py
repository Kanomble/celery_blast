#this script is part of the snakemake worklfow it produces a csv file with all reciprocal hits.
#the script blast_tables_to_html.py has the same starting code due to clarity this script has its own snakemake rule
import pandas as pd

rec_prot=pd.read_table(snakemake.input['rec_res'])
fw_res=pd.read_table(snakemake.input['fw_res'],header=None)
fw_res.columns=["qseqid", "sseqid", "evalue", "bitscore", "qgi", "sgi", "sacc", "staxids", "sscinames", "scomnames",
                  "stitle"]

fw_res['qseqid'] = fw_res['qseqid'].map(lambda line: line.split('.')[0])
fw_res['sacc_transformed'] = fw_res['sacc'].map(lambda line: line.split('.')[0])
rec_prot = rec_prot.rename(columns={"forward_genome_id": "sacc_transformed"})
rec_prot = rec_prot.rename(columns={"backward_genome_id": "qseqid"})
result_data = rec_prot.merge(fw_res,how='inner', on=['sacc_transformed','qseqid'])
#the backward blast is currently limited to output only the best match, but the best match can contain several hsps,
#thus it is possible that there are multiple lines of one qseqid present, which gets loaded by reading the dictionary for
#filtering reciprocal best hits
result_data = result_data.drop_duplicates('sacc_transformed', keep='first')
result_data = result_data.reset_index(drop=True)

result_data.to_csv(snakemake.output['result_csv'],header=list(result_data.columns))