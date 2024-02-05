'''build folders with filtered blast results for each query
input: blast_results_with_tax.table
output: blast_results_with_tax.table, target_sequence_ids.txt
This function filters the blast result table with taxonomic information for the best hits(10%) of each query and builds
a folder that contains the filtered dataframe as a csv as well as a text document that lists the names of each hit,
which is later used to perform the msa and phylo.
'''
import pandas as pd
import sys
import math
from os.path import isdir
from os import mkdir

RETURNCODE=5
try:
    result_data=pd.read_table(snakemake.input['results_table'], header=0, index_col=0)
    result_data['qseqid'] = result_data['qseqid'].map(lambda line: line.split('.')[0])
    result_data['sacc_transformed'] = result_data['sacc'].map(lambda line: line.split('.')[0])

    qseqids = pd.unique(result_data["qseqid"])
    for qseq in qseqids:
        qseq_df = result_data[result_data["qseqid"]==qseq]
        qseq_df = qseq_df.drop_duplicates(subset='sacc', keep='first')
        if len(qseq_df) > 50:
            best_hits = math.ceil(len(qseq_df)*0.1)
            qseq_df = qseq_df.nlargest(best_hits, "bitscore")
        elif len(qseq_df) > 30:
            best_hits = math.ceil(len(qseq_df)*0.50)
            qseq_df = qseq_df.nlargest(best_hits, "bitscore")

        if isdir(qseq) == False:
            mkdir(qseq)

        qseq_df = qseq_df.reset_index(drop=True)

        output_file_path = str(qseq) + '/' + 'blast_results_with_tax.table'
        qseq_df.to_csv(output_file_path,sep='\t', header=list(qseq_df.columns))

        try:
            sacc_list= list(qseq_df['sseqid'].unique())
            id_path= str(qseq) + '/' + 'target_sequence_ids.txt'
            output= open(id_path,'w')
            for sacc in sacc_list:
                output.write(sacc+'\n')
            output.close()

        except Exception as e:
            with open(snakemake.log['log'], 'w') as log_f:
                log_f.write('ERROR:' + "[-] Error during writing hit info as txt {} with Exception: {}".format(
                    snakemake.output['ids'], e))
            sys.exit(RETURNCODE)

        with open(snakemake.log['log'], 'w') as log_f:
            log_f.write('DONE: The results table has been filtered for the best hits, which where parsed to {}'.format(
                snakemake.output['ids']))

except Exception as e:
    with open(snakemake.log['log'], 'w') as log_f:
        log_f.write('ERROR:' + "[-] Error during parsing of csv: {} with Exception: {}".format(snakemake.output['result_csv'], e))
    sys.exit(RETURNCODE)
