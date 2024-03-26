'''build folders with filtered blast results for each query
input: blast_results_with_tax.table
output: blast_results_with_tax.table, target_sequence_ids.txt
This function filters the blast result table with taxonomic information for the best hits(10%) of each query and builds
a folder that contains the filtered dataframe as a csv as well as a text document that lists the names of each hit,
which is later used to perform the msa and phylo.
'''
#TODO rework limiting
import pandas as pd
import sys
import math
from os.path import isdir
from os import mkdir

RETURNCODE=5
with open(snakemake.log[0], "w") as log_f:

    try:
        log_f.write("INFO:starting to filter blast results ...\n")

        result_data=pd.read_table(snakemake.input['results_table'], delimiter="\t", header=None)
        result_data.columns = ["qseqid", "sseqid", "pident", "evalue", "bitscore", "slen", "qgi", "sgi", "sacc", "staxids", "sscinames",
                      "scomnames",
                      "stitle"]
        result_data['qseqid'] = result_data['qseqid'].map(lambda line: line.split('.')[0])
        result_data['sacc_transformed'] = result_data['sacc'].map(lambda line: line.split('.')[0])

        no_hits = False
        log_f.write("INFO:checking hits ...\n")

        if len(result_data["qseqid"]) == 1:
            if "no information" in list(result_data["qseqid"]):
                no_hits = True

        if no_hits == False:
            log_f.write("INFO:found {} hits in the result table.\n".format(len(result_data["qseqid"])))
        else:
            log_f.write("INFO:there are no hits for you query proteins.\n")

        queries = {}
        queryfile = open(snakemake.input['query_file'], "r")
        for line in queryfile.readlines():
            if ">" in line:
                prot_id = line.split(" ")[0].split(">")[1].strip()
                if "|" in prot_id:
                    prot_id = prot_id.split("|")[1]
                prot_id = prot_id.split(' ')[0].split(".")[0]

                queries[prot_id] = ""
            queries[prot_id] += line.rstrip()

        queryfile.close()

        qseqids = queries.keys()
        for qseq in qseqids:
            qseq_df = result_data[result_data["qseqid"]==qseq]
            log_f.write("INFO:working with query: {}\n".format(qseq))
            log_f.write("\tINFO:found {} hits for query sequence\n".format(len(qseq_df)))
            qseq_df = qseq_df.drop_duplicates(subset='sacc', keep='first')
            if len(qseq_df) > 50:
                best_hits = math.ceil(len(qseq_df)*0.1)
                qseq_df = qseq_df.nlargest(best_hits, "bitscore")
            elif len(qseq_df) > 30:
                best_hits = math.ceil(len(qseq_df)*0.50)
                qseq_df = qseq_df.nlargest(best_hits, "bitscore")

            log_f.write("\tINFO:checking if directory for {} exists ...\n".format(qseq))
            if isdir(qseq) == False:
                mkdir(qseq)

            qseq_df = qseq_df.reset_index(drop=True)

            output_file_path = str(qseq) + '/' + 'blast_results.table'
            log_f.write("\tINFO:writing output csv: {} to disc ...\n".format(output_file_path))

            qseq_df.to_csv(output_file_path,sep='\t', header=list(qseq_df.columns))

            try:
                sacc_list= list(qseq_df['sacc'].unique())
                id_path= str(qseq) + '/' + 'target_sequence_ids.txt'
                output= open(id_path,'w')
                for sacc in sacc_list:
                    output.write(sacc+'\n')
                output.close()

            except Exception as e:
                log_f.write('ERROR:' + " Error during writing hit info as txt {} with Exception: {}\n".format(id_path,e))
                sys.exit(RETURNCODE)
            log_f.write('DONE: The results table has been filtered for the best hits, which where parsed to {}\n'.format(id_path))

    except Exception as e:
        log_f.write('ERROR:'+"[-] Error during parsing of csv: {} with Exception: {}\n".format(output_file_path,e))
        sys.exit(RETURNCODE)