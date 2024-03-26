'''extract_reciprocal_best_hits
This script extracts reciprocal best hits (RBHs).
Based on the blast_fw_df (forward blast) and blast_bw_df (backward blast) output files of the rules for the
forward_blast and backward_blast, RBHs are inferred.
The inference is based on the pandas merge function.
Before RBH inference the forward dataframe is filtered based on the user provided bitscore threshold,
the default threshold is 50.
The columns for the qseqid and targetid as well as the corresponding taxid column are used for the merge comparison.
'''
import pandas as pd
from sys import exit
ERRORCODE=5
try:
    with open(snakemake.log['log'],'w') as logfile:
        forward_df = pd.read_table(snakemake.input['fw_res'],header=None)
        forward_df[7] = forward_df[7].map(lambda line : line.split('.')[0])
        logfile.write("INFO:starting to extract reciprocal best hits (RBHs) from blast output tables\n")
        logfile.write("INFO:loaded forward dataframe into pandas with length {}\n".format(len(forward_df)))

        #apply bitscore filter
        if len(forward_df) <= 0:
            logfile.write("WARNING:there are no hits in the forward BLAST table...\n")
            logfile.write("WARNING:adding original query sequence as hit ...\n")
        else:
            forward_df = forward_df[forward_df[4] >= snakemake.params['bitscore_filter']]
            if len(forward_df) <= 0:
                logfile.write("WARNING:after applying the bitscore filter there are no hits in the forward BLAST table...\n")

        forward_df = pd.DataFrame([forward_df[0][:],forward_df[7][:], forward_df[8][:]]).transpose()
        backward_df = pd.read_table(snakemake.input['bw_res'],header=None)
        backward_df[7] = backward_df[7].map(lambda line : line.split('.')[0])
        backward_df = pd.DataFrame([backward_df[0][:],backward_df[7][:], backward_df[8][:]]).transpose()
        forward_df[0] = forward_df[0].map(lambda line : line.split('.')[0])
        backward_df[0] = backward_df[0].map(lambda line: '_'.join(line.split("_")[0:2]).split(".")[0])
        logfile.write("INFO:loaded backward dataframe into pandas with length {}\n".format(len(backward_df)))


        result_df = pd.DataFrame(columns=["targetid","qseqid","taxid_x","taxid_y"])

        for qseq in forward_df[0].unique():
            logfile.write("INFO:analysis of reciprocal best hits for {}\n".format(qseq))

            qseq_df_fw = forward_df[forward_df[0] == qseq]
            qseq_df_bw = backward_df[backward_df[7] == qseq]
            qseq_df_fw.columns = ['qseqid','targetid','taxid']
            qseq_df_bw.columns = ['targetid','qseqid','taxid']
            qseq_result_df = qseq_df_bw.merge(qseq_df_fw,how='inner',on=['targetid','qseqid']).drop_duplicates()
            logfile.write("\tINFO:found {} RBHs\n".format(len(qseq_result_df)))
            result_df = pd.concat([result_df,qseq_result_df],ignore_index=True)



        # returncode for no reciprocal hits
        if len(result_df['targetid']) == 0:
            logfile.write("ERROR:there are no reciprocal best hits for the provided query sequences\n")
            exit(ERRORCODE)
        logfile.write("INFO:generating RBH output table with following headers, separated by a tab: forward_genome_id\tbackward_genome_id\tstaxids\n")
        # generating output file
        with open(snakemake.output['rec_best_hits'],'w') as recfile:
            recfile.write("forward_genome_id\tbackward_genome_id\tstaxids\n")
            for targetid, qseqid,taxid in zip(result_df['targetid'],result_df['qseqid'],result_df['taxid_y']):
                recfile.write("{}\t{}\t{}\n".format(targetid,qseqid, taxid))
        logfile.write("DONE\n")
except Exception as e:
    logfile.write("ERROR:something unexpected happened - exception {}\n".format(e))
    exit(ERRORCODE)


