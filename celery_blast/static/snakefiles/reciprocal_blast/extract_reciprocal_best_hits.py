#this script is used to identify RBHs it converts the queries and subjects of the forward and backward BLAST output tables into dictionaries
#and compares both dictionaries in order to identify RBHs, a detailed description is given in the bachelor-thesis of this project
import pandas as pd
import numpy as np
import itertools as it


def get_seq_match_dict_and_flat_list(df):
    #extract protein identifier for matches: key=protein identifier for forward_input_sequences
    #and values=protein identifier for matches 
    seq_matches_dict = {}
    for i in df[0].unique():
        if "." in i:
            key_id = i.split(".")[0]
            seq_matches_dict[key_id] = list(set(df[df[0] == i][6]))
        else:
            seq_matches_dict[i] = list(set(df[df[0] == i][6]))

    return seq_matches_dict

def extract_reciprocal_best_hits_and_return_protein_ids(seq_matches_backward_dict,seq_matches_forward_dict):
    result_set = []
    for forward_key in seq_matches_forward_dict.keys():
        for forward_value in seq_matches_forward_dict[forward_key]:
            if forward_value in seq_matches_backward_dict.keys():
                if forward_key in seq_matches_backward_dict[forward_value]:
                    result_set.append([forward_value,seq_matches_backward_dict[forward_value]])
    return result_set


forward_df = pd.read_table(snakemake.input['fw_res'],header=None)
forward_df = pd.DataFrame([forward_df[0][:],forward_df[6][:]]).transpose()
fw_seqs = get_seq_match_dict_and_flat_list(forward_df)

backward_df = pd.read_table(snakemake.input['bw_res'],header=None)
backward_df = pd.DataFrame([backward_df[0][:],backward_df[6][:]]).transpose()
bw_seqs = get_seq_match_dict_and_flat_list(backward_df)
best_hits = extract_reciprocal_best_hits_and_return_protein_ids(bw_seqs,fw_seqs)

out = open(snakemake.output['rec_best_hits'],'w')
out.write('forward_genome_id\tbackward_genome_id\n')
for prot_id_pair in best_hits:
	out.write(str(prot_id_pair[0])+'\t'+str(prot_id_pair[1][0])+'\n')
out.close()
