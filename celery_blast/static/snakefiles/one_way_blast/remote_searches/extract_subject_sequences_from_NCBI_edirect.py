import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio import Entrez


try:
    col = ["searches"]
    target_df = pd.read_csv(snakemake.input['ids'], header=None, names=col)
    Entrez.email = snakemake.params['user_email']
    target_list = target_df["searches"].drop_duplicates() #deletes duplicates
    handle = Entrez.efetch(db="protein", id=target_list, retmode="xml")
    record = Entrez.read(handle)
    handle.close()
    output = open(snakemake.output['fasta_file'],'w') # could use the name of the file used as input
    for rec in record:
        output.write('>'+ rec['GBSeq_locus'] + ' ' + rec['GBSeq_definition']+"\n")
        output.write(rec['GBSeq_sequence']+"\n")
    output.close()
    with open(snakemake.log['log'], 'w') as log_f:
        log_f.write('DONE: The fasta file {} has been parsed.'.format(snakemake.output['fasta_file']))
except Exception as e:
    with open(snakemake.log['log'], 'w') as log_f:
        log_f.write("ERROR:[-] Error during parsing of fasta: {} with Exception: {}".format(snakemake.output['fasta_file'],e))
    raise Exception("[-] Error during parsing of fasta: {} with Exception: {}.".format(snakemake.output['fasta_file'],e))