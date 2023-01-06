import pandas as pd
from Bio import Entrez

try:
    col = ["searches"]
    target_df = pd.read_csv(snakemake.input['ids'], header=None, names=col)
    Entrez.email = snakemake.params['user_email']
    target_list = target_df["searches"].drop_duplicates() #deletes duplicates

    end = len(target_list)
    begin = 0
    step = 500
    steps = 500
    records=[]
    while begin < end:
        if step >= end:
            step = end
        splitted_ids = target_list[begin:step]
        for attempt in range(10):
            try:
                handle = Entrez.efetch(id=splitted_ids, db="protein", retmode="xml")
                record = Entrez.read(handle)
                handle.close()
            except Exception as e:
                if attempt == 9:
                    raise Exception

            else:
                records.append(record)
            break
        begin += steps
        step += steps

    output = open(snakemake.output['fasta_file'],'w') # could use the name of the file used as input
    for record in records:
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