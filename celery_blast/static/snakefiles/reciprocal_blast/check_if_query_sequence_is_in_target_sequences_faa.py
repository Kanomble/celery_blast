'''
This script checks if the query sequence resides in the target sequence fasta file.
If it is in the fasta file nothing will be changed except the name of the fasta file,
if the query sequence does not reside in the fasta file, the sequence gets added.
'''

from sys import exit
ERRORCODE=16

with open(snakemake.log['log'],'w') as logfile:
    try:
        logfile.write("INFO:starting to check if query sequence resides in target sequence file\n")
        qseqids=snakemake.params['qseqs']
        logfile.write("INFO:working with {} qseqids\n".format(len(qseqids)))
        logfile.write("INFO:parsing query file and building id to sequence dictionary\n")
        query_dict={}

        with open(snakemake.input['query_file'],'r') as queryfile:
            for line in queryfile.readlines():
                if line.startswith(">"):
                    qseq_id = line.split(" ")[0].split(">")[1].split(".")[0]
                    if qseq_id not in qseqids:
                        logfile.write("WARNING:{} does not reside in the snakemake QSEQS list\n".format(qseq_id))
                        logfile.write("WARNING:have you deleted appended a query sequence manually?\n")
                        logfile.write("ERROR")
                        exit(ERRORCODE)
                    query_dict[qseq_id] = line
                else:
                    query_dict[qseq_id] += line

        logfile.write("INFO:trying to find matching ids\n")

        for qseq_id in qseqids:
            with open(str(qseq_id)+'/target_sequences_raw.faa','r') as targetfile:
                lines=targetfile.readlines()
                logfile.write("\tINFO:working with target sequence file - total number of lines: {}\n".format(len(lines)))
                with open(str(qseq_id)+'/target_sequences.faa','w') as new_targetfile:
                    target_ids = []
                    for line in lines:
                        if line.startswith(">"):
                            target_id = line.split(" ")[0].split(">")[1].split(".")[0]
                            target_ids.append(target_id)
                            new_targetfile.write(line)
                        else:
                            new_targetfile.write(line)
                    if qseq_id not in target_ids:
                        logfile.write("INFO:query sequence {} does not reside in target fasta file\n".format(qseq_id))
                        logfile.write("INFO:trying to attach query sequence to target_sequences.faa\n")
                        new_targetfile.write(query_dict[qseq_id])

        logfile.write("DONE\n")
    except Exception as e:
        logfile.write("ERROR:{}\n".format(e))
        exit(ERRORCODE)