from Bio import Entrez
from sys import exit

ERRORCODE = 3

with open(snakemake.log['log'],"w") as logfile:
    try:
        logfile.write("INFO:starting to fetch subject sequences for the backward BLAST.\n")
        logfile.write("INFO:setting up BioPython.\n")
        Entrez.email = snakemake.params['user_email']
        logfile.write("INFO:opening gi_list file {}.\n".format(snakemake.input['gi_list']))
        target_list = []
        with open(snakemake.input['gi_list']) as idfile:
            for line in idfile.readlines():
                line = line.rstrip()
                target_list.append(line)

        logfile.write("INFO:amount of target sequences: {}.\n".format(len(target_list)))
        logfile.write("INFO:fetching sequences in steps of 500 ... \n")

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

        gi_list = []
        with open(snakemake.output['fasta_file'], 'w') as output:
            for record in records:
                for rec in record:
                    output.write('>' + rec['GBSeq_primary-accession'] + ' ' + rec['GBSeq_definition'] + "\n")
                    output.write(rec['GBSeq_sequence'] + "\n")
                    gi_list.append(rec['GBSeq_sequence'])

            # this section checks if the query sequences are among the gi lists extracted by ENTREZ ...
            logfile.write("INFO:checking if query sequences are present in backward BLAST fasta file output ...\n")
            logfile.write("INFO:open query file ...\n")
            queries = {}
            with open(snakemake.input["fw_queries"], "r") as queryfile:
                for line in queryfile.readlines():
                    if ">" in line:
                        prot_id = line.rstrip().split(">")[1].split(' ')[0]
                        queries[prot_id] = ""
                    else:
                        queries[prot_id] += line.rstrip()

            for key in queries.keys():
                if key not in gi_list:
                    header = ">" + key + "\n"
                    sequence = queries[key] + "\n"
                    output.write(header)
                    output.write(sequence)

        logfile.write('INFO: The fasta file {} has been created.\n'.format(snakemake.output['fasta_file']))
        logfile.write("DONE")
    except Exception as e:
        logfile.write("ERROR: exception during extraction of subject sequences with BioPython: {}.\n".format(e))
        exit(ERRORCODE)