from Bio import Entrez

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

        with open(snakemake.output['fasta_file'], 'w') as output:
            for record in records:
                for rec in record:
                    output.write('>' + rec['GBSeq_locus'] + ' ' + rec['GBSeq_definition'] + "\n")
                    output.write(rec['GBSeq_sequence'] + "\n")

        logfile.write('INFO: The fasta file {} has been created.\n'.format(snakemake.output['fasta_file']))
        logfile.write("DONE")
    except Exception as e:
        logfile.write("INFO: exception during extraction of subject sequences with BioPython: {}.\n".format(e))
        raise Exception("[-] ERROR during subject sequence download with exception: {}".format(e))