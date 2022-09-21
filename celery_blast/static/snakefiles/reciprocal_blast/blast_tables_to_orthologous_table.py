#this script is part of the snakemake worklfow it produces a csv file with all reciprocal hits.
#the script blast_tables_to_html.py has the same starting code due to clarity this script has its own snakemake rule
import pandas as pd
import sys
from Bio import Entrez

ERRORCODE=6
with open(snakemake.log['log'], 'w') as logfile:
    try:
        logfile.write("INFO:starting to construct result dataframe ...\n")
        Entrez.email = snakemake.config['user_email']
        rec_prot = pd.read_table(snakemake.input['rec_res'], index_col=False)
        fw_res = pd.read_table(snakemake.input['fw_res'], header=None)
        fw_res.columns = ["qseqid", "sseqid", "pident", "evalue", "bitscore", "qgi", "sgi", "sacc", "staxids",
                          "sscinames", "scomnames",
                          "stitle"]

        fw_res['qseqid'] = fw_res['qseqid'].map(lambda line: line.split('.')[0])
        fw_res['sacc_transformed'] = fw_res['sacc'].map(lambda line: line.split('.')[0])
        rec_prot = rec_prot.rename(columns={"forward_genome_id": "sacc_transformed"})
        rec_prot = rec_prot.rename(columns={"backward_genome_id": "qseqid"})

        result_data = rec_prot.merge(fw_res, how='inner', on=['sacc_transformed', 'qseqid', 'staxids'])
        # the backward blast is currently limited to output only the best match, but the best match can contain several hsps,
        # thus it is possible that there are multiple lines of one qseqid present, which gets loaded by reading the dictionary for
        # filtering reciprocal best hits
        result_data = result_data.drop_duplicates(['sacc_transformed', 'staxids'], keep='first')
        result_data = result_data.reset_index(drop=True)


        # Adds taxonomic information to pandas dataframe.
        # In order to retrieve the taxonomic information biopython calls are conducted with the query accession ids in the pandas dataframe.
        # The query file gets processed for extracting the fasta header line, this information is then loaded into a specific dataframe column.
        def add_taxonomic_information_to_result_dataframe(taxids, query_file):
            def read_query_file(query_file):
                queries = {}
                queryfile = open(query_file, "r")
                for line in queryfile.readlines():
                    if ">" in line:
                        prot_id = line.split(">")[1].split(' ')[0].split('.')[0]
                        line = ' '.join(line.split(">")[1].split(' ')[1:]).rstrip()
                        queries[prot_id] = line
                queryfile.close()
                return queries

            queries = read_query_file(query_file)

            logfile.write("INFO:working on {} taxids\n".format(len(taxids)))

            query_info = []
            taxid = []
            taxonomy = []
            genus = []
            superfamily = []
            family = []
            order = []
            classt = []
            phylum = []

            end = len(taxids)
            begin = 0
            step = 500
            steps = 500
            logfile.write("INFO:inference of taxonomic informations for all {} taxids\n".format(len(taxids)))
            while begin < end:
                if step >= end:
                    step = end
                splitted_ids = taxids[begin:step]
                # range(X) tries for biopython calls
                for attempt in range(10):
                    try:
                        handle = Entrez.efetch(id=splitted_ids, db="taxonomy", retmode="xml")
                        record = Entrez.read(handle)
                        handle.close()
                    except Exception as e:
                        if attempt == 9:
                            logfile.write(
                                "ERROR:inference of taxonomic informations failed with exception {}\n".format(e))
                            sys.exit(ERRORCODE)
                    else:
                        for i in range(len(record)):
                            taxonomy.append(record[i]['ScientificName'])
                            if 'AkaTaxIds' in record[i].keys():
                                for akaid in record[i]['AkaTaxIds']:
                                    if int(akaid) in splitted_ids:
                                        taxid.append(akaid)
                                        logfile.write("\tINFO: AkaTaxIds detected: {}\n".format(akaid))
                                        break
                                else:
                                    taxid.append(record[i]['TaxId'])
                            else:
                                taxid.append(record[i]['TaxId'])
                            for j in record[i]['LineageEx']:
                                if j['Rank'] == 'genus':
                                    genus.append(j['ScientificName'])
                                if j['Rank'] == 'superfamily':
                                    superfamily.append(j['ScientificName'])
                                if j['Rank'] == 'family':
                                    family.append(j['ScientificName'])
                                if j['Rank'] == 'order':
                                    order.append(j['ScientificName'])
                                if j['Rank'] == 'class':
                                    classt.append(j['ScientificName'])
                                if j['Rank'] == 'phylum':
                                    phylum.append(j['ScientificName'])

                            if (len(taxonomy) != len(genus)):
                                genus.append('unknown')
                            if (len(taxonomy) != len(superfamily)):
                                superfamily.append('unknown')
                            if (len(taxonomy) != len(family)):
                                family.append('unknown')
                            if (len(taxonomy) != len(order)):
                                order.append('unknown')
                            if (len(taxonomy) != len(classt)):
                                classt.append('unknown')
                            if (len(taxonomy) != len(phylum)):
                                phylum.append('unknown')

                            if len(record) != len(splitted_ids):
                                missing_ids = [m_taxid for m_taxid in splitted_ids if m_taxid not in taxid]
                                for m_taxid in missing_ids:
                                    logfile.write(
                                        "WARNING: problem during fetching of taxonomic information for: {}\n".format(
                                            m_taxid))
                                    taxid.append(m_taxid)
                                    taxonomy.append('unknown')
                                    genus.append('unknown')
                                    superfamily.append('unknown')
                                    family.append('unknown')
                                    order.append('unknown')
                                    classt.append('unknown')
                                    phylum.append('unknown')
                        break
                logfile.write("INFO: Done with chunk: {} - {}\n".format(begin, step))

                begin += steps
                step += steps

            # compare length of all informations - they should have the same length
            if (len(genus) == len(taxids) and len(family) == len(taxids) and len(superfamily) == len(taxids)):
                columns = ["staxids", 'organism_name_taxdb', 'genus', 'family', 'superfamily', 'order', 'class',
                           'phylum']
                tax_db = pd.DataFrame([taxid, taxonomy, genus, family, superfamily, order, classt, phylum])
                tax_db = tax_db.transpose()
                tax_db.columns = columns
                return tax_db
            else:
                logfile.write(
                    "ERROR:Informations from biopython calls have not the same length as the dataframe,\
                                dataframe cant get combined with taxonomic information'\n")
                sys.exit(ERRORCODE)


        taxids_from_df = result_data['staxids'].unique()
        logfile.write("INFO:constructing result dataframe\n")
        result_data.to_csv(snakemake.output['result_csv'], header=list(result_data.columns))
        taxonomy_dataframe = add_taxonomic_information_to_result_dataframe(taxids_from_df, snakemake.input['query_file'])
        taxonomy_dataframe['staxids'] = taxonomy_dataframe['staxids'].astype('int64')
        result_data = result_data.merge(taxonomy_dataframe, on='staxids')
        result_data['query_info'] = result_data['stitle']

        logfile.write("INFO:writing result dataframe\n")
        logfile.write("DONE\n")
        logfile.write("INFO:constructing result dataframe\n")

        logfile.write("INFO:writing result dataframe\n")
        result_data.to_csv(snakemake.output['taxonomy_result_csv'], header=list(result_data.columns))
        logfile.write("DONE\n")

    except Exception as e:
        logfile.write("ERROR:{}\n".format(e))
        sys.exit(ERRORCODE)