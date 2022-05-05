#this script is part of the snakemake worklfow it produces a csv file with all reciprocal hits.
#the script blast_tables_to_html.py has the same starting code due to clarity this script has its own snakemake rule
import pandas as pd
import sys
from Bio import Entrez

rec_prot=pd.read_table(snakemake.input['rec_res'],index_col=False)
fw_res=pd.read_table(snakemake.input['fw_res'],header=None)
fw_res.columns=["qseqid", "sseqid", "evalue", "bitscore", "qgi", "sgi", "sacc", "staxids", "sscinames", "scomnames",
                  "stitle"]

fw_res['qseqid'] = fw_res['qseqid'].map(lambda line: line.split('.')[0])
fw_res['sacc_transformed'] = fw_res['sacc'].map(lambda line: line.split('.')[0])
rec_prot = rec_prot.rename(columns={"forward_genome_id": "sacc_transformed"})
rec_prot = rec_prot.rename(columns={"backward_genome_id": "qseqid"})


try:
    result_data = rec_prot.merge(fw_res, how='inner', on=['sacc_transformed', 'qseqid', 'staxids'])
    # the backward blast is currently limited to output only the best match, but the best match can contain several hsps,
    # thus it is possible that there are multiple lines of one qseqid present, which gets loaded by reading the dictionary for
    # filtering reciprocal best hits
    result_data = result_data.drop_duplicates(['sacc_transformed', 'staxids'], keep='first')
    result_data = result_data.reset_index(drop=True)
except:
    sys.exit(123)


# Adds taxonomic information to pandas dataframe.
# In order to retrieve the taxonomic information biopython calls are conducted with the query accession ids in the pandas dataframe.
# The query file gets processed for extracting the fasta header line, this information is then loaded into a specific dataframe column.
def add_taxonomic_information_to_result_dataframe(result_data, query_file):
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


    df = result_data.copy()
    '''
    queries = read_query_file(query_file)

    dataframes = []
    unique_queries = list(df['qseqid'].unique())
    for query in unique_queries:
        # print("processing : {}".format(query))
        dataframe = df.loc[df['qseqid'] == query].copy()
        dataframe['sacc'] = dataframe['sacc'].map(lambda protid: protid.split(".")[0])
        dataframe = dataframe.drop_duplicates(subset=['sacc'], keep="first")
        # print("Current length: {}".format(len(dataframe)))
        if len(dataframe) >= 4999:
            dataframe = dataframe[0:4999]
            # print("New Length : {}".format(len(dataframe)))
        staxids = []
        for ids in list(dataframe['staxids']):
            if type(ids) == str:
                staxids.append(ids.split(";")[0])
            elif type(ids) == int:
                staxids.append(ids)

        result_record = []

        end = len(dataframe[dataframe['qseqid'] == query])
        begin = 0
        step = 500
        steps = 500
        while begin < end:
            if step >= end:
                step = end
            # print("\t {} to {}".format(begin, step))
            splitted_ids = staxids[begin:step]
            # range(X) tries for biopython calls
            for attempt in range(10):
                try:
                    # print("length ids : {}".format(len(splitted_ids)))
                    handle = Entrez.efetch(id=splitted_ids, db="taxonomy", retmode="xml")
                    record = Entrez.read(handle)
                    handle.close()
                except Exception as e:
                    # print("attempt : {} queries : {}".format(attempt, query))
                    if attempt == 9:
                        raise Exception('[-] Couldnt connect to NCBI server with Exception : {} after 10 tries ...'.format(e))

                else:
                    for rec in record:
                        result_record.append(rec)
                    # print("result record length : {}".format(len(result_record)))
                    break
            begin += steps
            step += steps

        query_info = []
        taxonomy = []
        genus = []
        superfamily = []
        family = []
        order = []
        classt = []
        phylum = []
        for i in range(len(result_record)):
            query_info.append(queries[query])
            taxonomy.append(result_record[i]['ScientificName'])
            # lineageEx = 'linex,'
            for j in result_record[i]['LineageEx']:
                # lineageEx += ','+str(j)
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

            # lineageExList.append(lineageEx)
            if (len(taxonomy) != len(genus)):
                genus.append('unknown')
            if (len(taxonomy) != len(superfamily)):
                superfamily.append('unknown')
            if (len(taxonomy) != len(family)):
                family.append('unknown')
            if (len(taxonomy) != len(order)):
                order.append('unknown')
            if (len(taxonomy) != len(phylum)):
                phylum.append('unknown')
            if (len(taxonomy) != len(classt)):
                classt.append('unknown')
        #print(len(genus) == len(taxonomy))

        # dataframe['taxonomic_name'] = taxonomy
        # compare length of all informations - they should have the same length
        if (len(genus) == len(dataframe) and len(family) == len(dataframe) and len(superfamily) == len(
                dataframe) and len(
                query_info) == len(dataframe)):
            dataframe['genus'] = genus
            dataframe['superfamily'] = superfamily
            dataframe['family'] = family
            dataframe['order'] = order
            dataframe['phylum'] = phylum
            dataframe['class'] = classt
            dataframe['query_info'] = query_info
            # dataframe['lineageEx'] = lineageExList
        else:
            raise Exception('[-] Informations from biopython calls have not the same length as the dataframe,\
                            dataframe cant get combined with taxonomic information')
        dataframes.append(dataframe)

    result_df = pd.concat(dataframes)
    '''
    return df #result_df

result_data.to_csv(snakemake.output['result_csv'],header=list(result_data.columns))
try:
    big_result_data = add_taxonomic_information_to_result_dataframe(result_data, snakemake.input['query_file'])
except:
    #returncode for error in taxonomic information adding
    sys.exit(124)
big_result_data.to_csv(snakemake.output['taxonomy_result_csv'],header=list(big_result_data.columns))