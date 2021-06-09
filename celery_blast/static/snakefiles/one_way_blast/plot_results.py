from Bio import Entrez
import pandas as pd
import math
import altair as alt

Entrez.email = "lukas.becker@hhu.de"

queries = {}
queryfile = open(snakemake.input['query_file'], "r")
for line in queryfile.readlines():
    if ">" in line:
        prot_id = line.split(">")[1].split(' ')[0]
        line = ' '.join(line.split(">")[1].split(' ')[1:]).rstrip()
        queries[prot_id] = line
queryfile.close()

df = pd.read_csv(snakemake.input['blast_results'], delimiter="\t", header=None)
df.columns = ["qseqid", "sseqid", "evalue", "bitscore", "qgi", "sgi", "sacc", "staxids", "sscinames", "scomnames",
              "stitle"]
unique_queries = list(df["qseqid"].unique())
dataframes = []

for query in unique_queries:
    dataframe = df.loc[df['qseqid'] == query].copy()
    dataframe['sseqid'] = dataframe['sseqid'].map(lambda protid: protid.split(".")[0])
    dataframe = dataframe.drop_duplicates(subset=['sseqid'], keep="first")
    protein_ids = list(dataframe['sseqid'])


    handle = Entrez.elink(dbfrom="protein", db="taxonomy", id=protein_ids)
    record = Entrez.read(handle)
    handle.close()

    linked = []
    for i in range(len(record)):
        linked.append(record[i]["LinkSetDb"][0]["Link"][0]['Id'])

    handle = Entrez.efetch(id=linked, db="taxonomy", retmode="xml")
    record = Entrez.read(handle)
    handle.close()

    query_info = []
    taxonomy = []
    genus = []
    superfamily = []
    family = []
    for i in range(len(record)):
        query_info.append(queries[query])
        taxonomy.append(record[i]['ScientificName'])
        for j in record[i]['LineageEx']:
            if j['Rank'] == 'genus':
                genus.append(j['ScientificName'])
            if j['Rank'] == 'superfamily':
                superfamily.append(j['ScientificName'])
            if j['Rank'] == 'family':
                family.append(j['ScientificName'])
        if (len(taxonomy) != len(genus)):
            genus.append('unknown')
        if (len(taxonomy) != len(superfamily)):
            superfamily.append('unknown')
        if (len(taxonomy) != len(family)):
            superfamily.append('unknown')
    # print(len(genus) == len(taxonomy))
    # dataframe['taxonomic_name'] = taxonomy
    if (len(genus) == len(dataframe) and len(family) == len(dataframe) and len(superfamily) == len(dataframe) and len(
            query_info) == len(dataframe)):
        dataframe['genus'] = genus
        dataframe['superfamily'] = superfamily
        dataframe['family'] = family
        dataframe['query_info'] = query_info

    dataframes.append(dataframe)

result_df = pd.concat(dataframes)

cols = math.ceil(math.sqrt(len(result_df['qseqid'].unique())))
bar = alt.Chart(result_df).mark_bar().encode(
    y="count()",
    x="genus",
    color="genus",
    tooltip=["count()"]
).facet(facet='query_info', columns=cols)
bar.save(snakemake.output['genus_bars'])
bar.save(snakemake.params['genus_bars_static'])
