import pandas as pd
import re
'''seq_id_to_taxinfo
    this function writes a new newick file with informations from the old newick file and reciprocal
    result dataframe. the re package is used to identify sequence ids in the old newick file, the ids are
    used to parse the reciprocal result dataframe for identification of related scientific names
    :param newick_treefile
        :type str
            :description path to newick treefile
    :param result_df
        :type str
            :description path to reciprocal blast result dataframe
'''
newick_treefile= snakemake.input['nw']
result_df= snakemake.input['df']
try:
    #reading newick treefile
    with open(newick_treefile,'r') as treefile:
        treefile_content = treefile.readlines()
    #searching for scientific names
    hits = re.findall('[(,]{1}([a-zA-Z_0-9]+[a-zA-Z_0-9._]+):{1}',treefile_content[0])
    hits = [hit.split(':')[0] for hit in hits]
    result_df = pd.read_csv(result_df, index_col=0)
    hit_to_name = {}
    scientific_names = result_df[result_df['sacc'].isin(hits)][['sscinames', 'sacc']]
    for name, sacc in zip(scientific_names['sscinames'], scientific_names['sacc']):
        if type(name) == str:
            hit_to_name[sacc] = '_'.join(name.split(' ')) + '_' + sacc
        else:
            hit_to_name[sacc] = 'NaN'

    # print(scientific_names)
    # scientific_names = ['_'.join(name.split(' '))+ '_'+ sacc if type(name) == str else 'NaN' for name, sacc in zip(scientific_names['sscinames'],scientific_names['sacc']) ]

    # check for duplicate scientific names and add counter to distinguish between them,
    # scientific_names includes the sscienames and sacc for more precise tree labeling

    # replacing sequence identifier with scientific names
    for hit in hit_to_name.keys():
        treefile_content[0] = re.sub(hit, hit_to_name[hit], treefile_content[0])
    #writing new newick file with changes header
    new_treefile = newick_treefile.split("/")
    new_treefile = new_treefile[0:len(new_treefile)-1]
    new_treefile = '/'.join(new_treefile)
    new_treefile = new_treefile + '/transformed_treefile.nwk'
    with open(new_treefile,'w') as treefile:
        treefile.write(treefile_content[0])
    #return treefile_content[0]
except Exception as e:
    raise Exception("[-] Error during parsing of newick file: {} with Exception: {}".format(newick_treefile,e))