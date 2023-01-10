#TODO refactor
import pandas as pd
import re
import sys
'''seq_id_to_taxinfo
    this function writes a new newick file with information from the old newick file and reciprocal
    result dataframe. the re package is used to identify sequence ids in the old newick file, the ids are
    used to parse the reciprocal result dataframe for identification of related scientific names
    :param newick_treefile
        :type str
            :description path to newick treefile
    :param result_df
        :type str
            :description path to reciprocal blast result dataframe
'''

RETURNCODE=6
newick_treefile= snakemake.input['nw']
result_df= snakemake.input['df']
log_file = snakemake.log['log']
try:
    #reading newick treefile
    with open(newick_treefile,'r') as treefile:
        treefile_content = treefile.readlines()
    # no content in file logging and parsing
    if len(treefile_content) < 1:
        with open(log_file, 'w') as log_f:
            log_f.write('ERROR: The Newick file used in this command has no content.')
        with open(snakemake.output['tree_file'], 'w') as new_file:
            new_file.write('')
    else:
        #searching for scientific names
        hits = re.findall('[(,]{1}([a-zA-Z_0-9]+[a-zA-Z_0-9._-]+):{1}', treefile_content[0])
        hits = [hit.split(':')[0] for hit in hits]
        result_df = pd.read_table(result_df, index_col=0)
        hit_to_name = {}
        scientific_names = result_df[result_df['sseqid'].isin(hits)][['sscinames', 'sacc', 'stitle', 'sseqid']]
        # it is important that scientific names does not have ',' or ':' in its entries because the treefile uses these to distinguish between different entries
        scientific_names['stitle'] = scientific_names['stitle'].replace(np.nan, 'unknown')
        for title in scientific_names['stitle']:
            if ',' in title:
                new_title = title.replace(',', '-')
                scientific_names['stitle'] = scientific_names['stitle'].replace(title, new_title)
        for title in scientific_names['stitle']:
            if '[' in title:
                new_title = title.replace('[', '_')
                scientific_names['stitle'] = scientific_names['stitle'].replace(title, new_title)
        for title in scientific_names['stitle']:
            if ']' in title:
                new_title = title.replace(']', '_')
                scientific_names['stitle'] = scientific_names['stitle'].replace(title, new_title)
        for title in scientific_names['stitle']:
            if '(' in title:
                new_title = title.replace('(', '_')
                scientific_names['stitle'] = scientific_names['stitle'].replace(title, new_title)
        for title in scientific_names['stitle']:
            if ')' in title:
                new_title = title.replace(')', '_')
                scientific_names['stitle'] = scientific_names['stitle'].replace(title, new_title)
        for title in scientific_names['stitle']:
            if ':' in title:
                new_title = title.replace(':', '_')
                scientific_names['stitle'] = scientific_names['stitle'].replace(title, new_title)
        for title in scientific_names['stitle']:
            if ';' in title:
                new_title = title.replace(']', '_')
                scientific_names['stitle'] = scientific_names['stitle'].replace(title, new_title)

        scientific_names['sscinames'] = scientific_names['sscinames'].replace(np.nan, 'unknown')
        for name in scientific_names['sscinames']:
            print(name)
            if ',' in name:
                new_name = name.replace(',', '-')
                scientific_names['sscinames'] = scientific_names['sscinames'].replace(name, new_name)
        for name in scientific_names['sscinames']:
            if ':' in name:
                new_name = name.replace(':', '_')
                scientific_names['sscinames'] = scientific_names['sscinames'].replace(name, new_name)
        for name in scientific_names['sscinames']:
            if ';' in name:
                new_name = name.replace(';', '_')
                scientific_names['sscinames'] = scientific_names['sscinames'].replace(name, new_name)
        for name in scientific_names['sscinames']:
            if '(' in name:
                new_name = name.replace('(', '_')
                scientific_names['sscinames'] = scientific_names['sscinames'].replace(name, new_name)
        for name in scientific_names['sscinames']:
            if ')' in name:
                new_name = name.replace(')', '_')
                scientific_names['sscinames'] = scientific_names['sscinames'].replace(name, new_name)
        for name in scientific_names['sscinames']:
            if '[' in name:
                new_name = name.replace('[', '_')
                scientific_names['sscinames'] = scientific_names['sscinames'].replace(name, new_name)
        for name in scientific_names['sscinames']:
            if ']' in name:
                new_name = name.replace(']', '_')
                scientific_names['sscinames'] = scientific_names['sscinames'].replace(name, new_name)
        # prepare entries
        for name, sacc, stitle, sseqid in zip(scientific_names['sscinames'], scientific_names['sacc'],
                                              scientific_names['stitle'], scientific_names['sseqid']):
            if type(name) == str:
                hit_to_name[sseqid] = '_'.join(name.split(' ')) + '_' + sacc + '_' + stitle
            else:
                hit_to_name[sseqid] = 'NaN'

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
        with open(snakemake.log['log'], 'w') as log_f:
            log_f.write('DONE: The Newick file has been updated with taxonomic information.')
except Exception as e:
    with open(snakemake.log['log'], 'w') as log_f:
        log_f.write('ERROR:'+"[-] Error during updating of sequence ids with taxonomic information: {} with Exception: {}".format(t,e))
    raise Exception("[-] Error during parsing of newick file: {} with Exception: {}".format(newick_treefile,e))
