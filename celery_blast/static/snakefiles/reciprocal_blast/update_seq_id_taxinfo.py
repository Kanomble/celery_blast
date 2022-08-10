import pandas as pd
import re
from sys import exit
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
ERRORCODE=14
newick_treefile= snakemake.input['nw']
result_df= snakemake.input['df']

with open(snakemake.log['log'],'w') as logfile:
    try:
        logfile.write("INFO:starting to change treefile headers ...\n")
        #reading newick treefile
        logfile.write("INFO:loading treefile\n")
        with open(newick_treefile,'r') as treefile:
            treefile_content = treefile.readlines()
        #searching for scientific names

        if len(treefile_content) > 0:
            logfile.write("INFO:parsing treefile headers\n")
            hits = re.findall('[(,]{1}([a-zA-Z_0-9]+[a-zA-Z_0-9._]+):{1}',treefile_content[0])
            hits = [hit.split(':')[0] for hit in hits]

            result_df = pd.read_csv(result_df, index_col=0)
            hit_to_name = {}
            scientific_names = result_df[result_df['sacc'].isin(hits)][['sscinames', 'sacc']]
            logfile.write("INFO:filling dictionary with new headers\n")
            for name, sacc in zip(scientific_names['sscinames'], scientific_names['sacc']):
                if type(name) == str:
                    hit_to_name[sacc] = '_'.join(name.split(' ')) + '_' + sacc
                else:
                    hit_to_name[sacc] = 'NaN'

            # check for duplicate scientific names and add counter to distinguish between them,
            # scientific_names includes the sscienames and sacc for more precise tree labeling

            # replacing sequence identifier with scientific names
            logfile.write("INFO:replacing headers\n")
            for hit in hit_to_name.keys():
                treefile_content[0] = re.sub(hit, hit_to_name[hit], treefile_content[0])
            #writing new newick file with changes header
            logfile.write("INFO:writing new treefile\n")
            new_treefile = newick_treefile.split("/")
            new_treefile = new_treefile[0:len(new_treefile)-1]
            new_treefile = '/'.join(new_treefile)
            new_treefile = new_treefile + '/transformed_treefile.nwk'
            with open(new_treefile,'w') as treefile:
                treefile.write(treefile_content[0])
            #return treefile_content[0]
        else:
            logfile.write("WARNING:no content in treefile\n")
            logfile.write("INFO:writing new empty treefile\n")
            new_treefile = newick_treefile.split("/")
            new_treefile = new_treefile[0:len(new_treefile)-1]
            new_treefile = '/'.join(new_treefile)
            new_treefile = new_treefile + '/transformed_treefile.nwk'
            with open(new_treefile,'w') as treefile:
                pass
        logfile.write("DONE\n")
    except Exception as e:
        logfile.write("ERROR:exception in treefile header transformation - {}\n".format(e))
        exit(ERRORCODE)