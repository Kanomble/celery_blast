'''ete3_tree_to_png.py
input:newick tree file
output: png showing phylogenetic tree
builds phylogenetic tree based on newick file and saves it as a png
'''
#TODO add logging; try to modify labeling of the tree branches


import pandas as pd
import ete3
import shutil
from os import environ
import re


environ['QT_QPA_PLATFORM'] = 'offscreen'
environ['XDG_RUNTIME_DIR'] = '../../../tmp'
#ete3_tree_to_png.py
try:
    filename= snakemake.input['tree']
    query= filename.split('/')[0]
    with open(snakemake.input['tree'], 'r') as t:
        tree = t.readlines()
    if len(tree) < 1:
        with open(snakemake.log['log'], 'w') as log_f:
            log_f.write('ERROR: The Newick file {} used in this command has no content.'.format(filename))
        shutil.copyfile('../../../static/images/no_results.svg', snakemake.params['static_pic'])
        shutil.copyfile('../../../static/images/no_results.svg', snakemake.output['pic'])
    else:
        tree = tree[0]
        tree = ete3.Tree(tree)
        ts = ete3.TreeStyle()
        ts.show_branch_length = True
        ts.title.add_face(ete3.TextFace('Phylogenetic tree of '+ str(query),fsize=20),column=0)
        tree.render(snakemake.params['static_pic'], tree_style=ts)
        tree.render(snakemake.output['pic'], tree_style=ts)
        with open(snakemake.log['log'], 'w') as log_f:
            log_f.write('DONE: The Newick tree file {} has been rendered to svg {}.'.format(filename,snakemake.output['pic']))
except Exception as e:
    with open(snakemake.log['log'], 'w') as log_f:
        log_f.write('ERROR:'+"[-] Error during parsing of svg: {} with Exception: {}".format(snakemake.output['pic'],e))
    raise Exception("[-] Error during parsing of svg: {} with Exception: {}.".format(snakemake.output['pic'],e))

