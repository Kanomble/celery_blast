import ete3
from os import environ

environ['QT_QPA_PLATFORM'] = 'offscreen'
environ['XDG_RUNTIME_DIR'] = '../../../tmp'
#ete3_tree_to_png.py
try:
    filename= snakemake.input['tree']
    query= filename.split('/')[0]
    with open(snakemake.input['tree'], 'r') as t:
        tree = t.readlines()
    tree = tree[0]
    tree = ete3.Tree(tree)
    ts = ete3.TreeStyle()
    ts.show_branch_length = True
    ts.title.add_face(ete3.TextFace('Phylogenetic tree of '+ str(query),fsize=20),column=0)
    tree.render(snakemake.params['static_pic'], tree_style=ts)
    tree.render(snakemake.output['pic'], tree_style=ts)
except Exception as e:
    raise Exception("[-] Error during parsing of png: {} with Exception: {}".format(t,e))