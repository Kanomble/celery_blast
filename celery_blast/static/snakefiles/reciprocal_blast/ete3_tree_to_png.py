import ete3
from os import environ
from sys import exit
import shutil
environ['QT_QPA_PLATFORM'] = 'offscreen'
environ['XDG_RUNTIME_DIR'] = '../../../tmp'

ERRORCODE=15
with open(snakemake.log['log'],'w') as logfile:
    try:
        logfile.write("INFO:starting to convert treefile to png ... \n")
        logfile.write("INFO:loading treefile\n")
        filename= snakemake.input['tree']
        query= filename.split('/')[0]
        with open(snakemake.input['tree'], 'r') as t:
            tree = t.readlines()
        if len(tree) < 1:
            with open(snakemake.log['log'], 'w') as log_f:
                log_f.write('ERROR: The Newick file used in this command has no content.')
            shutil.copyfile('../../../static/images/no_results.svg', snakemake.params['static_pic'])
            shutil.copyfile('../../../static/images/no_results.svg', snakemake.output['pic'])
        else:
            tree = tree[0].replace("'","")
            logfile.write("INFO:converting tree to ete3 object\n")
            tree = ete3.Tree(tree)
            ts = ete3.TreeStyle()
            ts.show_branch_length = True
            ts.title.add_face(ete3.TextFace('Phylogenetic tree of '+ str(query),fsize=20),column=0)
            logfile.write("INFO:generating png image\n")
            tree.render(snakemake.params['static_pic'], tree_style=ts)
            tree.render(snakemake.output['pic'], tree_style=ts)
        logfile.write("DONE\n")
    except Exception as e:
        try:
            with open(snakemake.input['nw'], 'r') as t:
                tree = t.readlines()
            if len(tree) < 1:
                with open(snakemake.log['log'], 'w') as log_f:
                    log_f.write('ERROR: The Newick file used in this command has no content.')
                shutil.copyfile('../../../static/images/no_results.svg', snakemake.params['static_pic'])
                shutil.copyfile('../../../static/images/no_results.svg', snakemake.output['pic'])
            else:
                tree = tree[0].replace("'", "")
                logfile.write("INFO:using normal treefile\n")
                tree = ete3.Tree(tree)
                ts = ete3.TreeStyle()
                ts.show_branch_length = True
                ts.title.add_face(ete3.TextFace('Phylogenetic tree of ' + str(query), fsize=20), column=0)
                logfile.write("INFO:generating png image\n")
                tree.render(snakemake.params['static_pic'], tree_style=ts)
                tree.render(snakemake.output['pic'], tree_style=ts)
        except Exception as e:
            logfile.write("ERROR:treefile couldn't get converted to png image - {}\n".format(e))
            exit(ERRORCODE)