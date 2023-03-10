import matplotlib.pyplot as plt
import matplotlib._color_data as mcd
import pandas as pd
import math
import random
from sys import exit
overlap = [name for name in mcd.CSS4_COLORS]
overlap.remove("lightgrey")

ERRORCODE=10
with open(snakemake.log[0],'w') as logfile:
    try:
        logfile.write("INFO:generating plots for the results of the reciprocal BLAST pipeline\n")
        logfile.write("INFO:loading reciprocal result dataframe into pandas\n")
        rec_prot=pd.read_table(snakemake.input['rec_res'])
        logfile.write("INFO:loading forward BLAST results dataframe into pandas\n")
        fw_res=pd.read_table(snakemake.input['fw_res'],header=None)
        fw_res.columns=["qseqid", "sseqid", "pident", "evalue", "bitscore", "qgi", "sgi", "sacc", "staxids", "sscinames", "scomnames",
                          "stitle", "slen"]

        fw_res['qseqid'] = fw_res['qseqid'].map(lambda line: line.split('.')[0])
        fw_res['sacc'] = fw_res['sacc'].map(lambda line: line.split('.')[0])
        rec_prot = rec_prot.rename(columns={"forward_genome_id": "sacc"})
        rec_prot = rec_prot.rename(columns={"backward_genome_id": "qseqid"})
        logfile.write("INFO:merging forward BLAST dataframe with reciprocal results dataframe\n")
        result_data = rec_prot.merge(fw_res,how='inner', on=['sacc','qseqid','staxids'])
        #the backward blast is currently limited to output only the best match, but the best match can contain several hsps,
        #thus it is possible that there are multiple lines of one qseqid present, which gets loaded by reading the dictionary for
        #filtering reciprocal best hits
        result_data = result_data.drop_duplicates(['sacc','staxids'], keep='first')
        result_data = result_data.reset_index(drop=True)

        logfile.write("INFO:reading query file\n")
        queries = open(snakemake.input['query_file'],'r')
        lines = queries.readlines()
        queries.close()

        queries = []
        for line in lines:
            if ">" in line:
                queries.append(line.split(" ")[0].split(".")[0].split(">")[1])

        logfile.write("INFO:parsing result dataframe\n")
        accid_taxids={}
        for query in queries:
            accid_taxids[query] = result_data[result_data['qseqid'] == query]


        rows = math.floor(math.sqrt(len(queries)))
        columns = math.ceil(math.sqrt(len(queries)))

        if rows * columns < len(queries):
            rows += 1

        logfile.write("INFO:generating e-value figure\n")
        fig, axs = plt.subplots(rows, columns, figsize=(15, 6),
                                facecolor='w', edgecolor='k', constrained_layout=True)
        # fig.subplots_adjust(hspace = .5, wspace=.001)
        axisticks_boolean = 0
        if len(queries) > 1:
            axs = axs.ravel()
            for ax in axs:
                ax.set_axis_off()

            for i in range(len(queries)):
                axs[i].set_axis_on()
                axs[i].grid()
                axs[i].set_facecolor("lightgrey")
                cl = overlap[random.randint(0, len(overlap) - 1)]
                axs[i].scatter(list(accid_taxids[queries[i]]['evalue']), range(len(accid_taxids[queries[i]]['evalue'])),
                               color=cl)
                axs[i].set_title(str(queries[i]))

                if i == axisticks_boolean:
                    axisticks_boolean += math.ceil(math.sqrt(len(queries)))
                    axs[i].set_xticks([0, 0.0005, 0.001])#0, 0.0005, 0.001
                    axs[i].set_xticklabels([0, 0.0005, 0.001])
                else:
                    axs[i].get_xaxis().set_visible(False)
                    axs[i].get_yaxis().set_visible(False)
                axs[i].set_xlim(0, 0.001)
            fig.suptitle("E-Value distribution in target sequences", fontsize=16)
            fig.supylabel("index of reciprocal hit")
            fig.supxlabel("evalue ranging from 0 to 0.001")
            plt.savefig(snakemake.output['evalue_plot'],dpi=300)
            # fig.tight_layout()
        else:
            axs.grid()
            axs.set_facecolor("lightgrey")
            cl = overlap[random.randint(0, len(overlap) - 1)]
            axs.scatter(list(accid_taxids[queries[0]]['evalue']), range(len(accid_taxids[queries[0]]['evalue'])), color=cl)
            axs.set_title(str(queries[0]))
            plt.savefig(snakemake.output['evalue_plot'],dpi=300)

        logfile.write("INFO:generating taxids to hit figure\n")
        fig, axs = plt.subplots(rows, columns, figsize=(15, 6),
                                facecolor='w', edgecolor='k', constrained_layout=True)

        if len(queries) > 1:

            axs = axs.ravel()
            for ax in axs:
                ax.set_axis_off()

            for i in range(len(accid_taxids.keys())):
                axs[i].set_axis_on()
                if (len(accid_taxids[queries[i]]) != 0):

                    hit_distribution = {}
                    for hit in accid_taxids[queries[i]]['staxids'].unique():
                        val = accid_taxids[queries[i]][accid_taxids[queries[i]]['staxids'] == hit]['qseqid'].count()
                        if val <= 20:
                            if val not in hit_distribution.keys():
                                hit_distribution[val] = 1
                            else:
                                hit_distribution[val] += 1

                    xvalues = []
                    yvalues = []
                    for key in sorted(hit_distribution.keys()):
                        xvalues.append(key)
                        yvalues.append(hit_distribution[key])

                    axs[i].grid()
                    axs[i].bar(x=xvalues, height=yvalues, color=overlap[random.randint(0, len(overlap) - 1)], width=0.8,
                               edgecolor="black")
                    axs[i].set_title(str(queries[i]))


                    if len(list(hit_distribution.keys())) == 1:
                        axs[i].set_xticks(range(1, max(list(hit_distribution.keys())) + 1, 1))
                    if len(list(hit_distribution.keys())) >= 10:
                        axs[i].tick_params('x', labelrotation=45)
                else:
                    axs[i].grid()
                    axs[i].set_title(str(queries[i]))

            fig.suptitle("Amount of hits in database organisms (limit 20)", fontsize=16)
            fig.supylabel("number of distinct taxids")
            fig.supxlabel("amount of hits in backward genome (limit 20)")
            plt.savefig(snakemake.output['taxids_hits_plot'],dpi=300)
        else:
            hit_distribution = {}
            for hit in accid_taxids[queries[0]]['staxids'].unique():
                val = accid_taxids[queries[0]][accid_taxids[queries[0]]['staxids'] == hit]['qseqid'].count()
                if val not in hit_distribution.keys():
                    hit_distribution[val] = 1
                else:
                    hit_distribution[val] += 1
            xvalues = []
            yvalues = []
            for key in sorted(hit_distribution.keys()):
                xvalues.append(key)
                yvalues.append(hit_distribution[key])

            axs.grid()
            axs.set_facecolor("lightgrey")
            cl = overlap[random.randint(0, len(overlap) - 1)]
            axs.bar(x=xvalues, height=yvalues, color=cl)
            axs.set_title(str(queries[0]))
            axs.set_xticks(range(1, max(list(hit_distribution.keys())) + 1, 1))
            plt.savefig(snakemake.output['taxids_hits_plot'],dpi=300)
        logfile.write("DONE\n")
    except Exception as e:
        logfile.write("ERROR:{}\n".format(e))
        exit(ERRORCODE)