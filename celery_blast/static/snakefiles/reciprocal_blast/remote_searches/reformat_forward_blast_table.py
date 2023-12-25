import pandas as pd

with open(snakemake.log["log"], "w") as logfile:
    try:

        logfile.write("INFO:starting to reformat forward BLAST dataframe.\n")
        logfile.write("INFO:loading forward BLAST dataframe with pandas.\n")

        fw_res = pd.read_table(snakemake.input["fw_res"], header=None)
        fw_res.columns = ["qseqid", "sseqid", "pident", "evalue", "bitscore", "qgi", "sgi", "sacc", "staxids",
                          "sscinames", "scomnames",
                          "stitle","slen"]
        dataframe_dict = {
            "qseqid":[], "sseqid":[], "pident":[], "evalue":[], "bitscore":[], "qgi":[], "sgi":[], "sacc":[], "staxids":[],
            "sscinames":[], "scomnames":[],"stitle":[],"slen":[]
        }
        counter = 0
        for taxids, sscinames, scomnames in zip(fw_res.staxids, fw_res.sscinames, fw_res.scomnames):
            try:
                taxid_list = taxids.split(";")
                scinames_list = sscinames.split(";")
                scomnames_list = scomnames.split(";")
                if len(taxid_list) == len(scinames_list) == len(scomnames_list):
                    temp_df = fw_res.iloc[counter]
                    for taxid, sciname, comname in zip(taxid_list, scinames_list, scomnames_list):
                        for index, val in zip(temp_df.index,temp_df.values):
                            if index == "staxids":
                                dataframe_dict['staxids'].append(int(taxid))
                            elif index == "sscinames":
                                dataframe_dict["sscinames"].append(sciname)
                            elif index == "scomnames":
                                dataframe_dict["scomnames"].append(comname)
                            else:
                                dataframe_dict[index].append(val)
                counter += 1
            except Exception as e:
                logfile.write("WARNING: error in row: {}\n".format(
                    counter
                ))
        logfile.write("INFO:done parsing original forward BLAST dataframe.\n")
        # rename original forward BLAST dataframe
        new_fw_dataframe = pd.DataFrame.from_dict(dataframe_dict)
        logfile.write("INFO:writing new forward BLAST dataframe to disc.\n")
        new_fw_dataframe.to_csv(snakemake.output["new_fw_res"], sep="\t", index=0, header=False)
        logfile.write("DONE\n")
    except Exception as e:
        logfile.write("ERROR: error during reformatting forward BLAST dataframe with exception: {}\n".format(e))
        raise Exception("[-] ERROR during reformatting forward BLAST dataframe with exception: {}".format(e))

