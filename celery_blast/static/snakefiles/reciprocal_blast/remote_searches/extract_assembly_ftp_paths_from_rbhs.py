import pandas as pd
from Bio import Entrez
from sys import exit
from typing import TextIO


def set_protein_assembly_file(ftp_path:str)->str:
    try:
        if type(ftp_path) == str:
            protein_genome = ftp_path.split('/')[-1:][0]
            protein_genome = ftp_path + '/' + str(protein_genome) + '_protein.faa.gz'
            return protein_genome
        else:
            return ftp_path
    except:
        raise Exception("[-] Problem during parsing the ftp_path column in the refseq assembly summary file")

'''read_assembly

    This function opens the assembly summary file (provided by NCBI and downloaded during CATHI setup or within the NCBI-DB
    transaction dashboard. Datatypes of the taxid, species taxid and excluded from refseq columns are transformed to str.
    The first two rows are skipped. Column headers are optained from the second comment line.

    :param summary_file_path - refseq or genbank 
        :type str
    :param assemblies
        :type list[str]

    :returns assembly_table
        :type pd.DataFrame

'''
def read_assembly(summary_file_path:str, assemblies:list)->pd.DataFrame:
    try:
        # skipping the first line with .readline()
        # --> second line resides the header information for the assembly summary file
        with open(summary_file_path, 'r') as rfile:
            line = rfile.readline()
            line = rfile.readline()
            header = line.replace('#', '').replace(" ", '').rstrip().split("\t")

        assembly_table = pd.read_table(summary_file_path, skiprows=[0, 1], header=None,
                                     dtype={20: str,  # 20 excluded from refseq
                                            5: str,  # 5 taxid
                                            6: str,  # 6 species taxid
                                            'ftp_path': str})
        assembly_table.columns = header
        assembly_table = assembly_table.astype({"taxid": str})

        # extract necessary data fields: assembly number, names, taxids and the correct ftp_filepath for downloading with gzip
        assembly_table = assembly_table[
            ['assembly_accession', 'organism_name', 'taxid', 'species_taxid', 'assembly_level', 'ftp_path']]
        # python lambda function applied to each row in the dataframe
        assembly_table['ftp_path'] = assembly_table['ftp_path'].apply(lambda row: set_protein_assembly_file(row))

        assembly_table = assembly_table[assembly_table['assembly_accession'].isin(assemblies)]
        return assembly_table
    except Exception as e:
        raise ValueError(
            "exception during pandas parsing of {} file ...\n\tException: {}".format(summary_file_path, e))


'''read_assembly_summary_files

    This function extracts entries based on assembly identifiers. Those identifiers are obtained from BLAST searches 
    against online databases with a subsequent entrez search against the protein database to fetch identical protein records, in which
    the assembly identifier are located. Those identifers are then used to parse the assembly_summary_files of the relevant databases.
    Databases are RefSeq (GCF) and GenBank (GCA). The function returns a pandas DataFrame that can be used in the synteny calculation 
    dashboard for remote BLAST projects.

    :param refseq_assemblies - assemblies starting with GCF
        :type list[str]
    :param genbank_assemblies- assemblies starting with GCA
        :type list[str]
    :param assembly_file_path
        :type str

    :returns combined_assemblies_dataframe
        :type pd.DataFrame

'''
def read_assembly_summary_files(refseq_assemblies: list, genbank_assemblies: list, assembly_file_path:str) -> pd.DataFrame:
    try:
        refseq_table = 0
        genbank_table = 0

        if len(refseq_assemblies) != 0:
            refseq_file_path = assembly_file_path + 'assembly_summary_refseq.txt'
            refseq_assembly = read_assembly(refseq_file_path, refseq_assemblies)
            refseq_table = 1

        if len(genbank_assemblies) != 0:
            genbank_file_path = assembly_file_path + 'assembly_summary_genbank.txt'
            genbank_assembly = read_assembly(genbank_file_path, genbank_assemblies)
            genbank_table = 1

        if refseq_table == 0 and genbank_table == 0:
            return pd.DataFrame()

        elif refseq_table == 1 and genbank_table == 1:
            return pd.concat([refseq_assembly, genbank_assembly])

        elif refseq_table == 1 and genbank_table == 0:
            return refseq_assembly
        elif refseq_table == 0 and genbank_table == 1:
            return genbank_assembly

    except Exception as e:
        raise Exception(
            "[-] ERROR during parsing of assembly summary files to obtain a dataframe with correct FTP paths with exception: {}".format(
                e))


'''read_result_dataframe

    This function reads the RBH result file and 
    returns a list of unique protein entries (RBHs).

    :param result_table_path - path to the rbh_table.tsf file
        :type str
    :retunrs sacc_list
        :type list[str]
'''
def read_result_dataframe(result_table_path: str) -> list:
    try:
        result_table = pd.read_csv(result_table_path, sep="\t")
        sacc_list = list(result_table.sacc.unique())
        return sacc_list
    except Exception as e:
        raise Exception("[-] ERROR reading reciprocal result dataframe with exception: {}".format(e))


'''write_ipg_table

    This function gets a list of protein identifiers and searches in the protein database for identical protein entries 
    to obtain a list of assembly identifers. This function implements the BioPython Entrez.efetch function.
    It writes a tab delimited table based on the output of the Entrez search.

    :param proteins - list of protein identifier
        :type list[str]
    :param ipg_table_path - path to project directory
        :type str
    :param logfile
        :type TextIO
    
    :returns returncode
        :type int
'''
def write_ipg_table(proteins: list, ipg_table_path: str, logfile:TextIO) -> int:
    try:
        with open(ipg_table_path, "w") as ipgtable:

            logfile.write("INFO:starting to fetch IPG records from protein database\n")
            end = len(proteins)
            logfile.write("INFO:fetching IPGs from {} proteins\n".format(end))
            logfile.write("INFO:fetching IPGs in steps of 500 protein identifier\n")
            begin = 0
            step = 500
            steps = 500
            while begin < end:
                if step >= end:
                    step = end
                splitted_ids = proteins[begin:step]
                for attempt in range(10):
                    try:
                        search = Entrez.efetch(id=splitted_ids, db="protein", rettype="ipg", retmode="text")
                        lines = search.readlines()
                        search.close()
                    except Exception as e:
                        if attempt == 9:
                            logfile.write("WARNING:ten attempts for step: {}\n".format(step))
                            raise Exception("INFO:ERROR during entrez.efetch with exception: {}".format(e))

                    else:
                        for line in lines:
                            ipgtable.write(line.decode())
                    break

                begin += steps
                step += steps
            logfile.write("INFO:DONE writing IPG table\n")
        return 0
    except Exception as e:
        logfile.write("ERROR:ERROR fetching and writing ipg table with exception: {}".format(e))
        raise Exception(
            "[-] ERROR fetching and writing ipg table with exception: {}".format(ipg_table_path, e))


'''return_ipg_table

    This function reads the identical protein groups table and returns a pandas dataframe.

    :param ipg_table_filepath
        :type str

    :returns ipg_table
        :type pd.DataFrame
'''
def return_ipg_table(ipg_table_filepath: str) -> pd.DataFrame:
    try:
        ipg_table = pd.read_csv(ipg_table_filepath, sep="\t")
        ipg_table.columns = ["entrez_id", "source", "nucleotide_accession", "start", "stop", "strand", "protein",
                             "protein_name", "organism", "strain", "assembly_accession"]
        return ipg_table
    except Exception as e:
        raise Exception("[-] ERROR reading ipg table with exception: {}".format(e))


'''extract_assembly_ftp_paths_from_remote_blast_search_rbhs

    This function is a wrapper for the functions that are used to create a merged datframe of the refseq/genbank assembly
    summary files and a dataframe of identical protein groups obtained from the reciprocal results (RBHs) of the remote BLAST
    project. An Entrez.efetch search is performed (write_ipg_table function).
    It outputs a dataframe that can be used for downloading GenBank files for synteny calculation.

    :param proteins
        :type list[str]       
    :param ipg_table_path - path to the identical protein groups table
        :type str
    :param path_to_synteny_dataframe - path to the result dataframe that can be used in the synteny dashboard
        :type str
    :param assembly_file_path
        :type str
    :param logfile
        :type TextIO
        
    :returns combined_result_dataframe - combination of the assembly summary tables and the ipg table
        :type pd.DataFrame
'''


def extract_assembly_ftp_paths_from_remote_blast_search_rbhs(proteins: list,
                                                             ipg_table_path: str,
                                                             path_to_synteny_dataframe: str,
                                                             assembly_file_path:str,
                                                             logfile:TextIO) -> pd.DataFrame:
    try:


        # loop over each
        if write_ipg_table(proteins, ipg_table_path, logfile) != 0:
            raise Exception("entrez search for identical protein groups resulted in an error")

        logfile.write("INFO:reading ipg_table.table file\n")
        ipg_table = return_ipg_table(ipg_table_path)

        logfile.write("INFO:filtering ipg_table ...\n")
        # filter for refseq and genbank entries
        refseq_assemblies_table = ipg_table[ipg_table['source'] == "RefSeq"]
        genbank_assemblies_table = ipg_table[ipg_table['source'] == "INSDC"]

        logfile.write("INFO:reading assembly summary files ...\n")
        # extract assembly summary information
        concat_assemblies_table = read_assembly_summary_files(list(refseq_assemblies_table.assembly_accession.unique()),
                                                              list(
                                                                  genbank_assemblies_table.assembly_accession.unique()),
                                                              assembly_file_path)
        logfile.write("INFO:merging ipg_table and assembly_summary_table on assembly_accession column\n")
        # merge dataframes
        ipg_table = ipg_table[['protein', 'assembly_accession', 'source']]
        ipg_table['protein'] = ipg_table.protein.apply(lambda x: x.split(".")[0])
        ipg_table = ipg_table[ipg_table.protein.isin(proteins)]

        concat_assemblies_table = concat_assemblies_table.merge(ipg_table, left_on="assembly_accession",
                                                                right_on="assembly_accession")

        logfile.write("INFO:filering assembly_table ...\n")
        concat_assemblies_table = concat_assemblies_table.dropna()
        concat_assemblies_table = concat_assemblies_table.drop_duplicates(subset=["assembly_accession", "protein"],
                                                                          keep="first")
        concat_assemblies_table = concat_assemblies_table[concat_assemblies_table.ftp_path != "na/na_protein.faa.gz"]

        # sort dataframe by taxids and reset the index
        columns = list(concat_assemblies_table.columns)
        concat_assemblies_table = concat_assemblies_table.sort_values(by="taxid").reset_index()
        concat_assemblies_table = concat_assemblies_table[columns]
        logfile.write("INFO:writing assembly file to disc ...\n")
        # save dataframe to disc
        concat_assemblies_table.to_csv(path_to_synteny_dataframe)

    except Exception as e:
        raise Exception(
            "[-] ERROR during extraction of FTP-paths for remote BLAST project RBHs with exception: {}".format(e))


try:
    Entrez.email = snakemake.params['user_email']
    ERRORCODE = 21

    with open(snakemake.log['log'],"w") as logfile:
        try:
            logfile.write("INFO:start to fetch identical protein groups (ipgs) from RBHs\n")
            logfile.write("INFO:setting up variables ...\n")
            assembly_file_path = snakemake.params['assembly_file_path']
            ipg_table_path = snakemake.output['ipg_table_path']
            result_dataframe_path = snakemake.input['result_dataframe']
            path_to_assembly_dataframe = snakemake.output['genome_assembly_table']
            logfile.write("INFO:starting procedure with following variables:\n")
            logfile.write("\t{}\n\t{}\n\t{}\n\t{}\n".format(
                "assembly_file_path:{}".format(assembly_file_path),
                "ipg_table_path:{}".format(ipg_table_path),
                "result_dataframe_path:{}".format(result_dataframe_path),
                "path_to_assembly_datafrane:{}".format(path_to_assembly_dataframe)
            ))

            logfile.write("INFO:reading RBH table file\n")
            proteins = read_result_dataframe(result_dataframe_path)

            extract_assembly_ftp_paths_from_remote_blast_search_rbhs(proteins, ipg_table_path, path_to_assembly_dataframe, assembly_file_path, logfile)
            logfile.write("DONE\n")
        except Exception as e:
            #DtypeWarning: Columns (34,35,36) have mixed types.Specify dtype option on import or set low_memory=False.
            logfile.write("ERROR: ERROR with exception: {}\n".format(e))
            raise Exception(e)
except Exception as e:

    exit(ERRORCODE)