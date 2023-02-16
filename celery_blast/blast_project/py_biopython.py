''' py_biopython.py

    Functions that use the biopython package, more informations: https://biopython.org/wiki/Documentation.
    The user_email parameter in the following functions is used for errorhandling by ncbi.

'''

from Bio import Entrez
from .py_django_db_services import get_project_by_id

''' get_species_taxid_by_name 

    Utilization in the clean_species_name form field of CreateTaxonomicFileForm.

    :param user_email
        :type str
        
    :param scientific_name
        :type str
    
    :returns taxonomic nodes defined in Entrez.esearch dictionary instance
        :type list  
'''


def get_species_taxid_by_name(user_email: str, scientific_name: str) -> list:
    try:
        Entrez.email = user_email
        search = Entrez.esearch(term=scientific_name, db="taxonomy", retmode="xml")
        record = Entrez.read(search)
        search.close()
        taxids = record['IdList']
        return taxids
    except Exception as e:
        raise Exception(
            "there is no taxonomic node defined by your specified scientific name: {} : {}".format(scientific_name, e))


'''check_if_protein_identifier_correspond_to_backward_taxid
    
    Used in the ProjectCreationForm form validation. Checks if the provided query sequences match to the specified backward organism.
    
    :param protein_identifier
        :type list
    :param taxonomic_identifier
        :type str
    :param user_email
        :type str
'''


def check_if_protein_identifier_correspond_to_backward_taxid(protein_identifier: list, taxonomic_identifier: str,
                                                             user_email: str) -> int:
    try:
        Entrez.email = user_email
        search = Entrez.elink(dbfrom='protein', id=protein_identifier, linkname="protein_taxonomy")
        record = Entrez.read(search)
        search.close()
        taxonomic_ids = []
        for rec in record:
            taxid = rec['LinkSetDb'][0]['Link'][0]['Id']
            if int(taxid) != int(taxonomic_identifier):
                if taxid not in taxonomic_ids:
                    taxonomic_ids.append(taxid)

        if len(taxonomic_ids) != 0:
            search = Entrez.efetch(id=taxonomic_ids, db='taxonomy', retmode='xml')
            record = Entrez.read(search)
            search.close()

            for rec in record:
                taxid_to_check = rec['TaxId']
                for lineage in rec['LineageEx']:
                    if lineage['Rank'] == 'genus':
                        species_level_tax_id = lineage['TaxId']
                        if int(species_level_tax_id) == int(taxonomic_identifier):
                            taxonomic_ids.remove(taxid_to_check)

            if len(taxonomic_ids) != 0:
                return 1

        return 0
    except Exception as e:
        raise Exception(
            "[-] Problem during validation of protein identifiers and taxid of backward organisms with exception {}".format(
                e))


'''get_list_of_species_taxid_by_name
    sometimes there are mutliple taxonomic nodes for one organism name (e.g. get_species_taxids.sh -n bacillus = 1386, 55087)
    therefore this function can be used to iterate over all available nodes. Those nodes will be written into one file that is 
    than processed for database parsing.
    
    :param user_email
        :type str
    :param scientific_name
        :type str
    :returns taxonomic nodes (list:int)
'''


def get_list_of_species_taxid_by_name(user_email: str, scientific_name: str) -> list:
    try:
        Entrez.email = user_email
        search = Entrez.esearch(term=scientific_name, db="taxonomy", retmode="xml")
        record = Entrez.read(search)
        taxid = record['IdList']
        return taxid
    except Exception as e:
        raise Exception(
            "there is no taxonomic node defined by your specified scientific name: {} : {}".format(scientific_name, e))


'''get_list_of_species_taxids_by_list_of_scientific_names

    this function iterates over a list of scientific/taxonomic names and converts those names into taxonomic identifier.
    Those identifier are stored in a list. Exception can occure if taxonomic names are not specified. Those strings
    will be appended to the errors list.
    
    :param user_email
        :type str
    :param scientific_names
        :type list[str]
    
    :returns taxonomic_nodes, errors
        :type list[int], list[str]

'''


def get_list_of_species_taxids_by_list_of_scientific_names(user_email: str, scientific_names: list) -> list:
    try:
        Entrez.email = user_email
        taxonomic_nodes = []
        errors = []
        for name in scientific_names:
            try:
                search = Entrez.esearch(term=name, db="taxonomy", retmode="xml")
                record = Entrez.read(search)
                taxids = record['IdList']
                if len(taxids) == 0:
                    errors.append(name)
                taxonomic_nodes.extend(taxids)
            except:
                errors.append(name)
                continue
        if len(taxonomic_nodes) == 0:
            raise Exception("There are no taxonomic nodes for the provided scientific names : {}".format(
                ' '.join(scientific_names)))
        return taxonomic_nodes, errors
    except Exception as e:
        raise Exception("[-] ERROR: {}".format(e))


'''check_given_taxonomic_node

    This functions checks if a provided integer corresponds to a real taxonomic node.
    If this is true, it returns the provided integer.
    
    :param user_email
        :type str
    :param taxid
        :type int
    
    :return taxid
        :type int
'''


def check_given_taxonomic_node(user_email: str, taxid: int) -> int:
    try:
        Entrez.email = user_email
        handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
        record = Entrez.read(handle)
        handle.close()
        if (len(record) != 0):
            return taxid
        else:
            raise Exception("there is no scientific name defined by your procided taxonomic node")
    except Exception as e:
        raise Exception("exception occured during validation of taxonomic node : {}".format(e))


''' calculate_pfam_and_protein_links_from_queries

    This function returns a dictionary 'prot_to_pfam', that contains several keys addressing relevant data for query sequences that reside in the coressponding file.
    The QUERIES key, which points to a list containing all query accession ids, the PFAM, TIGR, REFSEQ, CDD keys, which point to additional key (query_id) value (link to corresponding database) pairs.
    The Definition and Length keys contain key (query_id) value (sequence description, sequence length) as well.
    
    prot_to_pfam datastructure:
        
        1. QUERIES: accession ids defined in the corresponding fasta file
        2. PFAM: The Pfam database is a large collection of protein families, each represented by multiple sequence alignments and hidden Markov models (HMMs)
        3. TIGR: Protein Family Models are a hierarchical collection of curated Hidden Markov Model-based and BLAST-based protein families (HMMs and BlastRules),
                and Conserved Domain Database architectures used to assign names, gene symbols,
                publications and EC numbers to the prokaryotic RefSeq proteins that meet the criteria for inclusion in a family.
                HMMs and BlastRules also contribute to structural annotation by NCBI's Prokaryotic Genome Annotation Pipeline (PGAP).
                
        4. CDD: Protein Family Model search by SPARCLE.
'''


def calculate_pfam_and_protein_links_from_queries(user_email, project_id):
    try:
        # load BlastProject from database
        project = get_project_by_id(project_id)
        # setup filepath to BlastProject
        path_to_queries = project.get_project_dir()
        path_to_queries += '/' + project.project_query_sequences
        # open query fasta file and read protein identifier
        with open(path_to_queries, 'r') as queries:
            lines = queries.readlines()

        queries = []
        for line in lines:
            if ">" in line:
                queries.append(line.split(" ")[0].split(".")[0].split(">")[1])

        # retrieving information from NCBIs protein database - biopython
        Entrez.email = user_email
        handle = Entrez.efetch(db="protein", id=queries, retmode="xml")
        record = Entrez.read(handle)
        handle.close()

        # setup dictionary that is returned to ajax_wp_to_links
        prot_to_pfam = {'QUERIES': [], 'PFAM': {}, 'TIGR': {}, 'REFSEQ': {}, 'CDD': {}, 'Definition': {}, 'Length': {}}
        for i in range(len(record)):
            prot_to_pfam['QUERIES'].append(record[i]['GBSeq_locus'])
            prot_to_pfam['Definition'][record[i]['GBSeq_locus']] = record[i]['GBSeq_definition']
            prot_to_pfam['Length'][record[i]['GBSeq_locus']] = record[i]['GBSeq_length']

            # setting up standard urls
            pfam = 'http://pfam.xfam.org/family/'
            jvci = 'https://www.ncbi.nlm.nih.gov/genome/annotation_prok/evidence/'
            cdd = 'https://www.ncbi.nlm.nih.gov/Structure/sparcle/archview.html?archid='
            refseq = 'https://www.ncbi.nlm.nih.gov/protein/' + str(record[i]['GBSeq_locus'])
            prot_to_pfam['REFSEQ'][record[i]['GBSeq_locus']] = refseq

            try:
                if "Domain architecture ID" in record[i]['GBSeq_comment']:
                    cdd += record[i]['GBSeq_comment'].split('Domain architecture ID')[1].split(';')[0].strip()
                    prot_to_pfam['CDD'][record[i]['GBSeq_locus']] = cdd

                if "EMBL-EBI" in record[i]['GBSeq_comment']:
                    pfam += record[i]['GBSeq_comment'].split("EMBL-EBI")[1].split("::")[1].split(";")[
                        0].strip().rstrip()
                    prot_to_pfam['PFAM'][record[i]['GBSeq_locus']] = pfam

                if "TIGR" in record[i]['GBSeq_comment']:
                    jvci += 'TIGR' + record[i]['GBSeq_comment'].split("TIGR")[1].split(";")[0].split(".")[
                        0].strip().rstrip()
                    prot_to_pfam['TIGR'][record[i]['GBSeq_locus']] = jvci

            except Exception as e:
                continue

        return prot_to_pfam
    except Exception as e:
        raise Exception("couldn't parse entrez.efetch with query ids with exception : {}".format(e))


''' fetch_protein_records

    Function wrapper for Entrez.efetch function. 
    Takes as input a list with maximal 500 protein accession ids and an email
    address for the current user.
    
    Returns the Entrez.Parser.ListElement (from Bio import Entrez) that can be used by the
    parse_entrez_xml function. Failures contains a list with unfetchable protein identifiers.
    
    :param proteins
        :type list
    :param email
        :type str
    
    :returns records, failures
        :type tuple(Entrez.Parser.ListElement,list)
'''


def fetch_protein_records(proteins: list, email: str):
    try:
        Entrez.email = email
        handle = Entrez.efetch(db="protein", id=proteins, retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        failures = []
        fetched_protein_ids_with_version = []
        fetched_protein_ids_without_version = []
        for rec in records:
            fetched_protein_ids_with_version.append(rec['GBSeq_accession-version'])
            fetched_protein_ids_without_version.append(rec['GBSeq_locus'])
        for protein in proteins:
            if protein not in fetched_protein_ids_with_version and protein not in fetched_protein_ids_without_version:
                failures.append(protein)

        return records, failures
    except Exception as e:
        raise Exception("[-] ERROR fetching protein xml data from Biopython with exception: {}".format(e))
