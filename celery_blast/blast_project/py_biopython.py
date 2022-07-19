''' py_biopython.py

functions that use the biopython package, more informations: https://biopython.org/wiki/Documentation

'''
from Bio import Entrez
from .py_django_db_services import get_project_by_id
''' get_species_taxid_by_name 

utilization in the clean_species_name form field of CreateTaxonomicFileForm

:param user_email
    :type str
:param scientific_name
    :type str
:returns taxonomic node (int) defined in Entrez.esearch dictionary instance
'''
def get_species_taxid_by_name(user_email:str,scientific_name:str)->str:
    try:
        Entrez.email = user_email
        search = Entrez.esearch(term=scientific_name, db="taxonomy", retmode="xml")
        record = Entrez.read(search)
        taxid = record['IdList'][0]
        return taxid
    except Exception as e:
        raise Exception("there is no taxonomic node defined by your specified scientific name: {} : {}".format(scientific_name, e))

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
def get_list_of_species_taxid_by_name(user_email:str,scientific_name:str)->list:
    try:
        Entrez.email = user_email
        search = Entrez.esearch(term=scientific_name, db="taxonomy", retmode="xml")
        record = Entrez.read(search)
        taxid = record['IdList']
        return taxid
    except Exception as e:
        raise Exception("there is no taxonomic node defined by your specified scientific name: {} : {}".format(scientific_name, e))

'''get_list_of_species_taxids_by_list_of_scientific_names
this function iterates over a list of scientific/taxonomic names and converts those names into taxonomic identifier.
Those identifier are stored in a list. Exception can occure if taxonomic names are not specified. 

'''
def get_list_of_species_taxids_by_list_of_scientific_names(user_email:str,scientific_names:list)->list:
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
            raise Exception("There are no taxonomic nodes for the provided scientific names : {}".format(' '.join(scientific_names)))
        return taxonomic_nodes, errors
    except Exception as e:
        raise Exception("[-] ERROR: {}".format(e))

#TODO documentation
def check_given_taxonomic_node(user_email, taxid):
    try:
        Entrez.email = user_email
        handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
        record = Entrez.read(handle)
        handle.close()
        if(len(record) != 0):
            return taxid
        else:
            raise Exception("there is no scientific name defined by your procided taxonomic node")
    except Exception as e:
        raise Exception("exception occured during validation of taxonomic node : {}".format(e))

#TODO documentation
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
def calculate_pfam_and_protein_links_from_queries(user_email,project_id):
    try:
        project = get_project_by_id(project_id)
        path_to_queries = project.get_project_dir()

        #TODO outsource function
        path_to_queries += '/' + project.project_query_sequences
        queries = open(path_to_queries, 'r')
        lines = queries.readlines()
        queries.close()

        queries = []
        for line in lines:
            if ">" in line:
                queries.append(line.split(" ")[0].split(".")[0].split(">")[1])

        Entrez.email = user_email

        #print(path_to_queries,queries)
        handle = Entrez.efetch(db="protein", id=queries, retmode="xml")
        record = Entrez.read(handle)
        handle.close()

        prot_to_pfam = {'QUERIES': [], 'PFAM': {}, 'TIGR': {}, 'REFSEQ': {}, 'CDD': {}, 'Definition':{}, 'Length':{}}
        for i in range(len(record)):
            prot_to_pfam['QUERIES'].append(record[i]['GBSeq_locus'])
            prot_to_pfam['Definition'][record[i]['GBSeq_locus']] = record[i]['GBSeq_definition']
            prot_to_pfam['Length'][record[i]['GBSeq_locus']] = record[i]['GBSeq_length']
            # record[i]['GBSeq_locus']
            pfam = 'http://pfam.xfam.org/family/'
            jvci = 'https://www.ncbi.nlm.nih.gov/genome/annotation_prok/evidence/'
            cdd = 'https://www.ncbi.nlm.nih.gov/Structure/sparcle/archview.html?archid='
            refseq = 'https://www.ncbi.nlm.nih.gov/protein/' + str(record[i]['GBSeq_locus'])
            prot_to_pfam['REFSEQ'][record[i]['GBSeq_locus']] = refseq
            try:
                if "Domain architecture ID" in record[i]['GBSeq_comment']:
                    cdd += record[i]['GBSeq_comment'].split('Domain architecture ID')[1].split(';')[0].strip()
                    #print(cdd)
                    prot_to_pfam['CDD'][record[i]['GBSeq_locus']] = cdd
                if "EMBL-EBI" in record[i]['GBSeq_comment']:
                    pfam += record[i]['GBSeq_comment'].split("EMBL-EBI")[1].split("::")[1].split(";")[0].strip().rstrip()
                    # prot_to_pfam.append((record[i]['GBSeq_locus'],pfam))
                    prot_to_pfam['PFAM'][record[i]['GBSeq_locus']] = pfam
                if "TIGR" in record[i]['GBSeq_comment']:
                    jvci += 'TIGR' + record[i]['GBSeq_comment'].split("TIGR")[1].split(";")[0].split(".")[
                        0].strip().rstrip()
                    # jvci += record[i]['GBSeq_comment'].split("JCVI")[1].split("::")[1].split(";")[0].strip().rstrip()
                    # prot_to_pfam.append((record[i]['GBSeq_locus'],jvci))
                    prot_to_pfam['TIGR'][record[i]['GBSeq_locus']] = jvci
                    #prot_to_pfam.append((record[i]['GBSeq_locus'], 'no entry found'))
                    #prot_to_pfam[record[i]['GBSeq_locus']] = "no entry in pfam database"
                    #print(record[i]['GBSeq_locus'], pfam)
            except Exception as e:
                continue
            #print(record[i]['GBSeq_comment'])
            #print("#############")
            #json.dumps
        return prot_to_pfam
    except Exception as e:
        raise Exception("couldn't parse entrez.efetch with query ids with exception : {}".format(e))