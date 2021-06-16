''' py_biopython.py

functions that use the biopython package, more informations: https://biopython.org/wiki/Documentation

'''
from Bio import Entrez
import json
from .py_django_db_services import get_project_by_id
''' get_species_taxid_by_name

utilization in the clean_species_name form field of CreateTaxonomicFileForm

:param user_email
    :type str
:param scientific_name
    :type str
:returns taxonomic node (int) defined in Entrez.esearch dictionary instance
'''
def get_species_taxid_by_name(user_email,scientific_name):
    try:
        Entrez.email = user_email
        search = Entrez.esearch(term=scientific_name, db="taxonomy", retmode="xml")
        record = Entrez.read(search)
        taxid = record['IdList'][0]
        return taxid
    except Exception as e:
        raise Exception("there is no taxonomic node defined by your specified scientific name: {} : {}".format(scientific_name, e))

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
def calculate_pfam_and_protein_links_from_queries(user_email,project_id):
    try:
        project = get_project_by_id(project_id)
        path_to_queries = project.get_project_dir()
        path_to_queries += '/' + project.project_query_sequences
        queries = open(path_to_queries, 'r')
        lines = queries.readlines()
        queries.close()

        queries = []
        for line in lines:
            if ">" in line:
                queries.append(line.split(" ")[0].split(".")[0].split(">")[1])

        Entrez.email = user_email

        print(path_to_queries,queries)
        handle = Entrez.efetch(db="protein", id=queries, retmode="xml")
        record = Entrez.read(handle)
        handle.close()

        prot_to_pfam = {'QUERIES': [], 'PFAM': {}, 'TIGR': {}, 'REFSEQ': {}}
        for i in range(len(record)):
            prot_to_pfam['QUERIES'].append(record[i]['GBSeq_locus'])
            # record[i]['GBSeq_locus']
            pfam = 'http://pfam.xfam.org/family/'
            jvci = 'https://www.ncbi.nlm.nih.gov/genome/annotation_prok/evidence/'
            refseq = 'https://www.ncbi.nlm.nih.gov/protein/' + str(record[i]['GBSeq_locus'])
            prot_to_pfam['REFSEQ'][record[i]['GBSeq_locus']] = refseq

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
            #print(record[i]['GBSeq_comment'])
            #print("#############")
            #json.dumps
        return prot_to_pfam
    except Exception as e:
        raise Exception("couldn't parse entrez.efetch with query ids with exception : {}".format(e))