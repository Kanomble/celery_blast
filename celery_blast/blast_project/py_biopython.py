from Bio import Entrez

''' get_species_taxid_by_name

utilization in the clean_species_name form field of CreateTaxonomicFileForm

:param user_email
    :type str
:param scientific_name
    :type str
'''
def get_species_taxid_by_name(user_email,scientific_name):
    try:
        Entrez.email = user_email
        search = Entrez.esearch(term=scientific_name, db="taxonomy", retmode="xml")
        record = Entrez.read(search)
        taxid = record['IdList'][0]
        return taxid
    except Exception as e:
        raise Exception("[-] There is no taxonomic node defined by your specified scientific name: {} Exception: {}".format(scientific_name,e))
