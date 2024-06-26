from Bio import Entrez

from .py_django_db_services import get_one_way_project_by_id, get_one_way_remote_project_by_id


# TODO documentation - duplicate from blastproject module
def calculate_pfam_and_protein_links_from_one_way_queries(user_email, project_id, remote):
    try:
        if remote == 0:
            project = get_one_way_project_by_id(project_id)
            path_to_queries = project.get_project_dir()
            path_to_queries += '/' + project.project_query_sequences
            queries = open(path_to_queries, 'r')
            lines = queries.readlines()
            queries.close()
        elif remote == 1:
            project = get_one_way_remote_project_by_id(project_id)
            path_to_queries = project.get_project_dir()
            path_to_queries += '/' + project.r_project_query_sequences
            queries = open(path_to_queries, 'r')
            lines = queries.readlines()
            queries.close()
        else:
            raise Exception("[-] ERROR in calculate_pfam_and_protein_links_from_one_way_queries remote != 1 or 0")

        queries = []
        for line in lines:
            if ">" in line:
                queries.append(line.split(" ")[0].split(".")[0].split(">")[1])

        Entrez.email = user_email

        handle = Entrez.efetch(db="protein", id=queries, retmode="xml")
        record = Entrez.read(handle)
        handle.close()

        prot_to_pfam = {'QUERIES': [], 'PFAM': {}, 'TIGR': {}, 'REFSEQ': {}, 'CDD': {}, 'Definition': {}, 'Length': {}}
        for i in range(len(record)):
            prot_to_pfam['QUERIES'].append(record[i]['GBSeq_locus'])
            prot_to_pfam['Definition'][record[i]['GBSeq_locus']] = record[i]['GBSeq_definition']
            prot_to_pfam['Length'][record[i]['GBSeq_locus']] = record[i]['GBSeq_length']
            pfam = 'http://pfam.xfam.org/family/'
            jvci = 'https://www.ncbi.nlm.nih.gov/genome/annotation_prok/evidence/'
            cdd = 'https://www.ncbi.nlm.nih.gov/Structure/sparcle/archview.html?archid='
            refseq = 'https://www.ncbi.nlm.nih.gov/protein/' + str(record[i]['GBSeq_locus'])
            prot_to_pfam['REFSEQ'][record[i]['GBSeq_locus']] = refseq
            try:
                if "Domain architecture ID" in record[i]['GBSeq_comment']:
                    cdd += record[i]['GBSeq_comment'].split('Domain architecture ID')[1].split(';')[0].strip()
                    # print(cdd)
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
