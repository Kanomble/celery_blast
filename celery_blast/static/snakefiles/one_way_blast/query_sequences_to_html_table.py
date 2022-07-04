from Bio import Entrez
import pandas as pd
from sys import exit
ERRORCODE=3
''' get_target_header

'''
def get_target_header(target_file):
    try:
        with open(target_file,'r') as proteins:
            header = []
            for line in proteins.readlines():
                if line.startswith(">"):
                    if "|" in line:
                        line = line.replace("|"," ")
                    acc = line.split(" ")[0].split(">")[-1].split(".")[0]
                    header.append(acc)
        return header
    except FileNotFoundError:
        raise FileNotFoundError("[-] There is no target fasta file for: {}".format(target_file))

''' fetch_protein_records
Function wrapper for Entrez.efetch function. 
Takes as input a list with maximal 500 protein accession ids and an email
address for the current user.

Returns the Entrez.Parser.ListElement (from Bio import Entrez) that can be used by the
parse_entrez_xml function.
'''
def fetch_protein_records(proteins:list,email:str):
    try:
        Entrez.email = email
        handle = Entrez.efetch(db="protein", id=proteins, retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        return records
    except Exception as e:
        raise Exception("[-] ERROR fetching protein xml data from Biopython with exception: {}".format(e))


''' parse_entrez_xml
function that generates an html table for the provided input sequences. 
Input sequences can derive from different sources, e.g. from BlastProjects or OneWayBlastProjects.
Input sequences should be protein sequences. The function returns an clean dictionary with following fields:

1. queries: list of all available protein accessions (in NCBI databases)

 --> Note: for all of the three dictionaries, you have to provide the accession name as a key
     for a top-level dictionary:

     e.g. --> protein_informations['features'][accession]

2. features: dictionary for all available features of the relevant protein sequence 
    keys --> locations values: note (definition and db_xref (CDD)) protein_informations['features'][accession] = {'locations':[region_name,note,db_xref]}
3. links: dictionary for jvci, cdd, protfam and refseq links
    keys --> jvci, cdd, protfam, refseq values: external links
4. general_information: dictionary for all other informations
    keys --> accession:
'''
def parse_entrez_xml(records) -> dict:
    try:
        protein_informations = {'queries': [], 'features': {}, 'links': {}, 'general_information': {}}
        for record in records:
            accession = record['GBSeq_locus']
            length = record['GBSeq_length']
            definition = record['GBSeq_definition']
            organism = record['GBSeq_organism']
            protein_informations['queries'].append(accession)
            protein_informations['general_information'][accession] = (definition, length, organism)

            # parsing GBSeq_feature_table
            protein_informations['features'][accession] = {}
            for region in record['GBSeq_feature-table']:
                if region['GBFeature_key'] == 'Region':

                    location = region['GBFeature_location']
                    protein_informations['features'][accession][location] = []

                    for feature in region['GBFeature_quals']:
                        if feature['GBQualifier_name'] == 'region_name':
                            region_name = feature['GBQualifier_value']
                            protein_informations['features'][accession][location].append(region_name)
                        elif feature['GBQualifier_name'] == 'note':
                            note = feature['GBQualifier_value']
                            protein_informations['features'][accession][location].append(note)
                        elif feature['GBQualifier_name'] == 'db_xref':
                            CDD = feature['GBQualifier_value']
                            protein_informations['features'][accession][location].append(CDD)

            # links to databases
            pfam = 'http://pfam.xfam.org/family/'
            jvci = 'https://www.ncbi.nlm.nih.gov/genome/annotation_prok/evidence/'
            cdd = 'https://www.ncbi.nlm.nih.gov/Structure/sparcle/archview.html?archid='
            refseq = 'https://www.ncbi.nlm.nih.gov/protein/' + accession

            protein_informations['links'][accession] = {'pfam': '', 'jvci': '', 'cdd': '', 'refseq': ''}
            if "Domain architecture ID" in record['GBSeq_comment']:
                linker = record['GBSeq_comment'].split('Domain architecture ID')[1].split(';')[0].strip()
                cdd += linker
                protein_informations['links'][accession]['cdd'] = cdd
            if "EMBL-EBI" in record['GBSeq_comment']:
                linker = record['GBSeq_comment'].split("EMBL-EBI")[1].split("::")[1].split(";")[0].strip().rstrip()
                pfam += linker
                protein_informations['links'][accession]['pfam'] = pfam
            if "TIGR" in record['GBSeq_comment']:
                linker = 'TIGR' + record['GBSeq_comment'].split("TIGR")[1].split(";")[0].split(".")[0].strip().rstrip()
                jvci += linker

                protein_informations['links'][accession]['jvci'] = jvci
        return protein_informations
    except Exception as e:
        raise Exception("[-] ERROR during parsing entrez records with exception: {}".format(e))

''' create_pandas_df_and_html_table
Function that takes as input the protein_informations dictionary from parse_entrez_xml and
outputs an html and tsf table with all informations for the provided query sequences.

Returns an integer:
    0 : Function returns successfully
'''
def create_pandas_df_and_html_table(proteins: list, protein_informations: dict, path_to_html_output: str, path_to_csv_output:str) -> int:
    try:
        # 'queries':[],'features':{},'links':{},'general_information':{}
        df = pd.DataFrame(
            columns=['Accession ID', 'Organism', 'Sequence Length', 'Sequence Definition', 'Features', 'PFAM', 'CDD',
                     'TIGR'])

        for protein in proteins:
            if protein in protein_informations['queries']:
                jvci = protein_informations['links'][protein]['jvci']
                cdd = protein_informations['links'][protein]['cdd']
                pfam = protein_informations['links'][protein]['pfam']
                refseq = protein_informations['links'][protein]['refseq']

                features = []
                for feature in protein_informations['features'][protein]:
                    assoc = ' '.join(protein_informations['features'][protein][feature])
                    features.append(feature)
                    features.append(assoc)

                values_to_add = {'Accession ID': protein,
                                 'Organism': protein_informations['general_information'][protein][2],
                                 'Sequence Definition': protein_informations['general_information'][protein][0],
                                 'Sequence Length': protein_informations['general_information'][protein][1],
                                 'Features': features,
                                 'PFAM': pfam, 'CDD': cdd, 'TIGR': jvci}

            else:
                jvci, cdd, pfam, refseq = '', '', '', ''
                length = ''
                definition = ''
                features = ''
                organism = ''
                values_to_add = {'Accession ID': protein,
                                 'Organism': organism,
                                 'Sequence Length': length,
                                 'Sequence Definition': definition,
                                 'Features': features,
                                 'PFAM': pfam, 'CDD': cdd, 'TIGR': jvci}

            row_to_add = pd.Series(values_to_add)
            df = df.append(row_to_add, ignore_index=True)


        pd.set_option('colheader_justify', 'left')
        html_string = '''
        <html>
          <head>
            <title>Sequence Information</title>
            <!-- DataTables stylesheets-->
            <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.10.24/css/jquery.dataTables.css" crossorigin="anonymous">
            <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/select/1.3.2/css/select.dataTables.min.css" crossorigin="anonymous">
            <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/buttons/1.7.0/css/buttons.dataTables.min.css" crossorigin="anonymous">
          </head>

          <body>
            <div id="result_table" style="display:none">
                {table}
            </div>
          </body>

            <script src="https://code.jquery.com/jquery-3.6.0.js" integrity="sha256-H+K7U5CnXl1h5ywQfKtSj8PCmoN9aaq30gDh27Xc0jk=" crossorigin="anonymous"></script>
            <!-- input scripts for DataTables: https://datatables.net/ -->
            <script type="text/javascript" charset="utf8" src="https://cdn.datatables.net/1.10.24/js/jquery.dataTables.js"></script>
            <script type="text/javascript" charset="utf8" src="https://cdn.datatables.net/select/1.3.2/js/dataTables.select.min.js"></script>
            <script type="text/javascript" charset="utf8" src="https://cdn.datatables.net/buttons/1.7.0/js/dataTables.buttons.min.js"></script>
            <script type="text/javascript" charset="utf8" src="https://cdn.datatables.net/buttons/1.7.0/js/buttons.html5.min.js"></script>
            <script>
            $(document).ready(function(){{
                var table = document.getElementsByTagName('table');
                table[0].id='myTable'
                $('#myTable').DataTable(
                    {{
                        dom: 'Bfrtip',
                        "lengthMenu": [ 10 ],
                        buttons: [
                            'copy',
                            'csv',
                            'selectAll',
                            'selectNone',
                            'selectRows'
                        ],
                        select: true,
                        "deferRender": true,
                        {datatable}
                    }}
                );
                var result_table = document.getElementById('result_table');
                result_table.style.display = "block";
            }});
            </script>

        </html>
        '''
        datatables = """
        "columns":[
            {"data":"index"},
            {"data":"Accession ID"},
            {"data":"Organism"},
            {"data":"Sequence Length"},
            {"data":"Sequence Definition"},
            {"data":"Features"},
            {"data":"PFAM",fnCreatedCell: function (nTd, sData, oData, iRow, iCol) {
                    if(oData.PFAM) {
                        $(nTd).html("<a href='"+oData.PFAM+"'>"+oData.PFAM+"</a>");
                    }
                }},
            {"data":"CDD",fnCreatedCell: function (nTd, sData, oData, iRow, iCol) {
                    if(oData.CDD) {
                        $(nTd).html("<a href='"+oData.CDD+"'>"+oData.CDD+"</a>");
                    }
                }
            },
            {"data":"TIGR",fnCreatedCell: function (nTd, sData, oData, iRow, iCol) {
                    if(oData.TIGR) {
                        $(nTd).html("<a href='"+oData.TIGR+"'>"+oData.TIGR+"</a>");
                    }
                }
            }
        ]

        """
        # OUTPUT AN HTML FILE
        with open(path_to_html_output, 'w') as f:
            f.write(html_string.format(table=df.to_html(classes='mystyle'), datatable=datatables))

        df.to_csv(path_to_csv_output, sep="\t", header=df.columns)

        return 0
    except Exception as e:
        raise Exception("[-] ERROR during creation of pandas html tables with exception: {}".format(e))

with open(snakemake.log['log'],'w') as logfile:
    try:
        logfile.write("INFO:starting to create query sequence html table\n")
        logfile.write("INFO:reading query sequence file\n")
        proteins = get_target_header(snakemake.input['target_file'])
        logfile.write("INFO:fetching protein records with biopython\n")
        records = fetch_protein_records(proteins,snakemake.params['email'])
        logfile.write("INFO:parsing records and extracting informations\n")
        protein_informations = parse_entrez_xml(records)
        logfile.write("INFO:writing html file with pandas\n")
        create_pandas_df_and_html_table(proteins,protein_informations,snakemake.output['output_html'],snakemake.output['output_csv'])
        logfile.write("DONE\n")
    except Exception as e:
        logfile.write("ERROR:{}\n".format(e))
        exit(ERRORCODE)