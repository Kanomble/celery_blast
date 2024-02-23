import pandas as pd
from sys import exit

ERRORCODE = 22
'''read_rpsbproc_output

    This function reads the output generated by rpsbproc and creates three pandas dataframes.
    One dataframe containing all queries. One dataframe containing all domains and one dataframe 
    containing all identified sites.

    :param rpsbproc_output_path
        :type str
    :return (query_df, full_domain_df, full_sites_df)
        :type tuple(pd.DataFrame, pd.DataFrame, pd.DataFrame)
'''
def read_rpsbproc_output(rpsbproc_output_path: str) -> tuple:
    try:
        query_dict = {"query_id": [], "seq_type": [], "seq_length": [], "definition_line": []}
        domain_dict = {"session_ordinal": [], "query_id": [], "hit_type": [], "PSSM_id": [], "from": [], "to": [],
                       "evalue": [], "bitscore": [], "accession": [], "short_name": [], "incomplete": [],
                       "superfamily": []}
        sites_dict = {"session_ordinal": [], "query_id": [], "annot_type": [],
                      "title": [], "residue": [], "complete_size": [], "mapped_size": [], "source_domain": []}

        with open(rpsbproc_output_path, "r") as rpsbfile:
            lines = rpsbfile.readlines()

            query = False
            domain = False
            sites = False

            for line in lines:

                line = line.rstrip()

                if line.startswith("QUERY"):
                    query = True
                    domain = False
                    sites = False


                elif line.startswith("DOMAINS"):
                    query = False
                    domain = True
                    sites = False

                elif line.startswith("SITES"):
                    query = False
                    domain = False
                    sites = True

                elif line.startswith("ENDDOMAINS"):
                    domain = False

                elif line.startswith("ENDSITES"):
                    sites = False

                elif line.startswith("ENDQUERY"):
                    query = False

                # there should be just one line for each query - including the key QUERY
                if query == True:
                    query_line = line.split("\t")

                    query_dict["query_id"].append(query_line[1])
                    query_dict["seq_type"].append(query_line[2])
                    query_dict["seq_length"].append(query_line[3])
                    query_dict["definition_line"].append(query_line[4])

                elif domain == True:
                    # therte are multiple lines within the DOMAINS section for each query
                    if line.startswith("DOMAINS") == False:
                        domain_line = line.split("\t")

                        domain_dict["session_ordinal"].append(domain_line[0])
                        domain_dict["query_id"].append(domain_line[1])
                        domain_dict["hit_type"].append(domain_line[2])
                        domain_dict["PSSM_id"].append(
                            "https://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?uid=" + domain_line[3])
                        domain_dict["from"].append(domain_line[4])
                        domain_dict["to"].append(domain_line[5])
                        domain_dict["evalue"].append(domain_line[6])
                        domain_dict["bitscore"].append(domain_line[7])
                        domain_dict["accession"].append(domain_line[8])
                        domain_dict["short_name"].append(domain_line[9])
                        domain_dict["incomplete"].append(domain_line[10])
                        if domain_line[11] != "-":
                            domain_dict["superfamily"].append(
                                "https://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?uid=" + domain_line[11])
                        else:
                            domain_dict["superfamily"].append(domain_line[11])

                elif sites == True:
                    # therte are multiple lines within the SITES section for each query
                    if line.startswith("SITES") == False:
                        sites_line = line.split("\t")

                        sites_dict["session_ordinal"].append(sites_line[0])
                        sites_dict["query_id"].append(sites_line[1])
                        sites_dict["annot_type"].append(sites_line[2])
                        sites_dict["title"].append(sites_line[3])
                        sites_dict["residue"].append(sites_line[4])
                        sites_dict["complete_size"].append(sites_line[5])
                        sites_dict["mapped_size"].append(sites_line[6])
                        sites_dict["source_domain"].append(
                            "https://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?uid=" + sites_line[7])

        # creation of pandas dataframes
        query_df = pd.DataFrame.from_dict(query_dict)
        sites_df = pd.DataFrame.from_dict(sites_dict)
        domain_df = pd.DataFrame.from_dict(domain_dict)

        # fill sites and domains based on query dataframe to ensure there is at least one line for
        # queries without sites or domains
        full_sites_df = query_df.merge(sites_df, on="query_id", how="outer").fillna("no sites available")
        full_domain_df = query_df.merge(domain_df, on="query_id", how="outer").fillna("no sites available")

        return (query_df, full_domain_df, full_sites_df)
    except Exception as e:
        raise Exception("[-] ERROR during rpsbproc output parsing with exception: {}".format(e))


'''create_domains_html

    This function creates a html-datatable document with identified protein domain sites (rpsblast -> rpsbproc).

    :param full_domain_dataframe
        :type pd.DataFrame
    :param savep
        :type str

    :returns returncode
        :type int

'''
def create_domains_html(full_domain_dataframe: pd.DataFrame, savep: str) -> int:
    try:
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
            <div id="rpsbproc_domain_result_table" style="display:none">
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
                table[0].id='rpsbprocDomainTable'
                $('#rpsbprocDomainTable').DataTable(
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
                var result_table = document.getElementById('rpsbproc_domain_result_table');
                result_table.style.display = "block";
            }});
            </script>

        </html>
        '''
        datatables = """
                        "columns":[
                            {"data":"index"},
                            {"data":"query_id"},
                            {"data":"seq_length"},
                            {"data":"definition_line"},
                            {"data":"hit_type"},
                            {"data":"PSSM_id",fnCreatedCell: function (nTd, sData, oData, iRow, iCol) {
                                    var dat = "Link to CDD"
                                    if(oData.PSSM_id != "-") {
                                        console.log(oData.PSSM_id);
                                        $(nTd).html("<a href='"+oData.PSSM_id+"'>"+dat+"</a>");
                                    } else {
                                        $(nTd).html("-");
                                    }
                                }},
                            {"data":"from"},
                            {"data":"to"},
                            {"data":"evalue"},
                            {"data":"bitscore"},
                            {"data":"accession"},
                            {"data":"short_name"},

                            {"data":"superfamily",fnCreatedCell: function (nTd, sData, oData, iRow, iCol) {
                                    var dat = "Link to superfamily"
                                    if(oData.superfamily != "-") {
                                        $(nTd).html("<a href='"+oData.superfamily+"'>"+dat+"</a>");
                                    } else {
                                        $(nTd).html("-");
                                    }
                                }},
                        ]
        """
        with open(savep, 'w') as f:
            f.write(html_string.format(table=full_domain_dataframe.to_html(classes='mystyle'), datatable=datatables))
        return 0
    except Exception as e:
        raise Exception("[-] ERROR creating CDD domain dataframe with exception: {}".format(e))


'''create_sites_html

    This function creates a html-datatable document with identified protein domain sites (rpsblast -> rpsbproc).

    :param full_sites_dataframe
        :type pd.DataFrame
    :param savep
        :type str

    :returns returncode
        :type int

'''
def create_sites_html(full_sites_dataframe: pd.DataFrame, savep: str) -> int:
    try:
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
            <div id="rpsbproc_sites_result_table" style="display:none">
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
                table[0].id='rpsbprocSitesTable'
                $('#rpsbprocSitesTable').DataTable(
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
                var result_table = document.getElementById('rpsbproc_sites_result_table');
                result_table.style.display = "block";
            }});
            </script>

        </html>
        '''
        datatables = """
                        "columns":[
                            {"data":"index"},
                            {"data":"query_id"},
                            {"data":"seq_length"},
                            {"data":"definition_line"},
                            {"data":"annot_type"},
                            {"data":"title"},
                            {"data":"residue"},
                            {"data":"complete_size"},
                            {"data":"mapped_size"},
                            {"data":"source_domain",fnCreatedCell: function (nTd, sData, oData, iRow, iCol) {
                                    var dat = "Link to domain"
                                    console.log(oData.source_domain);
                                    if(oData.source_domain != "no sites available") {
                                        $(nTd).html("<a href='"+oData.source_domain+"'>"+dat+"</a>");
                                    } else {
                                        $(nTd).html("-");
                                    }
                                }},
                        ]
        """
        with open(savep, 'w') as f:
            f.write(html_string.format(table=full_sites_dataframe.to_html(classes='mystyle'), datatable=datatables))
        return 0
    except Exception as e:
        raise Exception("[-] ERROR creating CDD sites dataframe with exception: {}".format(e))


###### MAIN SCRIPT #######

with open(snakemake.log['log'],'w') as logfile:
    try:
        logfile.write("INFO:start parsing rpsbproc output ...\n")
        query_df, full_domain_df, full_sites_df = read_rpsbproc_output(snakemake.input["rpsbproc_output"])
        logfile.write("INFO:DONE parsing rpsbproc output\n")
        logfile.write("INFO:dataframe length:\n\tquery_df:{}\n\tfull_domain_df:{}\n\tfull_sites_df:{}\n".format(
            len(query_df), len(full_domain_df), len(full_sites_df)))
        logfile.write("INFO:starting to slice dataframes\n")
        full_domain_df = full_domain_df[["query_id", "seq_length", "definition_line",
                                         "hit_type", "PSSM_id", "from", "to", "evalue", "bitscore", "accession",
                                         "short_name", "superfamily"]]

        full_sites_df = full_sites_df[["query_id", "seq_length", "definition_line",
                                       "annot_type", "title", "residue", "complete_size", "mapped_size", "source_domain"
                                       ]]
        logfile.write("INFO:done slicing, writing HTML documents ...\n")
        create_domains_html(full_domain_df, savep=snakemake.output["domain_html"])
        create_sites_html(full_sites_df, savep=snakemake.output["sites_html"])

        logfile.write("INFO:writing csv files ...\n")

        full_domain_df.to_csv(snakemake.output["domain_csv"])
        full_sites_df.to_csv(snakemake.output["sites_csv"])

        logfile.write("DONE\n")
    except Exception as e:
        logfile.write("ERROR: Rpsbproc result parsing resulted in an error with exception: {}".format(e))
        exit(ERRORCODE)