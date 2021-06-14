from .models import OneWayBlastProject, OneWayRemoteBlastProject
from .py_django_db_services import get_one_way_remote_project_by_id, get_one_way_project_by_id
from os.path import isdir
from shutil import rmtree
from django.db import IntegrityError, transaction
import pandas as pd

#TODO documentation
def delete_one_way_blast_project_and_associated_directories_by_id(project_id):
    try:
        with transaction.atomic():
            project = OneWayBlastProject.objects.get(id=project_id)
            if isdir('media/one_way_blast/' + str(project_id)):
                rmtree('media/one_way_blast/' + str(project_id))
            if isdir('static/images/result_images/one_way_blast/'+str(project_id)):
                rmtree('static/images/result_images/one_way_blast/'+str(project_id))
            project.delete()
    except Exception as e:
        raise IntegrityError("couldnt delete one way blast project entry : {}".format(e))

#TODO documentation
def delete_one_way_remote_blast_project_and_associated_directories_by_id(project_id):
    try:
        with transaction.atomic():
            project = OneWayRemoteBlastProject.objects.get(id=project_id)
            if isdir('media/one_way_blast/remote_searches/' + str(project_id)):
                rmtree('media/one_way_blast/remote_searches/' + str(project_id))
            if isdir('static/images/result_images/one_way_blast/remote_searches/'+str(project_id)):
                rmtree('static/images/result_images/one_way_blast/remote_searches/'+str(project_id))
            project.delete()
    except Exception as e:
        raise IntegrityError("couldnt delete one way blast project entry : {}".format(e))

#TODO documentation
#loads the reciprocal results table that is written with one of the last rules in the snakefiles
def get_one_way_html_results(project_id,filename, remote):
    try:
        if remote == 0:
            with open("media/one_way_blast/"+str(project_id)+"/"+filename) as res:
                data = res.readlines()
        elif remote == 1:
            with open("media/one_way_blast/remote_searches/"+str(project_id)+"/"+filename) as res:
                data = res.readlines()
        return data
    except Exception as e:
        raise FileNotFoundError("Couldn't read file {} with Exception: {}".format(e))

def filter_blast_table_by_genus(path_to_blast_table, genus):
    result_data = pd.read_table(path_to_blast_table, header=None)
    result_data.columns = ["qseqid", "sseqid", "evalue", "bitscore", "qgi", "sgi", "sacc", "staxids", "sscinames",
                           "scomnames",
                           "stitle"]

    result_data = result_data[result_data['genus'] == genus].drop_duplicates(subset=['qseqid'], keep="first")
    return result_data

def create_html_table_from_pandas_dataframe(dataframe):
    pd.set_option('colheader_justify', 'left')
    html_string = '''
    <html>
      <head>
        <title>BLAST Result Table</title>
        <!-- DataTables stylesheets-->
        <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.10.24/css/jquery.dataTables.css" crossorigin="anonymous">
        <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/select/1.3.2/css/select.dataTables.min.css" crossorigin="anonymous">
        <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/buttons/1.7.0/css/buttons.dataTables.min.css" crossorigin="anonymous">
      </head>

      <body>
        <div id="blast_results_table" style="display:none">
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
                    "lengthMenu": [ 100 ],
                    buttons: [
                        'copy',
                        'csv',
                        'selectAll',
                        'selectNone',
                        'selectRows'
                    ],
                    select: true
                }}
            );
            var result_table = document.getElementById('blast_results_table');
            result_table.style.display = "block";
        }});
        </script>

    </html>
    '''
    return html_string.format(table=dataframe.to_html(classes='mystyle'))

def view_best_blast_results_for_genus(project_id, remote, genus):
    if remote == 0:
        blast_project = get_one_way_remote_project_by_id(project_id)

    elif remote == 1:
        blast_project = get_one_way_project_by_id(project_id)

    blast_table_path = blast_project.get_project_dir() + '/blast_results.table'
    dataframe = filter_blast_table_by_genus(blast_table_path, genus)
    html_dataframe = create_html_table_from_pandas_dataframe(dataframe)
    return html_dataframe