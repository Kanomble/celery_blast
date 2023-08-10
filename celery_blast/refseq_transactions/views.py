# MAIN VIEWS REFSEQ TRANSACTIONS FOR BLAT DATABASES
from blast_project.py_django_db_services import get_database_by_id
from blast_project.py_services import delete_blastdb_and_associated_directories_by_id
from blast_project.views import failure_view
from django.contrib.auth.decorators import login_required
from django.http import JsonResponse
from django.shortcuts import render, redirect
from django.views.decorators.csrf import csrf_exempt

from .forms import RefseqDatabaseForm
from .py_refseq_transactions import get_downloaded_databases, get_databases_in_progress, \
    get_databases_without_tasks, create_blastdatabase_table_and_directory, \
    read_database_table_by_database_id_and_return_json, get_failed_tasks
from .py_services import refseq_file_exists, get_database_download_and_formatting_task_result_progress, \
    genbank_file_exists
from .tasks import download_refseq_assembly_summary, download_blast_databases_based_on_summary_file

''' dashboard
    
    dashboard for refseq databases, lists all available BLAST databases, 
    download button for the refseq assembly summary file
    
    url: refseq_transactions_dashboard
    template: refseq_transactions_dashboard.html
    context: BlastDatabase instances, defined by their TaskResult status
'''


@login_required(login_url='login')
def dashboard(request):
    try:
        context = {}

        if (refseq_file_exists()):
            context['refseq_exists'] = True

        if (genbank_file_exists()):
            context['genbank_exists'] = True

        refseq_database_form = RefseqDatabaseForm(request.user)
        executed_databases = get_downloaded_databases()
        not_executed_databases = get_databases_without_tasks()
        download_in_progress_databases = get_databases_in_progress()
        failed_databases = get_failed_tasks()

        context['RefseqDatabaseForm'] = refseq_database_form
        context['ActiveBlastDatabases'] = executed_databases
        context['DownloadInProgressBlastDatabases'] = download_in_progress_databases
        context['UnactiveBlastDatabases'] = not_executed_databases
        context['FailedDatabases'] = failed_databases

        return render(request, 'refseq_transactions/refseq_transactions_dashboard.html', context)
    except Exception as e:
        return failure_view(request, e)


''' download_refseq_assembly_summary_view
    
    triggers a @shared_task function in order to download the current summary file
    from the ncbi ftp directory
    redirects to the refseq_transactions_dashboard url e.g. the dashboard view above
    
    downloadable summary files: genbank/refseq
    
    url: download_refseq_assembly_summary
    redirect: refseq_transactions_dashboard.html
    
'''


@login_required(login_url='login')
def download_refseq_assembly_summary_view(request, summary_file:str):
    try:
        download_refseq_assembly_summary(summary_file)
        return redirect('refseq_transactions_dashboard')
    except Exception as e:
        return failure_view(request, e)


''' create_blast_database_model_and_directory
    
    creates the blast database model, saves it into the database and 
    simultaneously creates a directory with a summary file of the blast database.
    
    this summary file is used as input for the celery task execution, it lists the assembly_accession, 
    organism_name, taxid, species_taxid, assembly_level and ftp_path.
    
    if the form is valid the user gets redirected to the refseq_transactions_dashboard
    if not the form validation errors are rendered and parsed into the refseq_transactions_dashboard.html template
    
    url: create_refseq_database_metadata
    redirect: refseq_transactions_dashboard.html
    template: refseq_transactions_dashboard.html
    method: POST
'''


@login_required(login_url='login')
def create_blast_database_model_and_directory(request):
    try:
        if request.method == 'POST':
            context = {}
            refseq_database_form = RefseqDatabaseForm(request.user, request.POST, request.FILES)
            # validate form
            if refseq_database_form.is_valid():
                create_blastdatabase_table_and_directory(refseq_database_form)
                return redirect('refseq_transactions_dashboard')

            # validation error
            else:
                print(refseq_database_form.errors)
                if (refseq_file_exists()):
                    context['refseq_exists'] = True

                # user stays at the page because of validation errors
                executed_databases = get_downloaded_databases()
                not_executed_databases = get_databases_without_tasks()
                download_in_progress_databases = get_databases_in_progress()

                context['ActiveBlastDatabases'] = executed_databases
                context['DownloadInProgressBlastDatabases'] = download_in_progress_databases
                context['UnactiveBlastDatabases'] = not_executed_databases
                context['RefseqDatabaseForm'] = refseq_database_form

            return render(request, 'refseq_transactions/refseq_transactions_dashboard.html', context)
        # should never been executed
        else:
            return failure_view(request,
                                "error creating blast database model and directory, there is no GET request for this view ...")
    except Exception as e:
        return failure_view(request, e)


''' delete_blast_database_model_and_directory
    
    triggered by button form submit in the refseq_transactions_dashboard.html template
    
    url: delete_blast_database database_id
    redirect: refseq_transactions_dashboard.html
    method: POST
'''


@login_required(login_url='login')
def delete_blast_database_model_and_directory(request, database_id):
    try:
        if request.method == "POST":
            delete_blastdb_and_associated_directories_by_id(database_id)
            return redirect('refseq_transactions_dashboard')
        else:
            return failure_view(request, "error deleting this database, request method is not a POST method.")
    except Exception as e:
        return failure_view(request, e)


''' display_blast_database_details_view
    
    this view displays informations about the underlying database. 
    the template uses the datatable package https://datatables.net/ in combination with ajax 
    to asynchronously load the summary file assembly entries of the associated blast database into 
    a html table, that is rendered by datatables on the client side.
    
    url: database_details database_id
    template: datatable_blast_database_details.html
    context: One BlastDatabase instance, lists the summary file in a Datatable
    method: GET

'''


@login_required(login_url='login')
def display_blast_database_details_view(request, database_id):
    try:
        context = {}
        blastdb = get_database_by_id(database_id)
        context['Database'] = blastdb
        return render(request, 'refseq_transactions/datatable_blast_database_details.html', context)
    except Exception as e:
        return failure_view(request, e)


''' ajax_call_for_database_details
    
    ajax call of the datatable "ajax" function, take a look at the template.
    reads the blast database summary file with pandas, 
    transforms the dataframe to a json object and send it back to the client
    
    url: ajax_call database_id
    template: datatable_blast_database_details.html
    method: GET
'''


@csrf_exempt
def ajax_call_for_database_details(request, database_id):
    try:
        if request.is_ajax and request.method == "GET":
            table_data = read_database_table_by_database_id_and_return_json(database_id)
            return JsonResponse({"data": table_data}, status=200)
    except Exception as e:
        return JsonResponse({"error": "{}".format(e)}, status=400)


'''ajax_call_for_database_download_progress

    This function is used for displaying the progress of the BLAST database download and formatting procedure.
    The progress is displayed within a HTML DataTables instance. 
    The progress is saved within a celery task database entry.
    
    :params request
        :type request
        
    :params database_id
        :type int
    
    :returns progress, status
        :type JsonResponse
'''


def ajax_call_for_database_download_progress(request, database_id):
    try:
        if request.is_ajax and request.method == "GET":
            # progress = read_database_download_and_format_logfile(database_id)
            database = get_database_by_id(database_id)
            if database.database_download_and_format_task.status == 'SUCCESS':
                return JsonResponse({"progress": 100}, status=200)
            else:
                progress = get_database_download_and_formatting_task_result_progress(database_id)
                return JsonResponse({"progress": progress}, status=200)
    except Exception as e:
        return JsonResponse({"error": "{}".format(e)}, status=400)


'''download_and_format_blast_database
    
    This function triggers the download and formatting procedure for BLAST databases.
    First the user has to crate a database table. This view function spawns a (long running) celery task.
    The task progress can be viewed in the celery_worker container.
    
    :params request
        :type request
        
    :params database_id
        :type int
        
    :returns redirection to refseq_transactions_dashboard
        :type redirection

'''


def download_and_format_blast_database(request, database_id):
    try:
        task = download_blast_databases_based_on_summary_file.delay(database_id)
        return redirect('refseq_transactions_dashboard')
    except Exception as e:
        return failure_view(request, e)
