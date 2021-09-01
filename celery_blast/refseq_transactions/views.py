from django.shortcuts import render

# Create your views here.
from django.shortcuts import render, redirect
from django.contrib.auth.decorators import login_required
from django.http import JsonResponse
from django.core import serializers


from .py_refseq_transactions import get_downloaded_databases, get_failed_tasks, get_databases_in_progress,\
                                    get_databases_without_tasks, create_blastdatabase_table_and_directory, \
                                    read_database_table_by_database_id_and_return_json
from .py_services import refseq_file_exists, read_database_download_and_format_logfile, get_database_download_and_formatting_task_result
from .forms import RefseqDatabaseForm
from .tasks import download_refseq_assembly_summary_file, download_blast_databases, download_blast_databases_based_on_summary_file

from blast_project.py_django_db_services import get_database_by_id
from blast_project.views import failure_view
from blast_project.py_services import delete_blastdb_and_associated_directories_by_id
from django.views.decorators.csrf import csrf_exempt


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

        if(refseq_file_exists()):
            context['refseq_exists'] = True

        refseq_database_form = RefseqDatabaseForm()
        executed_databases = get_downloaded_databases()
        not_executed_databases = get_databases_without_tasks()
        download_in_progress_databases = get_databases_in_progress()

        context['RefseqDatabaseForm'] = refseq_database_form
        context['ActiveBlastDatabases'] = executed_databases
        context['DownloadInProgressBlastDatabases'] = download_in_progress_databases
        context['UnactiveBlastDatabases'] = not_executed_databases

        return render(request, 'refseq_transactions/refseq_transactions_dashboard.html', context)
    except Exception as e:
        return failure_view(request,e)

''' download_refseq_assembly_summary_view
    
    triggers a @shared_task function in order to download the current refseq summary file
    from the ncbi ftp directory
    redirects to the refseq_transactions_dashboard url e.g. the dashboard view above
    
    url: download_refseq_assembly_summary
    redirect: refseq_transactions_dashboard.html
    
'''
@login_required(login_url='login')
def download_refseq_assembly_summary_view(request):
    try:
        download_refseq_assembly_summary_file.delay()
        return redirect('refseq_transactions_dashboard')
    except Exception as e:
        return failure_view(request,e)

''' create_blast_database_model_and_directory
    
    creates the blast database model, saves it into the database and 
    simultaneously creates a directory with a summary file of the blast database.
    
    this summary file is used as input for the celery task execution, it lists the assembly_accession, 
    organism_name, taxid, species_taxid, assembly_level and ftp_path.
    
    if the form is valid the user gets redirected to the refseq_transactions_dashboard
    if not the form validation errors are rendered and parsed into the refseq_transactions_dashboard.html template
    
    url: download_refseq_assemblies
    redirect: refseq_transactions_dashboard.html
    template: refseq_transactions_dashboard.html
    method: POST
'''
@login_required(login_url='login')
def create_blast_database_model_and_directory(request):
    try:
        if request.method == 'POST':
            context = {}
            refseq_database_form = RefseqDatabaseForm(request.POST,request.FILES)
            #validate form
            if refseq_database_form.is_valid():
                #do something here if validation succeeds
                create_blastdatabase_table_and_directory(refseq_database_form)
                return redirect('refseq_transactions_dashboard')
            #validation error
            else:
                #user stays at the page because of validation errors
                context['RefseqDatabaseForm'] = refseq_database_form
                return render(request,'refseq_transactions/refseq_transactions_dashboard.html',context)
        # should never been executed
        else:
            return redirect('refseq_transactions_dashboard')
    except Exception as e:
        return failure_view(request,e)

''' delete_blast_database_model_and_directory
    
    triggered by button form submit in the refseq_transactions_dashboard.html template
    
    url: delete_blast_database database_id
    redirect: refseq_transactions_dashboard.html
    method: POST
'''
@login_required(login_url='login')
def delete_blast_database_model_and_directory(request,database_id):
    try:
        if request.method == "POST":
            delete_blastdb_and_associated_directories_by_id(database_id)
            return redirect('refseq_transactions_dashboard')
        # should never been executed
        else:
            return redirect('refseq_transactions_dashboard')
    except Exception as e:
        return failure_view(request,e)

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
def display_blast_database_details_view(request,database_id):
    try:
        context={}
        blastdb = get_database_by_id(database_id)
        context['Database'] = blastdb
        return render(request,'refseq_transactions/datatable_blast_database_details.html',context)
    except Exception as e:
        return failure_view(request,e)

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
            return JsonResponse({"data":table_data}, status=200)
    except Exception as e:
        return JsonResponse({"error": "{}".format(e)}, status=400)


#TODO documentation
def ajax_call_for_database_download_progress(request, database_id):
    try:
        if request.is_ajax and request.method == "GET":
            #progress = read_database_download_and_format_logfile(database_id)
            progress = get_database_download_and_formatting_task_result(database_id)
            return JsonResponse({"progress":progress},status=200)
    except Exception as e:
        return JsonResponse({"error": "{}".format(e)}, status=400)

#TODO documentation
def download_and_format_blast_database(request, database_id):
    try:
        task = download_blast_databases_based_on_summary_file.delay(database_id)
        #task = download_blast_databases.delay(database_id)
        return redirect('refseq_transactions_dashboard')
    except Exception as e:
        return failure_view(request,e)