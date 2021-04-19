from django.shortcuts import render

# Create your views here.
from django.shortcuts import render, redirect
from django.contrib.auth.decorators import login_required
from .py_refseq_transactions import get_databases_with_tasks, create_blastdatabase_table_and_directory, \
                                    get_databases_without_tasks
from .py_services import refseq_file_exists
from .forms import RefseqDatabaseForm
from .tasks import download_refseq_assembly_summary_file
from blast_project.views import failure_view
from blast_project.py_services import delete_blastdb_and_associated_directories_by_id

''' dashboard
    
    dashboard for refseq databases, lists all available BLAST databases, 
    download button for the refseq assembly summary file
    
    url: refseq_transactions_dashboard

'''
@login_required(login_url='login')
def dashboard(request):
    try:
        context = {}

        if(refseq_file_exists()):
            context['refseq_exists'] = True

        refseq_database_form = RefseqDatabaseForm()
        executed_databases = get_databases_with_tasks()
        not_executed_databases = get_databases_without_tasks()

        context['RefseqDatabaseForm'] = refseq_database_form
        context['ActiveBlastDatabases'] = executed_databases
        context['UnactiveBlastDatabases'] = not_executed_databases

        return render(request, 'refseq_transactions/refseq_transactions_dashboard.html', context)
    except Exception as e:
        return failure_view(request,e)

''' download_refseq_assembly_summary_view
    
    triggers a @shared_task function in order to download the current refseq summary file
    from the ncbi ftp directory
    redirects to the refseq_transactions_dashboard url e.g. the dashboard view above
    
    url: download_refseq_assembly_summary
    
'''
@login_required(login_url='login')
def download_refseq_assembly_summary_view(request):
    try:
        download_refseq_assembly_summary_file()
        return redirect('refseq_transactions_dashboard')
    except Exception as e:
        return failure_view(request,e)

''' download_refseq_assemblies_view
    
    triggers a @shared function in order to start a snakemake task
    if the form is valid the user gets redirected to the refseq_transactions_dashboard
    if not the form validation errors are rendered
    
    url: download_refseq_assemblies
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

    except Exception as e:
        return failure_view(request,e)

@login_required(login_url='login')
def delete_blast_database_model_and_directory(request,database_id):
    try:
        delete_blastdb_and_associated_directories_by_id(database_id)
        return redirect('refseq_transactions_dashboard')
    except Exception as e:
        return failure_view(request,e)