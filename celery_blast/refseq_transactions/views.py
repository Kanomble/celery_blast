# MAIN VIEWS REFSEQ TRANSACTIONS FOR BLAT DATABASES
import json

import pandas as pd
from blast_project.py_django_db_services import get_database_by_id
from blast_project.py_services import delete_blastdb_and_associated_directories_by_id, upload_file
from blast_project.setup_state import cathi_setup_is_running
from blast_project.views import failure_view
from django.contrib import messages
from django.contrib.auth.decorators import login_required
from django.http import Http404, JsonResponse
from django.shortcuts import render, redirect
from django.urls import reverse
from django.views.decorators.csrf import csrf_exempt
from django_celery_results.models import TaskResult
from time import sleep
from .forms import RefseqDatabaseForm
from .py_refseq_transactions import get_downloaded_databases, get_databases_in_progress, \
    get_databases_without_tasks, read_database_table_by_database_id_and_return_json, get_failed_tasks, \
    create_blastdb_dir_and_table_based_on_user_selection
from .py_services import refseq_file_exists, get_database_download_and_formatting_task_result_progress, \
    genbank_file_exists
from .tasks import create_blast_database_preview, download_refseq_assembly_summary, download_blast_databases_based_on_summary_file
from celery_blast.settings import TAXONOMIC_NODES


SUMMARY_DOWNLOAD_PROGRESS = {
    'PENDING': 5,
    'STARTED': 10,
    'PROGRESS': 50,
    'RETRY': 50,
    'SUCCESS': 100,
    'FAILURE': 'ERROR',
}

PREVIEW_CREATION_PROGRESS = {
    'PENDING': 5,
    'STARTED': 10,
    'PROGRESS': 50,
    'RETRY': 50,
    'SUCCESS': 100,
    'FAILURE': 'ERROR',
}


def _summary_download_progress_from_result(task_result):
    return _task_progress_from_result(task_result, SUMMARY_DOWNLOAD_PROGRESS)


def _task_progress_from_result(task_result, default_progress):
    progress = default_progress.get(task_result.status, 0)
    description = task_result.status

    if task_result.result:
        try:
            payload = json.loads(task_result.result)
            progress_payload = payload.get('progress', {})
            progress = progress_payload.get('percent', progress)
            description = progress_payload.get('description', description)
        except (TypeError, ValueError, AttributeError):
            pass

    return progress, description


def _serializable_refseq_form_data(form):
    cleaned_data = dict(form.cleaned_data)
    taxid_file = cleaned_data.get('taxid_file')
    if taxid_file:
        taxid_file_path = TAXONOMIC_NODES + taxid_file.name
        upload_file(taxid_file, taxid_file_path)
        cleaned_data['taxid_file'] = None
        cleaned_data['taxid_uploaded_file'] = taxid_file.name

    return cleaned_data


def build_dashboard_context(user, refseq_database_form=None):
    context = {}
    context['setup_running'] = cathi_setup_is_running()

    if refseq_file_exists():
        context['refseq_exists'] = True

    if genbank_file_exists():
        context['genbank_exists'] = True

    if refseq_database_form is None:
        refseq_database_form = RefseqDatabaseForm(user)

    context['RefseqDatabaseForm'] = refseq_database_form
    context['ActiveBlastDatabases'] = get_downloaded_databases()
    context['DownloadInProgressBlastDatabases'] = get_databases_in_progress()
    context['InactiveBlastDatabases'] = get_databases_without_tasks()
    context['FailedDatabases'] = get_failed_tasks()
    return context


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
        context = build_dashboard_context(request.user)
        return render(request, 'refseq_transactions/refseq_transactions_dashboard.html', context)
    except Exception as e:
        return failure_view(request, e)

'''update_ncbi_databases

    This view triggers a download process of the RefSeq and GenBank assembly summary files.

'''
@login_required(login_url='login')
def update_ncbi_databases_view(request):
    try:
        if cathi_setup_is_running():
            messages.warning(
                request,
                'CATHI setup is still running. RefSeq and GenBank summary downloads are available after setup finishes.',
            )
            return redirect('refseq_transactions_dashboard')

        download_refseq_assembly_summary.delay('RefSeq')
        download_refseq_assembly_summary.delay('GenBank')
        return redirect('refseq_transactions_dashboard')
    except Exception as e:
        return failure_view(request,e)

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
        if cathi_setup_is_running():
            messages.warning(
                request,
                'CATHI setup is still running. RefSeq and GenBank summary downloads are available after setup finishes.',
            )
            return redirect('refseq_transactions_dashboard')

        task = download_refseq_assembly_summary.delay(summary_file)
        return redirect(
            'assembly_summary_download_progress',
            summary_file=summary_file,
            task_id=task.id,
        )
    except Exception as e:
        return failure_view(request, e)


@login_required(login_url='login')
def assembly_summary_download_progress_view(request, summary_file: str, task_id: str):
    if summary_file not in {'RefSeq', 'GenBank'}:
        raise Http404("Unknown assembly summary file")

    context = {
        'summary_file': summary_file,
        'task_id': task_id,
        'progress_url': reverse(
            'assembly_summary_download_progress_ajax',
            kwargs={'task_id': task_id},
        ),
        'redirect_url': reverse('refseq_transactions_dashboard'),
    }
    return render(request, 'refseq_transactions/assembly_summary_download_progress.html', context)


@login_required(login_url='login')
def assembly_summary_download_progress_ajax(request, task_id: str):
    try:
        task_result = TaskResult.objects.get(task_id=task_id)
        progress, description = _summary_download_progress_from_result(task_result)
        return JsonResponse({
            'progress': progress,
            'status': task_result.status,
            'description': description,
            'complete': task_result.status == 'SUCCESS',
            'failed': task_result.status == 'FAILURE',
        })
    except TaskResult.DoesNotExist:
        return JsonResponse({
            'progress': SUMMARY_DOWNLOAD_PROGRESS['PENDING'],
            'status': 'PENDING',
            'description': 'waiting for worker',
            'complete': False,
            'failed': False,
        })


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
                cleaned_form_data = _serializable_refseq_form_data(refseq_database_form)
                task = create_blast_database_preview.delay(cleaned_form_data)
                return redirect(
                    'database_preview_creation_progress',
                    task_id=task.id,
                )

            # validation error
            else:
                context = build_dashboard_context(request.user, refseq_database_form)

            sleep(2)
            return render(request, 'refseq_transactions/refseq_transactions_dashboard.html', context)
        # should never been executed
        else:
            return failure_view(request,
                                "error creating blast database model and directory, there is no GET request for this view ...")
    except Exception as e:
        return failure_view(request, e)


@login_required(login_url='login')
def database_preview_creation_progress_view(request, task_id: str):
    context = {
        'task_id': task_id,
        'progress_url': reverse(
            'database_preview_creation_progress_ajax',
            kwargs={'task_id': task_id},
        ),
        'redirect_url': reverse('refseq_transactions_dashboard'),
    }
    return render(request, 'refseq_transactions/database_preview_creation_progress.html', context)


@login_required(login_url='login')
def database_preview_creation_progress_ajax(request, task_id: str):
    try:
        task_result = TaskResult.objects.get(task_id=task_id)
        progress, description = _task_progress_from_result(task_result, PREVIEW_CREATION_PROGRESS)
        return JsonResponse({
            'progress': progress,
            'status': task_result.status,
            'description': description,
            'complete': task_result.status == 'SUCCESS',
            'failed': task_result.status == 'FAILURE',
        })
    except TaskResult.DoesNotExist:
        return JsonResponse({
            'progress': PREVIEW_CREATION_PROGRESS['PENDING'],
            'status': 'PENDING',
            'description': 'waiting for worker',
            'complete': False,
            'failed': False,
        })


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
            sleep(1)
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
        is_ajax = request.headers.get('X-Requested-With') == 'XMLHttpRequest'
        if is_ajax:
            if request.method == "GET":
                table_data = read_database_table_by_database_id_and_return_json(database_id)
                return JsonResponse({"data": table_data}, status=200)
        return JsonResponse({"data":"No ajax request!"}, status=400)
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
        is_ajax = request.headers.get('X-Requested-With') == 'XMLHttpRequest'
        if is_ajax:
            if request.method == "GET":
                # progress = read_database_download_and_format_logfile(database_id)
                database = get_database_by_id(database_id)
                if database.database_download_and_format_task.status == 'SUCCESS':
                    return JsonResponse({"progress": 100}, status=200)
                else:
                    progress = get_database_download_and_formatting_task_result_progress(database_id)
                    return JsonResponse({"progress": progress}, status=200)
        return JsonResponse({"progress","No ajax request!"}, status=400)
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
        sleep(2)
        return redirect('refseq_transactions_dashboard')
    except Exception as e:
        return failure_view(request, e)


'''download_and_format_selected_proteomes
    
    This function should trigger a download and format process of the selected proteome entries.
    
    :param database_id
        :type int
        
'''
@csrf_exempt
@login_required(login_url='login')
def download_and_format_selected_proteomes(request, database_id:int):
    try:
        if request.method == "POST":
            form_data = request.POST
            data = form_data.dict()
            try:
                database_id = create_blastdb_dir_and_table_based_on_user_selection(data)
            except Exception as e:
                raise Exception("[-] ERROR creating new BLAST database directory and table with exception: {}".format(e))

            # celery task for downloading and formatting proteomes
            download_blast_databases_based_on_summary_file.delay(database_id)

            sleep(1)
            return redirect('refseq_transactions_dashboard')
        else:
            return failure_view("There is no GET request for this view function ...")

    except Exception as e:
        return failure_view(request, e)
