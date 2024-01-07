import pandas as pd
from django.shortcuts import render, redirect, HttpResponse
from django.contrib import messages
from django.http import JsonResponse
from django.utils import timezone
from django.db import IntegrityError, transaction
from django.contrib.auth import authenticate, login, logout
from django.contrib.auth.models import User
from django.contrib.auth.decorators import login_required
from .view_access_decorators import unauthenticated_user
from .forms import CreateUserForm, CreateTaxonomicFileForm, UploadMultipleFilesGenomeForm, \
    ProjectCreationForm, RemoteProjectCreationForm, BlastSettingsFormBackward, \
    BlastSettingsFormForward, UploadGenomeForm, \
    CreateTaxonomicFileForMultipleScientificNames, SymBLASTProjectSettingsForm
from .tasks import write_species_taxids_into_file, execute_reciprocal_blast_project, \
    execute_makeblastdb_with_uploaded_genomes, download_and_format_taxdb, \
    calculate_database_statistics_task, execute_remote_reciprocal_blast_project
from .py_services import list_taxonomic_files, upload_file, check_if_file_exists, get_remote_html_results, \
    delete_project_and_associated_directories_by_id, get_html_results, check_if_taxdb_exists, \
    read_task_logs_summary_table, download_project_directory, delete_domain_database, delete_remote_project_and_associated_directories_by_id
from .py_project_creation import create_blast_project, create_remote_blast_project
from .py_database_statistics import get_database_statistics_task_status, delete_database_statistics_task_and_output, \
    transform_normalized_database_table_to_json, get_database_selection_task_status
from .py_django_db_services import get_users_blast_projects, get_users_remote_blast_projects, get_project_by_id, \
    save_uploaded_genomes_into_database, get_remote_project_by_id, \
    save_uploaded_multiple_file_genomes_into_database, delete_failed_or_unknown_databases, get_domain_database_model
from one_way_blast.py_django_db_services import get_users_one_way_blast_projects, \
    get_users_one_way_remote_blast_projects
from .py_biopython import calculate_pfam_and_protein_links_from_queries
from refseq_transactions.py_refseq_transactions import get_downloaded_databases
from external_tools.tasks import setup_cathi_download_cdd_refseq_genbank_assembly_files
from Bio import Entrez
from os.path import isfile
# BLAST_PROJECT_DIR DEFAULT = 'media/blast_projects/'
# BLAST_DATABASE_DIR DEFAULT = 'media/databases/'
from celery_blast.settings import BLAST_PROJECT_DIR, BLAST_DATABASE_DIR, REMOTE_BLAST_PROJECT_DIR

'''setup_cathi_view

    This function executes the celery_task setup_cathi_download_cdd_refseq_genbank_assembly_files.
    The function will download and decompress the CDD database and the refseq and genbank assembly
    summary files.
    
'''
@login_required(login_url='login')
def setup_cathi_view(request):
    try:
        if request.method == "POST":
            setup_cathi_download_cdd_refseq_genbank_assembly_files.delay()
        else:
            raise Exception("[-] ERROR. There is no GET method for this view")
        return redirect("blast_project_dashboard")
    except Exception as e:
        return failure_view(request, e)

'''delete_domain_database_view
    
    This function deletes the domain database.
    
'''
@login_required(login_url="login")
def delete_domain_database_view(request):
    try:
        if request.method == "POST":
            delete_domain_database()
            return redirect("blast_project_dashboard")
        else:
            raise Exception("There is no GET method for this view")
    except Exception as e:
        return failure_view(request, e)

'''dashboard_view

    View for the first dashboard page, this page enables monitoring of blast_projects,
    created by the currently logged in user. All context functions do return QuerySets.
    
    :GET
        Uses the get_users_blast_projects utility function which returns a BlastProject Query-Set.
        The Query-Set inherits all projects from the currently logged in user.
        Uses the get_all_blast_databases utility functions to load all BlastDatabase db entries.
        display blast_projects and links to other view functions

'''


@login_required(login_url='login')
def dashboard_view(request):
    try:
        context = {}
        context['domain_database'] = get_domain_database_model()

        if request.method == 'GET':
            users_blast_projects = get_users_blast_projects(request.user.id)
            available_blast_databases = get_downloaded_databases()
            one_way_blast_projects = get_users_one_way_blast_projects(request.user.id)
            one_way_remote_blast_projects = get_users_one_way_remote_blast_projects(request.user.id)
            context['blast_projects'] = users_blast_projects
            context['ActiveBlastDatabases'] = available_blast_databases
            context['OneWayBlastProjects'] = one_way_blast_projects
            context['OneWayRemoteBlastProjects'] = one_way_remote_blast_projects
        return render(request, 'blast_project/blast_project_dashboard.html', context)
    except Exception as e:
        return failure_view(request, e)

'''active_table_view
    
    This view is part of the navigation bar.
'''
@login_required(login_url='login')
def active_table_view(request, selected_table:str):
    try:
        context = {}
        if request.method == 'GET':
            if selected_table == "reciprocal_blast_projects":
                users_blast_projects = get_users_blast_projects(request.user.id)
                context['blast_projects'] = users_blast_projects
                return render(request, 'blast_project/table_blast_projects.html', context)
            elif selected_table == "one_way_projects":
                one_way_blast_projects = get_users_one_way_blast_projects(request.user.id)
                context['OneWayBlastProjects'] = one_way_blast_projects
                return render(request, 'blast_project/table_one_way_blast_projects.html', context)
            elif selected_table == "one_way_remote_projects":
                one_way_remote_blast_projects = get_users_one_way_remote_blast_projects(request.user.id)
                context['OneWayRemoteBlastProjects'] = one_way_remote_blast_projects
                return render(request, 'blast_project/table_one_way_remote_blast_projects.html', context)
            elif selected_table == "databases":
                available_blast_databases = get_downloaded_databases()
                context['ActiveBlastDatabases'] = available_blast_databases
                return render(request, 'blast_project/table_downloaded_and_formatted_blast_databases.html', context)
            elif selected_table == "reciprocal_remote_blast_projects":
                users_remote_blast_projects = get_users_remote_blast_projects(request.user.id)
                context['blast_projects'] = users_remote_blast_projects
                return render(request, 'blast_project/table_remote_blast_projects.html', context)
            else:
                raise Exception("There is no table view with an selected_table value of: {}".format(selected_table))
        else:
            raise Exception("There is no POST request for this view function.")
    except Exception as e:
        return failure_view(request, e)

'''project_creation_view

    This view receives user provided form data and creates a BlastProject model object.
    During BlastProject creation a project specific sub-directory in media/blast_projects/ 
    and in static/images/result_images/ is created. If those sub-directories exists,
    they wont get touched by project creation. The view accepts POST requests.
    It does also check if the taxonomy database is loaded, if not the download_and_format_taxdb is triggered
    (check_if_taxdb_exists, download_and_format_taxdb).
    
    
    :POST
        :FORMS - based on form field data that can be found in blast_project/forms.py
            :ProjectCreationForm
            :BlastSettingsFormForward
            :BlastSettingsFormBackward
            :SymBLASTProjectSettingsForm - additional settings for the MSA and phylogeny procedures
    
    If the POST data is valid the project gets created and the function redirects to the
    associated project_details page. If the POST data is not valid, the function returns the 
    /blast_project/project_creation_dashboard template with validation errors.
'''


@login_required(login_url='login')
def project_creation_view(request):
    try:
        if request.method == 'POST':
            if request.POST['project_type'] == 'local':

                project_creation_form = ProjectCreationForm(request.user, request.POST, request.FILES)
                blast_settings_forward_form = BlastSettingsFormForward(request.POST)
                blast_settings_backward_form = BlastSettingsFormBackward(request.POST)
                symblast_project_settings_form = SymBLASTProjectSettingsForm(request.POST)
                # RETURN PROJECT DETAILS VIEW
                if project_creation_form.is_valid() and blast_settings_forward_form.is_valid() \
                        and blast_settings_backward_form.is_valid() and symblast_project_settings_form.is_valid():
                    query_sequence_file = project_creation_form.cleaned_data['query_sequence_file']
                    query_sequences = project_creation_form.cleaned_data['query_sequence_text']

                    try:
                        with transaction.atomic():
                            # user uploaded a fasta file
                            if query_sequence_file != None:
                                blast_project = create_blast_project(
                                    user=request.user,
                                    query_file_name=query_sequence_file.name,
                                    project_form=project_creation_form,
                                    fw_settings_form=blast_settings_forward_form,
                                    bw_settings_form=blast_settings_backward_form,
                                    symblast_settings_form=symblast_project_settings_form,
                                    filepath=BLAST_PROJECT_DIR)
                                path_to_query_file = BLAST_PROJECT_DIR + str(
                                    blast_project.id) + '/' + query_sequence_file.name
                                upload_file(query_sequence_file, path_to_query_file)

                            # user provided sequence identifier
                            elif query_sequences != '':
                                query_file_name = "target_sequences.faa"
                                blast_project = create_blast_project(
                                    user=request.user,
                                    query_file_name=query_file_name,
                                    project_form=project_creation_form,
                                    fw_settings_form=blast_settings_forward_form,
                                    bw_settings_form=blast_settings_backward_form,
                                    symblast_settings_form=symblast_project_settings_form,
                                    filepath=BLAST_PROJECT_DIR)

                                path_to_query_file = BLAST_PROJECT_DIR + str(blast_project.id) + '/' + query_file_name
                                if type(query_sequences) != Entrez.Parser.ListElement:
                                    raise Exception("wrong protein data for form field query_sequence_text")

                                with open(path_to_query_file, 'w') as qfile:
                                    for rec in query_sequences:
                                        txt = ">" + rec['GBSeq_primary-accession'] + " " + rec['GBSeq_definition'] + "\n" + \
                                              rec[
                                                  'GBSeq_sequence'] + "\n"
                                        qfile.write(txt)

                    except IntegrityError as e:
                        return failure_view(request, e)

                    return redirect('project_details', project_id=blast_project.id)



                else:  # RETURN PROJECT CREATION VIEW WITH VALIDATION ERROR

                    # gets returned at the end of this function if validation errors occurred
                    context = {'ProjectCreationForm': project_creation_form,
                               'BlastSettingsForwardForm': blast_settings_forward_form,
                               'BlastSettingsBackwardForm': blast_settings_backward_form,
                               'SymBLASTProjectSettingsForm':symblast_project_settings_form,
                               'taxdb': check_if_taxdb_exists()}
            else:
                # REMOTE BLAST PROJECT CREATION
                remote_project_creation_form = RemoteProjectCreationForm(request.user, request.POST, request.FILES)
                blast_settings_forward_form = BlastSettingsFormForward(request.POST)
                blast_settings_backward_form = BlastSettingsFormBackward(request.POST)
                symblast_project_settings_form = SymBLASTProjectSettingsForm(request.POST)
                # RETURN PROJECT DETAILS VIEW
                if remote_project_creation_form.is_valid() and blast_settings_forward_form.is_valid() \
                        and blast_settings_backward_form.is_valid() and symblast_project_settings_form.is_valid():
                    query_sequence_file = remote_project_creation_form.cleaned_data['r_query_sequence_file']
                    query_sequences = remote_project_creation_form.cleaned_data['r_query_sequence_text']
                    try:
                        with transaction.atomic():
                            # user uploaded a fasta file
                            if query_sequence_file != None:
                                blast_project = create_remote_blast_project(
                                    user=request.user,
                                    query_file_name=query_sequence_file.name,
                                    project_form=remote_project_creation_form,
                                    fw_settings_form=blast_settings_forward_form,
                                    bw_settings_form=blast_settings_backward_form,
                                    symblast_settings_form=symblast_project_settings_form,
                                    filepath=REMOTE_BLAST_PROJECT_DIR)
                                path_to_query_file = REMOTE_BLAST_PROJECT_DIR + str(
                                    blast_project.id) + '/' + query_sequence_file.name
                                upload_file(query_sequence_file, path_to_query_file)

                            # user provided sequence identifier
                            elif query_sequences != '':
                                query_file_name = "target_sequences.faa"
                                blast_project = create_remote_blast_project(
                                    user=request.user,
                                    query_file_name=query_file_name,
                                    project_form=remote_project_creation_form,
                                    fw_settings_form=blast_settings_forward_form,
                                    bw_settings_form=blast_settings_backward_form,
                                    symblast_settings_form=symblast_project_settings_form,
                                    filepath=REMOTE_BLAST_PROJECT_DIR)

                                path_to_query_file = REMOTE_BLAST_PROJECT_DIR + str(blast_project.id) + '/' + query_file_name
                                if type(query_sequences) != Entrez.Parser.ListElement:
                                    raise Exception("wrong protein data for form field query_sequence_text")

                                with open(path_to_query_file, 'w') as qfile:
                                    for rec in query_sequences:
                                        txt = ">" + rec['GBSeq_primary-accession'] + " " + rec['GBSeq_definition'] + "\n" + \
                                              rec[
                                                  'GBSeq_sequence'] + "\n"
                                        qfile.write(txt)

                    except IntegrityError as e:
                        return failure_view(request, e)

                    return redirect('remote_project_details', project_id=blast_project.id)

                else:  # RETURN PROJECT CREATION VIEW WITH VALIDATION ERRORS
                    # gets returned at the end of this function if validation errors occurred
                    context = {'RemoteProjectCreationForm': remote_project_creation_form,
                               'BlastSettingsForwardForm': blast_settings_forward_form,
                               'BlastSettingsBackwardForm': blast_settings_backward_form,
                               'SymBLASTProjectSettingsForm': symblast_project_settings_form,
                               'taxdb': check_if_taxdb_exists()}

        else:
            project_creation_form = ProjectCreationForm(request.user)
            remote_project_creation_form = RemoteProjectCreationForm(request.user)

            blast_settings_forward_form = BlastSettingsFormForward()
            blast_settings_backward_form = BlastSettingsFormBackward()
            symblast_project_settings_form = SymBLASTProjectSettingsForm()
            context = {'ProjectCreationForm': project_creation_form,
                       'BlastSettingsForwardForm': blast_settings_forward_form,
                       'BlastSettingsBackwardForm': blast_settings_backward_form,
                       'SymBLASTProjectSettingsForm': symblast_project_settings_form,
                       'RemoteProjectCreationForm': remote_project_creation_form,
                       'taxdb': check_if_taxdb_exists()}

        context['domain_database'] = get_domain_database_model()
        return render(request, 'blast_project/project_creation_dashboard.html', context)
    except Exception as e:
        return failure_view(request, e)


'''project_details_view

    View for project details. All result graphs and additional tasks, that can 
    be performed with the results (RBHs) are accessible via this page.
    
    :GET
    Loads the BlastProject model object associated to the given project_id.
    Returns the project_details_dashboard.html page, which renders the results according
    to the current status of the associated snakemake task. For more details consider
    to take a look at the respective template.
    
    :param project_id
        :type int 
'''


@login_required(login_url='login')
def project_details_view(request, project_id: int):
    try:
        blast_project = get_project_by_id(project_id)
        context = {'BlastProject': blast_project,
                   'Database': blast_project.project_forward_database}
        return render(request, 'blast_project/project_details_dashboard.html', context)
    except Exception as e:
        return failure_view(request, e)

@login_required(login_url='login')
def remote_project_details_view(request, project_id: int):
    try:
        blast_project = get_remote_project_by_id(project_id)
        selection_task = get_database_selection_task_status(project_id, 'remote')

        context = {'BlastProject': blast_project,
                   'Database': blast_project.r_project_forward_database,
                   'project_id':project_id,
                   'selection_task':selection_task
                   }

        if blast_project.r_project_execution_snakemake_task:
            if blast_project.r_project_execution_snakemake_task.status == 'SUCCESS':
                bokeh_plot_template = "blast_projects/remote_projects/" + str(project_id) + '/interactive_bokeh_plot.html'
                context['BokehPlot'] = bokeh_plot_template


        return render(request, 'blast_project/remote_project_details_dashboard.html', context)
    except Exception as e:
        return failure_view(request, e)

'''project_delete_view

    This view deletes the associated BlastProject model and its files 
    in media/blast_projects/BlastProject.id and static/images/result_images/BlastProject.id.
    
    :POST
    Deletes the associated BlastProject and all files.
    
    :GET
    Returns to project_details_view.
    
    :param project_id
        :type int
'''


@login_required(login_url='login')
def project_delete_view(request, project_id: int):
    try:
        if request.method == "POST":
            delete_project_and_associated_directories_by_id(project_id)
            return success_view(request)
        else:
            return project_details_view(request, project_id)
    except Exception as e:
        return failure_view(request, e)


@login_required(login_url='login')
def remote_project_delete_view(request, project_id: int):
    try:
        if request.method == "POST":
            delete_remote_project_and_associated_directories_by_id(project_id)
            return success_view(request)
        else:
            return project_details_view(request, project_id)
    except Exception as e:
        return failure_view(request, e)


'''ajax_wp_to_links --> OBSOLETE
    
    Asynchronous call to NCBI-Server for fetching information of the provided query sequences.
    This function is now obsolete as retrieving information is now done via the snakemake pipeline.
    
'''


def ajax_wp_to_links(request, project_id: int):
    try:
        is_ajax = request.headers.get('X-Requested-With') == 'XMLHttpRequest'
        if is_ajax:
            if request.method == "GET":
                # progress = read_database_download_and_format_logfile(database_id)
                prot_to_pfam = calculate_pfam_and_protein_links_from_queries(request.user.email, project_id)
                return JsonResponse(prot_to_pfam, status=200)
        return JsonResponse({"ERROR": "NOT OK"}, status=200)
    except Exception as e:
        return JsonResponse({"error": "{}".format(e)}, status=400)


'''execute_reciprocal_blast_project_view

    This view is executed by pressing the "execute snakemake" button in the project_details page.
    It triggers a celery_task, which in turn executes snakemake by a Popen command of the subprocess module.
    
    :POST
        Executes the asynchronous celery function: execute_reciprocal_blast_project which resides in tasks.py.

'''


@login_required(login_url='login')
def start_reciprocal_blast_project_view(request, project_id: int):
    try:
        if request.method == 'POST':
            blast_project = get_project_by_id(project_id)
            if blast_project.project_execution_snakemake_task:
                if blast_project.project_execution_snakemake_task.status != 'FAILURE':
                    execute_reciprocal_blast_project.delay(project_id)
                    return redirect('project_details', project_id=project_id)
                else:
                    execute_reciprocal_blast_project.delay(project_id)
            else:
                execute_reciprocal_blast_project.delay(project_id)
        return redirect('project_details', project_id=project_id)
    except Exception as e:
        return failure_view(request, e)


def start_remote_reciprocal_blast_project_view(request, project_id: int):
    try:
        if request.method == 'POST':
            blast_project = get_remote_project_by_id(project_id)
            if blast_project.r_project_execution_snakemake_task:
                if blast_project.r_project_execution_snakemake_task.status != 'FAILURE':
                    execute_remote_reciprocal_blast_project.delay(project_id)
                    return redirect('remote_project_details', project_id=project_id)
                else:
                    execute_remote_reciprocal_blast_project.delay(project_id)
            else:
                execute_remote_reciprocal_blast_project.delay(project_id)

        return redirect('remote_project_details', project_id=project_id)
    except Exception as e:
        return failure_view(request, e)

'''load_reciprocal_result_html_table_view

    This function is triggered by pressing the "Reciprocal Results Table" button within the 
    project details page. It loads the result table with all reciprocal best hits as HTML code
    and returns a HttpResponse, which will redirect the user to a new page in which the table data
    is displayed.
    
'''


@login_required(login_url='login')
def load_reciprocal_result_html_table_view(request, project_id):
    try:
        html_data = get_html_results(project_id, "reciprocal_results.html")
        return HttpResponse(html_data)
    except Exception as e:
        return failure_view(request, e)

@login_required(login_url='login')
def load_remote_reciprocal_result_html_table_view(request, project_id):
    try:
        html_data = get_remote_html_results(project_id, "reciprocal_results.html")
        return HttpResponse(html_data)
    except Exception as e:
        return failure_view(request, e)

''' create_taxonomic_file_view

    view for creation of taxonomic files, produced by the get_species_taxids.sh script.
    
    :GET
        display available taxonomic files
    :POST
        create taxonomic files
        synchronous call of write_species_taxids_into_file
'''


@login_required(login_url='login')
def create_taxonomic_file_view_old(request):
    try:
        taxform = CreateTaxonomicFileForm(request.user)
        if request.method == 'POST':
            taxform = CreateTaxonomicFileForm(request.user, request.POST)
            if taxform.is_valid():
                species_name, taxonomic_nodes = taxform.cleaned_data['species_name']
                task = write_species_taxids_into_file(taxonomic_nodes, species_name + '.taxids')
        taxid_files = list_taxonomic_files()
        taxid_files = zip(taxid_files[0], taxid_files[1])
        context = {'taxform': taxform, 'taxid_files': taxid_files}
        return render(request, 'blast_project/create_taxonomic_file.html', context)
    except Exception as e:
        return failure_view(request, e)


'''create_taxonomic_file_view
    
    view for creation of taxonomic files, produced by the get_species_taxids.sh script.
    
    :GET
        display available taxonomic files
    :POST
        create taxonomic files
        synchronous call of write_species_taxids_into_file

'''


@login_required(login_url='login')
def create_taxonomic_file_view(request):
    try:
        taxform = CreateTaxonomicFileForMultipleScientificNames(request.user)
        if request.method == 'POST':
            taxform = CreateTaxonomicFileForMultipleScientificNames(request.user, request.POST)
            if taxform.is_valid():
                filename = taxform.cleaned_data['filename']
                filename = filename + '.taxids'
                taxonomic_nodes = taxform.cleaned_data['species_names']
                task = write_species_taxids_into_file(taxonomic_nodes, filename)
        taxid_files = list_taxonomic_files()
        taxid_files = zip(taxid_files[0], taxid_files[1])
        context = {'taxform': taxform, 'taxid_files': taxid_files}
        return render(request, 'blast_project/create_taxonomic_file.html', context)

    except Exception as e:
        return failure_view(request, e)


# two upload genome options one for single and one for multiple files
# first view function - upload_genome_view for the single files
# second view function - upload_multiple_genomes_post view
@login_required(login_url='login')
def upload_genome_view(request):
    try:
        if request.method == "POST":
            # the UploadMultipleFilesGenomeForm instance has to reside here because both post requests
            # belong to one template, the upload_genome_files_dashboard.html
            multiple_files_genome_form = UploadMultipleFilesGenomeForm(request.user)
            upload_genome_form = UploadGenomeForm(request.user, request.POST, request.FILES)
            if upload_genome_form.is_valid():
                with transaction.atomic():
                    new_db = save_uploaded_genomes_into_database(
                        database_title=upload_genome_form.cleaned_data['database_title'],
                        database_description=upload_genome_form.cleaned_data['database_description'],
                        genome_file=upload_genome_form.cleaned_data['genome_fasta_file'],
                        assembly_entries=upload_genome_form.cleaned_data['assembly_entries'],
                        assembly_level=upload_genome_form.cleaned_data['assembly_level'],
                        taxonomic_node=upload_genome_form.cleaned_data['taxonomic_node'],
                        assembly_accession=upload_genome_form.cleaned_data['assembly_accession'],
                        user_email=request.user.email,
                        organism_name=upload_genome_form.cleaned_data['organism_name'],
                        taxmap_file=upload_genome_form.cleaned_data['taxmap_file'],
                        organism_file=upload_genome_form.cleaned_data['organism_name_file'],
                        assembly_accession_file=upload_genome_form.cleaned_data['assembly_accessions_file'],
                        assembly_level_file=upload_genome_form.cleaned_data['assembly_level_file']
                    )
                    genome_file_name = upload_genome_form.cleaned_data['database_title'].replace(' ',
                                                                                                 '_').upper() + '.database'

                    # inside transaction atomic blog
                    if upload_genome_form.cleaned_data['taxmap_file'] != None:
                        execute_makeblastdb_with_uploaded_genomes.delay(
                            new_db.id,
                            new_db.path_to_database_file + '/' + genome_file_name,
                            taxmap_file=True)
                    elif upload_genome_form.cleaned_data['taxonomic_node'] != None:
                        execute_makeblastdb_with_uploaded_genomes.delay(
                            new_db.id,
                            new_db.path_to_database_file + '/' + genome_file_name,
                            taxonomic_node=upload_genome_form.cleaned_data['taxonomic_node'])
                    else:
                        raise IntegrityError('couldnt trigger makeblastdb execution ...')
                return success_view(request)
            else:
                context = {'UploadGenomeForm': upload_genome_form,
                           'MultipleFileUploadGenomeForm': multiple_files_genome_form, }
        else:
            multiple_files_genome_form = UploadMultipleFilesGenomeForm(request.user)
            upload_genome_form = UploadGenomeForm(request.user)
            context = {'UploadGenomeForm': upload_genome_form,
                       'MultipleFileUploadGenomeForm': multiple_files_genome_form, }

        return render(request, 'blast_project/upload_genome_files_dashboard.html', context)
    except Exception as e:
        # deletes all failed or unknown subdirectories within the database directories
        # all directories without a corresponding database id
        returncode = delete_failed_or_unknown_databases()
        if returncode == 0:
            return failure_view(request, e)
        else:
            return failure_view(request,
                                "Error deleting unknown or failed databases, clean up your database and directories manually. {}".format(
                                    e))


# TODO implement view for disentangling big genome upload view ...
@login_required(login_url='login')
def upload_multiple_genomes_view(request):
    try:
        if request.method == 'POST':
            # the UploadGenomeForm instance has to reside here because both post requests
            # belong to one template, the upload_genome_files_dashboard.html
            upload_genome_form = UploadGenomeForm(request.user)
            extra_field_count = request.POST.get('extra_field_count')
            multiple_files_genome_form = UploadMultipleFilesGenomeForm(request.user, request.POST, request.FILES,
                                                                       extra=extra_field_count)
            if multiple_files_genome_form.is_valid():
                with transaction.atomic():

                    new_db = save_uploaded_multiple_file_genomes_into_database(multiple_files_genome_form.cleaned_data,
                                                                               int(extra_field_count) + 1,
                                                                               request.user.email)
                    genome_file_name = new_db.database_name.replace(' ', '_').upper() + '.database'

                    execute_makeblastdb_with_uploaded_genomes.delay(
                        new_db.id,
                        new_db.path_to_database_file + '/' + genome_file_name,
                        taxmap_file=True)

                    return success_view(request)

            else:
                context = {'UploadGenomeForm': upload_genome_form,
                           'MultipleFileUploadGenomeForm': multiple_files_genome_form, }
                return render(request, 'blast_project/upload_genome_files_dashboard.html', context)

        else:
            return failure_view(request, exception="This view is just for POST requests!")
    except Exception as e:
        return failure_view(request, e)


''' registration, login and logout views

    register an account with email and password, email can be used inside biopython functions
    user needs to authenticate otherwise they will get redirected to this login page

'''


# you dont need an account to trigger this view
def login_user(request):
    # login with django default authenticate method
    if request.method == 'POST':
        username = request.POST.get('username')
        password = request.POST.get('password')
        user = authenticate(request, username=username, password=password)
        if user is not None:
            login(request, user)
            return redirect('blast_project_dashboard')
        else:
            messages.info(request, 'Username OR password is incorrect')
            return render(request, 'blast_project/login.html')
    return render(request, 'blast_project/login.html')


@login_required(login_url='login')
def logout_user(request):
    logout(request)
    return redirect('login')


# registration view
def registration_view(request):
    user_form = CreateUserForm()
    if request.method == 'POST':
        user_form = CreateUserForm(request.POST, initial={'last_login': timezone.now()})
        if user_form.is_valid():
            try:
                username = user_form.cleaned_data.get('username')
                password = user_form.cleaned_data.get('password1')
                email = user_form.cleaned_data.get('email')
                user = User.objects.create_user(username, email, password=password, last_login=timezone.now())

                # group = Group.objects.get(name='customer')
                # user.groups.add(group)
                messages.success(request, 'Account was created for ' + username)
                return redirect('login')
            except Exception as e:
                return failure_view(request, e)
    context = {'form': user_form, }
    return render(request, 'blast_project/register.html', context)


''' failure view

    Standard view for exceptions raised by other view functions.
    
    :param exception
        :type str
'''


# if an exception occurres this page is rendered in order to evaluate the exception context
def failure_view(request, exception):
    context = {'exception': exception}
    return render(request, 'blast_project/failure.html', context)


@login_required(login_url='login')
def success_view(request):
    return render(request, 'blast_project/success.html')


'''database_statistics
    
    Function triggers execution of the database statistics optional postprocessing.

    Status not executed: NOTEXEC
    Status ongoing: PROGRESS
    Status finished: SUCCESS
    Status failed: FAILURE
    
    :param project_id
        :type int
'''


@login_required(login_url='login')
def database_statistics_dashboard(request, project_id):
    try:
        task_status = get_database_statistics_task_status(project_id)
        selection_task_status = get_database_selection_task_status(project_id, 'local')

        context = {'project_id': project_id, 'task_status': task_status}

        if task_status == 'SUCCESS':
            context['DatabaseStatisticsBokehPlot'] = str(project_id) + "/" + "interactive_bokeh_plot.html"
            taxonomic_units = ['genus', 'family', 'superfamily', 'order', 'class', 'phylum']
            for unit in taxonomic_units:
                project_path = BLAST_PROJECT_DIR + str(project_id) + "/" + unit + "_database_statistics_normalized.csv"
                if check_if_file_exists(project_path):
                    table = pd.read_csv(project_path, index_col=0, header=0)
                    number = len(table.columns)
                    key = unit + "_number"
                    context[key] = number
                    key_norm = unit + "_normalized_number"
                    number_queries = len(table.index)
                    table = table.round(2)
                    table = table.transpose()[(table == 0.0).sum() != number_queries]
                    number = len(table.transpose().columns)
                    context[key_norm] = number
                else:
                    key = unit + "_number"
                    error_phrase = "table does not exist, please recompute the database statistics by pressing the button"
                    context[key] = error_phrase
        context['selection_task'] = selection_task_status
        # calculate_database_statistics(project_id)
        return render(request, 'blast_project/database_statistics_dashboard.html', context)
    except Exception as e:
        return failure_view(request, e)


'''database_selection_phylogeny_task_status

    Asynchronous ajax call for the status of the selection constrained phylogenetic inference within the database
    statistics interactive bokeh plot.
    
    Status not executed: NOTEXEC
    Status ongoing: PROGRESS
    Status finished: SUCCESS
    Status failed: FAILURE
    
    :param project id
        :type int
    :param remote_or_local
        :type str
'''
@login_required(login_url='login')
def database_selection_phylogeny_task_status(request, project_id, remote_or_local:str):
    try:
        is_ajax = request.headers.get('X-Requested-With') == 'XMLHttpRequest'
        if is_ajax:
            if request.method == "GET":
                data = get_database_selection_task_status(project_id, remote_or_local)
                return JsonResponse({"data": data}, status=200)
        return JsonResponse({"ERROR": "NOT OK"}, status=200)
    except Exception as e:
        return JsonResponse({"error": "{}".format(e)}, status=400)


'''load_database_statistics_for_class_ajax
    
    Function returns database statistics json dataframe of the specified taxonomic unit.
    
'''


@login_required(login_url='login')
def load_database_statistics_for_taxonomic_unit_ajax(request, project_id, taxonomic_unit: str):
    try:
        is_ajax = request.headers.get('X-Requested-With') == 'XMLHttpRequest'
        if is_ajax:
            if request.method == "GET":
                data = transform_normalized_database_table_to_json(project_id, taxonomic_unit)
                return JsonResponse({"data": data}, status=200)
        return JsonResponse({"ERROR": "NOT OK"}, status=200)
    except Exception as e:
        return JsonResponse({"error": "{}".format(e)}, status=400)


'''execute_database_statistics_task

    Function triggers execution of the database statistics optional postprocessing.
    Executes the task function that includes py_database_statistic function database_statistics.

    :param project_id
        :type int
'''


@login_required(login_url='login')
def execute_database_statistics_task(request, project_id: int):
    try:
        taxonomic_units = ['genus', 'family', 'superfamily', 'order', 'class', 'phylum']
        calculate_database_statistics_task.delay(project_id, request.user.email, taxonomic_units)
        return redirect('database_statistics', project_id=project_id)
    except Exception as e:
        return failure_view(request, e)


'''delete_database_statistics

    Function for deletion of files and task object of the database statistics task.
    
    :param project_id
        :type int
 
'''


@login_required(login_url='login')
def delete_database_statistics(request, project_id):
    try:
        logfile = BLAST_PROJECT_DIR + str(project_id) + '/log/delete_database_statistics_task_and_output.log'
        delete_database_statistics_task_and_output(project_id, logfile=logfile)
        return redirect('database_statistics', project_id=project_id)
    except Exception as e:
        return failure_view(request, e)


'''ajax_call_to_logfiles

    This function sends task progress data based on available logfiles to the template.
    
    :param project_id
        :type int

'''
def ajax_call_to_logfiles(request, project_id: int):
    try:
        is_ajax = request.headers.get('X-Requested-With') == 'XMLHttpRequest'
        if is_ajax:
            if request.method == "GET":
                blast_project = get_project_by_id(project_id)
                if blast_project.project_execution_snakemake_task.status != "SUCCESS":
                    logfiles = blast_project.return_list_of_all_logfiles()
                    logfile_table = read_task_logs_summary_table("local")
                    logfile_table = logfile_table.loc[0:17, :]
                    queries = []
                    progress_without_subtasks = []
                    for logfile in logfiles:
                        if len(logfile.split("/")) == 2:
                            query = logfile.split("/")[0]
                            if query not in queries:
                                queries.append(query)
                            progress = logfile_table[logfile_table['logfile'] == query + "/" + logfile]['progress'].values

                            if len(progress) == 1:
                                progress_without_subtasks.append(progress[0])
                        else:
                            progress = logfile_table[logfile_table['logfile'] == logfile]['progress'].values

                            if len(progress) == 1:
                                progress_without_subtasks.append(progress[0])
                    progress_without_subtasks.sort()
                    return JsonResponse({"progress": max(progress_without_subtasks)}, status=200)
                elif blast_project.project_execution_snakemake_task.status == "SUCCESS":
                    return JsonResponse({"progress": 100}, status=200)
        return JsonResponse({"error": "POST not allowed"}, status=200)
    except Exception as e:
        return JsonResponse({"error": "{}".format(e)}, status=400)

def ajax_call_to_remote_logfiles(request, project_id: int):
    try:
        is_ajax = request.headers.get('X-Requested-With') == 'XMLHttpRequest'
        if is_ajax:
            if request.method == "GET":
                blast_project = get_remote_project_by_id(project_id)

                if blast_project.r_project_execution_snakemake_task.status != 'SUCCESS':

                    logfiles = blast_project.return_list_of_all_logfiles()
                    logfile_table = read_task_logs_summary_table("remote")
                    queries = []
                    progress_without_subtasks = []
                    for logfile in logfiles:
                        if len(logfile.split("/")) == 2:
                            query = logfile.split("/")[0]
                            if query not in queries:
                                queries.append(query)
                            progress = logfile_table[logfile_table['logfile'] == query + "/" + logfile]['progress'].values

                            if len(progress) == 1:
                                progress_without_subtasks.append(progress[0])
                        else:
                            progress = logfile_table[logfile_table['logfile'] == logfile]['progress'].values

                            if len(progress) == 1:
                                progress_without_subtasks.append(progress[0])
                    progress_without_subtasks.sort()
                    return JsonResponse({"progress": max(progress_without_subtasks)}, status=200)
                elif blast_project.r_project_execution_snakemake_task.status == 'SUCCESS':
                    return  JsonResponse({"progress":100},status=200)
        return JsonResponse({"error": "POST not allowed"}, status=200)
    except Exception as e:
        return JsonResponse({"error": "{}".format(e)}, status=400)

'''get_domain_database_download_task_status
    
    This ajax call tracks the status of the domain database download task.

'''
def get_domain_database_download_task_status(request):
    try:
        is_ajax = request.headers.get('X-Requested-With') == 'XMLHttpRequest'
        if is_ajax:
            if request.method == "GET":
                domain_database = get_domain_database_model()
                if domain_database.domain_database_download_task_result:
                    return JsonResponse({"progress_status":domain_database.domain_database_download_task_result.status})
                else:
                    return JsonResponse({"progress_status":"NOT_EXEC"})
            else:
                return JsonResponse({"error","POST not allowed"})
        else:
            raise Exception("No ajax request!")
    except Exception as e:
        return JsonResponse({"error":"{}".format(e)}, status=400)

'''send_logfile_content_view
    
    This view function is part of the pipeline monitoring. It is executed
    if the user presses an image within the progress bar of the project.
    
    :param project_id
        :type int
    :param logfile
        :type str
    
'''
@login_required(login_url='login')
def send_logfile_content_view(request, project_id: int, logfile: str) -> HttpResponse:
    try:
        logfile_path = BLAST_PROJECT_DIR + str(project_id) + '/log/' + logfile + ".log"
        if isfile(logfile_path):
            with open(logfile_path, 'r') as lfile:
                lines = lfile.readlines()

            return HttpResponse(lines, content_type="text/plain")
        else:
            return HttpResponse("couldnt find logfile: {} ...".format(logfile_path), content_type="text/plain")
    except Exception as e:
        return failure_view(request, e)


@login_required(login_url='login')
def send_query_sequence_information_view(request, project_id: int) -> HttpResponse:
    try:
        # BLAST_PROJECT_DIR+
        return render(request, str(project_id) + '/query_sequence_information.html')
    except Exception as e:
        return failure_view(request, e)

'''download_project_as_zip_archive_view
    
    This view function compresses the specified directory and creates a HttpResponse with the zipped data directory.
    The directory is simply specified by the project_id.
    
    :param project_id
        :type int
        
    :return response
        :type HttpResponse - with the zipped directory as attachment
'''
@login_required(login_url='login')
def download_project_as_zip_archive_view(request, project_id:int) -> HttpResponse:
    try:
        project_dir = BLAST_PROJECT_DIR + str(project_id)
        blast_project = get_project_by_id(project_id)
        file_wrapper = download_project_directory(project_dir)
        response = HttpResponse(file_wrapper, content_type='application/zip')
        response['Content-Disposition'] = 'attachment; filename="{filename}.zip"'.format(
            filename = blast_project.project_title.replace(" ", "_")
        )
        return response
    except Exception as e:
        return failure_view(request, e)

@login_required(login_url='login')
def download_remote_project_as_zip_archive_view(request, project_id:int) -> HttpResponse:
    try:
        project_dir = REMOTE_BLAST_PROJECT_DIR + str(project_id)
        blast_project = get_remote_project_by_id(project_id)
        file_wrapper = download_project_directory(project_dir)
        response = HttpResponse(file_wrapper, content_type='application/zip')
        response['Content-Disposition'] = 'attachment; filename="{filename}.zip"'.format(
            filename = blast_project.r_project_title.replace(" ", "_")
        )
        return response
    except Exception as e:
        return failure_view(request, e)

'''view_selection_phylogeny

    This function is part of the database statistics bokeh plot dashboard. It is similar to the load_reciprocal_result_view view function
    in blast_project/views.py. It loads the standalone HTML page for the phylogeny.

    :param request
        :type WSGIRequest
    :param project_id
        :type int
    :param query_id
        :type str
    :param remote_or_local
        :type str
        
'''
@login_required(login_url='login')
def view_selection_phylogeny(request, project_id: int, remote_or_local:str):
    try:
        if remote_or_local == 'local':
            html_path = BLAST_PROJECT_DIR
        elif remote_or_local == 'remote':
            html_path = REMOTE_BLAST_PROJECT_DIR
        else:
            raise Exception("[-] ERROR project neither local nor remote")
        html_data = get_html_results(project_id, "selection_sliced_phylogeny.html", html_result_path=html_path)
        return HttpResponse(html_data)
    except Exception as e:
        return failure_view(request, e)