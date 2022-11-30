import pandas as pd
from django.shortcuts import render, redirect, HttpResponse
from django.contrib import messages
from django.http import JsonResponse
from django.utils import timezone
from django.contrib.auth.models import User
from .view_access_decorators import unauthenticated_user
from django.contrib.auth.decorators import login_required
from django.contrib.auth import authenticate, login, logout
from .forms import CreateUserForm, CreateTaxonomicFileForm, UploadMultipleFilesGenomeForm, \
    ProjectCreationForm, BlastSettingsFormBackward, BlastSettingsFormForward, UploadGenomeForm, CreateTaxonomicFileForMultipleScientificNames
from .tasks import write_species_taxids_into_file, execute_reciprocal_blast_project, execute_makeblastdb_with_uploaded_genomes, download_and_format_taxdb, \
    calculate_database_statistics_task
from .py_services import list_taxonomic_files, upload_file, check_if_file_exists, \
    delete_project_and_associated_directories_by_id, get_html_results, check_if_taxdb_exists
from .py_project_creation import create_blast_project
from .py_database_statistics import get_database_statistics_task_status, delete_database_statistics_task_and_output,\
    transform_normalized_database_table_to_json
from django.db import IntegrityError, transaction

from .py_django_db_services import get_users_blast_projects, get_project_by_id, save_uploaded_genomes_into_database, \
    save_uploaded_multiple_file_genomes_into_database
from one_way_blast.py_django_db_services import  get_users_one_way_blast_projects, get_users_one_way_remote_blast_projects
from .py_biopython import calculate_pfam_and_protein_links_from_queries
from refseq_transactions.py_refseq_transactions import get_downloaded_databases

''' dashboard

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
        if request.method == 'GET':
            users_blast_projects = get_users_blast_projects(request.user.id)
            available_blast_databases = get_downloaded_databases()
            one_way_blast_projects = get_users_one_way_blast_projects(request.user.id)
            one_way_remote_blast_projects = get_users_one_way_remote_blast_projects(request.user.id)
            context['blast_projects'] = users_blast_projects
            context['ActiveBlastDatabases'] = available_blast_databases
            context['OneWayBlastProjects'] = one_way_blast_projects
            context['OneWayRemoteBlastProjects'] = one_way_remote_blast_projects
        return render(request,'blast_project/blast_project_dashboard.html',context)
    except Exception as e:
        return failure_view(request,e)

'''project_creation_view

    This view receives user provided form data and creates a BlastProject model object.
    During BlastProject creation a project specific sub-directory in media/blast_projects/ 
    and in static/images/result_images/ is created. The view accepts POST requests.
    It does also check if the taxonomy database is loaded, if not the download_and_format_taxdb is triggered
    (check_if_taxdb_exists, download_and_format_taxdb).
    
    
    :POST
        :FORMS - based on form field data that can be found in blast_project/forms.py
            :ProjectCreationForm
            :BlastSettingsFormForward
            :BlastSettingsFormBackward
    
    If the POST data is valid the project gets created and the function redirects to the
    associated project_details page. If the POST data is not valid, the function returns the 
    /blast_project/project_creation_dashboard template with validation errors.
'''
@login_required(login_url='login')
def project_creation_view(request):
    try:
        if request.method == 'POST':
            project_creation_form = ProjectCreationForm(request.user, request.POST, request.FILES)
            blast_settings_forward_form = BlastSettingsFormForward(request.POST)
            blast_settings_backward_form = BlastSettingsFormBackward(request.POST)

            if project_creation_form.is_valid() and blast_settings_forward_form.is_valid() and blast_settings_backward_form.is_valid():
                query_sequences = request.FILES['query_sequence_file']
                try:
                    with transaction.atomic():

                        blast_project = create_blast_project(
                            user=request.user,
                            query_file_name=query_sequences.name,
                            project_form=project_creation_form,
                            fw_settings_form=blast_settings_forward_form,
                            bw_settings_form=blast_settings_backward_form)
                        path_to_query_file = 'media/blast_projects/' + str(
                            blast_project.id) + '/' + query_sequences.name
                        upload_file(query_sequences, path_to_query_file)
                except IntegrityError as e:
                    return failure_view(request,e)

                return redirect('project_details',project_id=blast_project.id)

            else:
                #TODO what happens if taxdb is not there - downloading database
                if check_if_taxdb_exists():
                    taxdb=True
                else:
                    #what happens if task runs into any error?
                    taxdb=False
                    task = download_and_format_taxdb.delay()

                context = {'ProjectCreationForm':project_creation_form,
                           'BlastSettingsForwardForm':blast_settings_forward_form,
                           'BlastSettingsBackwardForm':blast_settings_backward_form,
                           'taxdb':taxdb}

        else:
            if check_if_taxdb_exists():
                project_creation_form = ProjectCreationForm(request.user)
                blast_settings_forward_form = BlastSettingsFormForward()
                blast_settings_backward_form = BlastSettingsFormBackward()

                context = {'ProjectCreationForm':project_creation_form,
                           'BlastSettingsForwardForm':blast_settings_forward_form,
                           'BlastSettingsBackwardForm':blast_settings_backward_form,
                           'taxdb':True}

            else:
                # what happens if task runs into any error?
                context = {'taxdb':False}
                task = download_and_format_taxdb.delay()

        return render(request,'blast_project/project_creation_dashboard.html',context)
    except Exception as e:
        return failure_view(request,e)

'''project_details_view

    View for project details. All result graphs and additional tasks, that can 
    be performed with the results (RBHs) are accessible via this page.
    
    :GET
    Loads the BlastProject model object associated to the given project_id.
    Returns the project_details_dashboard.html page, which renders the results according
    to the current status of the associated snakemake task. For more details consider
    to take a look at the respective template.
    
'''
@login_required(login_url='login')
def project_details_view(request, project_id:int):
    try:
        blast_project = get_project_by_id(project_id)
        context = {'BlastProject':blast_project,
                   'Database':blast_project.project_forward_database}
        return render(request,'blast_project/project_details_dashboard.html',context)
    except Exception as e:
        return failure_view(request,e)

#TODO documentation
@login_required(login_url='login')
def project_delete_view(request, project_id):
    try:
        delete_project_and_associated_directories_by_id(project_id)
        return success_view(request)
    except Exception as e:
        return failure_view(request,e)

#TODO documentation
def ajax_wp_to_links(request, project_id):
    try:
        if request.is_ajax and request.method == "GET":
            #progress = read_database_download_and_format_logfile(database_id)
            prot_to_pfam = calculate_pfam_and_protein_links_from_queries(request.user.email,project_id)
            return JsonResponse(prot_to_pfam,status=200)
        return JsonResponse({"ERROR":"NOT OK"},status=200)
    except Exception as e:
        return JsonResponse({"error": "{}".format(e)}, status=400)

#TODO documentation
@login_required(login_url='login')
def execute_reciprocal_blast_project_view(request, project_id):
    try:
        if request.method == 'POST':
            execute_reciprocal_blast_project.delay(project_id)
        return redirect('project_details', project_id=project_id)
    except Exception as e:
        return failure_view(request,e)

#TODO documentation
@login_required(login_url='login')
def load_reciprocal_result_html_table_view(request, project_id):
    try:
        html_data = get_html_results(project_id, "reciprocal_results.html")
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
            taxform = CreateTaxonomicFileForm(request.user,request.POST)
            if taxform.is_valid():
                species_name,taxonomic_nodes = taxform.cleaned_data['species_name']
                task = write_species_taxids_into_file(taxonomic_nodes,species_name+'.taxids')
        taxid_files = list_taxonomic_files()
        taxid_files = zip(taxid_files[0],taxid_files[1])
        context = {'taxform': taxform, 'taxid_files': taxid_files}
        return render(request, 'blast_project/create_taxonomic_file.html', context)
    except Exception as e:
        return failure_view(request,e)

@login_required(login_url='login')
def create_taxonomic_file_view(request):
    try:
        taxform = CreateTaxonomicFileForMultipleScientificNames(request.user)
        if request.method == 'POST':
            taxform = CreateTaxonomicFileForMultipleScientificNames(request.user,request.POST)
            if taxform.is_valid():
                filename = taxform.cleaned_data['filename']
                filename = filename + '.taxids'
                taxonomic_nodes = taxform.cleaned_data['species_names']
                task = write_species_taxids_into_file(taxonomic_nodes,filename)
        taxid_files = list_taxonomic_files()
        taxid_files = zip(taxid_files[0], taxid_files[1])
        context = {'taxform': taxform, 'taxid_files': taxid_files}
        return render(request, 'blast_project/create_taxonomic_file.html', context)

    except Exception as e:
        return failure_view(request, e)


#TODO documentation rename function to upload_single_genome_ ...
#two upload genome options one for single and one for multiple files
#first view function - upload_genome_view for the single files
#second view function - upload_multiple_genomes_post view
@login_required(login_url='login')
def upload_genome_view(request):
    try:
        if request.method == "POST":
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

                    genome_file_name = upload_genome_form.cleaned_data['database_title'].replace(' ','_').upper()+'.database'
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
                        raise Exception('couldnt trigger makeblastdb execution ...')
                    return success_view(request)
            else:
                context = {'UploadGenomeForm': upload_genome_form,
                           'MultipleFileUploadGenomeForm': multiple_files_genome_form, }
        else:
            multiple_files_genome_form = UploadMultipleFilesGenomeForm(request.user)
            upload_genome_form = UploadGenomeForm(request.user)
            context = {'UploadGenomeForm': upload_genome_form,
                       'MultipleFileUploadGenomeForm': multiple_files_genome_form,}

        #files = request.FILES.getlist('genome_fasta_files')
        return render(request,'blast_project/upload_genome_files_dashboard.html',context)
    except Exception as e:
        return failure_view(request,e)

#TODO implement view for disentangling big genome upload view ...
@login_required(login_url='login')
def upload_multiple_genomes_view(request):
    try:
        if request.method == 'POST':
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
                return render(request,'blast_project/upload_genome_files_dashboard.html',context)

        else:
            return failure_view(request,exception="This view is just for POST requests!")
    except Exception as e:
        return failure_view(request,e)

''' registration, login and logout views
register an account with email and password, email can be used inside biopython functions
user needs to authenticate otherwise they will get redirected to this login page

'''
#login user
@unauthenticated_user #you dont need an account to trigger this view
def login_user(request):

    #login with django default authenticate method
    if request.method == 'POST':
        username = request.POST.get('username')
        password = request.POST.get('password')
        user = authenticate(request, username=username, password=password)
        if user is not None:
            login(request, user)
            return redirect('blast_project_dashboard')
        else:
            messages.info(request,'Username OR password is incorrect')
            return render(request, 'blast_project/login.html')
    return render(request,'blast_project/login.html')

#logout view
def logout_user(request):
    logout(request)
    return redirect('login')

#registration view
@unauthenticated_user
def registration_view(request):
    userForm = CreateUserForm()
    if request.method == 'POST':
        userForm = CreateUserForm(request.POST,initial={'last_login':timezone.now()})
        if userForm.is_valid():
            try:
                username = userForm.cleaned_data.get('username')
                password = userForm.cleaned_data.get('password1')
                email = userForm.cleaned_data.get('email')
                user = User.objects.create_user(username, email, password=password, last_login=timezone.now())

                #group = Group.objects.get(name='customer')
                #user.groups.add(group)
                messages.success(request,'Account was created for '+ username)
                return redirect('login')
            except Exception as e:
                return failure_view(request,e)
    context = {'form': userForm, }
    return render(request,'blast_project/register.html',context)

''' failure view

    Standard view for exceptions raised by other view functions.
    
    :param exception
        :type str
'''
#if an exception occurres this page is rendered in order to evaluate the exception context
def failure_view(request,exception):
    context={'exception':exception}
    return render(request,'blast_project/failure.html', context)

#TODO documentation
@login_required(login_url='login')
def success_view(request):
    return render(request,'blast_project/success.html')

'''database_statistics
    
    Function triggers execution of the database statistics optional postprocessing.
    Executes the task function that includes py_database_statistic function database_statistics.
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
        task_status=get_database_statistics_task_status(project_id)
        context={'project_id':project_id,'task_status':task_status}
        if task_status == 'SUCCESS':
            taxonomic_units = ['genus', 'family', 'superfamily', 'order', 'class', 'phylum']
            for unit in taxonomic_units:
                project_path = "media/blast_projects/" + str(project_id) + "/" + unit + "_database_statistics_normalized.csv"
                if check_if_file_exists(project_path):
                    table = pd.read_csv(project_path,index_col=0,header=0)
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

        #calculate_database_statistics(project_id)
        return render(request,'blast_project/database_statistics_dashboard.html',context)
    except Exception as e:
        return failure_view(request, e)

#TODO documentation
'''database_statistics_details
    
'''
def database_statistics_details(request,project_id:int,taxonomic_unit:str):
    try:
        task_status=get_database_statistics_task_status(project_id)
        context={'project_id':project_id,'task_status':task_status}
        if task_status == 'SUCCESS':
            context_key_altair = 'DatabaseStatisticsAltairPlot'
            context_key_bokeh = 'DatabaseStatisticsBokehPlot'
            context[context_key_altair] = str(project_id) + "/" + taxonomic_unit + "_altair_plot_normalized.html"
            context[context_key_bokeh] = str(project_id) + "/" + taxonomic_unit + "_bokeh_plot.html"
            context['taxonomic_unit'] = taxonomic_unit
        return render(request,'blast_project/database_statistics_details.html',context)
    except Exception as e:
        return failure_view(request,e)

'''load_database_statistics_for_class_ajax
    
    Function returns database statistics json dataframe of the specified taxonomic unit.
    
'''
@login_required(login_url='login')
def load_database_statistics_for_taxonomic_unit_ajax(request, project_id, taxonomic_unit:str):
    try:
        if request.is_ajax and request.method == "GET":
            data = transform_normalized_database_table_to_json(project_id,taxonomic_unit)
            return JsonResponse({"data":data}, status=200)
        return JsonResponse({"ERROR":"NOT OK"},status=200)
    except Exception as e:
        return JsonResponse({"error": "{}".format(e)}, status=400)

'''database_statistics

    Function triggers execution of the database statistics optional postprocessing.
    Executes the task function that includes py_database_statistic function database_statistics.

    :param project_id
        :type int
'''
@login_required(login_url='login')
def execute_database_statistics_task(request, project_id):
    try:
        taxonomic_units = ['genus', 'family', 'superfamily', 'order', 'class', 'phylum']
        calculate_database_statistics_task.delay(project_id,request.user.email,taxonomic_units)
        return redirect('database_statistics',project_id=project_id)
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
        logfile='media/blast_projects/'+str(project_id)+'/log/delete_database_statistics_task_and_output.log'
        delete_database_statistics_task_and_output(project_id,logfile=logfile)
        return redirect('database_statistics',project_id=project_id)
    except Exception as e:
        return failure_view(request,e)


