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
from .tasks import write_species_taxids_into_file, execute_reciprocal_blast_project, execute_makeblastdb_with_uploaded_genomes
from .py_services import list_taxonomic_files, upload_file, \
    delete_project_and_associated_directories_by_id, get_html_results
from .py_project_creation import create_blast_project
from django.db import IntegrityError, transaction

from .py_django_db_services import get_users_blast_projects, get_all_blast_databases, get_project_by_id, save_uploaded_genomes_into_database, \
    save_uploaded_multiple_file_genomes_into_database
from one_way_blast.py_django_db_services import  get_users_one_way_blast_projects, get_users_one_way_remote_blast_projects
from .py_biopython import calculate_pfam_and_protein_links_from_queries
from refseq_transactions.py_refseq_transactions import get_downloaded_databases

''' dashboard

view for the first dashboard page, this page enables monitoring of blast_projects,
created by the currently logged in user.

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

#TODO documentation
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
            project_creation_form = ProjectCreationForm(request.user)
            blast_settings_forward_form = BlastSettingsFormForward()
            blast_settings_backward_form = BlastSettingsFormBackward()

        context = {'ProjectCreationForm':project_creation_form,
                   'BlastSettingsForwardForm':blast_settings_forward_form,
                   'BlastSettingsBackwardForm':blast_settings_backward_form}

        return render(request,'blast_project/project_creation_dashboard.html',context)
    except Exception as e:
        return failure_view(request,e)

#TODO documentation
@login_required(login_url='login')
def project_details_view(request, project_id):
    try:
        blast_project = get_project_by_id(project_id)
        #prot_to_pfam = calculate_pfam_and_protein_links_from_queries(request.user.email,project_id)
        context = {'BlastProject':blast_project,
                   'Database':blast_project.project_forward_database, }
                   #'ProtPfam':prot_to_pfam}
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
        print("HELLO")
        taxform = CreateTaxonomicFileForMultipleScientificNames(request.user)
        print("HELLO")
        if request.method == 'POST':
            taxform = CreateTaxonomicFileForMultipleScientificNames(request.user,request.POST)
            if taxform.is_valid():
                filename = taxform.cleaned_data['filename']
                filename = filename + '.taxids'
                taxonomic_nodes = taxform.cleaned_data['species_names']
                task = write_species_taxids_into_file(taxonomic_nodes,filename)
        print("HELLO1")
        taxid_files = list_taxonomic_files()
        taxid_files = zip(taxid_files[0], taxid_files[1])
        context = {'taxform': taxform, 'taxid_files': taxid_files}
        return render(request, 'blast_project/create_taxonomic_file.html', context)

    except Exception as e:
        return failure_view(request, e)



#upload_types: standard for GET one_file and multiple_files for POST
#TODO documentation rename function to upload_single_genome_ ...
#two upload genome views one for single and one for multiple files
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
def upload_multiple_genomes_post_view(request):
    try:
        if request.method == 'POST':
            upload_genome_form = UploadGenomeForm(request.user)
            extra_field_count = request.POST.get('extra_field_count')
            multiple_files_genome_form = UploadMultipleFilesGenomeForm(request.user, request.POST, request.FILES,
                                                                       extra=extra_field_count)
            if multiple_files_genome_form.is_valid():
                with transaction.atomic():

                    # TODO: implementation of correct functionality for multiple file uploads
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

returned if an exception ocurred within execution of view functions.

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