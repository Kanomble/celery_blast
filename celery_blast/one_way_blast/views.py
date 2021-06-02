from django.shortcuts import render, redirect, HttpResponse
from django.contrib.auth.decorators import login_required
from blast_project.views import failure_view, success_view
from .forms import OneWayProjectCreationForm, BlastSettingsForm
from django.db import transaction, IntegrityError
from .py_project_creation import create_one_way_blast_project
from .py_services import delete_one_way_blast_project_and_associated_directories_by_id
from .py_django_db_services import get_one_way_project_by_id
from blast_project.py_services import upload_file


#TODO documentation
@login_required(login_url='login')
def one_way_blast_project_creation_view(request):
    try:
        if request.method == 'POST':
            project_creation_form = OneWayProjectCreationForm(request.user, request.POST, request.FILES)
            blast_settings_form = BlastSettingsForm(request.POST)

            if project_creation_form.is_valid() and blast_settings_form.is_valid():
                query_sequences = request.FILES['query_sequence_file']
                try:
                    with transaction.atomic():

                        blast_project = create_one_way_blast_project(
                            user=request.user,
                            query_file_name=query_sequences.name,
                            project_form=project_creation_form,
                            settings_form=blast_settings_form)

                        path_to_query_file = 'media/one_way_blast/' + str(
                            blast_project.id) + '/' + query_sequences.name

                        upload_file(query_sequences, path_to_query_file)

                except IntegrityError as e:
                    return failure_view(request,e)
                return success_view(request)
        else:
            project_creation_form = OneWayProjectCreationForm(request.user)
            blast_settings_form = BlastSettingsForm()

        context = {'OneWayProjectCreationForm':project_creation_form,
                   'BlastSettingsForm':blast_settings_form}

        return render(request,'one_way_blast/one_way_blast_creation_dashboard.html',context)
    except Exception as e:
        return failure_view(request,e)


#TODO documentation
@login_required(login_url='login')
def one_way_project_details_view(request, project_id):
    try:
        blast_project = get_one_way_project_by_id(project_id)
        #prot_to_pfam = calculate_pfam_and_protein_links_from_queries(request.user.email,project_id)
        context = {'BlastProject':blast_project,
                    'Database':blast_project.project_database,}
                   #'ProtPfam':prot_to_pfam}
        return render(request,'blast_project/project_details_dashboard.html',context)
    except Exception as e:
        return failure_view(request,e)

#TODO documentation
@login_required(login_url='login')
def one_way_project_delete_view(request, project_id):
    try:
        delete_one_way_blast_project_and_associated_directories_by_id(project_id)
        return success_view(request)
    except Exception as e:
        return failure_view(request,e)
