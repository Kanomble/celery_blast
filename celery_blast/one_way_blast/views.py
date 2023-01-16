import os
from django.shortcuts import render, redirect, HttpResponse
from django.http import JsonResponse
from django.contrib.auth.decorators import login_required
from django.conf import settings
from blast_project.views import failure_view, success_view
from .forms import OneWayProjectCreationForm, BlastSettingsForm, OneWayRemoteProjectCreationForm
from .py_project_creation import create_one_way_blast_project, create_one_way_remote_blast_project
from .py_services import delete_one_way_blast_project_and_associated_directories_by_id,\
    delete_one_way_remote_blast_project_and_associated_directories_by_id, get_one_way_html_results
from .py_django_db_services import get_one_way_project_by_id, get_one_way_remote_project_by_id
from .py_biopython import calculate_pfam_and_protein_links_from_one_way_queries
from .tasks import execute_one_way_blast_project, execute_one_way_remote_blast_project

#TODO documentation
@login_required(login_url='login')
def one_way_blast_project_creation_view(request):
    try:
        if request.method == 'POST':
            print("HELLO!")
            if request.POST['project_type'] == 'local':
                project_creation_form = OneWayProjectCreationForm(
                    request.user, request.POST, request.FILES)
                blast_settings_form = BlastSettingsForm(request.POST)
                remote_project_creation_form = OneWayRemoteProjectCreationForm(request.user)

                if project_creation_form.is_valid() and blast_settings_form.is_valid():
                    create_one_way_blast_project(
                        user=request.user,
                        project_form=project_creation_form,
                        settings_form=blast_settings_form)

                    return success_view(request)
            elif request.POST['project_type'] == 'remote':
                project_creation_form = OneWayProjectCreationForm(request.user)
                remote_project_creation_form = OneWayRemoteProjectCreationForm(
                    request.user, request.POST, request.FILES)
                blast_settings_form = BlastSettingsForm(request.POST)
                if remote_project_creation_form.is_valid() and blast_settings_form.is_valid():
                    create_one_way_remote_blast_project(
                        request.user,
                        remote_project_creation_form,
                        blast_settings_form)

                    return success_view(request)

        else:
            project_creation_form = OneWayProjectCreationForm(request.user)
            blast_settings_form = BlastSettingsForm()
            remote_project_creation_form = OneWayRemoteProjectCreationForm(request.user)

        context = {'OneWayProjectCreationForm':project_creation_form,
                   'BlastSettingsForm':blast_settings_form,
                   'OneWayRemoteProjectCreationForm':remote_project_creation_form,
                   'BlastRemoteSettingsForm':blast_settings_form
                  }

        return render(request,'one_way_blast/one_way_blast_creation_dashboard.html',context)
    except Exception as e:
        return failure_view(request,e)

#TODO documentation
@login_required(login_url='login')
def one_way_project_details_view(request, project_id:int):
    #try:
    blast_project = get_one_way_project_by_id(project_id)
    bokeh_plot_template = "one_way_blast/"+str(project_id)+'/bokeh_plot.html'
    context = {'OneWayBlastProject':blast_project,
                'Database':blast_project.project_database,
               'BokehPlot':bokeh_plot_template,
               }
    return render(request,'one_way_blast/one_way_blast_details.html',context)
    #except Exception as e:
    #    return failure_view(request,e)

@login_required(login_url='login')
def one_way_remote_project_details_view(request, project_id):
    try:

        blast_project = get_one_way_remote_project_by_id(project_id)
        bokeh_plot_template = "one_way_blast/remote_searches/"+str(project_id)+'/bokeh_plot.html'
        #prot_to_pfam = calculate_pfam_and_protein_links_from_queries(request.user.email,project_id)
        context = {'OneWayRemoteBlastProject':blast_project,
                   'BokehPlot':bokeh_plot_template}
                   #'ProtPfam':prot_to_pfam}
        return render(request,'one_way_blast/one_way_remote_blast_details.html',context)
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

#TODO documentation
@login_required(login_url='login')
def one_way_remote_project_delete_view(request, project_id):
    try:
        delete_one_way_remote_blast_project_and_associated_directories_by_id(project_id)
        return success_view(request)
    except Exception as e:
        return failure_view(request,e)

#TODO documentation
def ajax_one_way_wp_to_links(request, project_id, remote):
    try:
        if request.is_ajax and request.method == "GET":
            #progress = read_database_download_and_format_logfile(database_id)
            prot_to_pfam = calculate_pfam_and_protein_links_from_one_way_queries(request.user.email,project_id,remote)
            return JsonResponse(prot_to_pfam,status=200)
        return JsonResponse({"ERROR":"NOT OK"},status=200)
    except Exception as e:
        return JsonResponse({"error": "{}".format(e)}, status=400)

#TODO documentation
@login_required(login_url='login')
def execute_one_way_blast_project_view(request, project_id):
    try:
        if request.method == 'POST':
            execute_one_way_blast_project.delay(project_id)
        return redirect('one_way_project_details',project_id=project_id)
    except Exception as e:
        return failure_view(request,e)

#TODO documentation
@login_required(login_url='login')
def execute_one_way_remote_blast_project_view(request, project_id):
    try:
        if request.method == 'POST':
            execute_one_way_remote_blast_project.delay(project_id)
        return redirect('one_way_remote_project_details',project_id=project_id)
    except Exception as e:
        return failure_view(request,e)

#TODO documentation
@login_required(login_url='login')
def load_one_way_result_html_table_view(request, project_id, remote):
    try:
        html_data = get_one_way_html_results(project_id, "blast_results.html", remote)
        return HttpResponse(html_data)
    except Exception as e:
        return failure_view(request, e)

#TODO documentation
@login_required(login_url='login')
def one_way_download_target_sequences(request, project_id, project_type, filename='',):
    try:
        if filename != '':
            if project_type == "one_way_blast":
                filepath = '/blast/reciprocal_blast/media/' + project_type + '/' + str(project_id) + '/' + filename
                if os.path.isfile(filepath):
                    with open(filepath,'r') as download_file:
                        content = [str(line) for line in download_file.readlines()]
                        content = ''.join(content)
                    response = HttpResponse(content,content_type="text/plain")
                    response['Contnt-Disposition'] = "attachment; filename={}".format(filename)
                else:
                    raise FileNotFoundError

            elif project_type == 'remote_searches':
                filepath = '/blast/reciprocal_blast/media/one_way_blast/' + project_type + '/' + str(project_id) + '/' + filename
                if os.path.isfile(filepath):
                    with open(filepath, 'r') as download_file:
                        content = [str(line) for line in download_file.readlines()]
                        content = ''.join(content)
                    response = HttpResponse(content, content_type="text/plain")
                    response['Contnt-Disposition'] = "attachment; filename={}".format(filename)
                else:
                    raise FileNotFoundError
            else:
                raise Exception("There is no such project_type available!")
            return response

    except Exception as e:
        return failure_view(request,e)