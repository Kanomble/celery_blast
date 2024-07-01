from os.path import isfile

from blast_project.views import failure_view, success_view
# BLAST_PROJECT_DIR DEFAULT = 'media/one_way_blast/'
from celery_blast.settings import ONE_WAY_BLAST_PROJECT_DIR
from django.contrib.auth.decorators import login_required
from django.http import JsonResponse
from django.shortcuts import render, redirect, HttpResponse

from .forms import OneWayProjectCreationForm, BlastSettingsForm, OneWayRemoteProjectCreationForm
from .py_biopython import calculate_pfam_and_protein_links_from_one_way_queries
from .py_django_db_services import get_one_way_project_by_id, get_one_way_remote_project_by_id
from .py_project_creation import create_one_way_blast_project, create_one_way_remote_blast_project
from .py_services import delete_one_way_blast_project_and_associated_directories_by_id, \
    delete_one_way_remote_blast_project_and_associated_directories_by_id, get_one_way_html_results,\
    read_snakemake_logfile, check_amount_of_blast_hits
from .tasks import execute_one_way_blast_project, execute_one_way_remote_blast_project
from time import sleep

'''one_way_blast_project_creation_view
    
    View function for the creation of one-way BLAST projects. There are two types of projects.
    The first uses local BLAST databases that have been previously formatted using this tool.
    The second uses remote BLAST databases. 
    Form fields are described within the forms.py file of this module.

'''


@login_required(login_url='login')
def one_way_blast_project_creation_view(request):
    try:
        if request.method == 'POST':
            if request.POST['project_type'] == 'local':
                project_creation_form = OneWayProjectCreationForm(
                    request.user, request.POST, request.FILES)
                blast_settings_form = BlastSettingsForm(request.POST)
                remote_project_creation_form = OneWayRemoteProjectCreationForm(request.user)

                if project_creation_form.is_valid() and blast_settings_form.is_valid():
                    blast_project = create_one_way_blast_project(
                        user=request.user,
                        project_form=project_creation_form,
                        settings_form=blast_settings_form)

                    return redirect('one_way_project_details', project_id=blast_project.id)
                else:
                    context = {'OneWayProjectCreationForm': project_creation_form,
                               'BlastSettingsForm': blast_settings_form,
                               'OneWayRemoteProjectCreationForm': remote_project_creation_form,
                               'BlastRemoteSettingsForm': blast_settings_form
                               }
                    return render(request, 'one_way_blast/one_way_blast_creation_dashboard.html', context)

            elif request.POST['project_type'] == 'remote':
                project_creation_form = OneWayProjectCreationForm(request.user)
                remote_project_creation_form = OneWayRemoteProjectCreationForm(
                    request.user, request.POST, request.FILES)
                blast_settings_form = BlastSettingsForm(request.POST)
                if remote_project_creation_form.is_valid() and blast_settings_form.is_valid():
                    blast_project = create_one_way_remote_blast_project(
                        request.user,
                        remote_project_creation_form,
                        blast_settings_form)

                    return redirect('one_way_remote_project_details', project_id=blast_project.id)
                else:
                    context = {'OneWayProjectCreationForm': project_creation_form,
                               'BlastSettingsForm': blast_settings_form,
                               'OneWayRemoteProjectCreationForm': remote_project_creation_form,
                               'BlastRemoteSettingsForm': blast_settings_form
                               }
                    return render(request, 'one_way_blast/one_way_blast_creation_dashboard.html', context)
        else:
            project_creation_form = OneWayProjectCreationForm(request.user)
            blast_settings_form = BlastSettingsForm()
            remote_project_creation_form = OneWayRemoteProjectCreationForm(request.user)

            context = {'OneWayProjectCreationForm': project_creation_form,
                       'BlastSettingsForm': blast_settings_form,
                       'OneWayRemoteProjectCreationForm': remote_project_creation_form,
                       'BlastRemoteSettingsForm': blast_settings_form
                       }

            return render(request, 'one_way_blast/one_way_blast_creation_dashboard.html', context)
    except Exception as e:
        return failure_view(request, e)


'''one_way_project_details_view

    View function that returns the detail website for the selected project.
    The interactive bokeh plot is loaded by this function.
    
'''


@login_required(login_url='login')
def one_way_project_details_view(request, project_id: int):
    try:
        blast_project = get_one_way_project_by_id(project_id)
        context = {}
        if blast_project.project_execution_task_result:
            if blast_project.project_execution_task_result.status == 'FAILURE':
                hits = check_amount_of_blast_hits(project_id, 'local')
                if hits == 0:
                    context['no_hits'] = True
                else:
                    context['no_hits'] = False
            else:
                bokeh_plot_template = "one_way_blast/" + str(project_id) + '/bokeh_plot.html'
                context['BokehPlot'] = bokeh_plot_template
        context['OneWayBlastProject'] = blast_project
        context['Database'] = blast_project.project_database
        return render(request, 'one_way_blast/one_way_blast_details.html', context)
    except Exception as e:
        return failure_view(request, e)


@login_required(login_url='login')
def one_way_remote_project_details_view(request, project_id):
    try:

        blast_project = get_one_way_remote_project_by_id(project_id)
        context = {}
        if blast_project.r_project_execution_task_result:
            if blast_project.r_project_execution_task_result.status == 'FAILURE':
                hits = check_amount_of_blast_hits(project_id, 'remote')
                if hits == 0:
                    context['no_hits'] = True
                else:
                    context['no_hits'] = False
            else:
                bokeh_plot_template = "one_way_blast/remote_searches/" + str(project_id) + '/bokeh_plot.html'
                context['BokehPlot'] = bokeh_plot_template
        context['OneWayRemoteBlastProject'] = blast_project
        return render(request, 'one_way_blast/one_way_remote_blast_details.html', context)
    except Exception as e:
        return failure_view(request, e)


@login_required(login_url='login')
def one_way_project_delete_view(request, project_id):
    try:
        delete_one_way_blast_project_and_associated_directories_by_id(project_id)
        return success_view(request)
    except Exception as e:
        return failure_view(request, e)


@login_required(login_url='login')
def one_way_remote_project_delete_view(request, project_id):
    try:
        delete_one_way_remote_blast_project_and_associated_directories_by_id(project_id)
        return success_view(request)
    except Exception as e:
        return failure_view(request, e)


# TODO documentation
def ajax_one_way_wp_to_links(request, project_id, remote):
    try:
        is_ajax = request.headers.get('X-Requested-With') == 'XMLHttpRequest'
        if is_ajax:
            if request.method == "GET":
                # progress = read_database_download_and_format_logfile(database_id)
                prot_to_pfam = calculate_pfam_and_protein_links_from_one_way_queries(request.user.email, project_id, remote)
                return JsonResponse(prot_to_pfam, status=200)
        return JsonResponse({"ERROR": "NOT OK"}, status=400)
    except Exception as e:
        return JsonResponse({"error": "{}".format(e)}, status=400)


# TODO documentation
@login_required(login_url='login')
def execute_one_way_blast_project_view(request, project_id):
    try:
        if request.method == 'POST':
            execute_one_way_blast_project.delay(project_id)
        sleep(1)
        return redirect('one_way_project_details', project_id=project_id)
    except Exception as e:
        return failure_view(request, e)


# TODO documentation
@login_required(login_url='login')
def execute_one_way_remote_blast_project_view(request, project_id):
    try:
        if request.method == 'POST':
            execute_one_way_remote_blast_project.delay(project_id)
        sleep(1)
        return redirect('one_way_remote_project_details', project_id=project_id)
    except Exception as e:
        return failure_view(request, e)


# TODO documentation
@login_required(login_url='login')
def load_one_way_result_html_table_view(request, project_id, remote):
    try:
        html_data = get_one_way_html_results(project_id, "blast_results.html", remote)
        return HttpResponse(html_data)
    except Exception as e:
        return failure_view(request, e)


'''one_way_download_target_sequences

    This function is executed by triggering the "Download Query Fasta Files" Button at the Local BLAST Project 
    details website.

'''
@login_required(login_url='login')
def one_way_download_target_sequences(request, project_id, project_type, filename=''):
    try:
        if filename != '':
            if project_type == "one_way_blast":
                filepath = '/blast/reciprocal_blast/' + ONE_WAY_BLAST_PROJECT_DIR + str(project_id) + '/' + filename
                if isfile(filepath):
                    with open(filepath, 'r') as download_file:
                        content = [str(line) for line in download_file.readlines()]
                        content = ''.join(content)
                    response = HttpResponse(content, content_type="text/plain")
                    response['Contnt-Disposition'] = "attachment; filename={}".format(filename)
                else:
                    raise FileNotFoundError

            elif project_type == 'remote_searches':
                filepath = '/blast/reciprocal_blast/' + ONE_WAY_BLAST_PROJECT_DIR + project_type + '/' + str(
                    project_id) + '/' + filename
                if isfile(filepath):
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
        return failure_view(request, e)


'''ajax_call_to_snakemake_logfiles

    This function sends task progress data based on available logfiles to the template.
    If the task has not written any logfile yet, the progress will be set to 0.
    
    :param project_id
        :type int
    :param remote - 0 = OneWayBlast 1 = OneWayRemoteBlast
        :type int
        
'''
def ajax_call_to_snakemake_logfiles(request, project_id: int, remote: int):
    try:
        is_ajax = request.headers.get('X-Requested-With') == 'XMLHttpRequest'
        if is_ajax:
            if request.method == "GET":
                if remote == 0:
                    blast_project = get_one_way_project_by_id(project_id)
                elif remote == 1:
                    blast_project = get_one_way_remote_project_by_id(project_id)
                else:
                    raise Exception("There is no one way blast project with this id")

            if remote == 1:
                if blast_project.r_project_execution_task_result:
                    if blast_project.r_project_execution_task_result.status == 'SUCCESS':
                        progress = 100
                    else:
                        progress = read_snakemake_logfile(project_id, remote)


                else:
                    progress = 0
                if blast_project.r_project_execution_task_result:
                    if progress == 100 and blast_project.r_project_execution_task_result.status == 'PROGRESS':
                        progress = 5

            else:
                if blast_project.project_execution_task_result:
                    if blast_project.project_execution_task_result.status == 'SUCCESS':
                        progress = 100
                    else:
                        progress = read_snakemake_logfile(project_id, remote)
                else:
                    progress = 0
                if blast_project.project_execution_task_result:
                    if progress == 100 and blast_project.project_execution_task_result.status == 'PROGRESS':
                        progress = 5


            return JsonResponse({"progress": progress}, status=200)
        return JsonResponse({"progress": "NOT OK"}, status=400)
    except Exception as e:
        return JsonResponse({"error": "{}".format(e)}, status=400)
