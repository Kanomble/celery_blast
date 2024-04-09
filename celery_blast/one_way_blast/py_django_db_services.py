from blast_project.models import BlastSettings
from django.db import IntegrityError
from django_celery_results.models import TaskResult

from .models import OneWayBlastProject, OneWayRemoteBlastProject


# TODO documentation
def get_one_way_project_by_id(project_id):
    return OneWayBlastProject.objects.get(id=project_id)


def get_one_way_remote_project_by_id(project_id):
    return OneWayRemoteBlastProject.objects.get(id=project_id)


def get_users_one_way_blast_projects(userid):
    return OneWayBlastProject.objects.get_blast_projects_by_userid(userid)


def get_users_one_way_remote_blast_projects(userid):
    return OneWayRemoteBlastProject.objects.get_blast_projects_by_userid(userid)

def check_if_one_way_project_title_exists(local_or_remote:str,new_project_title:str):
    try:
        if local_or_remote == "local":
            local_projects = OneWayBlastProject.objects.all()
            for project in local_projects:
                if project.project_title == new_project_title:
                    return True
        elif local_or_remote == "remote":
            remote_projects = OneWayRemoteBlastProject.objects.all()
            for project in remote_projects:
                if project.r_project_title == new_project_title:
                    return True
        else:
            raise Exception("[-] Project is neither remote nor local ...")
        return False
    except Exception as e:
        raise Exception("ERROR during project title form validation with excpetion: {}".format(e))


def create_one_way_remote_project_from_form(valid_project_form, user, settings, query_sequence_filename):
    try:
        blast_project = OneWayRemoteBlastProject.objects.create_one_way_remote_blast_project(
            r_project_title=valid_project_form.cleaned_data['r_project_title'],
            r_project_query_sequences=query_sequence_filename,
            r_project_user=user,
            r_project_settings=settings,
            r_project_database=valid_project_form.cleaned_data['r_project_database'],
            r_project_search_strategy=valid_project_form.cleaned_data['r_search_strategy'],
        )
        return blast_project
    except Exception as e:
        raise IntegrityError('couldnt create blast project with exception : {}'.format(e))


# TODO documentation
def create_one_way_project_from_form(valid_project_form, user, settings, query_sequence_filename):
    try:
        blast_project = OneWayBlastProject.objects.create_one_way_blast_project(
            project_title=valid_project_form.cleaned_data['project_title'],
            project_query_sequences=query_sequence_filename,
            project_user=user,
            project_settings=settings,
            project_database=valid_project_form.cleaned_data['project_database'],
        )
        return blast_project
    except Exception as e:
        raise IntegrityError('couldnt create blast project with exception : {}'.format(e))


# TODO documentation
def create_blast_settings_from_form(valid_settings_form):
    try:
        blast_settings = BlastSettings.objects.create(e_value=valid_settings_form.cleaned_data['e_value'],
                                                      word_size=valid_settings_form.cleaned_data['word_size'],
                                                      num_alignments=valid_settings_form.cleaned_data['num_alignments'],
                                                      num_threads=valid_settings_form.cleaned_data['num_threads'],
                                                      max_hsps=valid_settings_form.cleaned_data['max_hsps'],
                                                      max_target_seqs=10000)
        return blast_settings
    except Exception as e:
        raise IntegrityError('something went wrong during creation of blastsettings with Exception : {}'.format(e))


# TODO documentation
def update_one_way_blast_project_with_task_result_model(project_id, task_id):
    try:
        blast_project = OneWayBlastProject.objects.get(id=project_id)
        taskresult = TaskResult.objects.get(task_id=task_id)
        blast_project.project_execution_task_result = taskresult
        blast_project.save()
    except Exception as e:
        raise IntegrityError(
            'problem during updating of onewayblastproject model with task result instance exception : {}'.format(e))


# TODO documentation
def update_one_way_remote_blast_project_with_task_result_model(project_id, task_id):
    try:
        blast_project = OneWayRemoteBlastProject.objects.get(id=project_id)
        taskresult = TaskResult.objects.get(task_id=task_id)
        blast_project.r_project_execution_task_result = taskresult
        blast_project.save()
    except Exception as e:
        raise IntegrityError(
            'problem during updating of onewayblastproject model with task result instance exception : {}'.format(
                e))
