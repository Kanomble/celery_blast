from .py_django_db_services import create_blast_settings_from_form, create_project_from_form
from django.db import IntegrityError, transaction


def create_blast_project(query_file_name,user,project_form,fw_settings_form,bw_settings_form):
    try:

        fw_settings = create_blast_settings_from_form('fw',fw_settings_form)
        bw_settings = create_blast_settings_from_form('bw',bw_settings_form)
        blastproject = create_project_from_form(project_form,user,fw_settings,bw_settings,query_file_name)
        return blastproject
    except Exception as e:
        raise IntegrityError('couldnt create project with exception : {}'.format(e))