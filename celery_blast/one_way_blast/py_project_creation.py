from .py_django_db_services import create_one_way_project_from_form,create_blast_settings_from_form
from django.db import IntegrityError

def create_one_way_blast_project(user, query_file_name, project_form, settings_form):
    try:

        settings = create_blast_settings_from_form(settings_form)
        blastproject = create_one_way_project_from_form(project_form, user, settings, query_file_name)
        return blastproject
    except Exception as e:
        raise IntegrityError('couldnt create project with exception : {}'.format(e))