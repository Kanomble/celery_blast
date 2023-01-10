from .py_django_db_services import create_blast_settings_from_form, create_project_from_form
from django.db import IntegrityError

'''create_blast_project

wrapper function for utilization of functions that create new database entries for reciprocal blast projects.

    :param query_file_name
        :type string
    :param user
        :type django user model
    :param project_form
        :type django form (BlastProject)
    :param fw_settings_form
        :type django form (ForwardBlastSettings)
    :param bw_settings_form
        :type django form (BackwardBlastSettings)
    :return blastproject
        :type django model (BlastProject)

'''
def create_blast_project(query_file_name,user,project_form,fw_settings_form,bw_settings_form,filepath='media/blast_projects/'):
    try:
        fw_settings = create_blast_settings_from_form('fw',fw_settings_form)
        bw_settings = create_blast_settings_from_form('bw',bw_settings_form)
        blastproject = create_project_from_form(project_form,user,fw_settings,bw_settings,query_file_name,filepath=filepath)
        return blastproject
    except Exception as e:
        raise IntegrityError('couldnt create project with exception : {}'.format(e))