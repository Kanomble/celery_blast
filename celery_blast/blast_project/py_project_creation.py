from .py_django_db_services import create_blast_settings_from_form,\
    create_project_from_form, create_remote_project_from_form
from django.db import IntegrityError
from celery_blast.settings import BLAST_PROJECT_DIR, REMOTE_BLAST_PROJECT_DIR

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
    :param symblast_settings_form
        :type django form (SymBLASTProjectSettings)
    :return blastproject
        :type django model (BlastProject)

'''


def create_blast_project(query_file_name, user, project_form,
                         fw_settings_form, bw_settings_form,
                         symblast_settings_form,
                         filepath=BLAST_PROJECT_DIR):
    try:
        fw_settings = create_blast_settings_from_form('fw', fw_settings_form)
        bw_settings = create_blast_settings_from_form('bw', bw_settings_form)
        blastproject = create_project_from_form(project_form, user, fw_settings, bw_settings, query_file_name,
                                                symblast_settings_form,
                                                filepath=filepath)
        return blastproject
    except Exception as e:
        raise IntegrityError('couldnt create project with exception : {}'.format(e))


def create_remote_blast_project(query_file_name, user, project_form,
                         fw_settings_form, bw_settings_form,
                         symblast_settings_form,
                         filepath=REMOTE_BLAST_PROJECT_DIR):
    try:
        fw_settings = create_blast_settings_from_form('fw', fw_settings_form)
        bw_settings = create_blast_settings_from_form('bw', bw_settings_form)
        blastproject = create_remote_project_from_form(project_form, user, fw_settings, bw_settings, query_file_name,
                                                symblast_settings_form,
                                               filepath=filepath)

        if project_form.cleaned_data['r_entrez_query'] != None:
            entrez_query = project_form.cleaned_data['r_entrez_query']
            update_snakemake_remote_configuration_with_entrez_query(blastproject.id, entrez_query)
            blastproject.r_entrez_query = entrez_query
            blastproject.save()

        return blastproject
    except Exception as e:
        raise IntegrityError('couldnt create project with exception : {}'.format(e))


def update_snakemake_remote_configuration_with_entrez_query(project_id, entrez_query):
    try:
        project_dir = REMOTE_BLAST_PROJECT_DIR + str(project_id)
        with open(project_dir + '/snakefile_config', 'r') as config:
            lines = config.readlines()

        with open(project_dir + "/snakefile_config", "w") as config:
            for line in lines:
                if "entrez_query" in str(line):
                    line = "entrez_query: \"{}\"\n".format(entrez_query)
                config.write(line)

    except Exception as e:
        raise IntegrityError(
            'Couldnt alter entrez query field in snakemake configuration file with Exception : {}'.format(e))