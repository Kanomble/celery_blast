from .py_django_db_services import create_one_way_project_from_form,create_blast_settings_from_form, \
    create_one_way_remote_project_from_form, get_one_way_remote_project_by_id
from django.db import IntegrityError, transaction
from blast_project.py_services import upload_file
from .models import OneWayRemoteBlastProject

#TODO Documentation
def create_one_way_blast_project(user, query_file_name,query_file, project_form, settings_form):
    try:
        with transaction.atomic():
            settings = create_blast_settings_from_form(settings_form)
            blast_project = create_one_way_project_from_form(project_form, user, settings, query_file_name)

            path_to_query_file = blast_project.get_project_dir() + '/' + query_file_name

            upload_file(query_file, path_to_query_file)
            return blast_project
    except Exception as e:
        raise IntegrityError('couldnt create one way blast project with exception : {}'.format(e))


def create_one_way_remote_blast_project(user, project_form, settings_form, request):
    try:
        with transaction.atomic():
            query_sequences = request.FILES['r_query_sequence_file']
            settings = create_blast_settings_from_form(settings_form)
            #valid_project_form,user,settings,query_sequence_filename
            blast_project = create_one_way_remote_project_from_form(
                project_form,
                user,
                settings,
                query_sequences.name)

            path_to_query_file = blast_project.get_project_dir() + '/' + query_sequences.name

            upload_file(query_sequences, path_to_query_file)

            if project_form.cleaned_data['r_entrez_query'] != None:
                entrez_query = project_form.cleaned_data['r_entrez_query']
                update_snakemake_configuration_with_entrez_query(blast_project.id,entrez_query)
                blast_project.r_entrez_query = entrez_query
                blast_project.save()


            return blast_project
    except Exception as e:
        raise IntegrityError('couldnt create one way remote blast project with exception : {}'.format(e))

def update_snakemake_configuration_with_entrez_query(project_id,entrez_query):
    try:
        project = get_one_way_remote_project_by_id(project_id)
        project_dir = project.get_project_dir()
        with open(project_dir+'/snakefile_config','r') as config:
            lines = config.readlines()

        with open(project_dir+"/snakefile_config","w") as config:
             for line in lines:
                 if "entrez_query" in str(line):
                     line = "entrez_query: \"{}\"\n".format(entrez_query)
                 config.write(line)

    except Exception as e:
        raise IntegrityError('Couldnt alter entrez query field in snakemake configuration file with Exception : {}'.format(e))
