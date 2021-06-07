from .py_django_db_services import create_one_way_project_from_form,create_blast_settings_from_form, create_one_way_remote_project_from_form
from django.db import IntegrityError, transaction
from blast_project.py_services import upload_file

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

'''
taxid_file = valid_blastdatabase_form.cleaned_data['taxid_file']
taxid_file_path = 'media/' + 'taxonomic_node_files/' + taxid_file.name
upload_file(taxid_file, taxid_file_path)
'''
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

            if project_form.cleaned_data['r_taxid_file'] != None:
                taxid_file = request.FILES['r_taxid_file']
                taxid_file_path = 'media/' + 'taxonomic_node_files/' + taxid_file.name
                upload_file(taxid_file, taxid_file_path)
                blast_project.r_attached_taxonomic_node_file = taxid_file_path
                blast_project.save()
            elif project_form.cleaned_data['r_taxid_uploaded_file'] != '':
                taxid_file = project_form.cleaned_data['r_taxid_uploaded_file']
                taxid_file_path = 'media/' + 'taxonomic_node_files/' + taxid_file
                blast_project.r_attached_taxonomic_node_file = taxid_file_path
                blast_project.save()

            return blast_project
    except Exception as e:
        raise IntegrityError('couldnt create one way remote blast project with exception : {}'.format(e))
