from .models import OneWayBlastProject
from blast_project.models import BlastSettings
from django.db import IntegrityError

#TODO documentation
def create_one_way_project_from_form(valid_project_form,user,settings,query_sequence_filename):
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

#TODO documentation
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