from .models import BlastProject, BlastDatabase, AssemblyLevels, BlastSettings
from django_celery_results.models import TaskResult
from django.db import IntegrityError
from celery.result import AsyncResult

#following functions are utilized in the dashboard_view
''' get_users_blast_projects
    
    returns a query-set of blast_projects from the logged in user
    
    :param userid (request.user.id)
        :type int
    :returns django.db.models.query.QuerySet of BlastProjects
'''
def get_users_blast_projects(userid):
    return BlastProject.objects.get_blast_projects_by_userid(userid)

''' get_all_blast_databases

    :returns all blastdatabases as a query-set
'''
def get_all_blast_databases():
    return BlastDatabase.objects.all()



#TODO documentation
def get_project_by_id(project_id):
    return BlastProject.objects.get(id=project_id)

#TODO documentation
def create_project_from_form(valid_project_form,user,fw_settings,bw_settings,query_sequence_filename):
    try:
        blast_project = BlastProject.objects.create_blast_project(
            project_title=valid_project_form.cleaned_data['project_title'],
            search_strategy='blastp',
            project_query_sequences=query_sequence_filename,
            project_user=user,
            project_forward_settings=fw_settings,
            project_backward_settings=bw_settings,
            project_database=valid_project_form.cleaned_data['project_database'],
            species_name_for_backward_blast=valid_project_form.cleaned_data['species_name_for_backward_blast']
        )
        return blast_project
    except Exception as e:
        raise IntegrityError('couldnt create blast project with exception : {}'.format(e))
#TODO documentation
def create_blast_settings_from_form(fwOrBw,valid_settings_form):
    try:
        if(fwOrBw == 'fw'):
            blast_settings = BlastSettings.objects.create(e_value=valid_settings_form.cleaned_data['fw_e_value'],
                                            word_size=valid_settings_form.cleaned_data['fw_word_size'],
                                            num_alignments=valid_settings_form.cleaned_data['fw_num_alignments'],
                                            max_target_seqs=valid_settings_form.cleaned_data['fw_max_target_seqs'],
                                            num_threads=valid_settings_form.cleaned_data['fw_num_threads'],
                                            max_hsps=valid_settings_form.cleaned_data['fw_max_hsps'])
        elif(fwOrBw == 'bw'):
            blast_settings = BlastSettings.objects.create(e_value=valid_settings_form.cleaned_data['bw_e_value'],
                                            word_size=valid_settings_form.cleaned_data['bw_word_size'],
                                            num_alignments=valid_settings_form.cleaned_data['bw_num_alignments'],
                                            max_target_seqs=valid_settings_form.cleaned_data['bw_max_target_seqs'],
                                            num_threads=valid_settings_form.cleaned_data['bw_num_threads'],
                                            max_hsps=valid_settings_form.cleaned_data['bw_max_hsps'])
        else:
            raise IntegrityError('fwOrBw is wrong ...')

        return blast_settings
    except Exception as e:
        raise IntegrityError('something went wrong during creation of blastsettings with Exception : {}'.format(e))

#TODO documentation
def update_blast_project_with_task_result_model(project_id,task_id):
    try:
        blast_project = BlastProject.objects.get(id=project_id)
        taskresult = TaskResult.objects.get(task_id=task_id)
        blast_project.project_execution_snakemake_task = taskresult
        blast_project.save()
    except Exception as e:
        raise IntegrityError('problem during updating of blastproject model with task result instance exception : {}'.format(e))

#TODO documentation
def update_blast_database_with_task_result_model(database_id,task_id):
    try:
        blastdb = BlastDatabase.objects.get(id=database_id)
        taskresult = TaskResult.objects.get(task_id=task_id)
        blastdb.database_download_and_format_task = taskresult
        blastdb.save()
    except Exception as e:
        raise IntegrityError('problem during updating of database model with task result instance exception : {}'.format(e))

#TODO documentation
def get_database_by_id(database_id):
    try:
        blastdb=BlastDatabase.objects.get(id=database_id)
        return blastdb
    except Exception as e:
        raise IntegrityError('there is no database with this {} id : {}'.format(database_id,e))

#TODO documentation
def get_all_succeeded_databases():
    return BlastDatabase.objects.get_databases_with_succeeded_tasks()

#TODO documentation
def create_and_save_refseq_database_model(database_name,database_description,assembly_levels,assembly_entries,attached_taxonomic_file=None):
    try:

        #create model refseq genome objects (s. models.py file)
        #path_to_database_file = 'media/' + 'databases/' + 'refseq_databases/' + database_description.replace(' ','_').upper() + '.database.faa'
        if attached_taxonomic_file != None:
            blast_database = BlastDatabase.objects.create(
                database_name=database_name,
                database_description=database_description,
                assembly_entries=assembly_entries,
                attached_taxonomic_node_file=attached_taxonomic_file)
        else:
            blast_database = BlastDatabase.objects.create(
                database_name=database_name,
                database_description=database_description,
                assembly_entries=assembly_entries)


        #get all associated assembly levels (max. 4)
        assembly_levels_models = AssemblyLevels.objects.filter(assembly_level__in=assembly_levels)

        for assembly_level in assembly_levels_models:
            blast_database.assembly_levels.add(assembly_level)

        blast_database.path_to_database_file = 'media/databases/' + str(blast_database.id)
        blast_database.save()
        return blast_database
    except Exception as e:
        raise IntegrityError('couldnt save refseq genome model into database with exception : {}'.format(e))