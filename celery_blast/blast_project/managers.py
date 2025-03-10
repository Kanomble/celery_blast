from django.db import models
from celery_blast.settings import BLAST_PROJECT_DIR, BLAST_DATABASE_DIR, REMOTE_BLAST_PROJECT_DIR

''' BlastProjectManager
    
    This manager class provides functionality for creating BLAST projects. It is a hub for the creation of the 
    BlastProject model and all important necessary files. 
    The create_blast_project function is used in py_django_db_services.create_project_from_form.
    
'''


class BlastProjectManager(models.Manager):
    def create_blast_project(
            self, project_title,
            search_strategy,
            project_query_sequences,
            project_user,
            project_forward_settings, project_backward_settings,
            project_settings,
            project_forward_database, project_backward_database,
            species_name_for_backward_blast,
            filepath=BLAST_PROJECT_DIR):
        # overwriting the create method
        blast_project = self.create(
            project_title=project_title, search_strategy=search_strategy,
            project_query_sequences=project_query_sequences,
            project_user=project_user,
            project_forward_settings=project_forward_settings, project_backward_settings=project_backward_settings,
            project_forward_database=project_forward_database,
            project_backward_database=project_backward_database,
            species_name_for_backward_blast=species_name_for_backward_blast)

        blast_project.initialize_project_directory(filepath=filepath)
        blast_project.write_snakemake_configuration_file(project_settings, filepath=filepath)

        return blast_project

    '''
    Functions returning Query-Sets
    '''

    # returns all executed projects
    def get_executed_projects(self):
        return self.filter(project_execution_snakemake_task__isnull=False)

    # return projects from username
    def get_blast_projects_by_userid(self, userid):
        return self.filter(project_user__id=userid)


''' RemoteBlastProjectManager

'''
class RemoteBlastProjectManager(models.Manager):
    def create_blast_project(
            self, r_project_title,
            r_search_strategy,
            r_project_query_sequences,
            r_project_user,
            r_project_forward_settings, r_project_backward_settings,
            r_project_settings,
            r_project_forward_database, r_project_backward_database,
            r_species_name_for_backward_blast,
            filepath=REMOTE_BLAST_PROJECT_DIR):
        # overwriting the create method
        blast_project = self.create(
            r_project_title=r_project_title, r_search_strategy=r_search_strategy,
            r_project_query_sequences=r_project_query_sequences,
            r_project_user=r_project_user,
            r_project_forward_settings=r_project_forward_settings, r_project_backward_settings=r_project_backward_settings,
            r_project_forward_database=r_project_forward_database,
            r_project_backward_database=r_project_backward_database,
            r_species_name_for_backward_blast=r_species_name_for_backward_blast)

        blast_project.initialize_project_directory(filepath=filepath)
        blast_project.write_snakemake_configuration_file(r_project_settings, filepath=filepath)

        return blast_project

    '''
    Functions returning Query-Sets
    '''

    # returns all executed projects
    def get_executed_projects(self):
        return self.filter(r_project_execution_snakemake_task__isnull=False)

    # return projects from username
    def get_blast_projects_by_userid(self, userid):
        return self.filter(r_project_user__id=userid)
