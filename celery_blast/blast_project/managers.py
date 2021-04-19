from django.db import models

# example for model manager that is invoked as __init__ for the blastproject model
# allows customization of queries for the database
''' BlastProjectManager
    
'''
class BlastProjectManager(models.Manager):
    # functions
    def create_blast_project(
            self, project_title,
            search_strategy,
            project_query_sequences,
            timestamp,
            project_user,
            project_forward_settings, project_backward_settings,
            project_forward_database, project_backward_database):
        # calling the create method (objects.create) ..
        blast_project = self.create(
            project_title=project_title, search_strategy=search_strategy,
            project_query_sequences=project_query_sequences,
            timestamp=timestamp, project_user=project_user,
            project_forward_settings=project_forward_settings, project_backward_settings=project_backward_settings,
            project_forward_database=project_forward_database, project_backward_database=project_backward_database)

        blast_project.initialize_project_directory()

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

class BlastDatabaseManager(models.Manager):
    '''
    Functions returning Query-Sets
    '''
    def get_databases_with_executed_tasks(self):
        return self.filter(database_download_and_format_task__isnull=False)
