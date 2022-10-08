from django.db import models

# allows customization of queries for the database
''' BlastProjectManager
    
'''
class BlastProjectManager(models.Manager):
    # functions
    def create_blast_project(
            self, project_title,
            search_strategy,
            project_query_sequences,
            project_user,
            project_forward_settings, project_backward_settings,
            project_forward_database, project_backward_database,
            species_name_for_backward_blast,
            filepath='media/blast_projects/'):
        # calling the create method (objects.create) ..
        blast_project = self.create(
            project_title=project_title, search_strategy=search_strategy,
            project_query_sequences=project_query_sequences,
            project_user=project_user,
            project_forward_settings=project_forward_settings, project_backward_settings=project_backward_settings,
            project_forward_database=project_forward_database,
            project_backward_database=project_backward_database,
            species_name_for_backward_blast=species_name_for_backward_blast)

        blast_project.initialize_project_directory(filepath=filepath)
        blast_project.write_snakemake_configuration_file(filepath=filepath)

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

