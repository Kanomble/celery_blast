from django.db import models

class OneWayBlastProjectManager(models.Manager):
    # functions
    def create_one_way_blast_project(
            self, project_title,
            project_query_sequences,
            project_user,
            project_settings,
            project_database):
        # calling the create method (objects.create) ..
        blast_project = self.create(
            project_title=project_title,
            project_query_sequences=project_query_sequences,
            project_user=project_user,
            project_settings=project_settings,
            project_database=project_database)

        blast_project.initialize_project_directory()
        blast_project.write_snakemake_configuration_file()

        return blast_project

    '''
    Functions returning Query-Sets
    '''
    # returns all executed projects
    def get_executed_projects(self):
        return self.filter(project_execution_task_result__isnull=False)

    # return projects from username
    def get_blast_projects_by_userid(self, userid):
        return self.filter(project_user__id=userid)
