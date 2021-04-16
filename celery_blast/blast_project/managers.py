from django.db import models

# example for model manager that is invoked as __init__ for the blastproject model
# allows customization of queries for the database
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

    # returns all executed projects
    def get_executed_projects(self):
        return self.filter(project_execution_snakemake_task__isnull=False)