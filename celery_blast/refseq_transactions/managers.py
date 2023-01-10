from django.db import models
class BlastDatabaseManager(models.Manager):
    '''
    Functions returning Query-Sets
    '''
    def get_databases_with_executed_tasks(self):
        return self.filter(database_download_and_format_task__isnull=False)

    def get_databases_without_executed_tasks(self):
        return self.filter(database_download_and_format_task__isnull=True)

    def get_databases_with_succeeded_tasks(self):
        return self.filter(database_download_and_format_task__status='SUCCESS')

    def get_databases_with_failed_tasks(self):
        return self.filter(database_download_and_format_task__status='FAILURE')

    def get_databases_with_task_on_progress(self):
        return self.filter(database_download_and_format_task__status='PROGRESS')