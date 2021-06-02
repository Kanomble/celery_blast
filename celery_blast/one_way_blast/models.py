from django.db import models
from django.contrib.auth.models import User
from django_celery_results.models import TaskResult
from blast_project.models import BlastSettings, BlastDatabase
# Create your models here.

class OneWayBlastProject(models.Model):
    project_title = models.CharField(
        max_length=200, blank=False, unique=True,
        verbose_name="title of this one way blast project"
    )

    project_query_sequences = models.CharField(
        max_length=200,
        verbose_name="query sequence filepath"
    )

    project_user = models.ForeignKey(
        User,
        on_delete=models.CASCADE,
        verbose_name="user who created this project"
    )

    project_settings = models.OneToOneField(
        BlastSettings,
        on_delete=models.CASCADE,
        verbose_name="BLAST settings"
    )

    project_database = models.ForeignKey(
        BlastDatabase,
        on_delete=models.CASCADE,
        verbose_name="associated BLAST database"
    )

    project_execution_task_result = models.OneToOneField(
        TaskResult,
        on_delete=models.SET_NULL,
        blank=True, null=True,
        verbose_name="django_celery_results taskresult model for this project"
    )
