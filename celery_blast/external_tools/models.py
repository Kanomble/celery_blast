from django.db import models
from blast_project.models import BlastProject
from django_celery_results.models import TaskResult

from .managers import ExternalToolsManager, QuerySequenceManager
# Create your models here.

class ExternalTools(models.Model):
    associated_project = models.OneToOneField(
        BlastProject,
        on_delete=models.CASCADE,
        verbose_name="Associated BlastProject"
    )

    objects = ExternalToolsManager()
    '''
    #many to many seems to be wrong we should consider to use a one to many relationship
    query_sequences = models.ManyToManyField(
        to=QuerySequences,
        verbose_name="query sequences of the associated BlastProject"
    )
    '''

class QuerySequences(models.Model):
    query_accession_id = models.CharField(
        max_length=200,
        blank=False,unique=False,
        verbose_name="query sequence identifier"
    )
    multiple_sequence_alignment_task = models.OneToOneField(
        TaskResult,
        on_delete=models.CASCADE,
        blank=True,null=True,
        verbose_name="celery task for multiple sequence alignment performed by mafft in the bioinformatic tools container",
        related_name="msa_task"
    )
    phylogenetic_tree_construction_task = models.OneToOneField(
        TaskResult,
        on_delete=models.CASCADE,
        blank=True,null=True,
        verbose_name="celery task for constructing a phylogenetic tree performed by fasttree in the bioinformatic tools container",
        related_name="tree_task"
    )
    #many to one relationship
    external_tool_for_query_sequence = models.ForeignKey(
        ExternalTools,
        verbose_name="query sequence for the external tools model",
        related_name="query_sequences",
        on_delete=models.CASCADE
    )

    objects = QuerySequenceManager()
