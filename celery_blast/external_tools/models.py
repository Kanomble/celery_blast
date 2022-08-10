from django.db import models
from blast_project.models import BlastProject
from django_celery_results.models import TaskResult
import pandas as pd
from .managers import ExternalToolsManager, QuerySequenceManager, EntrezSearchManager
from django.contrib.auth.models import User

#TODO documentation - explain why ExternalTools model is usefull (ManyToOne Relationship)
class ExternalTools(models.Model):
    associated_project = models.OneToOneField(
        BlastProject,
        on_delete=models.CASCADE,
        verbose_name="Associated BlastProject"
    )

    objects = ExternalToolsManager()

    def initialize_external_tools_project(self):
        try:
            blast_project = self.associated_project
            query_sequence_id_list = blast_project.get_list_of_query_sequences()
            for qseqid in query_sequence_id_list:
                    QuerySequences.objects.create_query_sequence(qseqid,external_tools=self)

        except Exception as e:
            raise Exception("[-] couldnt extract query sequence ids from associated project with exception : {}".format(e))

    def update_query_sequences_msa_task(self, query_sequence_id, msa_task_id):
        try:
            if self.query_sequences.filter(query_accession_id=query_sequence_id).exists() == True:
                query_sequence = self.query_sequences.get(query_accession_id=query_sequence_id)
                taskresult = TaskResult.objects.get(task_id=msa_task_id)
                query_sequence.multiple_sequence_alignment_task = taskresult
                query_sequence.save()
            else:
                raise Exception("[-] couldnt update query sequence with multiple sequence alignment taskresult object")
        except Exception as e:
            raise Exception("[-] couldnt update query sequence object with exceptipon : {}".format(e))

    def update_query_sequences_phylo_task(self,query_sequence_id,phylo_task_id):
        try:
            if self.query_sequences.filter(query_accession_id=query_sequence_id).exists() == True:
                query_sequence = self.query_sequences.get(query_accession_id=query_sequence_id)
                taskresult = TaskResult.objects.get(task_id=phylo_task_id)
                query_sequence.phylogenetic_tree_construction_task = taskresult
                query_sequence.save()
            else:
                raise Exception("[-] couldnt update query sequence with phylo taskresult object")
        except Exception as e:
            raise Exception("[-] couldnt update query sequence object with exceptipon : {}".format(e))


    def update_for_all_query_sequences_msa_task(self, msa_task_id):
        try:
            query_sequences = self.query_sequences.get_queryset()
            for qseq in query_sequences:
                qseq.update_multiple_sequence_alignment_task(msa_task_id)
        except Exception as e:
            raise Exception("[-] couldnt update query sequences with taskresult object by performing msa for all queries with exception : {}".format(e))

    def update_for_all_query_sequences_phylo_task(self, phylo_task_id):
        try:
            query_sequences = self.query_sequences.get_queryset()
            for qseq in query_sequences:
                qseq.update_phylogenetic_tree_task(phylo_task_id)
        except Exception as e:
            raise Exception(
                "[-] couldnt update query sequences with taskresult object by performing msa for all queries with exception : {}".format(
                    e))

    #TODO refactoring!
    def check_if_msa_task_is_completed(self,query_sequence_id):
        try:
            if self.query_sequences.filter(query_accession_id=query_sequence_id).exists() == True:
                query_sequence = self.query_sequences.get(query_accession_id=query_sequence_id)
                if query_sequence.multiple_sequence_alignment_task:
                    if query_sequence.multiple_sequence_alignment_task.status == 'SUCCESS':
                        return True
                    elif query_sequence.multiple_sequence_alignment_task.status == 'FAILURE':
                        return False
                    else:
                        raise Exception("[-] couldnt check msa taskresult status of query sequence : {}".format(query_sequence_id))
                else:
                    return False
            else:
                raise Exception("[-] couldnt check msa taskresult status of query sequence : {}".format(query_sequence_id))
        except Exception as e:
            raise Exception("[-] couldnt check msa taskresult status for query sequence object with exceptipon : {}".format(e))

#TODO documentation
'''
Query sequences of reciprocal BLAST projects. 
This model combines the results of the RecBLAST for each query sequence to
multiple sequence alignments and phylogenetic tree task result objects. 

It can be used as a hub for new tasks.
'''
class QuerySequences(models.Model):
    query_accession_id = models.CharField(
        max_length=200,
        blank=False,unique=False,
        verbose_name="query sequence identifier"
    )
    multiple_sequence_alignment_task = models.ForeignKey(
        TaskResult,
        on_delete=models.CASCADE,
        blank=True,null=True,
        verbose_name="celery task for multiple sequence alignment performed by mafft in the bioinformatic tools container",
        related_name="msa_task",
        unique=False
    )
    phylogenetic_tree_construction_task = models.ForeignKey(
        TaskResult,
        on_delete=models.CASCADE,
        blank=True,null=True,
        verbose_name="celery task for constructing a phylogenetic tree performed by fasttree in the bioinformatic tools container",
        related_name="tree_task",
        unique=False
    )
    #many to one relationship
    external_tool_for_query_sequence = models.ForeignKey(
        ExternalTools,
        verbose_name="query sequence for the external tools model",
        related_name="query_sequences",
        on_delete=models.CASCADE
    )

    objects = QuerySequenceManager()

    def update_multiple_sequence_alignment_task(self, msa_task_id):
        try:
            taskresult = TaskResult.objects.get(task_id=msa_task_id)
            self.multiple_sequence_alignment_task = taskresult
            self.save()
        except Exception as e:
            raise Exception("[-] couldnt update query sequences with taskresult object for msa with exception : {}".format(e))

    def update_phylogenetic_tree_task(self, phylo_task_id):
        try:
            taskresult = TaskResult.objects.get(task_id=phylo_task_id)
            self.phylogenetic_tree_construction_task = taskresult
            self.save()
        except Exception as e:
            raise Exception(
                "[-] couldnt update query sequences with taskresult object for msa with exception : {}".format(e))

class EntrezSearch(models.Model):

    database = models.CharField(
        max_length=200, unique=False,
        default="pubmed", verbose_name="Search Database"
    )

    entrez_user = models.ForeignKey(
        User,
        on_delete=models.CASCADE,
        verbose_name="user who created this project",
        default=1)

    entrez_query = models.CharField(
        max_length=600,
        verbose_name="Search query.", default="Lipopolysaccharides AND review [PT]")

    fasta_file_name = models.CharField(max_length=200,
                                       blank=True,
                                       null=True,
                                       verbose_name="search associated fasta file")

    file_name = models.CharField(
        max_length=200, unique=True,
        verbose_name="File name for search results.", default="paper"
    )

    paper_entries = models.IntegerField(
        verbose_name="Amount of paper in result file.",
        default=0
    )

    search_task_result = models.OneToOneField(
        TaskResult,
        on_delete=models.CASCADE,
        verbose_name="TaskResult model for entrez searches",
        related_name="search_task",
        null=True
    )

    download_task_result = models.OneToOneField(
        TaskResult,
        on_delete=models.CASCADE,
        verbose_name="TaskResult model for downloads",
        related_name="download_task",
        null=True
    )

    timestamp = models.DateTimeField(auto_now=True)

    objects = models.Manager()
    edirect_objects = EntrezSearchManager()

    def get_paper_content(self):
        pandas_header = {}
        pandas_header['pubmed'] = ['Id', 'PubDate', 'Source', 'Title', 'ElocationID']
        pandas_header['protein'] = ['Id','Caption','Title','Organism']

        paper = pd.read_table(self.file_name, header=None)
        paper.columns = pandas_header[self.database]

        paper = paper.to_html(classes='entrezsearch" id="searchResultTable')
        return paper

    def get_pandas_table(self):
        pandas_header = {}
        pandas_header['pubmed'] = ['Id', 'PubDate', 'Source', 'Title', 'ElocationID']
        pandas_header['protein'] = ['Id','Caption','Title','Organism']

        paper = pd.read_table(self.file_name, header=None)
        paper.columns = pandas_header[self.database]
        return paper

    def get_paper_number(self):
        return len(pd.read_table(self.file_name, header=None))

    def update_paper_entries(self):
        paper_entries = len(pd.read_table(self.file_name, header=None))
        self.paper_entries=paper_entries
        self.save()
        return paper_entries
