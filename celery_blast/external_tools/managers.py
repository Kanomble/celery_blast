from django.db import models, IntegrityError
from blast_project.models import BlastProject

class ExternalToolsManager(models.Manager):
    def create_external_tools(self,project_id):
        try:
            blast_project = BlastProject.objects.get(id=project_id)
            external_tools = self.create(
                associated_blast_project=blast_project)
            return external_tools
        except Exception as e:
            raise IntegrityError("[-] couldnt save external tools model into database with exception : {}".format(e))

class QuerySequenceManager(models.Manager):
    def create_query_sequence(self,query_sequence_id,external_tools):
        try:
            query_sequence = self.create(
                query_accession_id=query_sequence_id,
                external_tool_for_query_sequence=external_tools)
            return query_sequence
        except Exception as e:
            raise IntegrityError("[-] couldnt save query sequence model into database with exception : {}".format(e))

    def create_query_sequences_based_on_qseqlist(self,query_sequence_list,external_tools):
        try:
            query_sequences = []
            for qseqid in query_sequence_list:
                query_sequence = self.create(query_accession_id=qseqid,external_tool_for_query_sequence=external_tools)
                query_sequences.append(query_sequence)
            return query_sequences
        except Exception as e:
            raise IntegrityError("[-] couldnt save query sequence model into database with exception : {}".format(e))