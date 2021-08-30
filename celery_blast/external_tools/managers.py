from django.db import models, IntegrityError
from blast_project.models import BlastProject

class ExternalToolsManager(models.Manager):
    def create_external_tools(self,project_id):
        try:
            if self.filter(associated_project_id=project_id).exists() == False:
                blast_project = BlastProject.objects.get(id=project_id)
                external_tools = self.create(
                    associated_project=blast_project)

                external_tools.initialize_external_tools_project()
            else:
                external_tools = self.get(associated_project_id=project_id)

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