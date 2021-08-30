from django.db import models, IntegrityError
from .models import QuerySequences
from blast_project.py_django_db_services import get_project_by_id, get_list_of_query_sequences

class ExternalToolsManager(models.Manager):
    def create_external_tools(self,project_id):
        try:
            blast_project = get_project_by_id(project_id)
            qseqids = get_list_of_query_sequences(blast_project)
            external_tools = self.create(
                associated_blast_project=blast_project)
            for query_sequence in qseqids:
                QuerySequences.objects.create(
                    query_accession_id=query_sequence,
                    external_tool_for_query_sequence=external_tools)
            return external_tools
        except Exception as e:
            raise IntegrityError("[-] couldnt save external tools model into database with exception : {}".format(e))