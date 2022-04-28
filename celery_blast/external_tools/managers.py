from django.db import models, IntegrityError
from blast_project.models import BlastProject
from django_celery_results.models import TaskResult
import os
import pandas as pd

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

    def get_external_tools_based_on_project_id(self, project_id):
        try:
            if self.filter(associated_project_id=project_id).exists() == False:
                raise IntegrityError("[-] there is no external tools object with your specified project id : {}".format(project_id))
            else:
                return self.get(associated_project_id=project_id)
        except Exception as e:
            raise IntegrityError(
                "[-] there is no external tools object with your specified project id : {}".format(project_id))


class QuerySequenceManager(models.Manager):
    def create_query_sequence(self,query_sequence_id,external_tools):
        try:
            query_sequence = self.create(
                query_accession_id=query_sequence_id,
                external_tool_for_query_sequence=external_tools)
            return query_sequence
        except Exception as e:
            raise IntegrityError("[-] couldnt save query sequence model into database with exception : {}".format(e))

class EntrezSearchManager(models.Manager):
    def create_entrez_search(self, database, entrez_query, file_path, task_result_id, entrez_user):
        file_name = file_path
        if os.path.isfile(file_name) == False:
            paper_entries = 0
        else:
            paper_entries = len(pd.read_table(file_name, header=None))
        task_result = TaskResult.objects.get(task_id=task_result_id)
        edirect_paper = self.create(database=database,
                                    entrez_query=entrez_query,
                                    file_name=file_name,
                                    paper_entries=paper_entries,
                                    task_result=task_result,
                                    entrez_user=entrez_user)
        return edirect_paper

    def get_entrezsearch_on_query(self, search_query, database):
        return self.filter(search_query=search_query,database=database)

    def get_all_entrez_searches_from_current_user(self,user_id):
        return self.filter(entrez_user__id=user_id)

