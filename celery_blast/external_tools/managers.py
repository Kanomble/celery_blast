from os.path import isfile
import pandas as pd
from blast_project.models import BlastProject, RemoteBlastProject
from django.db import models, IntegrityError
from django_celery_results.models import TaskResult
from external_tools import models as mdl

from .py_services import check_if_cdd_search_can_get_executed


class ExternalToolsManager(models.Manager):
    def create_external_tools(self, project_id:int,remote_or_local:str):
        try:

            if remote_or_local == 'local':
                blast_project = BlastProject.objects.get(id=project_id)
                external_tools = self.create(
                    associated_project=blast_project)

                external_tools.initialize_external_tools_project()
            else:
                remote_blast_project = RemoteBlastProject.objects.get(id=project_id)
                external_tools = self.create(
                    remote_associated_project=remote_blast_project)

                external_tools.initialize_external_tools_project()
            return external_tools
        except Exception as e:
            raise IntegrityError("[-] couldnt save external tools model into database with exception : {}".format(e))


    def get_external_tools_based_on_project_id(self, project_id, remote_or_local):
        try:
            if remote_or_local == "local":
                if self.filter(associated_project_id=project_id).exists() == False:
                    raise IntegrityError(
                        "[-] there is no external tools object with your specified project id : {}".format(project_id))
                else:
                    return self.get(associated_project_id=project_id)
            elif remote_or_local == "remote":
                if self.filter(remote_associated_project_id=project_id).exists() == False:
                    raise IntegrityError(
                        "[-] there is no external tools object with your specified project id : {}".format(project_id))
                else:
                    return self.get(remote_associated_project_id=project_id)
            else:
                raise Exception("[-] associated external tools object is not local nor remote ...")
        except Exception as e:
            raise IntegrityError(
                "[-] there is no external tools object with your specified project id : {} with exception: {}".format(
                    project_id, e))

    def get_all_associated_query_sequences(self, project_id, remote_or_local):
        try:
            external_tools = self.get_external_tools_based_on_project_id(project_id, remote_or_local)
            query_sequence_set = mdl.QuerySequences.objects.filter(external_tool_for_query_sequence=external_tools)
            return query_sequence_set
        except Exception as e:
            raise IntegrityError(
                "[-] ERROR fetching associated query sequences for external tools with project id: {} and exception: {}".format(
                    project_id, e))

    '''get_associated_query_sequence
        
        Returns the associated QuerySequence model instance in an QuerySet.
        
        :param project_id
            :type int
        :param query_sequence
            :type str
        
        :returns query_sequence
            :type django.db.models.query.QuerySet
    '''

    def get_associated_query_sequence(self, project_id: int, query_sequence: str, remote_or_local):
        try:
            external_tools = self.get_external_tools_based_on_project_id(project_id, remote_or_local)
            query_sequence = mdl.QuerySequences.objects.filter(external_tool_for_query_sequence=external_tools,
                                                               query_accession_id=query_sequence)
            return query_sequence
        except Exception as e:
            raise IntegrityError(
                "[-] ERROR fetching associated query sequences for external tools with project id: {} and exception: {}".format(
                    project_id, e))

    '''get_associated_query_sequence_and_return_cdd_task

        Returns the CDD search TaskResult object of a QuerySequence model instance.

        :param project_id
            :type int
        :param query_sequence
            :type str

        :returns task_result
            :type TaskResult
    '''

    def get_associated_query_sequence_and_return_cdd_task(self, project_id: int, query_sequence: str, remote_or_local:str):
        try:
            external_tools = self.get_external_tools_based_on_project_id(project_id, remote_or_local)
            query_sequence = mdl.QuerySequences.objects.filter(external_tool_for_query_sequence=external_tools,
                                                               query_accession_id=query_sequence)
            task_id = query_sequence[0].cdd_domain_search_task_id
            task_result = TaskResult.objects.get(task_id=task_id)
            return task_result
        except Exception as e:
            raise IntegrityError(
                "[-] ERROR fetching associated query sequences for external tools with project id: {} and exception: {}".format(
                    project_id, e))

    '''check_cdd_domain_search_tasks

        This function checks the CDD domain search task status of the associated query sequence models.
        It also evaluates the practicability of a CDD search. If it wouldnt make sense to perform a search, the 
        value for the query_sequence is "not valid" if its practicable the value will be "valid".
        
        :param project_id
            :type int

        :returns query_sequence_cdd_tasks
            :type dict[str] = tuple(int,str)
    '''

    def check_cdd_domain_search_task_status(self, project_id: int) -> dict:
        try:
            query_sequence_cdd_tasks = {}
            # query sequence set
            query_sequences = self.get_all_associated_query_sequences(project_id)
            for query_sequence in query_sequences:
                returncode = check_if_cdd_search_can_get_executed(query_sequence.query_accession_id, project_id)
                if returncode == 1:
                    query_sequence_cdd_tasks[query_sequence.query_accession_id] = query_sequence. \
                                                                                      check_if_cdd_search_is_complete(), \
                                                                                  'not valid'
                else:
                    query_sequence_cdd_tasks[query_sequence.query_accession_id] = query_sequence. \
                                                                                      check_if_cdd_search_is_complete(), \
                                                                                  'valid'
            return query_sequence_cdd_tasks
        except Exception as e:
            raise Exception(
                "[-] ERROR fetching CDD domain search task status for project: {} with exception: {}".format(project_id,
                                                                                                             e))

    '''get_cdd_searchable_queries
        
        Returns a list with QuerySequence models for potential CDD serches.
        
        :param project_id
            :type int
        
        :returns query_sequences_rdy_for_cdd
            :type list[QuerySequence]
    '''

    def get_cdd_searchable_queries(self, project_id: int) -> list:
        try:
            query_sequences_rdy_for_cdd = []
            # query sequence set
            query_sequences = self.get_all_associated_query_sequences(project_id)
            for query_sequence in query_sequences:
                returncode = check_if_cdd_search_can_get_executed(query_sequence.query_accession_id, project_id)
                if returncode == 0:
                    if query_sequence.check_if_cdd_search_is_complete() == 'NOTEXEC':
                        query_sequences_rdy_for_cdd.append(query_sequence)
            return query_sequences_rdy_for_cdd
        except Exception as e:
            raise Exception(
                "[-] ERROR during creation of query sequence list with "
                "not executed cdd searches for project {} with exception: {}".format(project_id, e))


class QuerySequenceManager(models.Manager):
    def create_query_sequence(self, query_sequence_id, query_sequence_info, external_tools):
        try:
            query_sequence = self.create(
                query_accession_id=query_sequence_id,
                query_accession_information=query_sequence_info,
                external_tool_for_query_sequence=external_tools)
            return query_sequence
        except Exception as e:
            raise IntegrityError("[-] couldnt save query sequence model into database with exception : {}".format(e))


class EntrezSearchManager(models.Manager):
    def create_entrez_search(self, database, entrez_query, file_path, search_task_result, entrez_user):
        file_name = file_path
        if isfile(file_name) == False:
            paper_entries = 0
        else:
            paper_entries = len(pd.read_table(file_name, header=None))
        task_result = TaskResult.objects.get(task_id=search_task_result)
        edirect_paper = self.create(database=database,
                                    entrez_query=entrez_query,
                                    file_name=file_name,
                                    paper_entries=paper_entries,
                                    search_task_result=task_result,
                                    entrez_user=entrez_user)
        return edirect_paper

    def get_entrezsearch_on_query(self, search_query, database):
        return self.filter(search_query=search_query, database=database)

    def get_all_entrez_searches_from_current_user(self, user_id):
        return self.filter(entrez_user__id=user_id)
