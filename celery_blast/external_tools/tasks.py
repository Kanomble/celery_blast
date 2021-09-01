from celery import shared_task
from celery.utils.log import get_task_logger
import requests
from .models import ExternalTools
from celery_progress.backend import ProgressRecorder

#logger for celery worker instances
logger = get_task_logger(__name__)

@shared_task(bind=True)
def execute_multiple_sequence_alignment(self, project_id, query_sequence_id):
    try:
        progress_recorder = ProgressRecorder(self)
        progress_recorder.set_progress(0, 100, "started process")


        external_tools = ExternalTools.objects.get_external_tools_based_on_project_id(project_id)
        external_tools.update_query_sequences_msa_task(query_sequence_id, str(self.request.id))

        url = "http://10.125.46.80:5001/perform_simple_msa/" + str(project_id) +'/' + query_sequence_id
        response = requests.post(url, json={"folder_path":query_sequence_id, "project_id":project_id})
        if response.status_code == 200:
            progress_recorder.set_progress(100, 100, "finished process")
            return 0
        elif response.status_code == 400:
            raise Exception("[-] status code for http request to bioinformatic tools container with url : {} is 400".format(url))
        elif response.status_code == 500:
            raise Exception(
                "[-] status code for http request to bioinformatic tools container with url : {} is 500".format(url))
        else:
            raise Exception("[-] error in response from bioinformatic tools container with url : {} status code is : {}"
                            .format(url,response.status_code))
    except Exception as e:
        raise Exception("[-] couldnt perform simple msa task with exception : {}".format(e))