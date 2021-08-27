from celery import shared_task
from celery.utils.log import get_task_logger
import requests
#logger for celery worker instances
logger = get_task_logger(__name__)

@shared_task(bind=True)
def execute_multiple_sequence_alignment(self, project_id, folder_path):
    try:
        url = "http://192.168.1.131:5001/perform_simple_msa/"+str(project_id)+'/'+folder_path
        returncode = requests.get(url)
        if returncode != 0:
            return 1
        return 0
    except Exception as e:
        raise Exception("couldnt perform simple msa task with exception : {}".format(e))