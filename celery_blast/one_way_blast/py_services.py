from .models import OneWayBlastProject
from os.path import isdir
from shutil import rmtree
from django.db import IntegrityError, transaction

def delete_one_way_blast_project_and_associated_directories(one_way_project_id):
    try:
        with transaction.atomic():
            project = OneWayBlastProject.objects.get(id=one_way_project_id)
            if isdir('media/blast_projects/' + str(one_way_project_id)):
                rmtree('media/blast_projects/' + str(one_way_project_id))
            if isdir('static/images/result_images/' + str(one_way_project_id)):
                rmtree('static/images/result_images/' + str(one_way_project_id))
            project.delete()
    except Exception as e:
        raise IntegrityError("couldnt delete blast project entry : {}".format(e))

