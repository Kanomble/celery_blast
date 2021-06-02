from .models import OneWayBlastProject
from os.path import isdir
from shutil import rmtree
from django.db import IntegrityError, transaction


#TODO documentation
def delete_one_way_blast_project_and_associated_directories_by_id(project_id):
    try:
        with transaction.atomic():
            project = OneWayBlastProject.objects.get(id=project_id)
            if isdir('media/one_way_blast/' + str(project_id)):
                rmtree('media/one_way_blast/' + str(project_id))
            if isdir('static/images/result_images/one_way_blast/'+str(project_id)):
                rmtree('static/images/result_images/one_way_blast/'+str(project_id))
            project.delete()
    except Exception as e:
        raise IntegrityError("couldnt delete one way blast project entry : {}".format(e))
