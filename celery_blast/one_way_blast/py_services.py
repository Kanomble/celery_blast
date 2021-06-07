from .models import OneWayBlastProject, OneWayRemoteBlastProject
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

#TODO documentation
def delete_one_way_remote_blast_project_and_associated_directories_by_id(project_id):
    try:
        with transaction.atomic():
            project = OneWayRemoteBlastProject.objects.get(id=project_id)
            if isdir('media/one_way_blast/remote_searches/' + str(project_id)):
                rmtree('media/one_way_blast/remote_searches/' + str(project_id))
            if isdir('static/images/result_images/one_way_blast/remote_searches/'+str(project_id)):
                rmtree('static/images/result_images/one_way_blast/remote_searches/'+str(project_id))
            project.delete()
    except Exception as e:
        raise IntegrityError("couldnt delete one way blast project entry : {}".format(e))

#TODO documentation
#loads the reciprocal results table that is written with one of the last rules in the snakefiles
def get_one_way_html_results(project_id,filename, remote):
    try:
        if remote == 0:
            with open("media/one_way_blast/"+str(project_id)+"/"+filename) as res:
                data = res.readlines()
        elif remote == 1:
            with open("media/one_way_blast/remote_searches/"+str(project_id)+"/"+filename) as res:
                data = res.readlines()
        return data
    except Exception as e:
        raise FileNotFoundError("Couldn't read file {} with Exception: {}".format(e))