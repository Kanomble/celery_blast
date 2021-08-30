import subprocess

from flask import (
    Blueprint, redirect, request, url_for, render_template, Response
)
from http import HTTPStatus
fmafft = Blueprint('fmafft',__name__)

@fmafft.route('/mafft_info',methods=['GET'])
def index():
    #print("[+] Connected to flask_api_bioinformatic_tools container, welcome!")
    return "Connected to flask_api_bioinformatic_tools container, welcome!"

@fmafft.route('/project_information/<int:project_id>',methods=['GET'])
def list_informations(project_id):
    return "listing view of project : {}".format(project_id)

@fmafft.route('/perform_simple_msa/<int:project_id>/<folder_path>',methods=['GET','POST'])
def perform_simple_msa(project_id,folder_path):
    if request.method == 'POST':
        project_id = request.json['project_id']
        folder_path = request.json['folder_path']
        path_to_project = 'data/blast_projects/' + str(project_id) + '/'
        path_to_query_file = path_to_project + folder_path + '/target_sequences.faa'
        output = path_to_project + folder_path + '/target_sequences.msa'
        cmd = "mafft {} > {}".format(path_to_query_file,output)
        try:
            process = subprocess.Popen(cmd, shell=True)
            returncode = process.wait(timeout=5000)
            if returncode != 0:
                raise Exception
            return Response("0",status=HTTPStatus.OK,mimetype="str")
        except subprocess.SubprocessError:
            return Response("1",status=HTTPStatus.INTERNAL_SERVER_ERROR,mimetype="str")
    else:
        return Response("1",status=HTTPStatus.BAD_REQUEST,mimetype="str")