import subprocess
import os

from flask import (
    Blueprint, redirect, request, url_for, render_template
)

fmafft = Blueprint('fmafft',__name__)

@fmafft.route('/mafft_info',methods=['GET'])
def index():
    #print("[+] Connected to flask_api_bioinformatic_tools container, welcome!")
    return "Connected to flask_api_bioinformatic_tools container, welcome!"

@fmafft.route('/project_information/<int:project_id>',methods=['GET'])
def list_informations(project_id):
    return "listing view of project : {}".format(project_id)

@fmafft.route('/perform_simple_msa/<int:project_id>/<folder_path>',methods=['GET'])
def perform_simple_msa(project_id,folder_path):
    path_to_project = 'data/blast_projects/' + str(project_id) + '/'
    path_to_query_file = path_to_project + folder_path + '/target_sequences.faa'
    output = path_to_project + folder_path + '/target_sequences.msa'
    cmd = "mafft {} > {}".format(path_to_query_file,output)
    #print("[*] received task: {}".format(cmd))
    #print("\t[*] {}".format(os.getcwd()))
    try:
        process = subprocess.Popen(cmd, shell=True)
        returncode = process.wait(timeout=5000)
        return str(returncode)
    except subprocess.SubprocessError as e:
        return "1"