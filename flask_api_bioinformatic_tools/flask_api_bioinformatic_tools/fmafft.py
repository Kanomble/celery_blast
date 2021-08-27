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