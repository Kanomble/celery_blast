import subprocess
import os

from flask import (
    Blueprint, redirect, request, url_for, render_template
)

fmafft = Blueprint('fmafft',__name__)

@fmafft.route('/mafft_info',methods=['GET'])
def index():
    print("[+] Connectied to flask_api_bioinformatic_tools container, welcome!")
    return "Hello World!"
