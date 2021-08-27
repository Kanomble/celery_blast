from django.shortcuts import render
from django.contrib.auth.decorators import login_required
from django.http.response import HttpResponse
from .py_services import get_list_of_query_sequence_folder
from blast_project.views import failure_view, success_view

# Create your views here.
@login_required(login_url='login')
def project_informations(request, project_id):
    try:
        context = {}
        qseqid_folder = get_list_of_query_sequence_folder(project_id)
        context['qseqid_folder'] = qseqid_folder
        return render(request,"external_tools/external_tools_dashboard.html",context)
    except Exception as e:
        return failure_view(request,e)