from django.shortcuts import render, redirect
from django.contrib.auth.decorators import login_required
from django.http.response import HttpResponse
from .py_services import get_list_of_query_sequence_folder
from blast_project.views import failure_view, success_view
from .tasks import execute_multiple_sequence_alignment
# Create your views here.
@login_required(login_url='login')
def project_informations(request, project_id):
    try:
        context = {}
        qseqid_folder = get_list_of_query_sequence_folder(project_id)
        context['qseqid_folder'] = qseqid_folder
        context['project_id'] = project_id
        return render(request,"external_tools/external_tools_dashboard.html",context)
    except Exception as e:
        return failure_view(request,e)

@login_required(login_url='login')
def perform_simple_msa(request,project_id,folder_path):
    try:
        if request.method == 'POST':
            context = {}
            returncode = execute_multiple_sequence_alignment.delay(project_id,folder_path)
            return redirect('external_project_informations',project_id=project_id)
        else:
            e = "There is no GET method for this view function"
            return failure_view(request,e)
    except Exception as e:
        return failure_view(request,e)