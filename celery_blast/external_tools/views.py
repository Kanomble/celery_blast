from django.shortcuts import render, redirect
from django.contrib.auth.decorators import login_required
from blast_project.views import failure_view, success_view
from .tasks import execute_multiple_sequence_alignment
from .models import ExternalTools

#TODO documentation
@login_required(login_url='login')
def project_informations(request, project_id):
    try:
        context = {}
        qseqids = ExternalTools.objects.get_external_tools_based_on_project_id(project_id)
        context['qseqids'] = qseqids
        context['project_id'] = project_id
        return render(request,"external_tools/external_tools_dashboard.html",context)
    except Exception as e:
        return failure_view(request,e)

#TODO documentation
@login_required(login_url='login')
def perform_simple_msa(request,project_id,query_sequence_id):
    try:
        if request.method == 'POST':
            #celery asynchronous task
            returncode = execute_multiple_sequence_alignment.delay(project_id,query_sequence_id)
            return redirect('external_project_informations',project_id=project_id)
        else:
            e = "There is no GET method for this view function"
            return failure_view(request,e)
    except Exception as e:
        return failure_view(request,e)
