from django.shortcuts import render

# Create your views here.
from django.shortcuts import render, redirect
from django.contrib.auth.decorators import login_required
from .py_services import refseq_file_exists
from .forms import RefseqTableForm
from .tasks import download_refseq_assembly_summary_file
from blast_project.views import failure_view

@login_required(login_url='login')
def dashboard(request):
    try:
        context = {}

        if(refseq_file_exists()):
            print("TRUE")
            context['refseq_exists'] = True

        refseq_form = RefseqTableForm()
        context['RefseqTableForm'] = refseq_form

        return render(request, 'refseq_transactions/refseq_transactions_dashboard.html', context)
    except Exception as e:
        return failure_view(request,e)

@login_required(login_url='login')
def download_refseq_assembly_summary_view(request):
    try:
        download_refseq_assembly_summary_file()
        return redirect('refseq_transactions_dashboard')
    except Exception as e:
        return failure_view(request,e)