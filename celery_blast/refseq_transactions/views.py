from django.shortcuts import render

# Create your views here.
from django.shortcuts import render, redirect
from django.contrib.auth.decorators import login_required

@login_required(login_url='login')
def dashboard(request):
    ''' NOT NEEDED ANYMORE
    try:
        projects = BlastProject.objects.filter(project_username=request.user)
        #query set of project specific databases
        target_genomes = Genomes.objects.filter(associated_project__in=projects)
        #deletes folders that are not included as projects ids in the database
        #delete_files_without_projects()
    except Exception as e:
        return failure_view(request,e)
    context = {
        'projects':projects,
        'genomes':target_genomes,
        }
    '''
    context={}
    return render(request,'refseq_transactions/refseq_transactions_dashboard.html',context)