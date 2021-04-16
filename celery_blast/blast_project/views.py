from django.shortcuts import render, redirect
from django.contrib import messages
from .view_access_decorators import unauthenticated_user
from django.contrib.auth.decorators import login_required
from django.contrib.auth import authenticate, login, logout
from .forms import CreateUserForm, CreateTaxonomicFileForm
from .tasks import write_species_taxids_into_file
from .py_services import list_taxonomic_files

from .py_django_db_services import get_users_blast_projects, get_all_blast_databases
''' dashboard

view for the first dashboard page, this page enables monitoring of blast_projects,
created by the currently logged in user.

:GET
    Uses the get_users_blast_projects utility function which returns a BlastProject Query-Set.
    The Query-Set inherits all projects from the currently logged in user.
    Uses the get_all_blast_databases utility functions to load all BlastDatabase db entries.
    display blast_projects and links to other view functions

'''
@login_required(login_url='login')
def dashboard_view(request):
    try:
        context = {}
        if request.method == 'GET':
            users_blast_projects = get_users_blast_projects(request.user.id)
            available_blast_databases = get_all_blast_databases()
            context['blast_projects'] = users_blast_projects
            context['blast_databases'] = available_blast_databases

        return render(request,'blast_project/blast_project_dashboard.html',context)
    except Exception as e:
        return failure_view(request,e)

@login_required(login_url='login')
def project_creation_view(request):
    try:
        context = {}
        return render(request,'blast_project/project_creation_dashboard.html',context)
    except Exception as e:
        return failure_view(request,e)

''' create_taxonomic_file_view

view for creation of taxonomic files, produced by the get_species_taxids.sh script.

:GET
    display available taxonomic files
:POST
    create taxonomic files
    synchronous call of write_species_taxids_into_file
'''
@login_required(login_url='login')
def create_taxonomic_file_view(request):
    try:
        taxform = CreateTaxonomicFileForm(request.user)
        if request.method == 'POST':
            taxform = CreateTaxonomicFileForm(request.user,request.POST)
            if taxform.is_valid():
                species_name,taxonomic_node = taxform.cleaned_data['species_name']
                task = write_species_taxids_into_file(taxonomic_node,species_name+'.taxids')
        taxid_files = list_taxonomic_files()
        context = {'taxform': taxform, 'taxid_files': taxid_files}
        return render(request, 'blast_project/create_taxonomic_file.html', context)
    except Exception as e:
        return failure_view(request,e)

''' registration, login and logout views
register an account with email and password, email can be used inside biopython functions
user needs to authenticate otherwise they will get redirected to this login page

'''
#login user
@unauthenticated_user #you dont need an account to trigger this view
def login_user(request):

    #login with django default authenticate method
    if request.method == 'POST':
        username = request.POST.get('username')
        password = request.POST.get('password')
        user = authenticate(request, username=username, password=password)
        if user is not None:
            login(request, user)
            return redirect('blast_project_dashboard')
        else:
            messages.info(request,'Username OR password is incorrect')
            return render(request, 'blast_project/login.html')
    return render(request,'blast_project/login.html')

#logout view
def logout_user(request):
    logout(request)
    return redirect('login')

#registration view
@unauthenticated_user
def registration_view(request):
    userForm = CreateUserForm()
    if request.method == 'POST':
        userForm = CreateUserForm(request.POST)
        if userForm.is_valid():
            try:
                userForm.save()
                username = userForm.cleaned_data.get('username')
                #group = Group.objects.get(name='customer')
                #user.groups.add(group)
                messages.success(request,'Account was created for '+ username)
                return redirect('login')
            except Exception as e:
                return failure_view(request,e)
    context = {'form': userForm, }
    return render(request,'blast_project/register.html',context)

''' failure view

returned if an exception ocurred within execution of view functions.

:param exception
    :type str
'''
#if an exception occurres this page is rendered in order to evaluate the exception context
def failure_view(request,exception):
    context={'exception':exception}
    return render(request,'blast_project/failure.html', context)