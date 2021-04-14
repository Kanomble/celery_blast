from django.shortcuts import render, redirect
from django.contrib import messages
from .view_access_decorators import unauthenticated_user
from django.contrib.auth.decorators import login_required
from django.contrib.auth import authenticate, login, logout
from .forms import CreateUserForm

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
    return render(request,'blast_project/blast_project_dashboard.html',context)

#register an account with email and password, email can be used inside biopython functions
#user needs to authenticate otherwise they will get redirected to this login page

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
            return redirect('main')
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

#if an exception occurres this page is rendered in order to evaluate the exception context
def failure_view(request,exception):
    context={'exception':exception}
    return render(request,'blast_project/failure.html', context)