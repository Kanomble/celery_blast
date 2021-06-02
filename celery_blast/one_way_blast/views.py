from django.shortcuts import render, redirect, HttpResponse
from django.contrib.auth.decorators import login_required
from django.contrib.auth import authenticate, login, logout
from blast_project.views import failure_view
@login_required(login_url='login')
def one_way_blast_dashboard(request):
    try:
        context = {}
        return render(request,'one_way_blast/one_way_blast_dashboard.html',context)
    except Exception as e:
        return failure_view(request,e)