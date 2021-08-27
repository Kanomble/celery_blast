from django.shortcuts import render
from django.contrib.auth.decorators import login_required
from django.http.response import HttpResponse

# Create your views here.
@login_required(login_url='login')
def index(request):
        return HttpResponse("hello world")