from django.shortcuts import render

# Create your views here.
from django.shortcuts import render, redirect
from django.contrib.auth.decorators import login_required

@login_required(login_url='login')
def dashboard(request):
    context={}
    return render(request,'refseq_transactions/refseq_transactions_dashboard.html',context)