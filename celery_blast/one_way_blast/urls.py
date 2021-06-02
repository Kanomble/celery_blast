'''
urls for blast projects
'''

from django.urls import path, include
from . import views

urlpatterns = [
    path('', views.one_way_blast_dashboard, name='one_way_blast_dashboard'),
]