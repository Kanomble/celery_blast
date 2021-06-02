'''
urls for blast projects
'''

from django.urls import path, include
from . import views

urlpatterns = [
    path('blast_project_creation',views.create_one_way_blast_project,name='one_way_project_creation')
]