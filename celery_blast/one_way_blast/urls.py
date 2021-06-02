'''
urls for blast projects
'''

from django.urls import path, include
from . import views

urlpatterns = [
    path('blast_project_creation',views.one_way_blast_project_creation_view,name='one_way_project_creation')
]