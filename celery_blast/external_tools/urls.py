from django.urls import path, include
from . import views

urlpatterns = [
    path('<int:project_id>/external_project_information', views.project_informations, name='external_project_informations'),
    path('<int:project_id>/<str:folder_path>/view_query_folder_informations',
         views.view_query_folder_informations,
         name='view_query_folder_informations')
]