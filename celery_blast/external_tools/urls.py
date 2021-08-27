from django.urls import path, include
from . import views

urlpatterns = [
    path('<int:project_id>/external_project_information', views.project_informations, name='external_project_informations'),
    path('<int:project_id>/<str:folder_path>/perform_simple_msa',
         views.perform_simple_msa,
         name='perform_simple_msa')
]