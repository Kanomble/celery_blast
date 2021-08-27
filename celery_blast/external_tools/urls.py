from django.urls import path, include
from . import views

urlpatterns = [
    path('<int:project_id>/external_project_information', views.project_informations, name='external_project_informations'),
]