from django.urls import path, include
from . import views

urlpatterns = [
    path('<int:project_id>/external_project_information', views.project_informations, name='external_project_informations'),
    path('<int:project_id>/<str:query_sequence_id>/perform_simple_msa',
         views.perform_simple_msa,
         name='perform_simple_msa'),
    path('<int:project_id>/<str:query_sequence_id>/perform_fasttree_phylobuild',
         views.perform_fasttree_phylobuild,
         name='perform_fasttree_phylobuild')
]