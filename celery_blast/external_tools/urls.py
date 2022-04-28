from django.urls import path, include
from . import views

urlpatterns = [
    path('entrez_dashboard', views.entrez_dashboard_view, name="entrez_dashboard"),
    path('<int:search_id>/search_details', views.search_detail_view, name='search_details'),
    path('<int:search_id>/search_delete', views.delete_search_view, name='search_delete'),

    path('<int:project_id>/external_project_information', views.project_informations, name='external_project_informations'),
    path('<int:project_id>/<str:query_sequence_id>/perform_simple_msa',
         views.perform_simple_msa,
         name='perform_simple_msa'),
    path('<int:project_id>/<str:query_sequence_id>/perform_fasttree_phylobuild',
         views.perform_fasttree_phylobuild,
         name='perform_fasttree_phylobuild'),
    path('<int:project_id>/perform_simple_msa_for_all_query_sequences',
         views.perform_simple_msa_for_all_query_sequences,
         name='perform_simple_msa_for_all_query_sequences'),
    #perform_fasttree_phylobuild_for_all_query_sequences
    path('<int:project_id>/perform_fasttree_phylobuild_for_all_query_sequences',
         views.perform_fasttree_phylobuild_for_all_query_sequences,
         name='perform_fasttree_phylobuild_for_all_query_sequences'),

]