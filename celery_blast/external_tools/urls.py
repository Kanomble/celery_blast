from django.urls import path, include
from . import views

urlpatterns = [
    path('entrez_dashboard', views.entrez_dashboard_view, name="entrez_dashboard"),
    path('<int:search_id>/search_details', views.search_detail_view, name='search_details'),
    path('<int:search_id>/search_delete', views.delete_search_view, name='search_delete'),
    path('<int:search_id>/search_details/download_protein_accessions',
         views.download_proteins_from_entrez_search, name='download_protein_accessions'),
    path('<int:search_id>/search_details/view_downloaded_sequences',
         views.view_downloaded_sequences, name='view_downloaded_sequences'),

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
    path('<int:search_id>/ajax_call_progress_fasta_download',
         views.ajax_call_progress_entrezsearch_to_fasta,
         name="ajax_call_progress_entrezsearch_to_fasta"),

    path('<int:query_sequence_id>/ajax_call_progress_phylogeny_task',
         views.ajax_call_progress_phylo_task,
         name="ajax_call_progress_phylogeny_task"),

    path('<int:query_sequence_id>/ajax_call_progress_msa_task',
         views.ajax_call_progress_msa_task,
         name="ajax_call_progress_msa_task"),
]