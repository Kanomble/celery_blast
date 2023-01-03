from django.urls import path, include
from . import views

cdd_search_urls = [
    path('<int:project_id>/project_details/cdd_domain_search_dashboard',
         views.cdd_domain_search_dashboard,
         name='cdd_domain_search_dashboard'),
    path('<int:project_id>/project_details/cdd_domain_search_dashboard/cdd_search_task',
         views.execute_cdd_domain_search_for_target_query,
         name='execute_cdd_domain_search'),
    path('<int:project_id>/project_details/cdd_domain_search_dashboard/<str:query_id>/cdd_domain_search_details',
         views.cdd_domain_search_details_view,
         name='cdd_domain_search_details'),
    path('<int:project_id>/project_details/cdd_domain_search_dashboard/<str:query_id>/delete_cdd_domain_search',
         views.delete_cdd_domain_search_view,
         name='delete_cdd_domain_search')
]

entrez_search_urls = [
    path('entrez_dashboard', views.entrez_dashboard_view, name="entrez_dashboard"),
    path('<int:search_id>/search_details', views.search_detail_view,
         name='search_details'),
    path('<int:search_id>/search_delete', views.delete_search_view,
         name='search_delete'),
    path('<int:search_id>/search_details/download_protein_accessions',
         views.download_proteins_from_entrez_search,
         name='download_protein_accessions'),
    path('<int:search_id>/search_details/view_downloaded_sequences',
         views.view_downloaded_sequences, name='view_downloaded_sequences'),
]

ajax_calls = [
    path('<int:search_id>/ajax_call_progress_fasta_download',
         views.ajax_call_progress_entrezsearch_to_fasta,
         name="ajax_call_progress_entrezsearch_to_fasta"),

    path('<int:query_sequence_id>/ajax_call_progress_phylogeny_task',
         views.ajax_call_progress_phylo_task,
         name="ajax_call_progress_phylogeny_task"),
    path('<int:query_sequence_id>/ajax_call_progress_msa_task',
         views.ajax_call_progress_msa_task,
         name="ajax_call_progress_msa_task"),
    path('<int:project_id>/<str:query_id>/cdd_domain_search_task_status',
         views.get_cdd_task_status_ajax_call,
         name="get_cdd_task_status_ajax_call")
]

phylogenetic_analysis_urls = [
    path('<int:project_id>/<str:query_sequence_id>/perform_simple_msa',
         views.perform_simple_msa,
         name='perform_simple_msa'),
    path('<int:project_id>/<str:query_sequence_id>/perform_fasttree_phylobuild',
         views.perform_fasttree_phylobuild,
         name='perform_fasttree_phylobuild'),
    path('<int:project_id>/perform_simple_msa_for_all_query_sequences',
         views.perform_simple_msa_for_all_query_sequences,
         name='perform_simple_msa_for_all_query_sequences'),
    # perform_fasttree_phylobuild_for_all_query_sequences
    path('<int:project_id>/perform_fasttree_phylobuild_for_all_query_sequences',
         views.perform_fasttree_phylobuild_for_all_query_sequences,
         name='perform_fasttree_phylobuild_for_all_query_sequences'),
    path('<int:project_id>/<str:query_sequence_id>/phylogenetic_information',
         views.phylogenetic_information, name='phylogenetic_information'),
]

urlpatterns = [
    path('', include(cdd_search_urls)),
    path('', include(entrez_search_urls)),
    path('', include(ajax_calls)),
    path('', include(phylogenetic_analysis_urls)),
    path('<int:project_id>/external_project_information', views.project_informations,
         name='external_project_informations'),
]