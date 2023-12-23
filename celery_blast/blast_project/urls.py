'''
urls for blast projects
'''

from django.urls import path, include
from . import views

registration_urls = [
    path('login', views.login_user, name='login'),
    path('logout', views.logout_user, name='logout'),
    path('register', views.registration_view, name='register'),
]

service_urls = [
    path('setup_cathi_view', views.setup_cathi_view, name='setup_cathi_view'),
    path('delete_domain_database',views.delete_domain_database_view, name='delete_domain_database'),
    path('<str:selected_table>/datatable_view', views.active_table_view, name='active_table_view'),

    path('create_taxonomic_file', views.create_taxonomic_file_view, name='species_taxids'),
    path('upload_genomes/',views.upload_genome_view,name='upload_genomes'),
    path('upload_multiple_genomes/', views.upload_multiple_genomes_view, name='upload_multiple_genomes'),

    path('project_creation',views.project_creation_view,name='project_creation'),
    path('<int:project_id>/project_details',views.project_details_view,name='project_details'),
    path('<int:project_id>/remote_project_details', views.remote_project_details_view, name='remote_project_details'),
    path('<int:project_id>/project_deletion',views.project_delete_view,name='project_deletion'),
    path('<int:project_id>/project_execution',views.start_reciprocal_blast_project_view,name='project_execution'),
    path('<int:project_id>/remote_project_deletion', views.remote_project_delete_view, name='remote_project_deletion'),
    path('<int:project_id>/remote_project_execution', views.start_remote_reciprocal_blast_project_view, name='remote_project_execution'),
    path('<int:project_id>/project_resulttable',views.load_reciprocal_result_html_table_view,name='reciprocal_results'),
    path('<int:project_id>/download_archive', views.download_project_as_zip_archive_view,
         name='download_archive'),
]

# this ajax call is currently not used - it has been replaced by the snakemake pipeline script
# query_sequences_to_html_table
ajax_urls = [
    path('<int:project_id>/ajax_wp_to_links',views.ajax_wp_to_links,name='ajax_wp_to_links'),
    path('<int:project_id>/ajax_call_to_logfiles',views.ajax_call_to_logfiles,name='ajax_call_to_logfiles'),
    path('<int:project_id>/project_details/database_statistics/ajax_call_to_selection_task',
         views.database_selection_phylogeny_task_status,name='ajax_selection_constrained_phylogeny'),
    path('domain_database_download_progress', views.get_domain_database_download_task_status,
         name='domain_database_download_status')
]

progress_urls = [
    path('<int:project_id>/project_details/<str:logfile>',
         views.send_logfile_content_view,
         name='send_logfile_content'),
    path('<int:project_id>/project_details/progress/query_sequence_information',
         views.send_query_sequence_information_view,
         name='send_query_sequence_information')
]
py_optional_postprocessing = [
    path('<int:project_id>/project_details/database_statistics', views.database_statistics_dashboard, name='database_statistics'),
    path('<int:project_id>/project_details/execute_database_statistics_task', views.execute_database_statistics_task,name='execute_database_statistics_task'),
    path('<int:project_id>/project_details/delete_database_statistics_task_and_output', views.delete_database_statistics,name='delete_database_statistics'),
    path('<int:project_id>/project_details/<str:taxonomic_unit>/taxonomic_unit_ajax_call', views.load_database_statistics_for_taxonomic_unit_ajax, name='ajax_call_for_taxonomic_unit'),
    path('<int:project_id>/project_details/database_statistics/selection_contrained_phylogeny', views.view_selection_phylogeny, name="view_selection_phylogeny")
]

success_failure_urls = [
    path('failure', views.failure_view, name='failure_view'),
    path('success',views.success_view,name='success_view')
]

urlpatterns = [
    path('', views.dashboard_view, name='blast_project_dashboard'),
    path('', include(py_optional_postprocessing)),
    path('', include(ajax_urls)),
    path('', include(registration_urls)),
    path('', include(service_urls)),
    path('', include(success_failure_urls)),
    path('', include(progress_urls)),

]