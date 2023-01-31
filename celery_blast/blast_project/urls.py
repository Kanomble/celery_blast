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
    path('create_taxonomic_file', views.create_taxonomic_file_view, name='species_taxids'),
    path('upload_genomes/',views.upload_genome_view,name='upload_genomes'),
    path('upload_multiple_genomes/', views.upload_multiple_genomes_view, name='upload_multiple_genomes'),

    path('project_creation',views.project_creation_view,name='project_creation'),
    path('<int:project_id>/project_details',views.project_details_view,name='project_details'),
    path('<int:project_id>/project_deletion',views.project_delete_view,name='project_deletion'),
    path('<int:project_id>/project_execution',views.execute_reciprocal_blast_project_view,name='project_execution'),
    path('<int:project_id>/project_resulttable',views.load_reciprocal_result_html_table_view,name='reciprocal_results')
]

#this ajax call is currently not used - it has been replaced by the snakemake pipeline script query_sequences_to_html_table
ajax_urls = [
    path('<int:project_id>/ajax_wp_to_links',views.ajax_wp_to_links,name='ajax_wp_to_links'),
    path('<int:project_id>/ajax_call_to_logfiles',views.ajax_call_to_logfiles,name='ajax_call_to_logfiles')
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
    path('<int:project_id>/project_details/database_statisitcs_details/<str:taxonomic_unit>', views.database_statistics_details,name='database_statistics_details'),
    path('<int:project_id>/project_details/delete_database_statistics_task_and_output', views.delete_database_statistics,name='delete_database_statistics'),
    path('<int:project_id>/project_details/<str:taxonomic_unit>/taxonomic_unit_ajax_call', views.load_database_statistics_for_taxonomic_unit_ajax, name='ajax_call_for_taxonomic_unit'),
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