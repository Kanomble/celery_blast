from django.urls import path, include

from . import views

service_urls = [
    path('<str:summary_file>/download_refseq_assembly_summary', views.download_refseq_assembly_summary_view,
         name='download_refseq_assembly_summary'),
    path('<str:summary_file>/<str:task_id>/assembly_summary_download_progress',
         views.assembly_summary_download_progress_view,
         name='assembly_summary_download_progress'),
    path('create_refseq_database_metadata', views.create_blast_database_model_and_directory,
         name='create_refseq_database_metadata'),
    path('<str:task_id>/database_preview_creation_progress',
         views.database_preview_creation_progress_view,
         name='database_preview_creation_progress'),
    path('<int:database_id>/delete_database', views.delete_blast_database_model_and_directory,
         name='delete_blast_database'),
    path('<int:database_id>/database_details', views.display_blast_database_details_view, name='database_details'),
    path('<int:database_id>/download_and_format_blast_database', views.download_and_format_blast_database,
         name='download_and_format_task'),
    path('<int:database_id>/download_and_format_selected_proteomes', views.download_and_format_selected_proteomes,
         name="download_and_format_selected_proteomes")
]

ajax_calls = [
    path('<int:database_id>/ajax_call', views.ajax_call_for_database_details, name='ajax_call'),
    path('<int:database_id>/ajax_call_progress', views.ajax_call_for_database_download_progress,
         name='ajax_call_progress'),
    path('<str:task_id>/assembly_summary_download_progress_ajax',
         views.assembly_summary_download_progress_ajax,
         name='assembly_summary_download_progress_ajax'),
    path('<str:task_id>/database_preview_creation_progress_ajax',
         views.database_preview_creation_progress_ajax,
         name='database_preview_creation_progress_ajax'),

]

urlpatterns = [
    path('', views.dashboard, name='refseq_transactions_dashboard'),
    path('', include(service_urls)),
    path('', include(ajax_calls))
]
