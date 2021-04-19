from django.urls import path,include
from . import views

service_urls = [
    path('download_refseq_assembly_summary', views.download_refseq_assembly_summary_view,name='download_refseq_assembly_summary'),
    path('create_refseq_database_metadata', views.create_blast_database_model_and_directory, name='create_refseq_database_metadata'),
    path('<int:database_id>/delete_database',views.delete_blast_database_model_and_directory,name='delete_blast_database')
]

urlpatterns = [
    path('', views.dashboard, name='refseq_transactions_dashboard'),
    path('',include(service_urls))
]