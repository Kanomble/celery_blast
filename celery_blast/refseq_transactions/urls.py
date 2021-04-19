from django.urls import path,include
from . import views

service_urls = [
    path('download_refseq_assembly_summary', views.download_refseq_assembly_summary_view,name='download_refseq_assembly_summary'),
    path('download_refseq_assemblies', views.download_refseq_assemblies_view ,name='download_refseq_assemblies'),
]

urlpatterns = [
    path('', views.dashboard, name='refseq_transactions_dashboard'),
    path('',include(service_urls))
]