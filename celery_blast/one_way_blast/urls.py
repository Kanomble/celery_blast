'''
urls for blast projects
'''

from django.urls import path, include
from . import views

service_urls = [
    path('blast_project_creation',views.one_way_blast_project_creation_view,name='one_way_project_creation'),
    path('<int:project_id>/one_way_project_details',views.one_way_project_details_view,name='one_way_project_details'),
    path('<int:project_id>/one_way_project_deletion',views.one_way_project_delete_view,name='one_way_project_deletion'),
    #path('<int:project_id>/one_way_project_execution',views.execute_reciprocal_blast_project_view,name='project_execution'),
    #path('<int:project_id>/one_way_project_resulttable',views.load_reciprocal_result_html_table_view,name='reciprocal_results')
]


urlpatterns = [
    path('', include(service_urls)),
]