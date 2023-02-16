'''
urls for blast projects
'''

from django.urls import path, include

from . import views

service_urls = [
    path('blast_project_creation', views.one_way_blast_project_creation_view, name='one_way_project_creation'),
    path('<int:project_id>/one_way_project_details', views.one_way_project_details_view,
         name='one_way_project_details'),
    path('<int:project_id>/one_way_project_deletion', views.one_way_project_delete_view,
         name='one_way_project_deletion'),
    path('<int:project_id>/one_way_remote_project_details', views.one_way_remote_project_details_view,
         name='one_way_remote_project_details'),
    path('<int:project_id>/one_way_remote_project_deletion', views.one_way_remote_project_delete_view,
         name='one_way_remote_project_deletion'),
    path('<int:project_id>/one_way_project_execution', views.execute_one_way_blast_project_view,
         name='one_way_project_execution'),
    path('<int:project_id>/one_way_remote_project_execution', views.execute_one_way_remote_blast_project_view,
         name='one_way_remote_project_execution'),
    path('<int:remote>/<int:project_id>/BLAST_results', views.load_one_way_result_html_table_view,
         name='one_way_html_results'),
    path('<int:project_id>/<str:filename>/<str:project_type>/target_sequences', views.one_way_download_target_sequences,
         name='one_way_target_sequence_download')
    # path('<int:project_id>/one_way_project_execution',views.execute_reciprocal_blast_project_view,
    # name='project_execution'), path('<int:project_id>/one_way_project_resulttable',
    # views.load_reciprocal_result_html_table_view,name='reciprocal_results')
]

ajax_urls = [
    path('<int:remote>/<int:project_id>/ajax_one_way_wp_to_links', views.ajax_one_way_wp_to_links,
         name='ajax_one_way_wp_to_links')
]
urlpatterns = [
    path('', include(service_urls)),
    path('', include(ajax_urls))
]
