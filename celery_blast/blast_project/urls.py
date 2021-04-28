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
    path('project_creation',views.project_creation_view,name='project_creation'),
    path('<int:project_id>/project_details',views.project_details_view,name='project_details'),
    path('<int:project_id>/project_deletion',views.project_delete_view,name='project_deletion'),
    path('<int:project_id>/project_execution',views.execute_reciprocal_blast_project_view,name='project_execution')
]


success_failure_urls = [
    path('failure', views.failure_view, name='failure_view'),
    path('success',views.success_view,name='success_view')
]
urlpatterns = [
    path('', views.dashboard_view, name='blast_project_dashboard'),

    path('', include(registration_urls)),
    path('', include(service_urls)),
    path('', include(success_failure_urls))

]