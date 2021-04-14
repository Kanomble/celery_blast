'''
urls for blast projects
'''

from django.urls import path, include
from . import views

registration_urls = [
    path('login', views.login_user, name='login'),
    path('logout', views.logout_user, name='logout'),
    path('register', views.registration_view, name='register')
]

service_urls = [
    path('create_taxonomic_file', views.create_taxonomic_file_view, name='species_taxids')
]


success_failure_urls = [
    path('failure', views.failure_view, name='failure_view')
]
urlpatterns = [
    path('', views.dashboard_view, name='blast_project_dashboard'),

    path('', include(registration_urls)),
    path('', include(service_urls)),
    path('', include(success_failure_urls))
]