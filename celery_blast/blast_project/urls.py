'''
urls for blast projects
'''

from django.urls import path
from . import views

urlpatterns = [
    path('', views.dashboard, name='blast_project_dashboard'),
    path('login', views.login_user, name='login'),
    path('logout', views.logout_user, name='logout'),
    path('register', views.registration_view, name='register'),
    path('failure', views.failure_view, name='failure_view'),
]