from django.urls import path, include
from . import views

urlpatterns = [
    path('hello_world', views.index, name='world'),
]