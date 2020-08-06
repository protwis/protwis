from django.urls import path

from . import views

urlpatterns = [
    path('releasenotes/', views.releasenotes, name='releasenotes'),
]
