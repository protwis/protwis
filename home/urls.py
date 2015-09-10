from django.conf.urls import patterns, url
#from home import views
from django.views.generic import TemplateView
from django.conf import settings

urlpatterns = patterns('',
    url(r'^$', TemplateView.as_view(template_name='home/index_{}.html'.format(settings.SITE_NAME)), {'site_title' : settings.SITE_TITLE}, name='index'), 
)