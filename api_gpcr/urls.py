from django.conf import settings

from api.urls import *
# from api_SITE_NAME import views
# from common.alignment_SITE_NAME import Alignment
views = getattr(__import__('api_' + settings.SITE_NAME, fromlist=['views']), 'views')
# Alignment = getattr(__import__('common.alignment_' + settings.SITE_NAME, fromlist=['Alignment']), 'Alignment')


urlpatterns.append(url(r'^plot/helixbox/(?P<entry_name>[^/].+)/$', views.HelixBoxView.as_view(), name='helixbox'))
urlpatterns.append(url(r'^plot/snake/(?P<entry_name>[^/].+)/$', views.SnakePlotView.as_view(), name='snakeplot'))