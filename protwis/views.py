# Imports
from django.shortcuts import render
from django.http import HttpResponse
from django.template import Context, loader


##
# Handle 404 Errors
# @param request WSGIRequest list with all HTTP Request
def error404(request, *args, **kwargs):
    try:
        template = loader.get_template('home/404.html')
        context = {'message': 'All: %s' % request}
    except Exception as e: 
        print(e)
        context = {}

    return HttpResponse(content=template.render(context), content_type='text/html; charset=utf-8', status=404)

def error500(request):
    try:
        template = loader.get_template('home/500.html')
        context = {'message': 'All: %s' % request}
    except Exception as e: 
        print(e)
        context = {}

    return HttpResponse(content=template.render(context), content_type='text/html; charset=utf-8', status=500)