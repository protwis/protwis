from django.conf import settings
from django.http import HttpResponse
from django.db import connection

import time,datetime,os

import uuid

class StatsMiddleware:
    def __init__(self, get_response):
        self.get_response = get_response
        # One-time configuration and initialization.

    def __call__(self, request):
        # Code to be executed for each request before
        # the view (and later middleware) are called.
        start_time = time.time()

        request_id = uuid.uuid4().hex

        text_file = open(os.path.join(settings.BASE_DIR, "logs/stats_start_stop.log"), "a")
        text_file.write('%s %s %s START %s %s\n' % (datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"),request.META.get('REMOTE_ADDR'),request_id, request.method, request.path ))
        text_file.close()

        response = self.get_response(request)

        # Code to be executed for each request/response after
        # the view is called.
        total = time.time() - start_time
        if settings.DEBUG:
            print(request.path,"Time to execute", round(total,2), "SQL queries",len(connection.queries))

        text_file = open(os.path.join(settings.BASE_DIR, "logs/stats.log"), "a")
        text_file.write('%s %s %s %s %s\n' % (datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"), round(total,2),request.META.get('REMOTE_ADDR'), request.method, request.path ))
        text_file.close()

        text_file = open(os.path.join(settings.BASE_DIR, "logs/stats_start_stop.log"), "a")
        text_file.write('%s %s %s FINISH %s %s %s\n' % (datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"),request.META.get('REMOTE_ADDR'),request_id, round(total,2), request.method, request.path ))
        text_file.close()

        if total>5:
            # If slower than 5 seconds
            text_file = open(os.path.join(settings.BASE_DIR, "logs/stats_slow.log"), "a")
            text_file.write('%s %s %s %s %s\n' % (datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"), round(total,2),request.META.get('REMOTE_ADDR'), request.method, request.path ))
            text_file.close()

        return response

    def process_exception(self, request, exception):
        text_file = open(os.path.join(settings.BASE_DIR, "logs/errors.log"), "a")
        text_file.write('%s %s %s %s "%s"\n' % (datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"),request.META.get('REMOTE_ADDR'), request.method, request.path,str(exception) ))
        text_file.close()