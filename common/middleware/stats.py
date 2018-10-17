from django.conf import settings

import time,datetime

class StatsMiddleware:
    def __init__(self, get_response):
        self.get_response = get_response
        # One-time configuration and initialization.

    def __call__(self, request):
        # Code to be executed for each request before
        # the view (and later middleware) are called.
        start_time = time.time()

        response = self.get_response(request)

        # Code to be executed for each request/response after
        # the view is called.
        total = time.time() - start_time

        text_file = open("logs/stats.log", "a")
        text_file.write('%s %s %s %s %s\n' % (datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"), round(total,2),request.META.get('REMOTE_ADDR'), request.method, request.path ))
        text_file.close()

        if total>5:
            # If slower than 5 seconds
            text_file = open("logs/stats_slow.log", "a")
            text_file.write('%s %s %s %s %s\n' % (datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"), round(total,2),request.META.get('REMOTE_ADDR'), request.method, request.path ))
            text_file.close()

        return response

class ProcessExceptionMiddleware(StatsMiddleware):
    def process_exception(self, request, exception):
        text_file = open("logs/errors.log", "a")
        text_file.write('%s %s %s %s "%s"\n' % (datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"),request.META.get('REMOTE_ADDR'), request.method, request.path,str(exception) ))
        text_file.close()
        return HttpResponse('Exception caught')