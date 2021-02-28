from django.conf import settings
from django.http import HttpResponse
from django.db import connection

import time, datetime, os
#import uuid

class StatsMiddleware:
    def __init__(self, get_response):
        self.get_response = get_response
        # One-time configuration and initialization.

    def __call__(self, request):
        '''Handling protwis request logs.'''
#       # Start/top logger - not needed as we do slow query logging instead
#        request_id = uuid.uuid4().hex
#        text_file = open(os.path.join(settings.BASE_DIR, "logs/stats_start_stop.log"), "a")
#        text_file.write('%s %s %s START %s %s\n' % (datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"),request.META.get('REMOTE_ADDR'),request_id, request.method, request.path ))
#        text_file.close()

        # start timer
        start_time = time.time()

        # Handle request
        response = self.get_response(request)

        # Code below is executed after handling the request

        # End timer
        total = time.time() - start_time

        if settings.DEBUG:
            print(request.path,"Time to execute", round(total,2), "SQL queries",len(connection.queries))

        bot_user = self.bot_detection(request)
        log_file = os.path.join(settings.BASE_DIR, "logs/stats.log")
        if bot_user:
            log_file = os.path.join(settings.BASE_DIR, "logs/stats_bots.log")

        text_file = open(log_file, "a")
        if bot_user:
            text_file.write('%s %s %s %s %s %s\n' % (datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"), round(total,2), request.META.get('HTTP_X_FORWARDED_FOR'), request.method, request.path, request.META.get('HTTP_USER_AGENT') ))
        else:
            text_file.write('%s %s %s %s %s\n' % (datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"), round(total,2), request.META.get('HTTP_X_FORWARDED_FOR'), request.method, request.path))

        text_file.close()

#       # Start/top logger - not needed as we do slow query logging instead
#        text_file = open(os.path.join(settings.BASE_DIR, "logs/stats_start_stop.log"), "a")
#        text_file.write('%s %s %s FINISH %s %s %s\n' % (datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"),request.META.get('REMOTE_ADDR'),request_id, round(total,2), request.method, request.path ))
#        text_file.close()

        # Extended logging of queries that are slower than 5 seconds
        if total>5:
            log_file = os.path.join(settings.BASE_DIR, "logs/stats_slow.log")
            if bot_user:
                log_file = os.path.join(settings.BASE_DIR, "logs/stats_slow_bots.log")
            text_file = open(log_file, "a")
            text_file.write('%s %s %s %s %s %s\n' % (datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"), round(total,2), request.META.get('HTTP_X_FORWARDED_FOR'), request.method, request.path, request.META.get('HTTP_USER_AGENT') ))
            text_file.close()

        return response

    def process_exception(self, request, exception):
        log_file = os.path.join(settings.BASE_DIR, "logs/errors.log")
        if self.bot_detection(request):
            log_file = os.path.join(settings.BASE_DIR, "logs/errors_bots.log")

        text_file = open(log_file, "a")
        text_file.write('%s %s %s %s - %s %s - %s "%s" "%s"\n' % (datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"), request.method, request.path, request.META.get('HTTP_REFERER'), request.META.get('HTTP_X_FORWARDED_FOR'), request.META.get('HTTP_USER_AGENT'), str(exception) ))
        text_file.close()

    @staticmethod
    def bot_detection(request):
        '''Simplistic bot detection to split log files.

        This function matches the user agent string for specific words matching
        all main search engine agents and >60% of any bot agents.
        Analyzed Feb-2021 on listings at useragentstring.com.
        '''
        bot_IDs = ["bot", "slurp", "crawler", "spider", "curl"]
        user_agent = request.META.get("HTTP_USER_AGENT", "").lower()
        return any(bot_ID in user_agent for bot_ID in bot_IDs)
