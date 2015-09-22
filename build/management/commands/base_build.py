from django.core.management.base import BaseCommand, CommandError
from django.core.management import call_command
from django.conf import settings
from django.db import connection

import datetime
import logging
from optparse import make_option
from multiprocessing import Queue, Process


class Command(BaseCommand):
    help = 'Basic functions for build scrips'

    logger = logging.getLogger(__name__)

    def add_arguments(self, parser):
        parser.add_argument('--njobs', action='store', dest='njobs', help='Number of jobs to run')

    def prepare_input(self, njobs, items):
        q = Queue()
        procs = list()
        num_items = len(items)
        
        # make sure not to use more jobs than proteins (chunk size will be 0, which is not good)
        if njobs > num_items:
            njobs = num_items

        chunk_size = int(num_items / njobs)
        connection.close()
        for i in range(0, njobs):
            first = chunk_size * i
            if i == njobs - 1:
                last = False
            else:
                last = chunk_size * (i + 1)
    
            p = Process(target=self.main_func, args=([(first, last)]))
            procs.append(p)
            p.start()

        for p in procs:
            p.join()