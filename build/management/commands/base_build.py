from django.core.management.base import BaseCommand, CommandError
from django.core.management import call_command
from django.conf import settings
from django.db import connection

import datetime
import logging
from multiprocessing import Queue, Process


class Command(BaseCommand):
    help = 'Basic functions for build scrips'

    logger = logging.getLogger(__name__)

    def add_arguments(self, parser):
        parser.add_argument('-p', '--proc',
            type=int,
            action='store',
            dest='proc',
            default=1,
            help='Number of processes to run')
        parser.add_argument('-t', '--test',
            action='store_true',
            dest='test',
            default=False,
            help='Include only a subset of data for testing')

    def prepare_input(self, proc, items, iteration=1):
        q = Queue()
        procs = list()
        num_items = len(items)

        if not num_items:
            return False

        # make sure not to use more jobs than proteins (chunk size will be 0, which is not good)
        if proc > num_items:
            proc = num_items

        chunk_size = int(num_items / proc)
        connection.close()
        for i in range(0, proc):
            first = chunk_size * i
            if i == proc - 1:
                last = False
            else:
                last = chunk_size * (i + 1)
    
            p = Process(target=self.main_func, args=([(first, last), iteration]))
            procs.append(p)
            p.start()

        for p in procs:
            p.join()