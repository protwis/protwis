from django.core.management.base import BaseCommand, CommandError
from django.core.management import call_command
from django.conf import settings
from django.db import connection
from django.db.models import Q
from django.template.loader import render_to_string
from contactnetwork.distances import *
from structure.models import *

# 5NX2 Contact Network 14.856671571731567
# 5WB1 Contact Network 15.49610686302185
# 5L7D Contact Network 16.39739680290222
# 5X93 Contact Network 17.055518865585327
# 5WS3 Contact Network 11.262089014053345
# 5XRA Contact Network 10.246203184127808

import time


class Command(BaseCommand):

    help = "Test distances"


    def handle(self, *args, **options):

        strucs = Structure.objects.all().exclude(structure_type__slug__startswith='af-').prefetch_related('pdb_code')[:100]
        ss = []
        for s in strucs:
            ss.append(s.pdb_code.index)

        for i in [3,30,50,100,200,310]:
            d = Distances()
            d.load_pdbs(ss[:i])
            current = time.time()
            d.fetch_distances()
            d.calculate()
            print(i,'pdbs',time.time()-current, 'time spent with individual records')
            del d

            d = Distances()
            d.load_pdbs(ss[:i])
            current = time.time()
            d.fetch_and_calculate(with_arr=True)
            print(i,'pdbs',time.time()-current, 'time spent with individual records and calc (pure DB)')
            del d


            d = Distances()
            d.load_pdbs(ss[:i])
            current = time.time()
            d.fetch_and_calculate()
            print(i,'pdbs',time.time()-current, 'time spent with individual records and calc (without individual distances) (pure DB)')
            del d


            d = Distances()
            d.load_pdbs(ss[:i])
            current = time.time()
            d.fetch_agg()
            print(i,'pdbs',time.time()-current, 'time spent with individual records (semi DB)')
            del d


            # d = Distances()
            # d.load_pdbs(ss[:i])
            # current = time.time()
            # d.fetch_distances_set()
            # print(i,'pdbs',time.time()-current, 'time spent with msg records')
            # d.calculate()
            # del d

            # d = Distances()
            # d.load_pdbs(ss[:i])
            # current = time.time()
            # d.fetch_distances_set_pickled()
            # print(i,'pdbs',time.time()-current, 'time spent with pickle records')
            # d.calculate()
            # del d


            # d = Distances()
            # d.load_pdbs(ss[:i])
            # current = time.time()
            # d.fetch_distances_set_json()
            # print(i,'pdbs',time.time()-current, 'time spent with json records')
            # d.calculate()
            # del d

            # d = Distances()
            # d.load_pdbs(ss[:i])
            # current = time.time()
            # d.fetch_distances_set_ujson()
            # print(i,'pdbs',time.time()-current, 'time spent with ujson records')
            # d.calculate()
            # del d

            print("########")
            # current = time.time()
            # d.fetch_distances_set()
            # print(i,'pdbs',time.time()-current, 'time spent with pickled records')
            # d.calculate()
