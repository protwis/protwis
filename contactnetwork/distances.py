from django.db.models import Avg, Variance, Count
from django.contrib.postgres.aggregates import ArrayAgg

from structure.models import Structure
from contactnetwork.models import *

from collections import OrderedDict

import numpy as np

class Distances():
    """A class to do distances"""
    def __init__(self):
        self.structures = []
        self.generic_numbers = OrderedDict()
        self.segments = OrderedDict()
        self.distances = {}

    def load_pdbs(self, pdbs):
        """Load a list of pdbs objects"""
        structures = Structure.objects.filter(pdb_code__index__in=pdbs).all()
        for s in structures:
            self.structures.append(s)


    def fetch_agg(self):
        ## REQUIRES PSQL SETTINGS TO HAVE MORE MEMORY
        # sudo nano /etc/postgresql/9.3/main/postgresql.conf
        # shared_buffers = 2GB  
        # work_mem = 100MB  
        # temp_buffers = 500MB
        # sudo /etc/init.d/postgresql restart 

        ds = list(Distance.objects.filter(structure__in=self.structures) \
                            .values('gns_pair').annotate(arr=ArrayAgg('distance')).values_list('gns_pair','arr'))
        # #print(ds)
        # print(ds.query)
        # print(len(ds))
        # for i,d in enumerate(ds):
        #     print(d)
        #     if i>100:
        #         break
        self.data = ds
        self.stats = {}
        self.stats_list = []
        for label,d in ds:
            var = round(np.var(d),2)
            mean = round(np.mean(d),2)
            dispersion = round(var/mean,3)
            self.stats[label] = {'var': var , 'mean': mean, 'dispersion': dispersion, 'count':len(d)}
            self.stats_list.append({'label':label, 'var': var , 'mean': mean, 'dispersion': dispersion, 'count':len(d)})
            # continue
        print(len(self.stats))

    def fetch_and_calculate(self, with_arr = False):
        ## REQUIRES PSQL SETTINGS TO HAVE MORE MEMORY
        # sudo nano /etc/postgresql/9.3/main/postgresql.conf
        # shared_buffers = 2GB  
        # work_mem = 100MB  
        # temp_buffers = 500MB
        # sudo /etc/init.d/postgresql restart 
        if with_arr:
            ds = list(Distance.objects.filter(structure__in=self.structures) \
                            .values('gns_pair') \
                            .annotate(mean = Avg('distance'), var = Variance('distance'), c = Count('distance'),arr=ArrayAgg('distance'),arr2=ArrayAgg('structure')).values_list('gns_pair','mean','var','c','arr','arr2'))
            for i,d in enumerate(ds):
                ds[i] += (d[2]/d[1],)
        else:
            ds = list(Distance.objects.filter(structure__in=self.structures) \
                            .values('gns_pair') \
                            .annotate(mean = Avg('distance'), var = Variance('distance'), c = Count('distance')).values_list('gns_pair','mean','var','c'))
            for i,d in enumerate(ds):
                ds[i] += (d[2]/d[1],)
        # # print(ds.query)
        # print(ds[1])
        stats_sorted = sorted(ds, key=lambda k: -k[-1]) 
        print(stats_sorted[1])

        #print(ds[1])


        self.stats = ds


    def fetch_distances(self):
        ds = Distance.objects.filter(structure__in=self.structures).all().values_list('distance','gns_pair')
        self.data = {}
        for d in ds:
            label = d[1]
            if label not in self.data:
                self.data[label] = []
            self.data[label].append(d[0]/100)

    def fetch_distances_set(self):
        dss = DistanceSet.objects.filter(structure__in=self.structures).only("distances").all()
        self.data = {}
        for ds in dss:
            distances = ds.get_distances()
            for d in distances:
                label = '_'.join([d[1], d[2]])
                if label not in self.data:
                    self.data[label] = []
                self.data[label].append(d[0]/100)

    def fetch_distances_set_pickled(self):
        dss = DistanceSet.objects.filter(structure__in=self.structures).only("distances_pickle").all()
        self.data = {}
        for ds in dss:
            distances = ds.get_distances_pickle()
            for d in distances:
                label = '_'.join([d[1], d[2]])
                if label not in self.data:
                    self.data[label] = []
                self.data[label].append(d[0]/100)

    def fetch_distances_set_json(self):
        dss = DistanceSet.objects.filter(structure__in=self.structures).only("distances_json").all()
        self.data = {}
        for ds in dss:
            distances = ds.get_distances_json()
            for d in distances:
                label = '_'.join([d[1], d[2]])
                if label not in self.data:
                    self.data[label] = []
                self.data[label].append(d[0]/100)

    def fetch_distances_set_ujson(self):
        dss = DistanceSet.objects.filter(structure__in=self.structures).only("distances_ujson").all()
        self.data = {}
        for ds in dss:
            distances = ds.get_distances_ujson()
            for d in distances:
                label = '_'.join([d[1], d[2]])
                if label not in self.data:
                    self.data[label] = []
                self.data[label].append(d[0]/100)

    def calculate(self):
        self.stats = {}
        self.stats_list = []
        for label,d in self.data.items():
            var = round(np.var(d),2)
            mean = round(np.mean(d),2)
            dispersion = round(var/mean,3)
            self.stats[label] = {'var': var , 'mean': mean, 'dispersion': dispersion, 'count':len(d)}
            self.stats_list.append({'label':label, 'var': var , 'mean': mean, 'dispersion': dispersion, 'count':len(d)})
            # continue
        print(len(self.stats))
        # stats_sorted = self.stats_list.sort(key=lambda x: x['dispersion'])
        # stats_sorted = sorted(self.stats_list, key=lambda k: -k['dispersion']) 
        # print(stats_sorted[1])
        #print(data)
            #print(d.interacting_pair.res1.generic_number.label)
