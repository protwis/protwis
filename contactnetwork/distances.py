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
        ds_with_key = {}
        if with_arr:
            ds = list(Distance.objects.filter(structure__in=self.structures) \
                            .values('gns_pair') \
                            .annotate(mean = Avg('distance'), var = Variance('distance'), c = Count('distance'),arr=ArrayAgg('distance'),arr2=ArrayAgg('structure')).values_list('gns_pair','mean','var','c','arr','arr2').filter(c__gte=int(0.5*len(self.structures))))
            for i,d in enumerate(ds):
                ds[i] += (d[2]/d[1],)
                ds_with_key[d[0]] = ds[i]
        else:
            ds = list(Distance.objects.filter(structure__in=self.structures).exclude(gns_pair__contains='8x').exclude(gns_pair__contains='12x').exclude(gns_pair__contains='23x').exclude(gns_pair__contains='34x').exclude(gns_pair__contains='45x') \
                            .values('gns_pair') \
                            .annotate(mean = Avg('distance'), var = Variance('distance'), c = Count('distance')).values_list('gns_pair','mean','var','c').filter(c__gte=int(len(self.structures)*0.5)))
            for i,d in enumerate(ds):
                ds[i] += (d[2]/d[1],)
                ds_with_key[d[0]] = ds[i]
        # # print(ds.query)
        # print(ds[1])
        stats_sorted = sorted(ds, key=lambda k: -k[-1])
        #print(ds[1])

        self.stats_key = ds_with_key
        self.stats = stats_sorted

    def calculate_window(self, list_of_gns = None):

        lookup = {}
        bins = []
        n = 3 # Size of windows.

        # Get all possible residues
        res = set()
        if list_of_gns:
            for r in list_of_gns:
                res.add(r)
        else:
            for i,pair in enumerate(sorted(self.stats_key.items())):
                res1, res2 = pair[0].split("_")
                res.add(res1)
                res.add(res2)
        res = sorted(list(res))

        # Make windows
        for r in res:
            # See if it needs to make a new bin
            if r not in lookup:
                # This is a new position, first in line, so create a window from here.
                seg = r.split("x")[0]
                num = int(r.split("x")[1])

                if len(str(num))>2:
                    # This is a bulge, reduce to non-bulge for window, but add to lookup
                    # Assume window is made, since there is no gap before a bulge
                    num =r.split("x")[1][:2]
                    lookup[r] = lookup['{}x{}'.format(seg,num)]

                    # Continue since no need to make a window, as it's made.
                    continue

                # Get the number so it's possible capture all the numbers in the bin
                b = []
                for i in range(num,num+n):
                    gn = "{}x{}".format(seg,str(i))
                    if gn not in res:
                        # If GN not present, exclude it from the label name
                        continue

                    b += [str(i)]
                label = '{}x{}'.format(seg,'_'.join(b))

                for b in b:
                    lookup['{}x{}'.format(seg,b)] = label

                bins.append(label)

        #Create all possible window pairs
        bin_pairs = OrderedDict()
        for i, b1 in enumerate(bins):
            # Calculate average dispersion to all other windows.
            # Track with i, to only do calculations once
            # now to pick all other windows going forward
            for i, b2 in enumerate(bins[i+1:]):
                label = "{}-{}".format(b1,b2) 
                bin_pairs[label] = [[],[],[],[]]

        #Populate all bin_pairs
        for pair in sorted(self.stats_key.items()):
            r = pair[0].split("_")

            # Find out which window each residue belongs to
            try:
                label1 = lookup[r[0]]
                label2 = lookup[r[1]]
            except:
                print("something off with ",r)
                continue
            if label1 == label2:
                # disregard pairs that belong to same window
                continue

            label = "{}-{}".format(label1,label2)
            
            if label not in bin_pairs:
                # This should not be possible, sanity check.
                print("LABEL MISSING IN PAIRS",label)
                continue

            # Add dispersion value
            bin_pairs[label][0].append(pair[1][1])
            bin_pairs[label][1].append(pair[1][2])
            bin_pairs[label][2].append(pair[1][3])
            bin_pairs[label][3].append(pair[1][4])

        pairs_to_remove = set()
        for label, ls in bin_pairs.items():
            for i,l in enumerate(ls):
                c = len(l)
                if c==0:
                    # Sanity check
                    # print("No Data in ",label)
                    pairs_to_remove.add(label)
                    continue
                avg = sum(l) / float(len(l))
                bin_pairs[label][i] = avg

        for p in pairs_to_remove:
            del bin_pairs[p]

        bin_pairs = OrderedDict(sorted(bin_pairs.items(), key=lambda x: -x[1][3]))
        self.bin_pairs = bin_pairs

        stats_window_reduced = []
        for label, ls in bin_pairs.items():
            r1 = label.split("-")[0].split("_")[0]
            r2 = label.split("-")[1].split("_")[0]
            label = ["{}_{}".format(r1,r2)]
            stats_window_reduced.append(label + ls)

        stats_window = []
        stats_window_key = {}
        for i,r1 in enumerate(res):
            for r2 in res[i+1:]:

                label = ["{}_{}".format(r1,r2)]
                label_window = "{}-{}".format(lookup[r1],lookup[r2])
                if lookup[r1]!=lookup[r2] and label_window not in pairs_to_remove:
                    value = bin_pairs["{}-{}".format(lookup[r1],lookup[r2])]
                else:
                    value = [[],[],[],[]]
                    continue

                stats_window.append(label + value)
                stats_window_key[label[0]] = label + value

        stats_window = sorted(stats_window, key=lambda k: -k[-1])
        stats_window_reduced = sorted(stats_window_reduced, key=lambda k: -k[-1])
        self.stats_window = stats_window
        self.stats_window_reduced = stats_window_reduced
        self.stats_window_key = stats_window_key
        # print(stats_window_key)

    def fetch_distances(self):
        ds = Distance.objects.filter(structure__in=self.structures).all().values_list('distance','gns_pair')
        self.data = {}
        for d in ds:
            label = d[1]
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
