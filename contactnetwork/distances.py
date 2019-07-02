from django.db.models import Avg, Variance, Count, Value, StdDev
from django.contrib.postgres.aggregates import ArrayAgg
from django.core.cache import cache

from structure.models import Structure
from contactnetwork.models import *
from residue.models import Residue, ResidueGenericNumber

from collections import OrderedDict

import numpy as np

class Distances():
    """A class to do distances"""
    def __init__(self):
        self.structures = []
        self.pdbs = []
        self.pconfs = []
        self.generic_numbers = OrderedDict()
        self.segments = OrderedDict()
        self.distances = {}

    def load_pdbs(self, pdbs):
        """Load a list of pdbs objects"""
        structures = Structure.objects.filter(pdb_code__index__in=pdbs).prefetch_related('protein_conformation').all()
        for s in structures:
            self.structures.append(s)
            self.pdbs.append(s.pdb_code.index)
            self.pconfs.append(s.protein_conformation)


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
            ds = list(Distance.objects.filter(structure__in=self.structures).exclude(gns_pair__contains='8x').exclude(gns_pair__contains='12x').exclude(gns_pair__contains='23x').exclude(gns_pair__contains='34x').exclude(gns_pair__contains='45x') \
                            .values('gns_pair') \
                            .annotate(mean = Avg('distance'), std = StdDev('distance'), c = Count('distance'), dis = Count('distance'),arr=ArrayAgg('distance'),arr2=ArrayAgg('structure__pdb_code__index'),arr3=ArrayAgg('gns_pair')).values_list('gns_pair','mean','std','c','dis','arr','arr2','arr3').filter(c__gte=int(0.8*len(self.structures))))
            for i,d in enumerate(ds):
                ds[i] = list(ds[i])
                ds[i][3] = d[2]/d[1]
                ds_with_key[d[0]] = ds[i]
        else:
            ds = list(Distance.objects.filter(structure__in=self.structures).exclude(gns_pair__contains='8x').exclude(gns_pair__contains='12x').exclude(gns_pair__contains='23x').exclude(gns_pair__contains='34x').exclude(gns_pair__contains='45x') \
                            .values('gns_pair') \
                            .annotate(mean = Avg('distance'), std = StdDev('distance'), c = Count('distance')).values_list('gns_pair','mean','std','c').filter(c__gte=int(len(self.structures)*0.8)))
            for i,d in enumerate(ds):
                ds[i] += (d[2]/d[1],)
                ds_with_key[d[0]] = ds[i]
        # # print(ds.query)
        # print(ds[1])
        # Assume that dispersion is always 4
        if len(self.structures)>1:
            stats_sorted = sorted(ds, key=lambda k: -k[3])
        else:
            stats_sorted = sorted(ds, key=lambda k: -k[1])
        #print(ds[1])

        self.stats_key = ds_with_key
        self.stats = stats_sorted

    def fetch_common_gns_tm(self):

        # ds = list(Distance.objects.filter(structure__in=self.structures).exclude(gns_pair__contains='8x').exclude(gns_pair__contains='12x').exclude(gns_pair__contains='23x').exclude(gns_pair__contains='34x').exclude(gns_pair__contains='45x') \
        #                 .values('gns_pair') \
        #                 .annotate(c = Count('distance')).values_list('gns_pair').filter(c__gte=int(len(self.structures)*0.9)))
        # common_gn = []
        # for i,d in enumerate(ds):
        #     res1, res2 = d[0].split("_")
        #     if res1 not in common_gn:
        #         common_gn.append(res1)
        #     if res2 not in common_gn:
        #         common_gn.append(res2)
        # common_gn.sort()
        res = Residue.objects.filter(protein_conformation__in=self.pconfs) \
                        .exclude(generic_number=None) \
                        .exclude(generic_number__label__contains='8x') \
                        .exclude(generic_number__label__contains='12x') \
                        .exclude(generic_number__label__contains='23x') \
                        .exclude(generic_number__label__contains='34x') \
                        .exclude(generic_number__label__contains='45x') \
                        .values('generic_number__label') \
                        .annotate(c = Count('generic_number__label')).filter(c__gte=int(len(self.structures)*0.9)).values_list('generic_number__label',flat=True).order_by('c') #.values_list('generic_number__label') #.filter(c__gte=int(len(self.structures)*0.9))

        common_gn = sorted(res)

        return common_gn

    def calculate_window(self, list_of_gns = None):

        lookup = {}
        bins = []
        n = 1 # Size of windows.

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
            for i,v in enumerate(pair[1][1:]):
                if i>=len(bin_pairs[label]):
                    bin_pairs[label].append([])
                if type(v) is list:
                    bin_pairs[label][i] += v
                else:
                    bin_pairs[label][i].append(v)

        pairs_to_remove = set()
        for label, ls in bin_pairs.items():
            for i,l in enumerate(ls):
                # Only avg over first 3 items
                if i>3:
                    continue
                c = len(l)
                if c==0:
                    # Sanity check
                    # print("No Data in ",label)
                    pairs_to_remove.add(label)
                    continue
                avg = sum(l) / float(len(l))
                bin_pairs[label][i] = avg
            means = ls[-3]
            pdbs = ls[-2]
            labels = ls[-1]
            d = {}

            for i,l in enumerate(labels):
                if l not in d:
                    d[l] = {}
                d[l][pdbs[i]] = means[i]/100
            bin_pairs[label].append(d)

            pdbs_per_line = 8
            table = "<table class='table'><tr>"
            # for pdb in self.pdbs[:5]:
            #     table += "<th>{}</th>".format(pdb)
            table += "</tr>"
            for i in set(labels):
                table += "<tr>" #<th>{}</th>".format(i)
                avg = sum(d[i].values()) / float(len(d[i]))
                for ii, pdb in enumerate(self.pdbs):
                    if ii % pdbs_per_line == 0:
                        # Make a new line to prevent too big tables.
                        table += "</tr><tr>"
                        for pdb in self.pdbs[ii:ii+pdbs_per_line]:
                            table += "<th>{}</th>".format(pdb)
                        table += "</tr><tr>"

                    if pdb in d[i]:
                        fold = d[i][pdb]/avg
                        if fold>1:
                            fold = 1/fold
                        color = '#%02x%02x%02x' % (255, round(255*(fold)), round(255*(fold)))
                        table += "<td bgcolor='{}'>{}</td>".format(color,d[i][pdb])
                    else:
                        table += "<td>N/A</td>"

                table += "</tr>"
            table += "</table>"
            bin_pairs[label].append(table)
            #print(d[l])

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
        if len(self.structures)>1:
            stats_window = sorted(stats_window, key=lambda k: -k[3])
            stats_window_reduced = sorted(stats_window_reduced, key=lambda k: -k[3])
        else:
            stats_window = sorted(stats_window, key=lambda k: -k[1])
            stats_window_reduced = sorted(stats_window_reduced, key=lambda k: -k[1])
        self.stats_window = stats_window
        self.stats_window_reduced = stats_window_reduced
        self.stats_window_key = stats_window_key

    def fetch_distances(self):
        ds = Distance.objects.filter(structure__in=self.structures).all().values_list('distance','gns_pair')
        self.data = {}
        for d in ds:
            label = d[1]
            if label not in self.data:
                self.data[label] = []
            self.data[label].append(d[0]/100)

    def fetch_distances_tm(self):
        ds = Distance.objects.filter(structure__in=self.structures) \
                .exclude(gns_pair__contains='8x').exclude(gns_pair__contains='12x').exclude(gns_pair__contains='23x').exclude(gns_pair__contains='34x').exclude(gns_pair__contains='45x') \
                .values_list('distance','gns_pair')

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

    def get_distance_matrix(self):
        # common GNs
        common_gn = self.fetch_common_gns_tm()

        all_gns = sorted(list(ResidueGenericNumber.objects.filter(scheme__slug='gpcrdb')\
            .exclude(label__startswith='8x') \
            .exclude(label__startswith='12x') \
            .exclude(label__startswith='23x') \
            .exclude(label__startswith='34x') \
            .exclude(label__startswith='45x') \
            .values_list('label',flat=True)))

        pdb_distance_maps = {}
        pdb_gns = {}
        for pdb in self.pdbs:
            cache_key = "distanceMap-" + pdb

            # Cached?
            if cache.has_key(cache_key):
                cached_data = cache.get(cache_key)
                distance_map = cached_data["map"]
                structure_gn = cached_data["gns"]
            else:
                # grab raw distance data per structure
                temp = Distances()
                temp.load_pdbs([pdb])
                temp.fetch_distances_tm()

                structure_gn = list(Residue.objects.filter(protein_conformation__in=temp.pconfs) \
                    .exclude(generic_number=None) \
                    .exclude(generic_number__label__startswith='8x') \
                    .exclude(generic_number__label__startswith='12x') \
                    .exclude(generic_number__label__startswith='23x') \
                    .exclude(generic_number__label__startswith='34x') \
                    .exclude(generic_number__label__startswith='45x') \
                    .values_list('generic_number__label',flat=True))

                # create distance map
                distance_map = np.full((len(all_gns), len(all_gns)), 0.0)

                for i,res1 in enumerate(all_gns):
                    for j in range(i+1, len(all_gns)):
                        # grab average value
                        res2 = all_gns[j]
                        if res1+"_"+res2 in temp.data:
                            distance_map[i][j] = temp.data[res1+"_"+res2][0]

                # store in cache
                store = {
                    "map" : distance_map,
                    "gns" : structure_gn
                }
                cache.set(cache_key, store, 60*60*24*14)

            pdb_gns[pdb] = structure_gn

            # Filtering indices to map to common_gns
            gn_indices = np.array([ all_gns.index(residue) for residue in common_gn ])
            pdb_distance_maps[pdb] = distance_map[gn_indices,:][:, gn_indices]

            if "average" in pdb_distance_maps:
                pdb_distance_maps["average"] +=  pdb_distance_maps[pdb]/len(self.pdbs)
            else:
                pdb_distance_maps["average"] =  pdb_distance_maps[pdb]/len(self.pdbs)

        # store distance map
        for pdb in self.pdbs:
            pdb_distance_maps[pdb] = np.nan_to_num(pdb_distance_maps[pdb]/pdb_distance_maps["average"])

        # calculate distance matrix
        distance_matrix = np.full((len(self.pdbs), len(self.pdbs)), 0.0)
        for i, pdb1 in enumerate(self.pdbs):
            for j in range(i+1, len(self.pdbs)):
                pdb2 = self.pdbs[j]

                # Get common GNs between two PDBs
                common_between_pdbs = sorted(list(set(pdb_gns[pdb1]).intersection(pdb_gns[pdb2])))

                # Get common between above set and the overall common set of GNs
                common_with_query_gns = sorted(list(set(common_gn).intersection(common_between_pdbs)))

                # Get the indices of positions that are shared between two PDBs
                gn_indices = np.array([ common_gn.index(residue) for residue in common_with_query_gns ])

                # Get distance between cells that have both GNs.
                distance = np.sum(np.absolute(pdb_distance_maps[pdb1][gn_indices,:][:, gn_indices] - pdb_distance_maps[pdb2][gn_indices,:][:, gn_indices]))
                distance_matrix[i, j] = pow(distance,2)/(len(gn_indices)*len(gn_indices))
                distance_matrix[j, i] = distance_matrix[i, j]

        return distance_matrix
