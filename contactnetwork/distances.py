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
        self.filtered_gns = False

        # Lower half membrane GNs (all TMs)
        self.lower_membrane_full = ["1x44", "1x45", "1x46", "1x47", "1x48", "1x49", "1x50", "1x51", "1x52", "1x53", "1x54", "1x55", "1x56", "1x57", "1x58", "1x59", "1x60", "1x61", "1x62", "1x63", "1x64", "2x34", "2x35", "2x36", "2x37", "2x38", "2x39", "2x40", "2x41", "2x42", "2x43", "2x44", "2x45", "2x46", "2x47", "2x48", "2x49", "2x50", "2x51", "2x52", "3x36", "3x37", "3x38", "3x39", "3x40", "3x41", "3x42", "3x43", "3x44", "3x45", "3x46", "3x47", "3x48", "3x49", "3x50", "3x51", "3x52", "3x53", "3x54", "3x55", "3x56", "3x57", "3x58", "3x59", "3x60", "4x33", "4x34", "4x35", "4x36", "4x37", "4x38", "4x39", "4x40", "4x41", "4x42", "4x43", "4x44", "4x45", "4x46", "4x47", "4x48", "4x49", "4x491", "4x50", "4x51", "4x52", "4x53", "4x54", "5x46", "5x461", "5x47", "5x48", "5x49", "5x50", "5x51", "5x52", "5x53", "5x54", "5x55", "5x56", "5x57", "5x58", "5x59", "5x60", "5x61", "5x62", "5x63", "5x64", "5x65", "5x66", "5x67", "5x68", "5x69", "5x70", "5x71", "5x72", "5x73", "5x74", "5x75", "5x76", "5x77", "6x18", "6x19", "6x20", "6x21", "6x22", "6x23", "6x24", "6x25", "6x26", "6x27", "6x28", "6x29", "6x30", "6x31", "6x32", "6x33", "6x34", "6x35", "6x36", "6x37", "6x38", "6x39", "6x40", "6x41", "6x42", "6x43", "6x44", "6x45", "6x46", "6x461", "6x47", "6x48", "7x43", "7x44", "7x45", "7x46", "7x47", "7x48", "7x49", "7x50", "7x51", "7x52", "7x521", "7x53", "7x54", "7x55", "7x56", "7x57", "7x58", "7x59", "7x60", "7x61", "7x62", "7x63"]

        # Lower half of membrane GNs without TM1 and TM4 (only G-prot interfacing TMs)
        self.lower_membrane_gprot = ["2x34", "2x35", "2x36", "2x37", "2x38", "2x39", "2x40", "2x41", "2x42", "2x43", "2x44", "2x45", "2x46", "2x47", "2x48", "2x49", "2x50", "2x51", "2x52", "3x36", "3x37", "3x38", "3x39", "3x40", "3x41", "3x42", "3x43", "3x44", "3x45", "3x46", "3x47", "3x48", "3x49", "3x50", "3x51", "3x52", "3x53", "3x54", "3x55", "3x56", "3x57", "3x58", "3x59", "3x60", "5x46", "5x461", "5x47", "5x48", "5x49", "5x50", "5x51", "5x52", "5x53", "5x54", "5x55", "5x56", "5x57", "5x58", "5x59", "5x60", "5x61", "5x62", "5x63", "5x64", "5x65", "5x66", "5x67", "5x68", "5x69", "5x70", "5x71", "5x72", "5x73", "5x74", "5x75", "5x76", "5x77", "6x18", "6x19", "6x20", "6x21", "6x22", "6x23", "6x24", "6x25", "6x26", "6x27", "6x28", "6x29", "6x30", "6x31", "6x32", "6x33", "6x34", "6x35", "6x36", "6x37", "6x38", "6x39", "6x40", "6x41", "6x42", "6x43", "6x44", "6x45", "6x46", "6x461", "6x47", "6x48", "7x43", "7x44", "7x45", "7x46", "7x47", "7x48", "7x49", "7x50", "7x51", "7x52", "7x521", "7x53", "7x54", "7x55", "7x56", "7x57", "7x58", "7x59", "7x60", "7x61", "7x62", "7x63"]

        self.intracell = ["1x60", "2x40", "3x56", "4x39", "5x68", "6x29", "7x55"]
        self.extracell = ["1x30", "2x66", "3x23", "4x64", "5x37", "6x61", "7x31"]

        # To be replaced with automated selection of lower half TM residues
        self.filter_gns = self.lower_membrane_gprot

    def load_pdbs(self, pdbs):
        """Load a list of pdbs objects"""
        structures = Structure.objects.filter(pdb_code__index__in=pdbs).filter(refined=False).prefetch_related('protein_conformation').all()
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
                ds[i][1] = ds[i][1] / distance_scaling_factor
                ds[i][2] = ds[i][2] / distance_scaling_factor
                ds[i][5] = [x / distance_scaling_factor for x in ds[i][5]]
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

        if self.filtered_gns:
            res = res.filter(generic_number__label__in=self.filter_gns)

        common_gn = sorted(res)

        return common_gn

    def fetch_conserved_gns_tm(self):
        res = Residue.objects.filter(protein_conformation__in=self.pconfs) \
                        .exclude(generic_number=None) \
                        .exclude(generic_number__label__contains='8x') \
                        .exclude(generic_number__label__contains='12x') \
                        .exclude(generic_number__label__contains='23x') \
                        .exclude(generic_number__label__contains='34x') \
                        .exclude(generic_number__label__contains='45x') \
                        .values('generic_number__label') \
                        .annotate(c = Count('generic_number__label')).filter(c=len(self.structures)).values_list('generic_number__label',flat=True).order_by('c')

        conserved_gn = sorted(res)

        return conserved_gn

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
                d[l][pdbs[i]] = means[i]/distance_scaling_factor
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
            self.data[label].append(d[0]/distance_scaling_factor)

    def fetch_distances_tm(self, distance_type = "CA"):
#                .filter(gn1__in=self.filter_gns).filter(gn2__in=self.filter_gns) \
        # ds = Distance.objects.filter(structure__in=self.structures) \
        #         .exclude(gns_pair__contains='8x').exclude(gns_pair__contains='12x').exclude(gns_pair__contains='23x').exclude(gns_pair__contains='34x').exclude(gns_pair__contains='45x')

        if distance_type == "HC":
            ds = Distance.objects.filter(structure__in=self.structures) \
                .exclude(gns_pair__contains='8x').exclude(gns_pair__contains='12x').exclude(gns_pair__contains='23x').exclude(gns_pair__contains='34x').exclude(gns_pair__contains='45x') \
                .values_list('distance_helix_center','gns_pair')
        elif distance_type == "CB":
            ds = Distance.objects.filter(structure__in=self.structures) \
                .exclude(gns_pair__contains='8x').exclude(gns_pair__contains='12x').exclude(gns_pair__contains='23x').exclude(gns_pair__contains='34x').exclude(gns_pair__contains='45x') \
                .values_list('distance_cb','gns_pair')
        else:
            ds = Distance.objects.filter(structure__in=self.structures) \
                .exclude(gns_pair__contains='8x').exclude(gns_pair__contains='12x').exclude(gns_pair__contains='23x').exclude(gns_pair__contains='34x').exclude(gns_pair__contains='45x') \
                .values_list('distance','gns_pair')

        if self.filtered_gns:
            ds = ds.filter(gn1__in=self.filter_gns).filter(gn2__in=self.filter_gns) \

        self.data = {}
        for d in ds:
            label = d[1]
            if label not in self.data:
                self.data[label] = []
            self.data[label].append(d[0]/distance_scaling_factor)

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

    def get_distance_matrix(self, normalize = True, cache_enabled = True):
        # common GNs
        common_gn = self.fetch_common_gns_tm()

#            .filter(label__in=self.filter_gns) \
        all_gns = ResidueGenericNumber.objects.filter(scheme__slug='gpcrdb')\
            .exclude(label__startswith='8x') \
            .exclude(label__startswith='12x') \
            .exclude(label__startswith='23x') \
            .exclude(label__startswith='34x') \
            .exclude(label__startswith='45x') \
            .values_list('label',flat=True)

        if self.filtered_gns:
            all_gns = all_gns.filter(label__in=self.filter_gns)

        all_gns = sorted(list(all_gns))

        pdb_distance_maps = {}
        pdb_gns = {}
        for pdb in self.pdbs:
            cache_key = "distanceMap-" + pdb

            # Cached?
            if cache.has_key(cache_key) and cache_enabled:
                cached_data = cache.get(cache_key)
                distance_map = cached_data["map"]
                structure_gn = cached_data["gns"]
            else:
                # grab raw distance data per structure
                temp = Distances()
                temp.load_pdbs([pdb])
                temp.fetch_distances_tm()

#                    .filter(generic_number__label__in=self.filter_gns) \
                structure_gn = Residue.objects.filter(protein_conformation__in=temp.pconfs) \
                    .exclude(generic_number=None) \
                    .exclude(generic_number__label__startswith='8x') \
                    .exclude(generic_number__label__startswith='12x') \
                    .exclude(generic_number__label__startswith='23x') \
                    .exclude(generic_number__label__startswith='34x') \
                    .exclude(generic_number__label__startswith='45x') \
                    .values_list('generic_number__label',flat=True)

                if self.filtered_gns:
                    structure_gn = structure_gn.filter(generic_number__label__in=self.filter_gns)

                structure_gn = list(structure_gn)

                # create distance map
                distance_map = np.full((len(all_gns), len(all_gns)), 0.0)

                for i,res1 in enumerate(all_gns):
                    for j in range(i+1, len(all_gns)):
                        # grab average value
                        res2 = all_gns[j]
                        if res1+"_"+res2 in temp.data:
                            distance_map[i][j] = temp.data[res1+"_"+res2][0]

                # store in cache
                if cache_enabled:
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
        if normalize:
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
                distance_matrix[i, j] = distance_matrix[j, i] = distance * distance/(len(gn_indices)*len(gn_indices))

        return distance_matrix
