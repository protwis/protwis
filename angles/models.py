from django.db import models
from django.db.models import Avg
import math, cmath
from django.contrib.postgres.aggregates import ArrayAgg
from structure.models import Structure
import time
from scipy.stats import circmean, circstd
import numpy as np

class ResidueAngle(models.Model):
    residue             = models.ForeignKey('residue.Residue', on_delete=models.CASCADE)
    structure           = models.ForeignKey('structure.Structure', on_delete=models.CASCADE)
    a_angle             = models.FloatField(default=0, null=True)
    b_angle             = models.FloatField(default=0, null=True)
    outer_angle         = models.FloatField(default=0, null=True)
    hse                 = models.IntegerField(default=0, null=True)
    sasa                = models.FloatField(default=0, null=True)
    rsa                 = models.FloatField(default=0, null=True)
    phi                 = models.FloatField(default=0, null=True)
    psi                 = models.FloatField(default=0, null=True)
    theta               = models.FloatField(default=0, null=True)
    tau                 = models.FloatField(default=0, null=True)
    core_distance       = models.FloatField(default=0, null=True)
    midplane_distance   = models.FloatField(default=0, null=True)
    mid_distance        = models.FloatField(default=0, null=True)
    ss_dssp             = models.CharField(max_length=1, null=True)
    ss_stride           = models.CharField(max_length=1, null=True)

    class Meta():
        db_table = 'residue_angles'
        unique_together = ("residue", "structure")

def get_angle_averages(pdbs,s_lookup,normalized = False, standard_deviation = False, split_by_amino_acid = False):
    start_time = time.time()
    pdbs_upper = [pdb.upper() for pdb in pdbs]

    if len(pdbs)==1:
        # Never get SD when only looking at a single pdb...
        standard_deviation = False

    if not s_lookup:
        # Get the list of unique protein_families among pdbs
        structures = Structure.objects.filter(pdb_code__index__in=pdbs_upper
             ).select_related('protein_conformation__protein'
             ).values('pk','pdb_code__index',
                    'protein_conformation__protein__parent__entry_name',
                    'protein_conformation__protein__parent__name',
                    'protein_conformation__protein__entry_name')
        pfs = set()
        s_lookup = {}
        for s in structures:
            protein, pdb_name,pf  = [s['protein_conformation__protein__parent__entry_name'],s['protein_conformation__protein__entry_name'],s['protein_conformation__protein__parent__name']]
            s_lookup[s['pk']] = [protein, pdb_name,pf]
            pfs |= {pf}
    group_angles = {}

    ds = ResidueAngle.objects.filter(structure__pdb_code__index__in=pdbs_upper) \
                        .exclude(residue__generic_number=None) \
                        .values_list('residue__generic_number__label','structure__pk','residue__amino_acid','core_distance','a_angle','outer_angle','tau','phi','psi', 'sasa', 'rsa','theta','hse') \
                        .order_by('residue__generic_number__label','residue__amino_acid')
                        #'core_distance','a_angle','outer_angle','tau','phi','psi', 'sasa', 'rsa','theta','hse'
    custom_angles = ['a_angle', 'outer_angle', 'phi', 'psi', 'theta', 'tau']
    index_names = {0:'core_distance',1:'a_angle',2:'outer_angle',3:'tau',4:'phi',5:'psi',6: 'sasa',7: 'rsa',8:'theta',9:'hse'}
    # First bin all those belonging to same receptor
    matrix = {}
    matrix_normalized = {}
    prev_key = ''
    for d in ds:
        d = list(d)
        gn = d[0]
        aa = d[2]

        if split_by_amino_acid:
            key = "{},{}".format(gn,aa)
        else:
            key = gn

        vals = d[3:]

        if normalized:
            pf = s_lookup[d[1]][2] # get the "receptor" level of the structure to group these regardless of species
            if key != prev_key:
                matrix_normalized[key] = {}
                matrix[key] = []
                prev_key = key
            if pf not in matrix_normalized[key]:
                matrix_normalized[key][pf] = []
            matrix_normalized[key][pf].append(vals)
        else:
            if key != prev_key:
                matrix[key] = []
                prev_key = key
            matrix[key].append(vals)

    # Calculate the average of averages
    if normalized:
        for key,pfs in matrix_normalized.items():
            means = []
            for pf,dists in pfs.items():
                if len(dists)==1:
                    means.append(dists[0])
                else:
                    grouped = list(zip(*dists))
                    mean_dists = []
                    for i,L in enumerate(grouped):
                        l = list(filter(None.__ne__, list(L)))
                        if len(l)>1:
                            if index_names[i] in custom_angles:
                                mean_dists.append(radial_average(l))
                            else:
                                mean_dists.append(round(sum(l)/len(l),2))
                        elif len(l)==1:
                            mean_dists.append(round(l[0],2))
                        else:
                            mean_dists.append(None)
                    means.append(mean_dists)
            matrix[key] = means

    for key,vals in matrix.items():
        grouped = list(zip(*vals))
        group_angles[key] = []
        for i,L in enumerate(grouped):
            #remove none
            l = list(filter(None, list(L)))
            if standard_deviation:
                if len(l)>1:
                    if index_names[i] in custom_angles:
                        group_angles[key].append(radial_stddev(l))
                        # print(key,index_names[i],l,radial_stddev(l),radial_average(l))
                    else:
                        group_angles[key].append(round(np.std(l, ddof=1),2))
                elif len(l)==1:
                    group_angles[key].append(0)
                else:
                    group_angles[key].append('')
            else:
                if len(l)>1:
                    if index_names[i] in custom_angles:
                        group_angles[key].append(radial_average(l))
                    else:
                        group_angles[key].append(round(sum(l)/len(l),2))
                elif len(l)==1:
                    group_angles[key].append(round(l[0],2))
                else:
                    group_angles[key].append('')

    return group_angles

def radial_average(L):
    # r = round(math.degrees(cmath.phase(sum(cmath.rect(1, math.radians(float(d))) for d in L)/len(L))),2)
    scipy = circmean(L, -180, 180)
    return scipy

def radial_stddev(L):
    scipy = circstd(L, 360,0)
    return scipy