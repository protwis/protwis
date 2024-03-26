from structure.models import Structure

from django.db import models
import numpy as np
import statistics


distance_scaling_factor = 10000

class InteractingResiduePair(models.Model):
    referenced_structure = models.ForeignKey('structure.Structure', on_delete=models.CASCADE)
    res1 = models.ForeignKey('residue.Residue', related_name='residue1', on_delete=models.CASCADE)
    res2 = models.ForeignKey('residue.Residue', related_name='residue2', on_delete=models.CASCADE)


    @classmethod
    def truncate(cls):
        from django.db import connection
        with connection.cursor() as cursor:
            cursor.execute('TRUNCATE TABLE "{0}" RESTART IDENTITY CASCADE'.format(cls._meta.db_table))

    def __str__(self):
        return '<{}-{}-{}>'.format(self.res1, self.res2, self.referenced_structure)

    class Meta():
        db_table = 'interacting_residue_pair'


class Interaction(models.Model):
    interacting_pair = models.ForeignKey('contactnetwork.InteractingResiduePair', on_delete=models.CASCADE)
    interaction_type = models.CharField(max_length=100)
    specific_type = models.CharField(max_length=100, null=True)

    # interaction_level -> 0 - normal definition, 1 - loosened definition
    interaction_level = models.IntegerField(null=False, default=0)
    atomname_residue1 = models.CharField(max_length=10, null=True, db_index=True)
    atomname_residue2 = models.CharField(max_length=10, null=True, db_index=True)

    @classmethod
    def truncate(cls):
        from django.db import connection
        with connection.cursor() as cursor:
            cursor.execute('TRUNCATE TABLE "{0}" RESTART IDENTITY CASCADE'.format(cls._meta.db_table))

    class Meta():
        db_table = 'interaction'


class InteractingPeptideResiduePair(models.Model):
    peptide_amino_acid_three_letter = models.CharField(max_length=3)
    peptide_amino_acid = models.CharField(max_length=1)
    peptide_sequence_number = models.IntegerField(null=False)
    peptide = models.ForeignKey('ligand.LigandPeptideStructure', related_name='peptide', on_delete=models.CASCADE)
    receptor_residue = models.ForeignKey('residue.Residue', related_name='receptor_residue', on_delete=models.CASCADE)

    def __str__(self):
        return '<{}{}-{}>'.format(self.peptide_amino_acid_three_letter, self.peptide_sequence_number, self.receptor_residue)

    class Meta():
        db_table = 'interacting_peptide_residue_pair'


class InteractionPeptide(models.Model):
    interacting_peptide_pair = models.ForeignKey('contactnetwork.InteractingPeptideResiduePair', related_name='pair', on_delete=models.CASCADE)
    peptide_atom = models.CharField(max_length=10)
    receptor_atom = models.CharField(max_length=10)
    interaction_type = models.CharField(max_length=100)
    specific_type = models.CharField(max_length=100, null=True)
    interaction_level = models.IntegerField()

    def __str__(self):
        return '<{}-{}>'

    class Meta():
        db_table = 'interaction_peptide'


class ConsensusInteraction(models.Model):
    gn1 = models.ForeignKey('residue.ResidueGenericNumber', on_delete=models.CASCADE, related_name='GN_1')
    gn2 = models.ForeignKey('residue.ResidueGenericNumber', on_delete=models.CASCADE, related_name='GN_2')
    protein_class = models.ForeignKey('protein.ProteinFamily', on_delete=models.CASCADE)
    state = models.ForeignKey('protein.ProteinState', on_delete=models.CASCADE)
    frequency = models.DecimalField(max_digits=5, decimal_places=2)
    structures = models.ManyToManyField('structure.Structure')
    proteins = models.ManyToManyField('protein.Protein')
    state_specific = models.BooleanField(default=False)

    class Meta():
        unique_together = ('gn1', 'gn2','protein_class','state')

class Distance(models.Model):
    structure = models.ForeignKey('structure.Structure', related_name='distances', on_delete=models.CASCADE, null=True)
    res1 = models.ForeignKey('residue.Residue', related_name='distance_residue1', on_delete=models.CASCADE, null=True)
    res2 = models.ForeignKey('residue.Residue', related_name='distance_residue2', on_delete=models.CASCADE, null=True)
    gn1 = models.CharField(max_length=100, null=True)
    gn2 = models.CharField(max_length=100, null=True)
    gns_pair = models.CharField(db_index=True, max_length=100, null=True)
    distance = models.IntegerField(null=True)
    distance_cb = models.IntegerField(null=True)
    distance_helix_center = models.IntegerField(null=True)

    @classmethod
    def truncate(cls):
        from django.db import connection
        with connection.cursor() as cursor:
            cursor.execute('TRUNCATE TABLE "{0}" RESTART IDENTITY CASCADE'.format(cls._meta.db_table))

    class Meta():
        db_table = 'distance'

def get_distance_averages(pdbs,s_lookup, interaction_keys,normalized = False, standard_deviation = False, split_by_amino_acid = False):
    ## Returned dataset is in ClassA GNs...
    matrix = {}
    group_distances = {}

    if len(pdbs)==1:
        # Never get SD when only looking at a single pdb...
        standard_deviation = False

    ds = list(Distance.objects.filter(structure__pdb_code__index__in=[ pdb.upper() for pdb in pdbs], gns_pair__in=interaction_keys) \
                         .values('gns_pair','distance','res1__amino_acid','res2__amino_acid','structure__pk'))
    if not normalized:
        for d in ds:
            if split_by_amino_acid:
                key = '{}{}{}'.format(d['gns_pair'],d['res1__amino_acid'],d['res2__amino_acid']).replace("_",",")
            else:
                key = d['gns_pair']
            if key not in matrix:
                matrix[key] = []
            matrix[key].append(d['distance'])

        # Calculate the average of averages
        for key,dists in matrix.items():
            if standard_deviation and len(pdbs)>1:
                if len(dists)==1:
                    stdevdists = 0
                else:
                    stdevdists = statistics.stdev(dists)
                group_distances[key] = stdevdists/distance_scaling_factor
            else:
                group_distances[key] = sum(dists)/len(dists)/distance_scaling_factor
    else:
        # NORMALIZE CODE
        for d in ds:
            if split_by_amino_acid:
                key = '{}{}{}'.format(d['gns_pair'],d['res1__amino_acid'],d['res2__amino_acid']).replace("_",",")
            else:
                key = d['gns_pair']
            dist = d['distance']
            pf = s_lookup[d['structure__pk']][2] # get the "receptor" level of the structure to group these regardless of species

            if key not in matrix:
                matrix[key] = {}
            if pf not in matrix[key]:
                matrix[key][pf] = []
            matrix[key][pf].append(dist)

        # Calculate the average of averages
        for key,pfs in matrix.items():
            # if key=='3x50_5x50':
            #     print(pfs)
            means = [sum(dists)/len(dists) for pf,dists in pfs.items() ]
            if standard_deviation and len(pdbs)>1:
                if len(means)==1:
                    stdevofmeans = 0
                else:
                    # stdevofmeans = np.std(means, ddof=1)
                    stdevofmeans = statistics.stdev(means)

                group_distances[key] = stdevofmeans/distance_scaling_factor
            else:
                meanofmeans = sum(means)/len(means)
                group_distances[key] = meanofmeans/distance_scaling_factor

    return group_distances
