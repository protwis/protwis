from django.core.management.base import BaseCommand, CommandError
from django.core.management import call_command
from django.conf import settings
from django.db import connection
from django.db.models import Q, F
from contactnetwork.distances import *
from protein.models import ProteinFamily

import time
import scipy

class Command(BaseCommand):

    help = "Build  distance representatives"


    def handle(self, *args, **options):
        self.receptor_representatives()

    def receptor_representatives(self):
        print('Script to decide distance representative for a state/receptor combination. Lowest average distance to all other structures for the same receptor/state')

        structures = Structure.objects.all().exclude(structure_type__slug__startswith='af-').prefetch_related(
            "pdb_code",
            "state",
            "protein_conformation__protein__parent__family")

        distinct_proteins = {}

        resolution_lookup = {}
        for s in structures:
            pdb = s.pdb_code.index
            resolution_lookup[pdb] = s.resolution
            state = s.state.slug
            slug = s.protein_conformation.protein.parent.family.slug
            name = s.protein_conformation.protein.parent.family.name

            key = '{}_{}'.format(name,state)

            if key not in distinct_proteins:
                distinct_proteins[key] = []

            distinct_proteins[key].append(pdb)

        for conformation, pdbs in distinct_proteins.items():
            print(conformation, "PDBS:",pdbs)
            number_of_pdbs = len(pdbs)
            if (number_of_pdbs==1):
                # Do not care when only one PDB for a conformation rep
                print("REPRESENTATIVE:", pdbs[0])
                s = Structure.objects.get(pdb_code__index=pdbs[0])
                s.distance_representative = True
                s.save()
            else:
                # Distances
                dis = Distances()
                dis.load_pdbs(pdbs)
                distance_matrix = dis.get_distance_matrix()

                # Calculate structures with lowest average distance (rank-by-vote fusion)
                ranking = np.zeros(len(distance_matrix))
                average = np.zeros(len(distance_matrix))
                for i in range(0,len(distance_matrix)):
                    ranking = ranking + scipy.stats.rankdata(distance_matrix[i,:], method='min')
                    average = average + distance_matrix[i,:]

                # check if single minimum
                lowest = np.where(ranking==min(ranking))[0]

                if len(lowest)>1:
                    lowest = lowest[np.where(average[lowest]==min(average))[0][0]]

                for i in range(0,len(distance_matrix)):
                    if i==lowest:
                        print("REPRESENTATIVE:",pdbs[i])
                    s = Structure.objects.get(pdb_code__index=pdbs[i])
                    s.distance_representative = (i==lowest)
                    s.save()
