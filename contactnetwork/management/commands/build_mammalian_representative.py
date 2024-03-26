from django.core.management.base import BaseCommand, CommandError
from django.core.management import call_command
from django.conf import settings
from django.db import connection
from django.db.models import Q, F
from contactnetwork.models import *
from protein.models import Protein, ProteinFamily, ProteinSegment
from common.alignment import Alignment

import time
from urllib.request import urlopen


class Command(BaseCommand):

    help = "Update structures with Mammalian status and if most 'human'"


    def handle(self, *args, **options):
        self.receptor_mammal_representatives()

    def receptor_mammal_representatives(self):
        # print('Script to label structures if they are mammal, and which are the closest structure to human')

        structures = Structure.objects.all().exclude(structure_type__slug__startswith='af-').prefetch_related(
            "pdb_code",
            "state",
            "protein_conformation__protein__parent__family",
            "protein_conformation__protein__species")

        distinct_proteins = {}
        is_mammal = {}

        ## Go through all structures and deduce if mammal and prepare receptor/state sets to find most "human"
        for s in structures:
            pdb = s.pdb_code.index
            state = s.state.slug
            slug = s.protein_conformation.protein.parent.family.slug
            name = s.protein_conformation.protein.parent.family.name
            species = s.protein_conformation.protein.species.common_name
            protein = s.protein_conformation.protein.parent

            if species not in is_mammal:
                mammal = self.check_uniprot_if_mammal(protein)
                is_mammal[species] = mammal
            else:
                mammal = is_mammal[species]

            # print(species, mammal)
            s.mammal = mammal
            s.save()

            key = '{}-{}'.format(slug,state)

            if key not in distinct_proteins:
                distinct_proteins[key] = []

            distinct_proteins[key].append([pdb, species, protein, s])

        print("DEBUG", is_mammal)



        for conformation, pdbs in distinct_proteins.items():
            p_slug, state = conformation.split("-")
            number_of_pdbs = len(pdbs)
            distinct_species = set(list(x[1] for x in pdbs))
            distinct_proteins = set(list(x[2] for x in pdbs))

            if 'Human' in distinct_species:
                # Human always best..
                best_species = 'Human'
            elif len(distinct_species) == 1:
                # If only one type.. then it most be the best match
                best_species = list(distinct_species)[0]
            else:
                # There are more than 1 species, and human is not in it.. do similarity
                a = Alignment()
                ref_p = Protein.objects.get(family__slug = p_slug, species__common_name = 'Human', sequence_type__slug = 'wt')
                a.load_reference_protein(Protein.objects.get(family__slug=p_slug, species__common_name='Human', sequence_type__slug='wt'))
                a.load_proteins(distinct_proteins)
                a.load_segments(ProteinSegment.objects.filter(slug__in=['TM1', 'TM2', 'TM3', 'TM4','TM5','TM6', 'TM7']))
                a.build_alignment()
                a.calculate_similarity()
                best_species = a.proteins[1].protein.species.common_name

            ## Now that we know which species to label as most "human" go through structures and label
            for pdb, species, protein, structure in pdbs:
                most_human = False
                if species == best_species:
                    most_human = True
                structure.closest_to_human = most_human
                structure.save()

    def check_uniprot_if_mammal(self,protein):
        # Small function to see if "Mammalia" is in a UP record

        remote_file_path = 'http://www.uniprot.org/uniprot/' + protein.accession + '.txt'
        found_mammal = False
        uf = urlopen(remote_file_path)
        for raw_line in uf:
            line = raw_line.decode('UTF-8')
            # end of file
            if line.startswith('//'):
                break
            if line.startswith('OC'):
                if "Mammalia" in line:
                    found_mammal = True

        return found_mammal
