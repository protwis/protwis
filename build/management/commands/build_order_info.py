from django.core.management.base import BaseCommand, CommandError
import logging

class Command(BaseCommand):
    
    table_population_order = {
        'web_resources': 'build_basic_webresources',
        'publication': 'build_publications',
        'ligand': 'build_ligands',
        'protein': 'build_proteins',
        'residue': 'build_residues',
        'structure': 'build_structures',
        'fragment': 'build_fragments',
        }

    def handle(self, *args, **options):
        print("The order of populating tables:\ntable_name ==> command")
        print('\n'.join(['\t\t'.join([x,y]) for [x,y] in self.table_population_order.items()]))


