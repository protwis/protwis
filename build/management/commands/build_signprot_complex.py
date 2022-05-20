from django.core.management.base import BaseCommand, CommandError
from django.conf import settings

from protein.models import Protein
from structure.models import Structure
from signprot.models import SignprotComplex
from structure.functions import ParseStructureCSV

import os
import yaml


class Command(BaseCommand):
    help = 'Build signalling protein complex data'


    def add_arguments(self, parser):
        parser.add_argument('-u', '--purge',
            action='store_true',
            dest='purge',
            help='Purge signprot complex entries')

    def handle(self, *args, **options):
        if options['purge']:
            SignprotComplex.objects.all().delete()
        self.create_signprot_complex()

    def create_signprot_complex(self):
        psc = ParseStructureCSV()
        psc.parse_g_proteins()
        psc.parse_arrestins()
        for pdb, data in psc.structures.items():
            if 'arrestin' in data or 'g_protein' in data:
                structure = Structure.objects.get(pdb_code__index=pdb)

                if 'g_protein' in data:
                    if data['g_protein']['beta_uniprot']=='':
                        b_protein = None
                        b_chain = None
                    else:
                        b_protein = Protein.objects.get(entry_name=data['g_protein']['beta_uniprot'])
                        b_chain = data['g_protein']['beta_chain']
                    if data['g_protein']['gamma_uniprot']=='':
                        g_protein = None
                        g_chain = None
                    else:
                        g_protein = Protein.objects.get(entry_name=data['g_protein']['gamma_uniprot'])
                        g_chain = data['g_protein']['gamma_chain']
                    
                    
                    signprot_complex, created = SignprotComplex.objects.get_or_create(protein=Protein.objects.get(entry_name=data['g_protein']['alpha_uniprot']),
                                                                                      structure=structure,
                                                                                      alpha=data['g_protein']['alpha_chain'], beta_chain=b_chain, gamma_chain=g_chain,
                                                                                      beta_protein=b_protein, gamma_protein=g_protein)
                if 'arrestin' in data:
                    signprot_complex, created = SignprotComplex.objects.get_or_create(protein=Protein.objects.get(entry_name=data['arrestin']['protein']), structure=structure,
                                                                                      alpha=data['arrestin']['chain'])
                structure.signprot_complex = signprot_complex
                structure.save()
