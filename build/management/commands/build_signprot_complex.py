from django.core.management.base import BaseCommand, CommandError
from django.conf import settings

from protein.models import Protein
from structure.models import Structure
from signprot.models import SignprotComplex
from structure.functions import ParseStructureCSV
from common.tools import test_model_updates

import os
import yaml
import django.apps

class Command(BaseCommand):
    help = 'Build signalling protein complex data'
    #Setting the variables for the test tracking of the model upadates
    tracker = {}
    all_models = django.apps.apps.get_models()[6:]
    test_model_updates(all_models, tracker, initialize=True)

    def add_arguments(self, parser):
        parser.add_argument('-u', '--purge',
            action='store_true',
            dest='purge',
            help='Purge signprot complex entries')
        parser.add_argument('--debug',
            default=False,
            action='store_true',
            help='Debug mode')

    def handle(self, *args, **options):
        self.options = options
        if self.options['purge']:
            print('Purging SignprotComplex model')
            SignprotComplex.objects.all().delete()
            self.tracker = {}
            test_model_updates(self.all_models, self.tracker, initialize=True)
            
        self.create_signprot_complex()
        test_model_updates(self.all_models, self.tracker, check=True)
        print(self.tracker)

    def create_signprot_complex(self):
        psc = ParseStructureCSV()
        psc.parse_g_proteins()
        psc.parse_arrestins()
        for pdb, data in psc.structures.items():
            if self.options['debug']:
                print(pdb,data)
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
