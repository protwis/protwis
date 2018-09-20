from django.core.management.base import BaseCommand, CommandError
from django.conf import settings

from protein.models import Protein
from structure.models import Structure
from signprot.models import SignprotComplex

import os
import yaml


class Command(BaseCommand):
    help = 'Build signalling protein complex data'

    signprot_complex_data_file = os.sep.join([settings.DATA_DIR, 'g_protein_data', 'complex_model_templates.yaml'])

    def handle(self, *args, **options):
        self.create_signprot_complex()

    def create_signprot_complex(self):
        with open(self.signprot_complex_data_file, 'r') as f:
            signprot_complex_data = yaml.load(f)
        for protein, data in signprot_complex_data.items():
            if type(data)==type([]):
                for i in data:
                    b_protein = Protein.objects.get(entry_name=i['beta']['protein'])
                    g_protein = Protein.objects.get(entry_name=i['gamma']['protein'])
                    structure = Structure.objects.get(pdb_code__index=i['pdb'])
                    signprot_complex, created = SignprotComplex.objects.get_or_create(protein=Protein.objects.get(entry_name=protein), 
                                                                                      structure=structure,
                                                                                      alpha=i['alpha'], beta_chain=i['beta']['chain'], gamma_chain=i['gamma']['chain'],
                                                                                      beta_protein=b_protein, gamma_protein=g_protein)
                    
                    structure.signprot_complex = signprot_complex
                    structure.save()
                    print(signprot_complex)
                    print(structure.signprot_complex)
                    
