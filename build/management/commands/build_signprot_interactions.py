from django.core.management.base import BaseCommand, CommandError
from django.conf import settings

from protein.models import Protein
from residue.models import Residue
from structure.models import Structure
from signprot.models import SignprotInteractions
from common.definitions import AMINO_ACID_GROUP_NAMES_OLD

import logging
import os
import yaml

class Command(BaseCommand):
    help = 'Create Signal Protein interactions with GPCRs'

    logger = logging.getLogger(__name__)


    def handle(self, *args, **options):

        try:
            self.create_manual_interactions()
        except Exception as msg:
            print(msg)
            self.logger.error(msg)

    def create_manual_interactions(self):
        self.logger.info('CREATING SIGNAL PROTEIN INTERACTIONS')
        # gprotein_residues = Residue.objects.filter(protein_conformation__protein__entry_name='gnaz_human').prefetch_related('protein_segment','display_generic_number','generic_number')
        gprotein = Protein.objects.get(entry_name='gnas2_human')
        crystal = Structure.objects.get(pdb_code__index="3SN6")
        protein = crystal.protein_conformation.protein.parent


        print(gprotein,crystal,protein)

        aa_names = AMINO_ACID_GROUP_NAMES_OLD
        annotation = [  {'pos': 135, 'aa': 'I', 'gprotseg': "H5",'segment': 'TM3', 'ligand': 'Gs', 'type': aa_names['hp'], 'gpcrdb': '3.54x54', 'gpnum': 'G.H5.16', 'gpaa': 'Q384', 'availability': 'interacting'},
                        {'pos': 136, 'aa': 'T', 'gprotseg': "H5",'segment': 'TM3', 'ligand': 'Gs', 'type': 'Polar (S/T)', 'gpcrdb': '3.55x55', 'gpnum': 'G.H5.12', 'gpaa': 'R380', 'availability': 'interacting'},
                        {'pos': 139, 'aa': 'F', 'gprotseg': "H5",'segment': 'ICL2', 'ligand': 'Gs', 'type': 'Aromatic', 'gpcrdb': '34.51x51', 'gpnum': 'G.H5.08', 'gpaa': 'F376', 'availability': 'interacting'},
                        {'pos': 139, 'aa': 'F', 'gprotseg': "S1",'segment': 'ICL2', 'ligand': 'Gs', 'type': 'Aromatic', 'gpcrdb': '34.51x51', 'gpnum': 'G.S1.02', 'gpaa': 'H41', 'availability': 'interacting'},
                        {'pos': 141, 'aa': 'Y', 'gprotseg': "H5",'segment': 'ICL2', 'ligand': 'Gs', 'type': 'Aromatic', 'gpcrdb': '34.53x53', 'gpnum': 'G.H5.19', 'gpaa': 'H387', 'availability': 'interacting'},
                        {'pos': 225, 'aa': 'E', 'gprotseg': "H5",'segment': 'TM5', 'ligand': 'Gs', 'type': 'Negative charge', 'gpcrdb': '5.64x64', 'gpnum': 'G.H5.12', 'gpaa': 'R380', 'availability': 'interacting'},
                        {'pos': 225, 'aa': 'E', 'gprotseg': "H5",'segment': 'TM5', 'ligand': 'Gs', 'type': 'Negative charge', 'gpcrdb': '5.64x64', 'gpnum': 'G.H5.16', 'gpaa': 'Q384', 'availability': 'interacting'},
                        {'pos': 229, 'aa': 'Q', 'gprotseg': "H5",'segment': 'TM5', 'ligand': 'Gs', 'type': 'Polar (N/Q/H)', 'gpcrdb': '5.68x68', 'gpnum': 'G.H5.13', 'gpaa': 'D381', 'availability': 'interacting'},
                        {'pos': 229, 'aa': 'Q', 'gprotseg': "H5",'segment': 'TM5', 'ligand': 'Gs', 'type': 'Polar (N/Q/H)', 'gpcrdb': '5.68x68', 'gpnum': 'G.H5.16', 'gpaa': 'Q384', 'availability': 'interacting'},
                        {'pos': 229, 'aa': 'Q', 'gprotseg': "H5",'segment': 'TM5', 'ligand': 'Gs', 'type': 'Polar (N/Q/H)', 'gpcrdb': '5.68x68', 'gpnum': 'G.H5.17', 'gpaa': 'R385', 'availability': 'interacting'},
                        {'pos': 274, 'aa': 'T', 'gprotseg': "H5",'segment': 'TM6', 'ligand': 'Gs', 'type': 'Polar (S/T)', 'gpcrdb': '6.36x36', 'gpnum': 'G.H5.24', 'gpaa': 'E392', 'availability': 'interacting'},
                        {'pos': 328, 'aa': 'R', 'gprotseg': "H5",'segment': 'TM7', 'ligand': 'Gs', 'type': 'Positive charge', 'gpcrdb': '7.55x55', 'gpnum': 'G.H5.24', 'gpaa': 'E392', 'availability': 'interacting'}, 
                        {'pos': 232, 'aa': 'K', 'segment': 'TM5', 'ligand': 'Gs', 'type': 'Positive charge', 'gpcrdb': '5.71x71', 'gprotseg': "H5", 'gpnum': 'G.H5.13', 'gpaa': 'D381', 'availability': 'interacting'}]

        accessible_gn = ['3.50x50', '3.53x53', '3.54x54', '3.55x55', '3.56x56', '34.50x50', '34.51x51', '34.52x52', '34.53x53', '34.54x54', '34.55x55', '34.56x56', '34.57x57', '5.61x61', '5.64x64', '5.65x65', 
                         '5.66x66', '5.67x67', '5.68x68', '5.69x69', '5.71x71', '5.72x72', '5.74x74', '5.75x75', '6.25x25', '6.26x26', '6.28x28', '6.29x29', '6.32x32', '6.33x33', '6.36x36', '6.37x37', '6.40x40', 
                         '7.55x55', '7.56x56', '8.47x47', '8.48x48', '8.49x49', '8.51x51']





        SignprotInteractions.objects.filter(structure=crystal).delete()

        for a in annotation:
            p_res = Residue.objects.get(sequence_number=a['pos'],protein_conformation=crystal.protein_conformation)
            s_res = Residue.objects.get(display_generic_number__label=a['gpnum'],protein_conformation__protein=gprotein)
            SignprotInteractions.objects.create(gpcr_residue=p_res,signprot_residue=s_res,structure=crystal,interaction_type=a['type'])
    

        # aa_names = definitions.AMINO_ACID_GROUP_NAMES
        # names_aa = dict(zip(aa_names.values(),aa_names.keys()))
        # names_aa['Polar (S/T)'] = 'pol_short'
        # names_aa['Polar (N/Q/H)'] = 'pol_long'


        self.logger.info('COMPLETED CREATING SIGNAL PROTEIN INTERACTIONS')
