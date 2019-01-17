from django.core.management.base import BaseCommand, CommandError
from django.core.management import call_command

from django.conf import settings
from django.db import connection
from django.utils.text import slugify
from django.db import IntegrityError
from contactnetwork.cube import compute_interactions
from contactnetwork.models import *
import contactnetwork.interaction as ci

from residue.models import ResidueGenericNumber, ResidueNumberingScheme, Residue, ResidueGenericNumberEquivalent


from structure.models import (Structure, StructureType, StructureSegment, StructureStabilizingAgent,PdbData,
    Rotamer, StructureSegmentModeling, StructureCoordinates, StructureCoordinatesDescription, StructureEngineering,
    StructureEngineeringDescription, Fragment)


import os
import yaml
from interaction.views import runcalculation,parsecalculation
from multiprocessing import Queue, Process, Value, Lock

class Command(BaseCommand):

    help = "Output all uniprot mappings"

    def prepare_input(self, proc, items, iteration=1):
        q = Queue()
        procs = list()
        num_items = len(items)
        num = Value('i', 0)
        lock = Lock()

        if not num_items:
            return False

        # make sure not to use more jobs than proteins (chunk size will be 0, which is not good)
        if proc > num_items:
            proc = num_items

        chunk_size = int(num_items / proc)
        connection.close()
        for i in range(0, proc):
            first = chunk_size * i
            if i == proc - 1:
                last = False
            else:
                last = chunk_size * (i + 1)
    
            p = Process(target=self.main_func, args=([(first, last), iteration,num,lock]))
            procs.append(p)
            p.start()

        for p in procs:
            p.join()

    def purge_contact_network(self,s):

        ii = Interaction.objects.filter(
            interacting_pair__referenced_structure=s
        ).all().delete()

        # for i in ii:
        #     i.delete()


    def build_contact_network(self,s,pdb_code):
        interacting_pairs = compute_interactions(pdb_code)

        for p in interacting_pairs:

            p.save_into_database()
            # # Create the pair
            # res1_seq_num = p.get_residue_1().id[1]
            # res2_seq_num = p.get_residue_2().id[1]
            # conformation = s.protein_conformation

            # # Get the residues
            # try:
            #     res1 = Residue.objects.get(sequence_number=res1_seq_num, protein_conformation=conformation)
            #     res2 = Residue.objects.get(sequence_number=res2_seq_num, protein_conformation=conformation)
            # except:
            #     # print('Error with pair between %s and %s (%s)' % (res1_seq_num,res2_seq_num,conformation))
            #     # print('Error with pair between %s and %s (%s)' % (res1_seq_num,res2_seq_num,conformation))
            #     continue

            # # Save the pair
            # pair = InteractingResiduePair()
            # pair.res1 = res1
            # pair.res2 = res2
            # pair.referenced_structure = s
            # pair.save()

            # # Add the interactions to the pair
            # for i in p.get_interactions():
            #     if type(i) is ci.VanDerWaalsInteraction:
            #         ni = Interaction()
            #         ni.interaction_type = 'VanDerWaals'
            #         ni.interacting_pair = pair
            #         ni.save()
            #     elif type(i) is ci.HydrophobicInteraction:
            #         ni = Interaction()
            #         ni.interaction_type = 'Hydrophobic'
            #         ni.interacting_pair = pair
            #         ni.save()
            #     elif type(i) is ci.PolarSidechainSidechainInteraction:
            #         ni = Interaction()
            #         ni.interaction_type = 'Polar'
            #         ni.interaction_type = 'PolarSidechainSidechain'
            #         ni.interacting_pair = pair
            #         ni.save()

            #         # ni.is_charged_res1 = i.is_charged_res1
            #         # ni.is_charged_res2 = i.is_charged_res2
            #     elif type(i) is ci.PolarBackboneSidechainInteraction:
            #         ni = Interaction()
            #         ni.interaction_type = 'Polar'
            #         ni.specific_type = 'PolarBackboneSidechain'
            #         ni.interacting_pair = pair
            #         ni.save()

            #         # ni.is_charged_res1 = i.is_charged_res1
            #         # ni.is_charged_res2 = i.is_charged_res2
            #         # ni.res1_is_sidechain = False
            #     elif type(i) is ci.PolarSideChainBackboneInteraction:
            #         ni = Interaction()
            #         ni.interaction_type = 'Polar'
            #         ni.specific_type = 'PolarSideChainBackbone'
            #         ni.interacting_pair = pair
            #         ni.save()

            #         # ni.is_charged_res1 = i.is_charged_res1
            #         # ni.is_charged_res2 = i.is_charged_res2
            #         # ni.res1_is_sidechain = True
            #     elif type(i) is ci.FaceToFaceInteraction:
            #         ni = Interaction()
            #         ni.interaction_type = 'Aromatic'
            #         ni.specific_type = 'FaceToFace'
            #         ni.interacting_pair = pair
            #         ni.save()
            #     elif type(i) is ci.FaceToEdgeInteraction:
            #         ni = Interaction()
            #         ni.interaction_type = 'Aromatic'
            #         ni.specific_type = 'FaceToEdge'
            #         ni.interacting_pair = pair
            #         ni.save()

            #         # ni.res1_has_face = True

            #     elif type(i) is ci.EdgeToFaceInteraction:
            #         ni = Interaction()
            #         ni.interaction_type = 'Aromatic'
            #         ni.specific_type = 'EdgeToFace'
            #         ni.interacting_pair = pair
            #         ni.save()

            #         # ni.res1_has_face = False
            #     elif type(i) is ci.PiCationInteraction:
            #         ni = Interaction()
            #         ni.interaction_type = 'Aromatic'
            #         ni.specific_type = 'PiCation'
            #         ni.interacting_pair = pair
            #         ni.save()

            #         #ni.res1_has_pi = True
            #     elif type(i) is ci.CationPiInteraction:
            #         ni = Interaction()
            #         ni.interaction_type = 'Aromatic'
            #         ni.specific_type = 'PiCation'
            #         ni.interacting_pair = pair
            #         ni.save()

            #         #ni.res1_has_pi = False

    def handle(self, *args, **options):

        self.ss = Structure.objects.filter(refined=False).all()
        self.structure_data_dir = os.sep.join([settings.DATA_DIR, 'structure_data', 'structures'])
        self.prepare_input(1, self.ss)

        # for s in Structure.objects.filter(refined=False).all():
        #   self.purge_contact_network(s)
        #   self.build_contact_network(s,s.pdb_code.index)

    def main_func(self, positions, iteration,count,lock):
        # filenames
        # if not positions[1]:
        #     filenames = self.filenames[positions[0]:]
        # else:
        #     filenames = self.filenames[positions[0]:positions[1]]
        ss = self.ss
        while count.value<len(ss):
            with lock:
                if count.value<len(ss):
                    s = ss[count.value]
                    count.value +=1
                else:
                    break 

            source_file_path = os.sep.join([self.structure_data_dir, s.pdb_code.index.upper() + ".yaml"])
            if os.path.isfile(source_file_path):
                with open(source_file_path, 'r') as f:
                    sd = yaml.load(f)
                    
            peptide_chain = ""
            if 'ligand' in sd and sd['ligand'] and sd['ligand']!='None':
                if isinstance(sd['ligand'], list):
                    ligands = sd['ligand']
                else:
                    ligands = [sd['ligand']]
                for ligand in ligands:
                    peptide_chain = ""
                    if 'chain' in ligand:
                        peptide_chain = ligand['chain']
            print(s,"Contact Network")
            self.purge_contact_network(s)
            self.build_contact_network(s,s.pdb_code.index)
            print(s,"Ligand Interactions")
            runcalculation(s.pdb_code.index,peptide_chain)
            parsecalculation(s.pdb_code.index,False)
            break