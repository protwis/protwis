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
        ).all()

        for i in ii:
            i.delete()

    def build_contact_network(self,s,pdb_code):
        interacting_pairs = compute_interactions(pdb_code)

        for p in interacting_pairs:
            # Create the pair
            res1_seq_num = p.get_residue_1().id[1]
            res2_seq_num = p.get_residue_2().id[1]
            conformation = s.protein_conformation

            # Get the residues
            try:
                res1 = Residue.objects.get(sequence_number=res1_seq_num, protein_conformation=conformation)
                res2 = Residue.objects.get(sequence_number=res2_seq_num, protein_conformation=conformation)
            except:
                # print('Error with pair between %s and %s (%s)' % (res1_seq_num,res2_seq_num,conformation))
                # print('Error with pair between %s and %s (%s)' % (res1_seq_num,res2_seq_num,conformation))
                continue

            # Save the pair
            pair = InteractingResiduePair()
            pair.res1 = res1
            pair.res2 = res2
            pair.referenced_structure = s
            pair.save()

            # Add the interactions to the pair
            for i in p.get_interactions():
                if type(i) is ci.VanDerWaalsInteraction:
                    ni = VanDerWaalsInteraction()
                    ni.interacting_pair = pair
                    ni.save()
                elif type(i) is ci.HydrophobicInteraction:
                    ni = HydrophobicInteraction()
                    ni.interacting_pair = pair
                    ni.save()
                elif type(i) is ci.PolarSidechainSidechainInteraction:
                    ni = PolarSidechainSidechainInteraction()
                    ni.interacting_pair = pair
                    ni.is_charged_res1 = i.is_charged_res1
                    ni.is_charged_res2 = i.is_charged_res2
                    ni.save()
                elif type(i) is ci.PolarBackboneSidechainInteraction:
                    ni = PolarBackboneSidechainInteraction()
                    ni.interacting_pair = pair
                    ni.is_charged_res1 = i.is_charged_res1
                    ni.is_charged_res2 = i.is_charged_res2
                    ni.res1_is_sidechain = False
                    ni.save()
                elif type(i) is ci.PolarSideChainBackboneInteraction:
                    ni = PolarBackboneSidechainInteraction()
                    ni.interacting_pair = pair
                    ni.is_charged_res1 = i.is_charged_res1
                    ni.is_charged_res2 = i.is_charged_res2
                    ni.res1_is_sidechain = True
                    ni.save()
                elif type(i) is ci.FaceToFaceInteraction:
                    ni = FaceToFaceInteraction()
                    ni.interacting_pair = pair
                    ni.save()
                elif type(i) is ci.FaceToEdgeInteraction:
                    ni = FaceToEdgeInteraction()
                    ni.interacting_pair = pair
                    ni.res1_has_face = True
                    ni.save()
                elif type(i) is ci.EdgeToFaceInteraction:
                    ni = FaceToEdgeInteraction()
                    ni.interacting_pair = pair
                    ni.res1_has_face = False
                    ni.save()
                elif type(i) is ci.PiCationInteraction:
                    ni = PiCationInteraction()
                    ni.interacting_pair = pair
                    ni.res1_has_pi = True
                    ni.save()
                elif type(i) is ci.CationPiInteraction:
                    ni = PiCationInteraction()
                    ni.interacting_pair = pair
                    ni.res1_has_pi = False
                    ni.save()

    def handle(self, *args, **options):

        self.ss = Structure.objects.filter(refined=False).all()
        self.prepare_input(16, self.ss)

        # for s in Structure.objects.filter(refined=False).all():
        #     print(s,s.pdb_code.index)
        #     self.purge_contact_network(s)
        #     self.build_contact_network(s,s.pdb_code.index)

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
            print(s)
            self.purge_contact_network(s)
            self.build_contact_network(s,s.pdb_code.index)