from django.core.management.base import BaseCommand

from protein.models import Protein
from structure.models import Structure
from protein.models import ProteinSegment
from common.alignment import Alignment
import pprint
from collections import OrderedDict


class Command(BaseCommand):
    
    def handle(self, *args, **options):
        Homology_model = MainTemplateSelection('gp132_human','Agonist')
        alignment = Homology_model.select_main_template(Homology_model.pairwise_alignment())

        self.stdout.write(alignment, ending='')


class HomologyModeling(object):
    
    def __init__(self,receptor,role):
        self.receptor = receptor
        self.role = role
        
    def receptor_data(self, receptor):
        recep = Protein.objects.get(entry_name=self.receptor)
        return recep

    def structure_data(self, pdb):
        self.pdb = pdb
        slist = Structure.objects.get(pdb_code_id__index=pdb)
        return slist.resolution

    def target_structures_data(self, role):
        self.role = role
        slist = Structure.objects.filter(endogenous_ligand_id__role_id__role=self.role).order_by('protein_id','-resolution').distinct('protein_id')
        return slist        
        
    def targets_data(self, structures_data):
#        indeces = tlist.values('protein_id')
#        indeces_ready = []        
#        for i in indeces:
#            if i not in indeces_ready:
#                indeces_ready.append(i)
#            else:
#                indeces_ready.append(0)
#        tlist_grouped = []
#        for i in range(0,len(indeces_ready)):
#            if indeces_ready[i] != 0:
#                tlist_grouped.append(tlist[i])
        plist = []
        for target in structures_data:
            plist.append(Protein.objects.get(id=target.protein_id))

        return plist
        

        
class MainTemplateSelection(HomologyModeling):
    # Inherits __init__ from HomologyModeling class

    def pairwise_alignment(self):
        self.targets = self.targets_data(self.target_structures_data(self.role))

        self.segments7TM = ProteinSegment.objects.filter(slug__in=['TM1','TM2','TM3','TM4','TM5','TM6','TM7'])
        
        # core functions from alignment.py
        a = Alignment()
        
        a.load_reference_protein(self.receptor_data(self.receptor))
        a.load_proteins(self.targets)
        a.load_segments(self.segments7TM)
        a.build_alignment()
        a.calculate_similarity()
        
        return a
    
    def select_main_template(self,alignment):
        # main template selection
        self.similarity_table = OrderedDict()

        for i in range(1,len(alignment.proteins)):
            self.similarity_table[alignment.proteins[i]] = int(alignment.proteins[i].similarity)
        best = [0,0]
        for key, value in self.similarity_table.items():
            if best[1]<value:            
                best = [key,value]
            elif best[1]==value:
                best+=[key,value]
            
        return best[0]
        
    
        
        
        