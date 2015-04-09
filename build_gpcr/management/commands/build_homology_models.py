from django.core.management.base import BaseCommand

from protein.models import Protein
from structure.models import Structure
from protein.models import ProteinSegment
from common.alignment import Alignment

class Command(BaseCommand):
    
    def handle(self, *args, **options):
        Homology_model = HomologyModeling('5ht2b_human','Agonist')
        Homology_model.main_template_select()
        self.stdout.write('Done', ending='')


class HomologyModeling(object):
    
    def __init__(self,receptor,role):
        self.receptor = receptor
        self.role = role
        
    def receptor_data(self):
        recep = Protein.objects.get(entry_name=self.receptor)
        return recep
        
    def targets_data(self):
        tlist = Structure.objects.filter(endogenous_ligand_id__role_id__role=self.role)
        plist = []
        for target in tlist:
            plist.append(Protein.objects.get(id=target.protein_id))
        return plist
        
    def structure_data(self, pdb):
        self.pdb = pdb
        slist = Structure.objects.get(pdb_code_id__index=pdb)
        return slist.resolution
        

    def main_template_select(self):
        protein = self.targets_data()
        
        a = Alignment()
        segments = ProteinSegment.objects.filter(slug__in=['TM1','TM2','TM3','TM4','TM5','TM6','TM7'])
        a.load_proteins(protein)
        a.load_positions(segments)
        