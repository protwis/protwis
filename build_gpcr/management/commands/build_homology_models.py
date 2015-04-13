from django.core.management.base import BaseCommand

from protein.models import Protein
from structure.models import Structure
from protein.models import ProteinSegment
from common.alignment import Alignment
from common.models import WebLink
import pprint
from collections import OrderedDict


class Command(BaseCommand):
    
    def handle(self, *args, **options):
        Homology_model = HomologyModeling('5ht1a_human','Agonist')
        multi_alignment = Homology_model.run_pairwise_alignment()
        main_template = Homology_model.select_main_template(multi_alignment)
        main_alignment = Homology_model.run_main_alignment(Homology_model.reference_protein,main_template)
        Homology_model.run_non_conserved_switcher(main_alignment)        
        
        self.stdout.write(Homology_model.statistics, ending='')


class HomologyModeling(object):
    ''' Class to build homology models for GPCRs. 
    
        @param reference_id: str, protein id \n
        @param role: str, endogenous ligand role \n
        @param testing: boolean, test program with existing structure, default: False \n
    '''
    def __init__(self, reference_id, role, testing=False):
        self.reference_id = reference_id
        self.role = role
        self.testing = testing
        
    def __repr__(self):
        return "<{}, {}>".format(self.reference_id, self.role)
        
    def get_receptor_data(self, receptor):
        ''' Get Protein object of receptor.
        
            @param receptor: str, protein id
        '''
        self.reference_protein = Protein.objects.get(entry_name=self.reference_id)
        self.uniprot_id = self.reference_protein.accession
        self.reference_sequence = self.reference_protein.sequence
        self.statistics = CreateStatistics(self.uniprot_id)
        return self.reference_protein

    def get_structure_data(self, pdb):
        ''' Get Structure object based on pdb code.
        
            @param pdb: str, pdb code
        '''
        slist = Structure.objects.get(pdb_code_id__index=pdb)
        return slist.resolution

    def get_target_structures_data(self, role):
        ''' Get all target Structure objects based on endogenous ligand role. Returns QuerySet object.
        
            @param role: str, endogenous ligand role
        '''
        self.role = role
        if self.testing == False:
            self.structures_datatable = Structure.objects.filter(endogenous_ligand_id__role_id__role=self.role).order_by('protein_id','-resolution').distinct('protein_id')
###TO BE FIXED        
        elif self.testing == True:
            self.structures_datatable = Structure.objects.filter(endogenous_ligand_id__role_id__role=self.role).order_by('protein_id','-resolution').distinct('protein_id')
###        
        return self.structures_datatable        
        
    def get_targets_protein_data(self, structures_data):
        ''' Get all target Protein objects based on Structure objects. Returns a list of Protein objects.
        
            @param structures_data: QuerySet, query set of Structure objects. Output of get_target_structures_data function.
        '''
        plist = []
        for target in structures_data:
            plist.append(Protein.objects.get(id=target.protein_id))
        return plist
        
    def run_pairwise_alignment(self, segments='default', reference=True, calculate_similarity=True, targets=None):
        ''' Creates pairwise alignment between reference and target receptor(s).
            Returns Alignment object.
            
            @param segments: list, list of segments to use, e.g.: ['TM1','IL1','TM2','EL1'] \n
            @param reference: boolean, if True, reference receptor is used as reference, default: True.
            @param calculate_similarity: boolean, if True, call Alignment.calculate_similarity, default: True.
            @param targets: list, list of Protein objects to use as targets. By default it uses all targets with role
            specified when initializing the HomologyModeling() class.
        '''
        if not targets:
            targets = self.get_targets_protein_data(self.get_target_structures_data(self.role))
        segments = ProteinSegment.objects.filter(slug__in=segments)
        
        # core functions from alignment.py
        a = Alignment()
        if reference==True:
            a.load_reference_protein(self.get_receptor_data(self.reference_id))
        a.load_proteins(targets)
        if segments=='default':
            a.load_segments(['TM1','TM2','TM3','TM4','TM5','TM6','TM7'])
        else:
            a.load_segments(segments)
        a.build_alignment()
        if calculate_similarity==True:
            a.calculate_similarity()   
        return a
    
    def select_main_template(self, alignment):
        ''' Select main template for homology model based on highest sequence similarity. Returns Structure object of main template.
        
            @param alignment: Alignment, aligment object where the first protein is the reference. Output of pairwise_alignment function.
        '''
        self.similarity_table = OrderedDict()
        self.identity_table = OrderedDict()
        
        for i in range(1,len(alignment.proteins)):
            self.similarity_table[alignment.proteins[i]] = int(alignment.proteins[i].similarity)
            self.identity_table[alignment.proteins[i]] = int(alignment.proteins[i].identity)  
        best = [0,0]
        for key, value in self.similarity_table.items():
            if best[1]<value:            
                best = [key,value]
            elif value==0:
                pass
            elif best[1]==value:
                best+=[key,value]
        print(best)
        main_structure = self.structures_datatable.get(protein_id__entry_name=best[0])

        self.main_pdb_id = str(WebLink.objects.get(id=main_structure.pdb_code_id))[-4:]
        self.main_template_sequence = best[0].sequence
        if len(main_structure.preferred_chain)>1:
            self.main_template_preferred_chain = main_structure.preferred_chain[0]
        else:
            self.main_template_preferred_chain = main_structure.preferred_chain   
            
        self.statistics.add_info("main_template", self.main_pdb_id)
        self.statistics.add_info("preferred_chain", self.main_template_preferred_chain)
        return main_structure
        
    def run_main_alignment(self, reference, main_template):
        ''' Creates an alignment between reference (Protein object) and main_template (Structure object) 
            where matching residues are depicted with the one-letter residue code, mismatches with '.', 
            gaps with '-', gaps due to shorter sequences with 'x'. returns a AlignedReferenceAndTemplate class.
            
            @param reference: Protein object, reference receptor
            @param main_template: Structure object, main template
        '''
        main_template_protein = Protein.objects.get(id=main_template.protein_id)
        a = self.run_pairwise_alignment(['TM1','TM2','TM3','TM4','TM5','TM6','TM7',], reference=False, calculate_similarity=False, targets=[reference, main_template_protein])
        ref = a.proteins[0].alignment
        temp = a.proteins[1].alignment
        reference_string = ''
        template_string = ''
        matching_string = ''       
        reference_dict = OrderedDict()
        template_dict = OrderedDict()
        segment_count = 0

        for ref_segment, temp_segment in zip(ref,temp):
            segment_count+=1
            for ref_position, temp_position in zip(ref_segment,temp_segment):
                
                if ref_position[1]!=False and temp_position[1]!=False:
                    if ref_position[0]==temp_position[0]:
                        reference_dict[ref_position[0]]=ref_position[2]
                        template_dict[temp_position[0]]=temp_position[2]
                        reference_string+=ref_position[2]
                        template_string+=temp_position[2]
                        if ref_position[2]==temp_position[2]:
                            matching_string+=ref_position[2]
                        else:
                            matching_string+='.'
                    else:
                        print("Error: Generic numbers don't align")
                            
                elif ref_position[1]!=False and temp_position[1]==False:
                    reference_dict[ref_position[0]]=ref_position[2]                    
                    reference_string+=ref_position[2]
                    if temp_position[2]=='-':
                        template_dict[temp_position[0]]='-'
                        template_string+='-'
                        matching_string+='-'
                    elif temp_position[2]=='_':
                        template_dict[temp_position[0]]='x'
                        template_string+='x'
                        matching_string+='x'
                        
                elif ref_position[2]=='-' and temp_position[1]!=False:
                    reference_dict[ref_position[0]]='-'
                    template_dict[temp_position[0]]=temp_position[2]
                    reference_string+='-'
                    template_string+=temp_position[2]
                    matching_string+='-'
                    
                elif ref_position[2]=='-' and temp_position[2]=='-':
                    reference_dict[ref_position[0]]='-'
                    template_dict[temp_position[0]]='-'
                    reference_string+='-'
                    template_string+='-'
                    matching_string+='-'  
            reference_dict["TM"+str(segment_count)+"_end"]='/'                     
            template_dict["TM"+str(segment_count)+"_end"]='/'  
            reference_string+='/'
            template_string+='/'
            matching_string+='/'

        output = AlignedReferenceAndTemplate(self.reference_id, self.main_pdb_id, reference_dict, template_dict, matching_string)           
        return output
        
    def run_non_conserved_switcher(aligned_reference_and_template, switch_bulges=True, switch_constrictions=True):
        '''
        '''
        ref_length = 0
        for ref_res, temp_res, aligned_res in zip(aligned_reference_and_template.reference_dict, aligned_reference_and_template.template_dict, aligned_reference_and_template.aligned_string):
            if ref_res!='-' and ref_res!='/':
                ref_length+=1
            
            # bulges and constrictions
            if aligned_res=='-':
                if switch_bulges==True:
                    if ref_res=='-':
                        print(aligned_reference_and_template.reference_dict[ref_res])
                
            
class BulgesAndConstrictions(object):
    def __init__(self):
        pass
    def bulge_in_reference():
        pass                    
 
class AlignedReferenceAndTemplate(object):
    ''' Representation class for HomologyModeling.run_main_alignment() function. 
    '''
    def __init__(self, reference_id, template_id, reference_dict, template_dict, aligned_string):
        self.reference_id = reference_id
        self.template_id = template_id
        self.reference_dict = reference_dict
        self.template_dict = template_dict
        self.aligned_string = aligned_string
        
    def __repr__(self):
        return "<{}, {}>".format(self.reference_id,self.template_id)
        
class CreateStatistics(object):
    ''' Statistics dictionary for HomologyModeling.
    '''
    def __init__(self, uniprot):
        self.uniprot = uniprot
        self.info_dict = OrderedDict()
    
    def __repr__(self):
        return "<{} \n {} \n>".format(self.uniprot, self.info_dict)
    
    def add_info(self, info_name, info):
        ''' Adds new information to the statistics dictionary.
        
            @param info_name: str, info name as dictionary key
            @param info: object, any object as value
        '''
        self.info_dict[info_name] = info
        
    
        
        
        