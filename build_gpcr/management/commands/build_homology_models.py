from django.core.management.base import BaseCommand

from protein.models import Protein, ProteinSegment, ProteinConformation
from residue.models import Residue
from structure.models import Structure
from common.alignment import Alignment
from django.conf import settings

import Bio.PDB as PDB
from collections import OrderedDict
import pprint

class Command(BaseCommand):
    
    def handle(self, *args, **options):
#        count=0
#        s = [struct.protein_conformation.protein.parent.entry_name for struct in Structure.objects.all()]
#        for protein in ProteinConformation.objects.all():
#            if protein.protein.entry_name not in s and count < 5:
#                Homology_model = HomologyModeling(protein.protein.entry_name, 'Inactive', ['Inactive'])
#                multi_alignment = Homology_model.run_pairwise_alignment()
#                main_template = Homology_model.select_main_template(multi_alignment)
#                main_alignment = Homology_model.run_main_alignment(alignment=multi_alignment)
#                non_conserved_alignment = Homology_model.run_non_conserved_switcher(main_alignment)
#    
#                self.stdout.write(Homology_model.statistics, ending='')
#                count+=1
                
        Homology_model = HomologyModeling('gp132_human', 'Inactive', ['Inactive'])
        multi_alignment = Homology_model.run_pairwise_alignment()
        main_template = Homology_model.select_main_template(multi_alignment)
        main_alignment = Homology_model.run_main_alignment(alignment=multi_alignment)
        non_conserved_alignment = Homology_model.run_non_conserved_switcher(main_alignment)
    
        self.stdout.write(Homology_model.statistics, ending='')

class HomologyModeling(object):
    ''' Class to build homology models for GPCRs. 
    
        @param reference_entry_name: str, protein entry name \n
        @param state: str, endogenous ligand state of reference \n
        @param query_states: list, list of endogenous ligand states to be applied for template search, 
        default: same as reference
    '''
    def __init__(self, reference_entry_name, state, query_states):
        self.reference_entry_name = reference_entry_name
        self.state = state
        self.query_states = query_states
        self.statistics = CreateStatistics(self.reference_entry_name)
        
    def __repr__(self):
        return "<{}, {}>".format(self.reference_entry_name, self.state)

    def get_structure_queryset(self, query_states):
        ''' Get all target Structure objects based on endogenous ligand state. Returns QuerySet object.
        
            @param query_states: list, list of endogenous ligand states to be applied for template search, 
            default: same as reference 
        '''
        return Structure.objects.filter(state__name__in=query_states).order_by('protein_conformation__protein__parent',
                                                                            'resolution').distinct('protein_conformation__protein__parent')
               
    def get_protein_objects(self, structures_data):
        ''' Get all target Protein objects based on Structure objects. Returns a list of Protein objects.
        
            @param structures_data: QuerySet, query set of Structure objects. Output of get_structure_queryset 
            function.
        '''
        return [Protein.objects.get(id=target.protein_conformation.protein.parent.id) for target in structures_data]
        
    def run_pairwise_alignment(self, segments=['TM1','TM2','TM3','TM4','TM5','TM6','TM7'], order_by='similarity', 
                               reference=True, calculate_similarity=True, targets=None):
        ''' Creates pairwise alignment between reference and target receptor(s).
            Returns Alignment object.
            
            @param segments: list, list of segments to use, e.g.: ['TM1','IL1','TM2','EL1'] \n
            @param reference: boolean, if True, reference receptor is used as reference, default: True.
            @param calculate_similarity: boolean, if True, call Alignment.calculate_similarity, default: True.
            @param targets: list, list of Protein objects to use as targets. By default it uses all targets with state
            specified when initializing the HomologyModeling() class.
        '''
        self.segments = segments
        if not targets:
            targets = self.get_protein_objects(self.get_structure_queryset(self.query_states))       
        
        # core functions from alignment.py
        
        a = Alignment()
        a.order_by = order_by
        
        if reference==True:
            self.reference_protein = Protein.objects.get(entry_name=self.reference_entry_name)
            self.uniprot_id = self.reference_protein.accession
            self.reference_sequence = self.reference_protein.sequence
            self.statistics.add_info('uniprot_id',self.uniprot_id)
            a.load_reference_protein(self.reference_protein)
        a.load_proteins(targets)
        self.segments = ProteinSegment.objects.filter(slug__in=self.segments)
        a.load_segments(self.segments)
        a.build_alignment()
        if calculate_similarity==True:
            a.calculate_similarity()
        
        return a
    
    def select_main_template(self, alignment):
        ''' Select main template for homology model based on highest sequence similarity. Returns Structure object of 
        main template.
        
            @param alignment: Alignment, aligment object where the first protein is the reference. Output of 
            pairwise_alignment function.
        '''
        self.similarity_table = self.create_similarity_table(alignment, self.get_structure_queryset(self.query_states))
        
        main_structure = list(self.similarity_table.items())[0][0]
        
        main_structure_protein = Protein.objects.get(id=main_structure.protein_conformation.protein.parent.id)
        
        self.main_pdb_id = main_structure.pdb_code.index 
        self.main_template_sequence = main_structure_protein.sequence
        if len(main_structure.preferred_chain)>1:
            self.main_template_preferred_chain = main_structure.preferred_chain[0]
        else:
            self.main_template_preferred_chain = main_structure.preferred_chain   
        
        self.statistics.add_info("main_template", self.main_pdb_id)
        self.statistics.add_info("preferred_chain", self.main_template_preferred_chain)
        
        return main_structure
        
    def run_main_alignment(self, alignment=None, reference=None, main_template=None, 
                           segments=['TM1','TM2','TM3','TM4','TM5','TM6','TM7']):
        ''' Creates an alignment between reference (Protein object) and main_template (Structure object) 
            where matching residues are depicted with the one-letter residue code, mismatches with '.', 
            gaps with '-', gaps due to shorter sequences with 'x'. returns a AlignedReferenceAndTemplate class.
            
            @param alignment: Alignment, output of run_pairwise_alignment. \n
            @param reference: Protein object, reference receptor. \n
            @param main_template: Structure object, main template. \n
            @param segments: list, list of segments to use, e.g.: ['TM1','IL1','TM2','EL1'].
        '''
        if alignment==None and reference!=None and main_template!=None:
            main_template_protein = Protein.objects.get(id=main_template.protein_conformation.protein.parent.id)        
            self.segments = segments
            a = self.run_pairwise_alignment(segments=self.segments, reference=False, 
                                                calculate_similarity=False, targets=[reference, main_template_protein])
        else:
            a = alignment
        ref = a.proteins[0].alignment
        temp = a.proteins[1].alignment
        reference_string, template_string, matching_string = '','',''   
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
                    
            reference_dict["TM"+str(segment_count)+"_end"]='/'                     
            template_dict["TM"+str(segment_count)+"_end"]='/'  
            reference_string+='/'
            template_string+='/'
            matching_string+='/'

        main_alignment = AlignedReferenceAndTemplate(self.reference_entry_name, self.main_pdb_id, 
                                                     reference_dict, template_dict, matching_string)           
        return main_alignment
        
    def run_non_conserved_switcher(self, alignment_input, switch_bulges=True, switch_constrictions=True):
        ''' Function to identify and switch non conserved residues in the alignment. Optionally,
            it can identify and switch bulge and constriction sites too. 
            
            @param alignment_input: AlignedReferenceAndTemplate, alignment of reference and main template with 
            alignment string. \n
            @param switch_bulges: boolean, identify and switch bulge sites. Default = True.
            @param switch_constrictions: boolean, identify and switch constriction sites. Default = True.
        '''
        ref_length = 0
        non_cons_count = 0
        matched_count = 0
        ref_bulge_list, temp_bulge_list, ref_const_list, temp_const_list = [],[],[],[]
        parse = GPCRDBParsingPDB()
        Bulge = Bulges()
        Const = Constrictions()
        main_pdb_array = parse.pdb_array_creator('./structure/PDB/{}_{}_GPCRDB.pdb'.format(self.main_pdb_id,self.main_template_preferred_chain))
        
        # bulges and constrictions
        if switch_bulges==True or switch_constrictions==True:
            structure_table_all = self.get_structure_queryset(['Inactive','Active'])
            alignment = self.run_pairwise_alignment(targets=self.get_protein_objects(structure_table_all))
            self.similarity_table_all = self.create_similarity_table(alignment, structure_table_all)
            for ref_res, temp_res, aligned_res in zip(alignment_input.reference_dict, alignment_input.template_dict, 
                                                      alignment_input.aligned_string):
                gn = ref_res
                gn_TM = parse.gn_num_extract(gn, 'x')[0]
                gn_num = parse.gn_num_extract(gn, 'x')[1]
                
                if aligned_res=='-':
                    if (alignment_input.reference_dict[ref_res]=='-' and 
                        alignment_input.reference_dict[parse.gn_indecer(gn,'x',-1)] not in 
                        ['-','/'] and alignment_input.reference_dict[parse.gn_indecer(gn,'x',+1)] not in ['-','/']): 
    
                        # bulge in template
                        if len(str(gn_num))==3:
                            if switch_bulges==True:
                                print('bulge_in_template')
                                bulge_template = Bulge.find_bulge_template(gn, self.similarity_table_all, 
                                                                           bulge_in_reference=False)
                                temp_bulge_list.append({gn:Bulge.template})
                                
                        # constriction in reference
                        else:
                            if switch_constrictions==True: 
                                print('const_in_ref')
                                constriction_template = Const.find_constriction_template(gn, self.similarity_table_all,
                                                                                    constriction_in_reference=True)
                                ref_const_list.append({parse.gn_indecer(gn, 'x', -1)+'-'+parse.gn_indecer(gn, 
                                                                                            'x', +1):Const.template})
                               
                    elif (alignment_input.template_dict[temp_res]=='-' and 
                          alignment_input.template_dict[parse.gn_indecer(gn,'x',-1)] not in 
                          ['-','/'] and alignment_input.template_dict[parse.gn_indecer(gn,'x',+1)] not in ['-','/']): 
                        
                        # bulge in reference
                        if len(str(gn_num))==3:
                            if switch_bulges==True:
                                print('bulge_in_ref')
                                bulge_template = Bulge.find_bulge_template(gn, self.similarity_table_all,
                                                                           bulge_in_reference=True)
                                ref_bulge_list.append({gn:Bulge.template})
                                               
                        # constriction in template
                        else:
                            if switch_constrictions==True:
                                print('const_in_temp')
                                constriction_template = Const.find_constriction_template(gn, self.similarity_table_all,
                                                                                    constriction_in_reference=False)
                                temp_const_list.append({parse.gn_indecer(gn, 'x', -1)+'-'+parse.gn_indecer(gn, 
                                                                                            'x', +1):Const.template})
                                
            self.statistics.add_info('reference_bulges', ref_bulge_list)
            self.statistics.add_info('template_bulges', temp_bulge_list)
            self.statistics.add_info('reference_constrictions', ref_const_list)
            self.statistics.add_info('template_constrictions', temp_const_list)
        
        # non-conserved residues
        non_cons_res_templates = []
        for ref_res, temp_res, aligned_res in zip(alignment_input.reference_dict, alignment_input.template_dict, 
                                                  alignment_input.aligned_string):
            if ref_res!='-' and ref_res!='/':
                ref_length+=1   
            
            gn = ref_res
            gn_TM = parse.gn_num_extract(gn, 'x')[0]
            gn_num = parse.gn_num_extract(gn, 'x')[1]
            
            if aligned_res=='.':
                non_cons_count+=1
                residues = Residue.objects.filter(generic_number__label=ref_res)
                proteins_w_this_gn = [res.protein_conformation.protein.parent for res in 
                                        residues if res.amino_acid==alignment_input.reference_dict[ref_res]]
                for struct in self.similarity_table:
                    if struct.protein_conformation.protein.parent in proteins_w_this_gn:
                        non_cons_res_templates.append({ref_res:struct})
                        matched_count+=1
                        break
                    else:
                        try:
                            main_pdb_array[ref_res.replace('x','.')] = OrderedDict([('N',main_pdb_array[ref_res.replace('x','.')]['N']),
                                                                        ('CA',main_pdb_array[ref_res.replace('x','.')]['CA']),
                                                                        ('C',main_pdb_array[ref_res.replace('x','.')]['C']),
                                                                        ('O',main_pdb_array[ref_res.replace('x','.')]['O'])])
                        except:
                            pass
        self.statistics.add_info('non_conserved_residue_templates', non_cons_res_templates)
        
    def create_similarity_table(self, alignment, structures_datatable):
        ''' Creates an ordered dictionary, where templates are sorted by similarity score.
        
            @param alignment: Alignment, alignment of sequences. Output of run_pairwise_alignment \n
            @param structures_datatable: QuerySet, involved template structures. Output of get_structure_queryset.
        '''
        similarity_table = OrderedDict()
        
        for protein in alignment.proteins:
            if protein.protein!=self.reference_protein:
                matches = structures_datatable.order_by('protein_conformation__protein__parent',
                    'resolution').distinct('protein_conformation__protein__parent').filter(protein_conformation__protein__parent__id=protein.protein.id)    
                similarity_table[list(matches)[0]] = int(protein.similarity)
        
        return similarity_table

class Bulges():
    ''' Class to handle bulges in GPCR structures.
    '''
    def __init__(self):
        pass
    
    def find_bulge_template(self, gn, similarity_table, bulge_in_reference):
        ''' Searches for bulge template, returns residues of template (5 residues if the bulge is in the reference, 4
            residues if the bulge is in the template). 
            
            @param gn: str, Generic number of bulge, e.g. 1x411 \n
            @param similarity_table: OrderedDict(), table of structures ordered by preference.
            Output of HomologyModeling().create_similarity_table(). \n
            @param bulge_in_reference: boolean, Set it to True if the bulge is in the reference, set it to False if the
            bulge is in the template.
        '''
        bulge_templates = []
        parse = GPCRDBParsingPDB()
        if bulge_in_reference==True:
            matches = Residue.objects.filter(generic_number__label=gn)
        elif bulge_in_reference==False:
            excludees = Residue.objects.filter(generic_number__label=gn)
            excludee_proteins = list(OrderedDict.fromkeys([res.protein_conformation.protein.parent.entry_name for res in excludees]))
            matches = Residue.objects.filter(generic_number__label=gn[:-1])
        for structure, value in similarity_table.items():  
            protein_object = Protein.objects.get(id=structure.protein_conformation.protein.parent.id)
            try:                            
                for match in matches:
                    if bulge_in_reference==True:
                        if match.protein_conformation.protein.parent==protein_object:
                            bulge_templates.append(structure)
                    elif bulge_in_reference==False:
                        if match.protein_conformation.protein.parent==protein_object and match.protein_conformation.protein.parent.entry_name not in excludee_proteins:
                            bulge_templates.append(structure)
            except:
                pass
        for temp in bulge_templates:
            try:
                if bulge_in_reference==True:
                    alt_bulge = parse.fetch_residues_from_pdb(temp, str(temp.preferred_chain)[0], 
                                                              [parse.gn_indecer(gn,'x',-2),
                                                               parse.gn_indecer(gn,'x',-1),gn,
                                                               parse.gn_indecer(gn,'x',+1),
                                                               parse.gn_indecer(gn,'x',+2)])
                elif bulge_in_reference==False:
                    alt_bulge = parse.fetch_residues_from_pdb(temp, str(temp.preferred_chain)[0], 
                                                              [parse.gn_indecer(gn,'x',-2),
                                                               parse.gn_indecer(gn,'x',-1),
                                                               parse.gn_indecer(gn,'x',+1),
                                                               parse.gn_indecer(gn,'x',+2)])
                self.template = temp              
                break
            except:
                self.template = None               
        try:
            return alt_bulge
        except:
            return None
            
class Constrictions():
    ''' Class to handle constrictions in GPCRs.
    '''
    def __init__(self):
        pass
    
    def find_constriction_template(self, gn, similarity_table, constriction_in_reference):
        ''' Searches for constriction template, returns residues of template (4 residues if the constriction is in the 
            reference, 5 residues if the constriction is in the template). 
            
            @param gn: str, Generic number of constriction, e.g. 7x44 \n
            @param similarity_table: OrderedDict(), table of structures ordered by preference.
            Output of HomologyModeling().create_similarity_table(). \n
            @param constriction_in_reference: boolean, Set it to True if the constriction is in the reference, set it 
            to False if the constriction is in the template.
        '''        
        constriction_templates = []
        parse = GPCRDBParsingPDB()
        if constriction_in_reference==True:
            excludees = Residue.objects.filter(generic_number__label=gn)
            excludee_proteins = list(OrderedDict.fromkeys([res.protein_conformation.protein.parent.entry_name for res in excludees]))
            matches = Residue.objects.filter(generic_number__label=parse.gn_indecer(gn,'x',-1))
        elif constriction_in_reference==False:
            matches = Residue.objects.filter(generic_number__label=gn)
        for structure, value in similarity_table.items():  
            protein_object = Protein.objects.get(id=structure.protein_conformation.protein.parent.id)
            try:                            
                for match in matches:
                    if constriction_in_reference==True:                        
                        if match.protein_conformation.protein.parent==protein_object and match.protein_conformation.protein.parent.entry_name not in excludee_proteins:
                            constriction_templates.append(structure)
                    elif constriction_in_reference==False:
                        if match.protein_conformation.protein.parent==protein_object:
                            constriction_templates.append(structure)
            except:
                pass
        for temp in constriction_templates:
            try:
                if constriction_in_reference==True:
                    alt_bulge = parse.fetch_residues_from_pdb(temp, str(temp.preferred_chain)[0], 
                                                              [parse.gn_indecer(gn,'x',-2),
                                                               parse.gn_indecer(gn,'x',-1),
                                                               parse.gn_indecer(gn,'x',+1),
                                                               parse.gn_indecer(gn,'x',+2)])
                elif constriction_in_reference==False:
                    alt_bulge = parse.fetch_residues_from_pdb(temp, str(temp.preferred_chain)[0], 
                                                              [parse.gn_indecer(gn,'x',-2),
                                                               parse.gn_indecer(gn,'x',-1),gn,
                                                               parse.gn_indecer(gn,'x',+1),
                                                               parse.gn_indecer(gn,'x',+2)])
                self.template = temp              
                break
            except:
                self.template = None               
        try:
            return alt_bulge
        except:
            return None

class AlignedReferenceAndTemplate(object):
    ''' Representation class for HomologyModeling.run_main_alignment() function. 
    '''
    def __init__(self, reference_entry_name, template_id, reference_dict, template_dict, aligned_string):
        self.reference_entry_name = reference_entry_name
        self.template_id = template_id
        self.reference_dict = reference_dict
        self.template_dict = template_dict
        self.aligned_string = aligned_string
        
    def __repr__(self):
        return "<{}, {}>".format(self.reference_entry_name,self.template_id)
        
class GPCRDBParsingPDB(object):
    ''' Class to manipulate cleaned pdb files of GPCRs.
    '''
    def __init__(self):
        pass
    
    def gn_num_extract(self, gn, delimiter):
        ''' Extract TM number and position for formatting.
        
            @param gn: str, Generic number \n
            @param delimiter: str, character between TM and position (usually 'x')
        '''
        try:
            split = gn.split(delimiter)
            return int(split[0]), int(split[1])
        except:
            return '/', '/'
            
    def gn_indecer(self, gn, delimiter, direction):
        ''' Get an upstream or downstream generic number from reference generic number.
        
            @param gn: str, Generic number \n
            @param delimiter: str, character between TM and position (usually 'x') \n 
            @param direction: int, n'th position from gn (+ or -)
        '''
        split = self.gn_num_extract(gn, delimiter)
        if split[0]!='/':
            if len(str(split[1]))==2:
                return str(split[0])+delimiter+str(split[1]+direction)
            elif len(str(split[1]))==3:
                if direction<0:
                    direction += 1
                return str(split[0])+delimiter+str(int(str(split[1])[:2])+direction)
        return '/'

    def fetch_residues_from_pdb(self, pdb, preferred_chain, generic_numbers, path='default'):
        ''' Fetches specific lines from pdb file by generic number (if generic number is
            not available then by residue number). Returns nested OrderedDict()
            with generic numbers as keys in the outer dictionary, and atom names as keys
            in the inner dictionary.
            
            @param pdb: Structure or str, 4 letter pdb code \n
            @param preferred_chain: str, preferred chain of structure \n
            @param generic_numbers: list, list of generic numbers to be fetched \n
            @param path: str, path to pdb file
        '''
        if path=='default':
            pdb_array = self.pdb_array_creator('./structure/PDB/'+str(pdb)+'_'+str(preferred_chain)+'_GPCRDB.pdb')
        else:
            pdb_array = self.pdb_array_creator(path+str(pdb)+'_'+str(preferred_chain)+'_GPCRDB.pdb')
        output = OrderedDict()
        for gn in generic_numbers:
            output[gn.replace('x','.')] = pdb_array[gn.replace('x','.')]
        return output

    def pdb_array_creator(self, filename):
        ''' Creates an OrderedDict() from a pdb file where residue numbers/generic numbers are 
            keys for the residues, and atom names are keys for the Bio.PDB.Residue objects.
            
            @param filename: str, file name with path.
        '''
        residue_array = OrderedDict()
        parser = PDB.PDBParser()
        pdb_struct = parser.get_structure('structure', filename)
        for model in pdb_struct:        
            for chain in model:
                for residue in chain:
                    if -8.1 < residue['CA'].get_bfactor() < 8.1:
                        gn = str(residue['CA'].get_bfactor())
                        if gn[0]=='-':
                            gn = gn[1:]+'1'
                        elif len(gn)==3:
                            gn = gn+'0'
                        residue_array[gn] = residue
                    else:                          
                        residue_array[str(residue.get_id()[1])] = residue
        return residue_array
        
class CreateStatistics(object):
    ''' Statistics dictionary for HomologyModeling.
    '''
    def __init__(self, reference):
        self.reference = reference
        self.info_dict = OrderedDict()
    
    def __repr__(self):
        return "<{} \n {} \n>".format(self.reference, self.info_dict)
    
    def add_info(self, info_name, info):
        ''' Adds new information to the statistics dictionary.
        
            @param info_name: str, info name as dictionary key
            @param info: object, any object as value
        '''
        self.info_dict[info_name] = info