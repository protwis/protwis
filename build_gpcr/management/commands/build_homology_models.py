from build.management.commands.base_build import Command as BaseBuild

from protein.models import Protein, ProteinConformation, ProteinAnomaly, ProteinState
from residue.models import Residue
from structure.models import *
from common.alignment import AlignedReferenceTemplate
import structure.structural_superposition as sp
import structure.assign_generic_numbers_gpcr as as_gn

import Bio.PDB as PDB
from modeller import *
from modeller.automodel import *
from collections import OrderedDict
import os
import logging
import pprint
from io import StringIO
import sys
import re
import zipfile
import shutil
import math
from datetime import datetime


startTime = datetime.now()
logger = logging.getLogger('homology_modeling')
hdlr = logging.FileHandler('./logs/homology_modeling.log')
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr) 
logger.setLevel(logging.INFO)

        
class Command(BaseBuild):  
    help = 'Build automated chimeric GPCR homology models'    
    
    def add_arguments(self, parser):
        super(Command, self).add_arguments(parser=parser)
        parser.add_argument('--update', help='Upload model to GPCRdb, overwrites existing entry', default=False, 
                            action='store_true')
        parser.add_argument('-r', help='Run program for specific receptor(s) by giving UniProt common name as argument', 
                            default=False, type=str, nargs='+')
        parser.add_argument('-z', help='Create zip file of model directory containing all built models', default=False,
                            action='store_true')
        parser.add_argument('--hmver', help='Homology modeling version', default=1.0, type=float)
        
    def handle(self, *args, **options):
        if not os.path.exists('./structure/homology_models/'):
            os.mkdir('./structure/homology_models')
        if not os.path.exists('./structure/PIR/'):
            os.mkdir('./structure/PIR')
        if options['update']==True:
            self.update = True
        else:
            self.update = False
        if options['hmver']!=1.0:
            self.version = options['hmver']
        else:
            self.version = 1.0
        if options['r']==False:
            structures = Structure.objects.all()
            struct_parent = [i.protein_conformation.protein.parent for i in structures]
            classA = Protein.objects.filter(parent__isnull=True, accession__isnull=False, species=1, 
                                            family__slug__istartswith='001')
            self.receptor_list = [i.entry_name for i in classA if i not in struct_parent]
            try:
                self.prepare_input(options['proc'], self.receptor_list)
            except Exception as msg:
                print(msg)
        elif len(options['r'])>1:
            self.receptor_list = options['r']
            try:
                self.prepare_input(options['proc'], self.receptor_list)
            except Exception as msg:
                print(msg)
        else:
            self.run_HomologyModeling(options['r'][0])
        
        os.chdir('./structure/')
        if options['z']==True:
            zipf = zipfile.ZipFile('../static/homology_models/homology_models_v{}.zip'.format(str(self.version)),'w',zipfile.ZIP_DEFLATED)
            for root, dirs, files in os.walk('homology_models'):
                for f in files:
                    zipf.write(os.path.join(root, f))
            zipf.close()
#        shutil.rmtree('homology_models')
        shutil.rmtree('PIR')

    def main_func(self, positions, iteration):
        if not positions[1]:
            receptor_list = self.receptor_list[positions[0]:]
        else:
            receptor_list = self.receptor_list[positions[0]:positions[1]]
        
        for receptor in receptor_list:
            self.run_HomologyModeling(receptor)
    
    def run_HomologyModeling(self, receptor):
#        try:
        state = 'Inactive'
        Homology_model = HomologyModeling(receptor, state, [state], update=self.update, version=self.version)
        alignment = Homology_model.run_alignment()
        Homology_model.build_homology_model(alignment)
        Homology_model.format_final_model()
        logger.info('Model built for {} {}'.format(receptor, state))
#        except Exception as msg:
#            print('Failed to build model {}\n{}'.format(receptor,msg))
#            logger.error('Failed to build model {}\n    {}'.format(receptor,msg))

        
class HomologyModeling(object):
    ''' Class to build homology models for GPCRs. 
    
        @param reference_entry_name: str, protein entry name \n
        @param state: str, endogenous ligand state of reference \n
        @param query_states: list, list of endogenous ligand states to be applied for template search, 
        default: same as reference
    '''
    segment_coding = {1:'TM1',2:'TM2',3:'TM3',4:'TM4',5:'TM5',6:'TM6',7:'TM7'}
    def __init__(self, reference_entry_name, state, query_states, update=False, version=1.0):
        self.reference_entry_name = reference_entry_name
        self.state = state
        self.query_states = query_states
        self.update = update
        self.version = version
        self.statistics = CreateStatistics(self.reference_entry_name)
        self.reference_protein = Protein.objects.get(entry_name=self.reference_entry_name)
        self.uniprot_id = self.reference_protein.accession
        self.reference_sequence = self.reference_protein.sequence
        self.reference_class = self.reference_protein.family.parent.parent.parent
        self.statistics.add_info('uniprot_id',self.uniprot_id)
        self.segments = []
        self.similarity_table = OrderedDict()
        self.similarity_table_all = OrderedDict()
        self.main_structure = None
        self.main_template_preferred_chain = ''
        self.loop_template_table = OrderedDict()
        self.loops = OrderedDict()
        self.prot_conf = ProteinConformation.objects.get(protein=self.reference_protein)
        self.template_source = OrderedDict()
        self.helix_end_mods = None
        for r in Residue.objects.filter(protein_conformation=self.prot_conf):
            if r.protein_segment.slug not in self.template_source:
                self.template_source[r.protein_segment.slug] = OrderedDict()
            try:
                self.template_source[r.protein_segment.slug][r.generic_number.label] = [None,None]
            except:
                self.template_source[r.protein_segment.slug][str(r.sequence_number)] = [None,None]
        
    def __repr__(self):
        return "<Hommod: {}, {}>".format(self.reference_entry_name, self.state)

    def upload_to_db(self, sections, rotamers):
        s_state=ProteinState.objects.get(name=self.state)
        formatted_model = self.format_final_model()
        new_entry = False
        try:
            hommod = StructureModel.objects.get(protein=self.reference_protein, state=s_state)
            hommod.main_template = self.main_structure
            hommod.pdb = formatted_model
            hommod.version = self.version
            hommod.save()
        except:
            hommod, created = StructureModel.objects.update_or_create(protein=self.reference_protein, state=s_state, 
                                                                      main_template=self.main_structure, 
                                                                      pdb=self.format_final_model(), 
                                                                      version=self.version)
            new_entry = True
        if not new_entry:
            StructureModelStatsSegment.objects.filter(homology_model=hommod).delete()
            StructureModelStatsRotamer.objects.filter(homology_model=hommod).delete()
        for s in sections:
            segs, created = StructureModelStatsSegment.objects.update_or_create(homology_model=hommod, segment=s[0],
                                                                                start=s[1],end=s[2],start_gn=s[3],
                                                                                end_gn=s[4],backbone_template=s[6])
        segs = StructureModelStatsSegment.objects.filter(homology_model=hommod)
        
        for r in rotamers:
            for s in segs:
                if s.start<=r[1]<=s.end:
                    m_seg = s
            rots, created = StructureModelStatsRotamer.objects.update_or_create(homology_model=hommod, segment=r[0],
                                                                                sequence_number=r[1],generic_number=r[2],
                                                                                rotamer_template=r[4],model_segment=m_seg)
                                   
    def right_rotamer_select(self, rotamer):
        if len(rotamer)>1:
            for i in rotamer:
                if i.pdbdata.pdb.startswith('COMPND')==False:
                    rotamer = i
                    break
        else:
            rotamer=rotamer[0]
        return rotamer
                                                            
    def format_final_model(self):
        pos_list = []
        for seg in self.template_source:
            for num in self.template_source[seg]:
                try:
                    num = str(Residue.objects.get(protein_conformation=self.prot_conf,generic_number__label=num).sequence_number)
                except:
                    pass
                pos_list.append(num)
        i = 0
        with open ('./structure/homology_models/{}_{}/{}_{}_model.pdb'.format(self.reference_entry_name, self.state,
                                                                              self.reference_entry_name, self.state), 
                                                                              'r+') as f:
            pdblines = f.readlines()
            out_list = []
            prev_num = None
            for line in pdblines:
                try:
                    if prev_num==None:
                        pdb_re = re.search('(ATOM[A-Z\s\d]{13}\S{3}\s+)(\d+)([A-Z\s\d.-]{49,53})',line)
                        prev_num = int(pdb_re.group(2))
                    pdb_re = re.search('(ATOM[A-Z\s\d]{13}\S{3}\s+)(\d+)([A-Z\s\d.-]{49,53})',line)
                    if int(pdb_re.group(2))>prev_num:
                        i+=1
                        prev_num = int(pdb_re.group(2))
                    whitespace = (len(pos_list[i])-len(pdb_re.group(2)))*-1
                    if whitespace==0:
                        out_line = pdb_re.group(1)+pos_list[i]+pdb_re.group(3)
                    else:
                        out_line = pdb_re.group(1)[:whitespace]+pos_list[i]+pdb_re.group(3)
                    out_list.append(out_line)
                except:
                    out_list.append(line)
        with open ('./structure/homology_models/{}_{}/{}_{}_model.pdb'.format(self.reference_entry_name, self.state,self.reference_entry_name,self.state), 'w') as f:
            f.write(''.join(out_list))
        return ''.join(out_list)


    def fetch_struct_helix_ends_from_db(self, structure):
        ''' Returns structure's helix end generic numbers after updating them with annotated data.
        '''
        raw = StructureSegment.objects.filter(structure=structure)
        annotated = StructureSegmentModeling.objects.filter(structure=structure)
        ends = OrderedDict()
        for i in raw:
            if i.protein_segment.slug[0]=='T' or i.protein_segment.slug=='H8':
                start_found = False
                while start_found==False:
                    try:
                        if Residue.objects.get(protein_conformation=structure.protein_conformation,sequence_number=i.start).generic_number==None:
                            i.start+=1
                        else:
                            start_found = True
                    except:
                        i.start+=1
                s = Residue.objects.get(protein_conformation=structure.protein_conformation,sequence_number=i.start)
                end_found = False                
                while end_found==False:
                    try:
                        if Residue.objects.get(protein_conformation=structure.protein_conformation,sequence_number=i.end).generic_number==None:
                            i.end-=1
                        else:
                            end_found = True
                    except:
                        i.end-=1
                e = Residue.objects.get(protein_conformation=structure.protein_conformation,sequence_number=i.end)
                ends[s.protein_segment.slug] = [s.generic_number.label,e.generic_number.label]
        for j in annotated:
            if j.protein_segment.slug[0]=='T' or j.protein_segment.slug=='H8':
                if j.start!=0:
                    while Residue.objects.get(protein_conformation=structure.protein_conformation,sequence_number=j.start).generic_number==None:
                        j.start+=1
                    sa = Residue.objects.get(protein_conformation=structure.protein_conformation,sequence_number=j.start)
                    ends[j.protein_segment.slug][0] = sa.generic_number.label
                if j.end!=0:
                    while Residue.objects.get(protein_conformation=structure.protein_conformation,sequence_number=j.end).generic_number==None:
                        j.end-=1
                    ea = Residue.objects.get(protein_conformation=structure.protein_conformation,sequence_number=j.end)
                    ends[j.protein_segment.slug][1] = ea.generic_number.label
        return ends      

    def fetch_struct_helix_ends_from_array(self, array):
        ''' Returns helix ends from structure array (GPCRDBParsingPDB.pdb_array_creator()).
        '''
        ends = OrderedDict()
        for seg_lab, seg in array.items():
            if seg_lab[0]=='T' or seg_lab=='H8':
                ends[seg_lab] = [list(seg.keys())[0].replace('.','x'),list(seg.keys())[-1].replace('.','x')]
        return ends
        
    def correct_helix_ends(self, main_structure, main_pdb_array, a):
        ''' Updates main template structure with annotated helix ends, if helix is too long, it removes residues, if it
            is too short, it superpositions residues from next closest template. Updates alignment with changes.
        '''
        raw_helix_ends = self.fetch_struct_helix_ends_from_array(main_pdb_array)
        anno_helix_ends = self.fetch_struct_helix_ends_from_db(main_structure)
        for lab,seg in a.template_dict.items():
            for gn,res in seg.items():
                try:
                    if lab[0] in ['H']:
                        if res!='-':
                            r = Residue.objects.get(protein_conformation=main_structure.protein_conformation,
                                                    generic_number__label=gn)
                            if len(Rotamer.objects.filter(structure=main_structure,residue=r))<1:
                                raise Exception()
                except:
                    a.template_dict[lab][gn] = 'x'       
        parser = GPCRDBParsingPDB()
        for raw_seg, anno_seg in zip(raw_helix_ends, anno_helix_ends):
            s_dif = parser.gn_comparer(raw_helix_ends[raw_seg][0],anno_helix_ends[anno_seg][0],main_structure.protein_conformation)
            e_dif = parser.gn_comparer(raw_helix_ends[raw_seg][1],anno_helix_ends[anno_seg][1],main_structure.protein_conformation)
            if s_dif<0:
                s_gn = Residue.objects.get(protein_conformation=main_structure.protein_conformation, generic_number__label=raw_helix_ends[raw_seg][0])
                seq_nums = [i for i in range(s_gn.sequence_number,s_gn.sequence_number-s_dif)]
                gns = [j.generic_number.label for j in list(Residue.objects.filter(protein_conformation=main_structure.protein_conformation, sequence_number__in=seq_nums))]
                for gn in gns:
                    a.template_dict[raw_seg][gn]='x'
                    a.alignment_dict[raw_seg][gn]='x'
            if e_dif>0:
                e_gn = Residue.objects.get(protein_conformation=main_structure.protein_conformation, generic_number__label=raw_helix_ends[raw_seg][1])
                seq_nums = [i for i in range(e_gn.sequence_number-e_dif+1,e_gn.sequence_number+1)]
                gns = [j.generic_number.label for j in list(Residue.objects.filter(protein_conformation=main_structure.protein_conformation, sequence_number__in=seq_nums))]
                for gn in gns:
                    a.template_dict[raw_seg][gn]='x'
                    a.alignment_dict[raw_seg][gn]='x'
        
        self.helix_ends = raw_helix_ends
        modifications = {'added':{'TM1':[[],[]],'TM2':[[],[]],'TM3':[[],[]],'TM4':[[],[]],'TM5':[[],[]],'TM6':[[],[]],'TM7':[[],[]], 'H8':[[],[]]},
                         'removed':{'TM1':[[],[]],'TM2':[[],[]],'TM3':[[],[]],'TM4':[[],[]],'TM5':[[],[]],'TM6':[[],[]],'TM7':[[],[]], 'H8':[[],[]]}}
        for ref_seg, temp_seg, align_seg in zip(a.reference_dict, a.template_dict, a.alignment_dict):
            mid = len(a.reference_dict[ref_seg])/2
            if ref_seg[0] not in ['T','H']:
                continue
            for ref_res, temp_res, align_res in zip(a.reference_dict[ref_seg],a.template_dict[temp_seg],a.alignment_dict[align_seg]):
                if a.template_dict[temp_seg][temp_res]=='-':
                    continue
                if a.reference_dict[ref_seg][ref_res]=='x':
                    if list(a.reference_dict[ref_seg].keys()).index(ref_res)<mid:    
                        modifications['removed'][ref_seg][0].append(ref_res)
                    else:
                        modifications['removed'][ref_seg][1].append(ref_res)
                    del a.reference_dict[ref_seg][ref_res]
                    del a.template_dict[temp_seg][temp_res]
                    del a.alignment_dict[align_seg][align_res]
                    try:
                        del main_pdb_array[ref_seg][ref_res.replace('x','.')]
                    except:
                        pass
                elif a.template_dict[temp_seg][temp_res]=='x' or (temp_seg[0]=='T' and temp_res.replace('x','.') not in list(main_pdb_array[temp_seg])):
                    if list(a.template_dict[temp_seg].keys()).index(temp_res)<mid:    
                        modifications['added'][temp_seg][0].append(temp_res)
                    else:
                        modifications['added'][temp_seg][1].append(temp_res)
            if ref_seg[0]=='T' or ref_seg=='H8':
                if len(modifications['added'][ref_seg][0])>0:
                    self.helix_ends[ref_seg][0] = modifications['added'][ref_seg][0][0]
                if len(modifications['added'][ref_seg][1])>0:
                    self.helix_ends[ref_seg][1] = modifications['added'][ref_seg][1][-1]               
                if len(modifications['removed'][ref_seg][0])>0:
                    self.helix_ends[ref_seg][0] = parser.gn_indecer(modifications['removed'][ref_seg][0][-1], 'x', 1)
                if len(modifications['removed'][ref_seg][1])>0:
                    self.helix_ends[ref_seg][1] = parser.gn_indecer(modifications['removed'][ref_seg][1][0], 'x', -1)
                if len(modifications['added'][ref_seg][0])>0:
                    for struct in self.similarity_table:
                        if struct!=main_structure:
                            alt_helix_ends = self.fetch_struct_helix_ends_from_db(struct)
                            try:
                                if parser.gn_comparer(alt_helix_ends[ref_seg][0],self.helix_ends[ref_seg][0],struct.protein_conformation)<=0:
                                    all_keys = list(a.reference_dict[ref_seg].keys())[:len(modifications['added'][ref_seg][0])+4]
                                    ref_keys = [i for i in all_keys if i not in modifications['added'][ref_seg][0]]
                                    reference = parser.fetch_residues_from_array(main_pdb_array[ref_seg],ref_keys)
                                    template = parser.fetch_residues_from_pdb(struct,all_keys)
                                    superpose = sp.OneSidedSuperpose(reference,template,4,0)
                                    sup_residues = superpose.run()
                                    new_residues = OrderedDict()
                                    for gn, atoms in sup_residues.items():
                                        if gn.replace('.','x') not in ref_keys:
                                            new_residues[gn] = atoms
                                            gn_ = gn.replace('.','x')
                                            a.template_dict[temp_seg][gn_] = PDB.Polypeptide.three_to_one(
                                                                             atoms[0].get_parent().get_resname())
                                            if a.template_dict[temp_seg][gn_]==a.reference_dict[ref_seg][gn_]:
                                                a.alignment_dict[ref_seg][gn_] = a.reference_dict[ref_seg][gn_]
                                            else:
                                                a.alignment_dict[ref_seg][gn_] = '.'
                                    for gn, atoms in main_pdb_array[ref_seg].items():
                                        if gn not in new_residues:
                                            new_residues[gn] = atoms
                                    main_pdb_array[ref_seg] = new_residues
                                    self.update_template_source(modifications['added'][ref_seg][0],struct,ref_seg)
                                    break
                            except:
                                pass
                if len(modifications['added'][ref_seg][1])>0:
                    for struct in self.similarity_table:
                        if struct!=main_structure:
                            alt_helix_ends = self.fetch_struct_helix_ends_from_db(struct)
                            try:
                                if parser.gn_comparer(alt_helix_ends[ref_seg][1],self.helix_ends[ref_seg][1],struct.protein_conformation)>=0:
                                    all_keys = list(a.reference_dict[ref_seg].keys())[-1*(len(modifications['added'][ref_seg][1])+4):]
                                    ref_keys = [i for i in all_keys if i not in modifications['added'][ref_seg][1]]
                                    reference = parser.fetch_residues_from_array(main_pdb_array[ref_seg],ref_keys)
                                    template = parser.fetch_residues_from_pdb(struct,all_keys)
                                    superpose = sp.OneSidedSuperpose(reference,template,4,1)
                                    sup_residues = superpose.run()
                                    new_residues = OrderedDict()
                                    for gn, atoms in sup_residues.items():
                                        if gn.replace('.','x') not in ref_keys:
                                            new_residues[gn]=atoms
                                    for gn, atoms in new_residues.items():
                                        gn_ = gn.replace('.','x')
                                        if gn_ in modifications['added'][ref_seg][1]:
                                            main_pdb_array[ref_seg][gn] = atoms
                                            a.template_dict[ref_seg][gn_] = PDB.Polypeptide.three_to_one(
                                                                            atoms[0].get_parent().get_resname())
                                            if a.template_dict[ref_seg][gn_]==a.reference_dict[ref_seg][gn_]:
                                                a.alignment_dict[ref_seg][gn_] = a.reference_dict[ref_seg][gn_]
                                            else:
                                                a.alignment_dict[ref_seg][gn_] = '.'
                                    self.update_template_source(modifications['added'][ref_seg][1],struct,segment=ref_seg)
                                    break
                            except:
                                pass
        self.helix_end_mods = modifications
        self.statistics.add_info('helix_end_mods',self.helix_end_mods)
        return main_pdb_array, a
        
    def update_template_source(self, keys, struct, segment, just_rot=False):
        for k in keys:
            if just_rot==True:
                try:
                    self.template_source[segment][k][1] = struct
                except:
                    pass
            else:
                try:
                    self.template_source[segment][k][0] = struct
                except:
                    pass
        
    def run_alignment(self, core_alignment=True, query_states='default', 
                      segments=['TM1','ICL1','TM2','ECL1','TM3','ICL2','TM4','ECL2','TM5','TM6','TM7','H8'], 
                      order_by='similarity'):
        ''' Creates pairwise alignment between reference and target receptor(s).
            Returns Alignment object.
            
            @param segments: list, list of segments to use, e.g.: ['TM1','ICL1','TM2','ECL1'] \n
            @param order_by: str, order results by identity, similarity or simscore
        '''
        if query_states=='default':
            query_states=self.query_states
        alignment = AlignedReferenceTemplate(self.reference_protein, segments, query_states, order_by)
        main_pdb_array = OrderedDict()
        if core_alignment==True:
#            print('Alignment: ',datetime.now() - startTime)
            enhanced_alignment = alignment.enhance_best_alignment()
#            print('Enhanced alignment: ',datetime.now() - startTime)
            if enhanced_alignment==None:
                return None
            self.segments = segments
            self.main_structure = alignment.main_template_structure           
            self.similarity_table = alignment.similarity_table
            self.similarity_table_all = self.run_alignment(core_alignment=False, 
                                                           query_states=['Inactive','Active'])[0].similarity_table
            self.main_template_preferred_chain = str(self.main_structure.preferred_chain)[0]
            self.statistics.add_info("main_template", self.main_structure)
            self.statistics.add_info("preferred_chain", self.main_template_preferred_chain)
            
            parse = GPCRDBParsingPDB()
            main_pdb_array = parse.pdb_array_creator(structure=self.main_structure)
            
            try:
                if len(alignment.reference_dict['H8'])==0:
                    del alignment.reference_dict['H8']
                    del alignment.template_dict['H8']
                    del alignment.alignment_dict['H8']
                    del main_pdb_array['H8']
            except:
                pass                
            for seg_l, seg in main_pdb_array.items():
                for gn, res in seg.items():
                    self.update_template_source([gn.replace('.','x')],self.main_structure,seg_l)
            end_correction = self.correct_helix_ends(self.main_structure, main_pdb_array, alignment)
            main_pdb_array = end_correction[0]
            alignment = end_correction[1]
            loops_in_ref = [i for i in list(self.template_source) if i[0] not in ['N','C','T','H']]
            for loop in loops_in_ref:
                loop_alignment = AlignedReferenceTemplate(self.reference_protein, [loop], ['Inactive','Active'], 
                                                          order_by='similarity', 
                                                          provide_main_template_structure=self.main_structure,
                                                          provide_similarity_table=self.similarity_table_all,
                                                          main_pdb_array=main_pdb_array)
                self.loop_template_table[loop] = loop_alignment.loop_table
                try:
                    if loop in list(alignment.alignment_dict.keys()) and self.main_structure in loop_alignment.loop_table:
                        temp_loop_table = OrderedDict([('aligned',100)])
                        try:
                            for lab, val in loop_alignment.loop_table.items():
                                temp_loop_table[lab] = val
                            self.loop_template_table[loop] = temp_loop_table
                        except:
                            pass
                except:
                    pass
            self.statistics.add_info('similarity_table', self.similarity_table)
            self.statistics.add_info('loops',self.loop_template_table)
#            print('Loop alignment: ',datetime.now() - startTime)
        return alignment, main_pdb_array
        
        
    def build_homology_model(self, ref_temp_alignment, switch_bulges=True, switch_constrictions=True, loops=True, 
                             switch_rotamers=True, N_and_C_termini=True):
        ''' Function to identify and switch non conserved residues in the alignment. Optionally,
            it can identify and switch bulge and constriction sites too. 
            
            @param ref_temp_alignment: AlignedReferenceAndTemplate, alignment of reference and main template with 
            alignment string. \n
            @param switch_bulges: boolean, identify and switch bulge sites. Default = True.
            @param switch_constrictions: boolean, identify and switch constriction sites. Default = True.
        '''
        a = ref_temp_alignment[0]
        main_pdb_array = ref_temp_alignment[1]
        ref_bulge_list, temp_bulge_list, ref_const_list, temp_const_list = [],[],[],[]
        parse = GPCRDBParsingPDB()
        try:
            if len(a.reference_dict['H8'])==0:
                del a.reference_dict['H8']
            if len(a.template_dict['H8'])==0:
                del a.template_dict['H8']
            if len(a.alignment_dict['H8'])==0:
                del a.alignment_dict['H8']
            if len(main_pdb_array['H8'])==0:
                del main_pdb_array['H8']
        except:
            pass
        
        # loops
        if loops==True:
            loop_stat = OrderedDict()
            for label, structures in self.loop_template_table.items():
                loop = Loops(self.reference_protein, label, structures, self.main_structure, self.helix_end_mods,
                             list(self.template_source))
                loop_template = loop.fetch_loop_residues(main_pdb_array)
                
                if type(loop.loop_output_structure)!=type([]):
                    loop_insertion = loop.insert_loop_to_arrays(loop.loop_output_structure, main_pdb_array, loop_template, 
                                                                a.reference_dict, a.template_dict, a.alignment_dict)
                else:
                    loop_insertion = loop.insert_ECL2_to_arrays(loop.loop_output_structure, main_pdb_array, loop_template,
                                                                a.reference_dict, a.template_dict, a.alignment_dict)
                main_pdb_array = loop_insertion.main_pdb_array
                a.reference_dict = loop_insertion.reference_dict
                a.template_dict = loop_insertion.template_dict
                a.alignment_dict = loop_insertion.alignment_dict
                if loop.new_label!=None:
                    loop_stat[loop.new_label] = loop.loop_output_structure
                    if 'dis' in loop.new_label:
                        self.update_template_source(list(self.template_source[label].keys())[1:-1],loop.loop_output_structure,label)
                    else:
                        self.update_template_source(list(self.template_source[label].keys()),loop.loop_output_structure,label)
                else:
                    loop_stat[label] = loop.loop_output_structure
                    if label=='ECL2' and loop.loop_output_structure!=None:
                        x50 = list(self.template_source[label].keys()).index('45x50')
                        self.update_template_source(list(self.template_source[label].keys())[:x50],loop.loop_output_structure[0],label)
                        self.update_template_source(list(self.template_source[label].keys())[x50:x50+3],loop.loop_output_structure[1],label)
                        self.update_template_source(list(self.template_source[label].keys())[x50+3:],loop.loop_output_structure[2],label)
                
            self.statistics.add_info('loops', loop_stat)
            self.loops = loop_stat
        
#        print('Integrate loops: ',datetime.now() - startTime)

        # bulges and constrictions
        if switch_bulges==True or switch_constrictions==True:
            for ref_seg, temp_seg, aligned_seg in zip(a.reference_dict, a.template_dict, a.alignment_dict):
                if ref_seg[0]=='T':
                    for ref_res, temp_res, aligned_res in zip(a.reference_dict[ref_seg], a.template_dict[temp_seg], 
                                                              a.alignment_dict[aligned_seg]):
                        gn = ref_res
                        gn_num = parse.gn_num_extract(gn, 'x')[1]
                        
                        if a.alignment_dict[aligned_seg][aligned_res]=='-':
                            if (a.reference_dict[ref_seg][ref_res]=='-' and 
                                a.reference_dict[ref_seg][parse.gn_indecer(gn,'x',-1)] not in 
                                ['-','/'] and a.reference_dict[ref_seg][parse.gn_indecer(gn,'x',+1)] not in ['-','/']): 
            
                                # bulge in template
                                if len(str(gn_num))==3:
                                    if switch_bulges==True:
                                        try:
                                            Bulge = Bulges(gn)
                                            bulge_template = Bulge.find_bulge_template(self.similarity_table_all, 
                                                                                       bulge_in_reference=False)
                                            bulge_site = OrderedDict([
                                                (parse.gn_indecer(gn,'x',-2).replace('x','.'), 
                                                 main_pdb_array[ref_seg][parse.gn_indecer(gn,'x',-2).replace('x','.')]),
                                                (parse.gn_indecer(gn,'x',-1).replace('x','.'), 
                                                 main_pdb_array[ref_seg][parse.gn_indecer(gn,'x',-1).replace('x','.')]),
                                                (gn.replace('x','.'), 
                                                 main_pdb_array[ref_seg][gn.replace('x','.')]),
                                                (parse.gn_indecer(gn,'x',+1).replace('x','.'), 
                                                 main_pdb_array[ref_seg][parse.gn_indecer(gn,'x',+1).replace('x','.')]),
                                                (parse.gn_indecer(gn,'x',+2).replace('x','.'), 
                                                 main_pdb_array[ref_seg][parse.gn_indecer(gn,'x',+2).replace('x','.')])]) 
                                            superpose = sp.BulgeConstrictionSuperpose(bulge_site, bulge_template)
                                            new_residues = superpose.run()
                                            switch_res = 0
                                            for gen_num, atoms in bulge_template.items():
                                                if switch_res!=0 and switch_res!=3:
                                                    gn__ = gen_num.replace('.','x')
                                                    self.update_template_source([gn__],Bulge.template,ref_seg)
                                                    main_pdb_array[ref_seg][gen_num] = new_residues[gen_num]
                                                    a.template_dict[temp_seg][gn__] = PDB.Polypeptide.three_to_one(
                                                                                       atoms[0].get_parent().get_resname())
                                                    if a.template_dict[temp_seg][gn__]==a.reference_dict[ref_seg][gn__]:
                                                        a.alignment_dict[aligned_seg][gn__]=a.template_dict[temp_seg][gn__]
                                                    else:
                                                        a.alignment_dict[aligned_seg][gn__]='.'
                                                switch_res+=1
                                            del main_pdb_array[ref_seg][gn.replace('x','.')]
                                            del a.reference_dict[ref_seg][gn]
                                            del a.template_dict[temp_seg][gn]
                                            del a.alignment_dict[aligned_seg][gn]
                                            temp_bulge_list.append({gn:Bulge.template})
                                        except:
                                            temp_bulge_list.append({gn:None})
                                        
                                # constriction in reference
                                else:
                                    if switch_constrictions==True:
                                        try:
                                            Const = Constrictions(gn)
                                            constriction_template = Const.find_constriction_template(
                                                                                            self.similarity_table_all,
                                                                                            constriction_in_reference=True)
                                            constriction_site = OrderedDict([
                                                (parse.gn_indecer(gn,'x',-2).replace('x','.'), 
                                                 main_pdb_array[ref_seg][parse.gn_indecer(gn,'x',-2).replace('x','.')]),
                                                (parse.gn_indecer(gn,'x',-1).replace('x','.'), 
                                                 main_pdb_array[ref_seg][parse.gn_indecer(gn,'x',-1).replace('x','.')]),
                                                (gn.replace('x','.'), 
                                                 main_pdb_array[ref_seg][gn.replace('x','.')]),
                                                (parse.gn_indecer(gn,'x',+1).replace('x','.'), 
                                                 main_pdb_array[ref_seg][parse.gn_indecer(gn,'x',+1).replace('x','.')]),
                                                (parse.gn_indecer(gn,'x',+2).replace('x','.'), 
                                                 main_pdb_array[ref_seg][parse.gn_indecer(gn,'x',+2).replace('x','.')])])                                      
                                            superpose = sp.BulgeConstrictionSuperpose(constriction_site, 
                                                                                      constriction_template)
                                            new_residues = superpose.run()                                  
                                            switch_res = 0
                                            for gen_num, atoms in constriction_template.items():
                                                if switch_res!=0 and switch_res!=3:
                                                    gn__ = gen_num.replace('.','x')
                                                    self.update_template_source([gn__],Const.template,ref_seg)
                                                    main_pdb_array[ref_seg][gen_num] = new_residues[gen_num]
                                                    a.template_dict[gn__] = PDB.Polypeptide.three_to_one(
                                                                                       atoms[0].get_parent().get_resname())
                                                    if a.template_dict[temp_seg][gn__]==a.reference_dict[ref_seg][gn__]:
                                                        a.alignment_dict[aligned_seg][gn__]=a.template_dict[temp_seg][gn__]
                                                switch_res+=1
                                            ref_const_list.append({gn:Const.template})
                                            del main_pdb_array[ref_seg][gn.replace('x','.')]
                                            del a.reference_dict[ref_seg][gn]
                                            del a.template_dict[temp_seg][gn]
                                            del a.alignment_dict[aligned_seg][gn]
                                        except:
                                            ref_const_list.append({gn:None})
                            elif (a.template_dict[ref_seg][temp_res]=='-' and 
                                  a.template_dict[temp_seg][parse.gn_indecer(gn,'x',-1)] not in 
                                  ['-','/'] and a.template_dict[temp_seg][parse.gn_indecer(gn,'x',+1)] not in ['-','/']): 
                                
                                # bulge in reference
                                if len(str(gn_num))==3:
                                    if switch_bulges==True:
                                        try:
                                            Bulge = Bulges(gn)
                                            bulge_template = Bulge.find_bulge_template(self.similarity_table_all,
                                                                                       bulge_in_reference=True)
                                            bulge_site = OrderedDict([
                                                (parse.gn_indecer(gn,'x',-2).replace('x','.'), 
                                                 main_pdb_array[ref_seg][parse.gn_indecer(gn,'x',-2).replace('x','.')]),
                                                (parse.gn_indecer(gn,'x',-1).replace('x','.'), 
                                                 main_pdb_array[ref_seg][parse.gn_indecer(gn,'x',-1).replace('x','.')]),
                                                (parse.gn_indecer(gn,'x',+1).replace('x','.'), 
                                                 main_pdb_array[ref_seg][parse.gn_indecer(gn,'x',+1).replace('x','.')]),
                                                (parse.gn_indecer(gn,'x',+2).replace('x','.'), 
                                                 main_pdb_array[ref_seg][parse.gn_indecer(gn,'x',+2).replace('x','.')])]) 
                                            superpose = sp.BulgeConstrictionSuperpose(bulge_site, bulge_template)
                                            new_residues = superpose.run()
                                            switch_res = 0
                                            for gen_num, atoms in bulge_template.items():
                                                if switch_res!=0 and switch_res!=4:
                                                    gn__ = gen_num.replace('.','x')
                                                    self.update_template_source([gn__],Bulge.template,ref_seg)
                                                    main_pdb_array[ref_seg][gen_num] = new_residues[gen_num]
                                                    a.template_dict[temp_seg][gn__] = PDB.Polypeptide.three_to_one(
                                                                                       atoms[0].get_parent().get_resname())
                                                    if a.template_dict[temp_seg][gn__]==a.reference_dict[ref_seg][gn__]:
                                                        a.alignment_dict[aligned_seg][gn__]=a.template_dict[temp_seg][gn__]
                                                switch_res+=1
                                            ref_bulge_list.append({gn:Bulge.template})
                                            if a.reference_dict[ref_seg][gn] == a.template_dict[temp_seg][gn]:
                                                a.alignment_dict[ref_seg][gn] = a.reference_dict[ref_seg][gn]
                                            else:
                                                a.alignment_dict[ref_seg][gn] = '.'
                                        except:
                                            ref_bulge_list.append({gn:None})
                                        
                                # constriction in template
                                else:
                                    if switch_constrictions==True:
                                        try:
                                            Const = Constrictions(gn)
                                            constriction_template = Const.find_constriction_template(
                                                                                           self.similarity_table_all,
                                                                                           constriction_in_reference=False)
                                            constriction_site = OrderedDict([
                                                (parse.gn_indecer(gn,'x',-2).replace('x','.'), 
                                                 main_pdb_array[ref_seg][parse.gn_indecer(gn,'x',-2).replace('x','.')]),
                                                (parse.gn_indecer(gn,'x',-1).replace('x','.'), 
                                                 main_pdb_array[ref_seg][parse.gn_indecer(gn,'x',-1).replace('x','.')]),
                                                (parse.gn_indecer(gn,'x',+1).replace('x','.'), 
                                                 main_pdb_array[ref_seg][parse.gn_indecer(gn,'x',+1).replace('x','.')]),
                                                (parse.gn_indecer(gn,'x',+2).replace('x','.'), 
                                                 main_pdb_array[ref_seg][parse.gn_indecer(gn,'x',+2).replace('x','.')])]) 
                                            superpose = sp.BulgeConstrictionSuperpose(constriction_site, 
                                                                                      constriction_template)
                                            new_residues = superpose.run()
                                            switch_res = 0
                                            for gen_num, atoms in constriction_template.items():
                                                if switch_res!=0 and switch_res!=4:
                                                    gn__ = gen_num.replace('.','x')
                                                    self.update_template_source([gn__],Const.template,ref_seg)
                                                    main_pdb_array[ref_seg][gen_num] = new_residues[gen_num]
                                                    a.template_dict[temp_seg][gn__] = PDB.Polypeptide.three_to_one(
                                                                                       atoms[0].get_parent().get_resname())
                                                    if a.template_dict[temp_seg][gn__]==a.reference_dict[ref_seg][gn__]:
                                                        a.alignment_dict[aligned_seg][gn__]=a.template_dict[temp_seg][gn__]
                                                switch_res+=1
                                            temp_const_list.append({gn:Const.template})
                                            if a.reference_dict[ref_seg][gn] == a.template_dict[temp_seg][gn]:
                                                a.alignment_dict[ref_seg][gn] = a.reference_dict[ref_seg][gn]
                                            else:
                                                a.alignment_dict[ref_seg][gn] = '.'
                                        except:
                                            temp_const_list.append({gn:None})
                                        
            self.statistics.add_info('reference_bulges', ref_bulge_list)
            self.statistics.add_info('template_bulges', temp_bulge_list)
            self.statistics.add_info('reference_constrictions', ref_const_list)
            self.statistics.add_info('template_constrictions', temp_const_list)
            
            # insert bulge to array in the right place
            if ref_bulge_list!=[]:
                out_pdb_array = OrderedDict()
                bulge_gns = []
                for bulge in ref_bulge_list:
                    if list(bulge.values())[0]!=None:
                        gn = list(bulge.keys())[0].replace('x','.')
                        bulge_gns.append(gn)
                for seg_id, residues in main_pdb_array.items():
                    seg = OrderedDict()
                    for key, value in residues.items():
                        seg[key] = value                
                        if str(key)+'1' in bulge_gns:
                            seg[str(key)+'1'] = main_pdb_array[seg_id][str(key)+'1']
                    out_pdb_array[seg_id] = seg
                main_pdb_array = out_pdb_array
            
            if temp_const_list!=[]:
                out_pdb_array = OrderedDict()
                const_gns = []
                for const in temp_const_list:
                    if list(const.values())[0]!=None:
                        gn = list(const.keys())[0].replace('x','.')
                        const_gns.append(gn)
                for seg_id, residues in main_pdb_array.items():
                    seg = OrderedDict()
                    for key, value in residues.items():
                        seg[key] = value
                        if parse.gn_indecer(key, '.', +1) in const_gns:
                            seg[gn] = main_pdb_array[seg_id][gn]
                    out_pdb_array[seg_id] = seg
                main_pdb_array = out_pdb_array
#        print('Integrate bulges/constrictions: ',datetime.now() - startTime)
                
        # check for inconsitencies with db
        pdb_db_inconsistencies = []
        for seg_label, segment in a.template_dict.items():
            try:
                for gn, res in segment.items():
                    try:
                        if res==PDB.Polypeptide.three_to_one(
                                            main_pdb_array[seg_label][gn.replace('x','.')][0].get_parent().get_resname()):
                            pass
                        elif 'x' in gn:
                            try:
                                Residue.objects.get(
                                        protein_conformation__protein=self.main_structure.protein_conformation.protein, 
                                        generic_number__label=gn)
                                pdb_db_inconsistencies.append({gn:a.template_dict[seg_label][gn]})
                            except:
                                pass
                    except:
                        pass
            except:
                pass
        
        if pdb_db_inconsistencies!=[]:
            for incons in pdb_db_inconsistencies:
                seg = self.segment_coding[int(list(incons.keys())[0][0])]
                seq_num = Residue.objects.get(
                                        protein_conformation__protein=self.main_structure.protein_conformation.protein, 
                                        generic_number__label=list(incons.keys())[0])
                temp_segment, temp_array = OrderedDict(), OrderedDict()
                for key, value in main_pdb_array[seg].items():
                    if key==str(seq_num.sequence_number):
                        temp_segment[list(incons.keys())[0].replace('x','.')] = value
                    else:
                        temp_segment[key] = value
                for seg_id, segment in main_pdb_array.items():
                    if seg_id==seg:
                        temp_array[seg_id] = temp_segment
                    else:
                        temp_array[seg_id] = segment
                main_pdb_array = temp_array
                a.template_dict[seg][list(incons.keys())[0]] = PDB.Polypeptide.three_to_one(
                            main_pdb_array[seg][list(incons.keys())[0].replace('x','.')][0].get_parent().get_resname())
                if a.reference_dict[seg][list(incons.keys())[0]]==a.template_dict[seg][list(incons.keys())[0]]:
                    a.alignment_dict[seg][list(incons.keys())[0]] = a.reference_dict[seg][list(incons.keys())[0]]
        
        self.statistics.add_info('pdb_db_inconsistencies', pdb_db_inconsistencies)
        path = "./structure/homology_models/{}_{}/".format(self.reference_entry_name,self.state)
        if not os.path.exists(path):
            os.mkdir(path)
#        self.write_homology_model_pdb(
#                                "./structure/homology_models/{}_{}/pre_switch.pdb".format(self.reference_entry_name, self.state), 
#                                main_pdb_array)        
#        print('Check inconsistencies: ',datetime.now() - startTime)
        # inserting loops for free modeling
        for label, template in loop_stat.items():
            if template==None:
                modeling_loops = Loops(self.reference_protein, label, self.similarity_table_all, self.main_structure, 
                                       self.helix_end_mods, list(self.template_source))
                modeling_loops.insert_gaps_for_loops_to_arrays(main_pdb_array, a.reference_dict, a.template_dict,
                                                               a.alignment_dict)
                main_pdb_array = modeling_loops.main_pdb_array
                a.reference_dict = modeling_loops.reference_dict
                a.template_dict = modeling_loops.template_dict
                a.alignment_dict = modeling_loops.alignment_dict
#        print('Free loops: ',datetime.now() - startTime)
        
        # Move H8 if needed
        trimmed_residues=[]
        if 'H8' in main_pdb_array and 'ICL4' not in main_pdb_array and len(self.helix_end_mods['removed']['TM7'][1])>0:
            unwind_num = math.ceil(len(self.helix_end_mods['removed']['TM7'][1])/2)
            trimmed_residues+=list(main_pdb_array['TM7'].keys())[(unwind_num*-1):]+list(main_pdb_array['H8'].keys())[:unwind_num]
#            reference = OrderedDict()
#            ref_nums = list(main_pdb_array['TM7'])[-2:]
#            for i in ref_nums:
#                reference[i] = main_pdb_array['TM7'][i]
#            h8_temp = self.template_source['H8'][list(self.template_source['H8'])[0]][0]
#            tm7 = list(Residue.objects.filter(protein_conformation=h8_temp.protein_conformation,protein_segment__slug='TM7'))
#            temp_nums = [i.generic_number.label for i in tm7[-2:]] 
#            tm7_template = parse.fetch_residues_from_pdb(h8_temp,temp_nums)
#            template = parse.add_two_ordereddict(tm7_template,main_pdb_array['H8'])
#            superpose = sp.OneSidedSuperpose(reference,template,2,1)
#            sup_residues = superpose.run()
#            new_h8 = OrderedDict()
#            for gn, atoms in sup_residues.items():
#                if gn[0]=='8':
#                    new_h8[gn] = atoms
#            main_pdb_array['H8'] = new_h8
#            moved_h8 = True
#        else:
#            moved_h8 = False
#        self.statistics.add_info('moved_H8', moved_h8)
        
        # N- and C-termini
        if N_and_C_termini==True:
            N_struct = self.template_source['TM1'][list(self.template_source['TM1'])[0]][0]
            N_term = Residue.objects.filter(protein_conformation=self.prot_conf, protein_segment__slug='N-term')
            N_term_temp = Residue.objects.filter(protein_conformation=N_struct.protein_conformation,
                                                 protein_segment__slug='N-term')
            last_five = [i.sequence_number for i in list(N_term_temp)[-5:]]
            if self.main_structure==N_struct and len(last_five)==5:
                try:
                    temp_coo = list(parse.fetch_residues_from_pdb(N_struct,last_five).values())
                except:
                    temp_coo = None
            elif len(last_five)==5:
                try:
                    temp_nums = last_five + [i for i in range(last_five[-1]+1,last_five[-1]+5)]
                    template = parse.fetch_residues_from_pdb(N_struct,temp_nums)
                    ref_nums = list(main_pdb_array['TM1'])[:4]
                    reference = OrderedDict()
                    for i in ref_nums:
                        reference[i] = main_pdb_array['TM1'][i]
                    superpose = sp.OneSidedSuperpose(reference,template,4,0)
                    sup_residues = superpose.run()
                    n_count2 = 0
                    temp_coo = []
                    for num, atoms in sup_residues.items():
                        if n_count2<5:
                            temp_coo.append(atoms)
                        n_count2+=1
                except:
                    temp_coo = None
            else:
                temp_coo = None
            r_i, t_i, a_i, arr_i = OrderedDict(),OrderedDict(),OrderedDict(),OrderedDict()
            N_r, N_t, N_a, N_arr = OrderedDict(),OrderedDict(),OrderedDict(),OrderedDict()
            n_count = 0
            for n in N_term:
                n_count+=1
                N_r[str(n.sequence_number)] = n.amino_acid
                N_a[str(n.sequence_number)] = '-'
                if n_count>len(N_term)-5:
                    try:
                        N_arr[str(n.sequence_number)] = temp_coo[-1*(len(N_term)-n_count+1)]
                        N_t[str(n.sequence_number)] = list(N_term_temp)[-1*(len(N_term)-n_count+1)].amino_acid
                        self.template_source['N-term'][str(n.sequence_number)][0] = N_struct
                    except:
                        N_t[str(n.sequence_number)] = '-'
                        N_arr[str(n.sequence_number)] = '-'
                else:
                    N_t[str(n.sequence_number)] = '-'
                    N_arr[str(n.sequence_number)] = '-'
            r_i['N-term'] = N_r
            t_i['N-term'] = N_t
            a_i['N-term'] = N_a
            arr_i['N-term'] = N_arr
            for r,t,al,arr in zip(a.reference_dict,a.template_dict,a.alignment_dict,main_pdb_array):
                r_i[r]=a.reference_dict[r]
                t_i[t]=a.template_dict[t]
                a_i[al]=a.alignment_dict[al]
                arr_i[arr]=main_pdb_array[arr]     
            a.reference_dict = r_i
            a.template_dict = t_i
            a.alignment_dict = a_i
            main_pdb_array = arr_i
            try:
                index = -1
                while self.template_source['H8'][list(self.template_source['H8'])[index]][0]==None:
                    index-=1
                C_struct = self.template_source['H8'][list(self.template_source['H8'])[index]][0]
                last_seg = 'H8'
            except:
                C_struct = self.template_source['TM7'][list(self.template_source['TM7'])[-1]][0]
                last_seg = 'TM7'
            C_term = Residue.objects.filter(protein_conformation=self.prot_conf, protein_segment__slug='C-term')
            C_term_temp = Residue.objects.filter(protein_conformation=C_struct.protein_conformation,
                                                 protein_segment__slug='C-term')
                                                 
            first_five = [i.sequence_number for i in list(C_term_temp)[:5]]
            if self.main_structure==N_struct and len(first_five)==5:
                try:
                    temp_coo2 = list(parse.fetch_residues_from_pdb(C_struct,first_five).values())
                except:
                    temp_coo2 = None
            elif len(first_five)==5:
                try:
                    temp_nums2 = [i for i in range(first_five[0]-4,first_five[0])] + first_five
                    template2 = parse.fetch_residues_from_array(C_struct,temp_nums2)
                    ref_nums2 = list(main_pdb_array[last_seg])[-4:]
                    reference2 = OrderedDict()
                    for i in ref_nums2:
                        reference2[i] = main_pdb_array[last_seg][i]
                    superpose2 = sp.OneSidedSuperpose(reference2,template2,4,1)
                    sup_residues2 = superpose2.run()
                    c_count2 = 0
                    temp_coo2 = []
                    for num, atoms in sup_residues2.items():
                        if c_count2<5:
                            temp_coo2.append(atoms)
                        c_count2+=1
                except:
                    temp_coo2 = None
            else:
                temp_coo2 = None
            
            a.reference_dict['C-term'],a.template_dict['C-term'] = OrderedDict(),OrderedDict()
            a.alignment_dict['C-term'],main_pdb_array['C-term'] = OrderedDict(),OrderedDict()
            c_count = -1
            for c in C_term:
                c_count+=1
                a.reference_dict['C-term'][str(c.sequence_number)] = c.amino_acid
                a.alignment_dict['C-term'][str(c.sequence_number)] = '-'
                if c_count<4:
                    try:
                        main_pdb_array['C-term'][str(c.sequence_number)] = temp_coo2[c_count]    
                        a.template_dict['C-term'][str(c.sequence_number)] = list(C_term_temp)[c_count].amino_acid
                        self.template_source['C-term'][str(c.sequence_number)][0] = C_struct 
                    except:
                        a.template_dict['C-term'][str(c.sequence_number)] = '-'
                        main_pdb_array['C-term'][str(c.sequence_number)] = '-'
                else:
                    a.template_dict['C-term'][str(c.sequence_number)] = '-'
                    main_pdb_array['C-term'][str(c.sequence_number)] = '-'
        
        # Shorten N- and C-termini 
        n_count=1
        orig_n_len = len(a.reference_dict['N-term'])
        if orig_n_len>5:
            for num in a.reference_dict['N-term']:
                if n_count<orig_n_len-4:
                    del a.reference_dict['N-term'][num]
                    del a.template_dict['N-term'][num]
                    del a.alignment_dict['N-term'][num]
                    del main_pdb_array['N-term'][num]
                    del self.template_source['N-term'][num]
                n_count+=1
        c_count=1
        orig_c_len = len(a.reference_dict['C-term'])
        if orig_c_len>5:
            for num in a.reference_dict['C-term']:
                if c_count>5:
                    del a.reference_dict['C-term'][num]
                    del a.template_dict['C-term'][num]
                    del a.alignment_dict['C-term'][num]
                    del main_pdb_array['C-term'][num]
                    del self.template_source['C-term'][num]
                c_count+=1
        
        # Shorten ICL3
        try:
            if len(a.reference_dict['ICL3_free'])>10:
                icl3_c = 0
                keys = list(self.template_source['ICL3'].keys())
                length = len(a.template_dict['ICL3_free'])
                for r_s,t_s,a_s,ar_s in zip(a.reference_dict['ICL3_free'],a.template_dict['ICL3_free'],
                                            a.alignment_dict['ICL3_free'],main_pdb_array['ICL3_free']):
                    icl3_c+=1
                    if 5<icl3_c<length-4:
                        del a.reference_dict['ICL3_free'][r_s]
                        del a.template_dict['ICL3_free'][t_s]
                        del a.alignment_dict['ICL3_free'][a_s]
                        del main_pdb_array['ICL3_free'][ar_s]
                        del self.template_source['ICL3'][keys[icl3_c-1]]
        except:
            pass

        if len(a.reference_dict['N-term'])==0:
            del a.reference_dict['N-term']
            del a.template_dict['N-term']
            del a.alignment_dict['N-term']
            del main_pdb_array['N-term']        
        
        # Correct template_source for incosistencies in position labeling
        temp_source = OrderedDict()
        for s_seg,t_seg in zip(self.template_source,a.template_dict):
            temp_source_seg = OrderedDict()
            for s, t in zip(self.template_source[s_seg],a.template_dict[t_seg]):
                if s!=t and 'x' in t and t not in temp_source_seg:
                    temp_source_seg[t] = self.template_source[s_seg][s]
                else:
                    temp_source_seg[s] = self.template_source[s_seg][s]
            temp_source[s_seg] = temp_source_seg      
        self.template_source = temp_source
        
        # non-conserved residue switching
        if switch_rotamers==True:
            non_cons_switch = self.run_non_conserved_switcher(main_pdb_array,a.reference_dict,a.template_dict,
                                                              a.alignment_dict)
            main_pdb_array = non_cons_switch[0]
            a.reference_dict = non_cons_switch[1]
            a.template_dict = non_cons_switch[2]
            a.alignment_dict = non_cons_switch[3]
            trimmed_residues+=non_cons_switch[4]
        else:
            for seg_id, seg in main_pdb_array.items():
                for key in seg:
                    if a.reference_dict[seg_id][str(key).replace('.','x')]!='-':
                        trimmed_residues.append(key)
        if 'ICL4_free' in main_pdb_array:
            freeICL4=True
        else:
            freeICL4=False
        if freeICL4==True:
            for i in list(main_pdb_array['H8']):
                if i not in trimmed_residues:
                    trimmed_residues.append(i)
        for i,j in self.helix_end_mods['added'].items():
            try:
                if j[0][-1].replace('x','.') not in trimmed_residues:
                    trimmed_residues.append(j[0][-1].replace('x','.'))
            except:
                pass
            try:
                if j[1][0].replace('x','.') not in trimmed_residues:
                    trimmed_residues.append(j[1][0].replace('x','.'))
            except:
                pass
        
        for i in ref_bulge_list+temp_bulge_list+ref_const_list+temp_const_list:
            i = list(i.keys())[0].replace('x','.')
            if parse.gn_indecer(i,'.',-2) not in trimmed_residues:
                trimmed_residues.append(parse.gn_indecer(i,'.',-2))
            if parse.gn_indecer(i,'.',-1) not in trimmed_residues:
                trimmed_residues.append(parse.gn_indecer(i,'.',-1))
            if parse.gn_indecer(i,'.',1) not in trimmed_residues:
                trimmed_residues.append(parse.gn_indecer(i,'.',1))
            if parse.gn_indecer(i,'.',2) not in trimmed_residues:
                trimmed_residues.append(parse.gn_indecer(i,'.',2))
                
#        print('Rotamer switching: ',datetime.now() - startTime) 
        
        # write to file
        path = "./structure/homology_models/{}_{}/".format(self.reference_entry_name,self.state)
        if not os.path.exists(path):
            os.mkdir(path)
        trimmed_res_nums = self.write_homology_model_pdb(path+self.uniprot_id+"_post.pdb", main_pdb_array, 
                                                         trimmed_residues=trimmed_residues)
        
        self.statistics.add_info('template_source',self.template_source)
        
        # Model with MODELLER
        self.create_PIR_file(a, path+self.uniprot_id+"_post.pdb")
        
        self.run_MODELLER("./structure/PIR/"+self.uniprot_id+"_"+self.state+".pir", path+self.uniprot_id+"_post.pdb", 
                          self.uniprot_id, 1, "{}_{}_model.pdb".format(self.reference_entry_name,self.state), atom_dict=trimmed_res_nums)
        os.remove(path+self.uniprot_id+"_post.pdb")
        
        # stat file
#        with open('./structure/homology_models/{}_{}/{}.stat.txt'.format(self.reference_entry_name, self.state, 
#                                                                         self.reference_entry_name), 'w') as stat_file:
#            for label, info in self.statistics.items():
#                stat_file.write('{} : {}\n'.format(label, info))

        with open(path+'{}_{}.stats.txt'.format(self.reference_entry_name,self.state),'w') as s_file:
            rot_table = []
            sections = []
                
            for seg, resis in self.template_source.items():
                list_keys = list(resis)
                first_gn = list_keys[0]
                first_temp = self.template_source[seg][first_gn][0]
                if 'x' in first_gn:
                    try:
                        first_seqnum = Residue.objects.get(protein_conformation=self.prot_conf,generic_number__label=list_keys[0]).sequence_number
                    except:
                        first_seqnum = int(list_keys[0])
                else:
                    first_seqnum = int(first_gn)
                    first_gn = None
                for gn, res in resis.items():
                    key = gn
                    if 'x' in gn:
                        seq_num = Residue.objects.get(protein_conformation=self.prot_conf,generic_number__label=gn).sequence_number
                        curr_seqnum = seq_num
                    else:
                        seq_num = int(gn)
                        curr_seqnum = seq_num
                        gn = None
                    rot_table.append([seg,seq_num,gn,self.reference_entry_name,res[1]])
                    seqnum_minus = False
                    if res[0]!=first_temp:
                        if seq_num==first_seqnum:# or key==list_keys[-1]:
                            if gn!=None:
                                prev_gn = gn
                            else:
                                prev_gn = None
                            seq_num = int(seq_num)
                        else:
                            if gn!=None:
                                prev_gn = list_keys[list_keys.index(key)-1]
                            else:
                                prev_gn = None
                            seq_num = int(seq_num)-1
                            seqnum_minus = True
                        sections.append([seg,first_seqnum,seq_num,first_gn,prev_gn,self.reference_entry_name,first_temp])
                        if prev_gn==None:
                            first_gn = None
                        else:
                            first_gn = key
                        first_seqnum = curr_seqnum
                        first_temp = res[0]
                    if key==list_keys[-1]:
                        prev_gn = gn
                        if seqnum_minus==True:
                            seq_num = int(seq_num)+1
                        else:
                            seq_num = int(seq_num)
                        sections.append([seg,first_seqnum,seq_num,first_gn,prev_gn,self.reference_entry_name,first_temp])
                    
            for sec in sections:
                s_file.write(str(sec)+"\n")
                for rot in rot_table:
                    if int(sec[1])<=int(rot[1])<=int(sec[2]):
                        s_file.write(str(rot)+"\n")
        
        if self.update==True:
            print('UPDATING DB')
            self.upload_to_db(sections, rot_table)

        print('MODELLER build: ',datetime.now() - startTime)
        pprint.pprint(self.statistics)
        print('################################')
        return self
    
    def run_non_conserved_switcher(self, main_pdb_array, reference_dict, template_dict, alignment_dict):
        ''' Switches non-conserved residues with best possible template. Returns refreshed main_pdb_array 
            (atom coordinates), reference_dict (reference generic numbers and residue ids), template_dict (template 
            generic numbers and residue ids) and alignment_dict (aligned reference and template dictionary). 
            
            @param main_pdb_array: nested OrderedDict(), output of GPCRDBParsingPDB().pdb_array_creator()
            @param reference_dict: reference dictionary of AlignedReferenceTemplate.
            @param template_dict: template dictionary of AlignedReferenceTemplate.
            @param alignment_dict: alignment dictionary of AlignedReferenceTemplate.
        '''
        parse = GPCRDBParsingPDB()
        ref_length = 0
        conserved_count = 0
        non_cons_count = 0
        trimmed_res_num = 0
        switched_count = 0
        non_cons_res_templates, conserved_residues = OrderedDict(), OrderedDict()
        trimmed_residues = []
        inconsistencies = []
        
        for incons in self.statistics.info_dict['pdb_db_inconsistencies']:
            inconsistencies.append(list(incons.keys())[0])
        for ref_seg, temp_seg, aligned_seg in zip(reference_dict, template_dict, alignment_dict):
            if len(ref_seg)>4:
                segment = ref_seg[:4]
            else:
                segment = ref_seg
            for ref_res, temp_res, aligned_res in zip(reference_dict[ref_seg], template_dict[temp_seg], 
                                                      alignment_dict[aligned_seg]):
                if reference_dict[ref_seg][ref_res]!='-':
                    ref_length+=1
                else:
                    trimmed_residues.append(ref_res)
                if '?' in temp_res:
                    trimmed_residues.append(ref_res)
                    trimmed_res_num+=1
                    non_cons_count+=1
                    continue
                if '-term' in ref_seg and template_dict[temp_seg][temp_res]=='-':
                    trimmed_residues.append(ref_res)
                    trimmed_res_num+=1
                    non_cons_count+=1
                    continue
                if (ref_res not in inconsistencies and
                    alignment_dict[aligned_seg][aligned_res]!='.' and
                    alignment_dict[aligned_seg][aligned_res]!='x' and 
                    alignment_dict[aligned_seg][aligned_res]!='-'):
                    conserved_residues[ref_res] = alignment_dict[aligned_seg][aligned_res]
                    conserved_count+=1
                    if 'x' not in ref_res:
                        num_in_loop = parse.gn_num_extract(ref_res,'|')[1]
                        this_res = list(Residue.objects.filter(protein_conformation=self.prot_conf,
                                                               protein_segment__slug=segment))[num_in_loop-1]
                        seq_num = str(this_res.sequence_number)
                        try:
                            self.update_template_source([seq_num],self.template_source[segment][seq_num][0],segment,
                                                        just_rot=True)
                            key_in_template_source = seq_num
                        except:
                            self.update_template_source([this_res.generic_number.label],
                                                        self.template_source[segment][this_res.generic_number.label][0],
                                                        segment,just_rot=True)
                            key_in_template_source = this_res.generic_number.label
                    else:
                        self.update_template_source([ref_res],self.template_source[segment][ref_res][0],segment,
                                                    just_rot=True)
                        key_in_template_source = ref_res
                    if '_dis' in ref_seg or (ref_seg=='ECL2' and self.template_source['ECL2'][key_in_template_source][0]!=self.main_structure 
                                             and '|' in ref_res):
                        trimmed_residues.append(ref_res)
                gn = ref_res    
                if (gn in inconsistencies or alignment_dict[aligned_seg][aligned_res]=='.' and 
                    reference_dict[ref_seg][gn]!=template_dict[temp_seg][gn]):
                    non_cons_count+=1
                    gn_ = str(ref_res).replace('x','.')
                    no_match = True
                    if '|' in gn_:
                        try:
                            list_num = int(gn.split('|')[1])-1                       
                            gn = list(Residue.objects.filter(protein_conformation__protein=self.reference_protein,
                                      protein_segment__slug=ref_seg.split('_')[0]))[list_num].generic_number.label
                            gn_ = gn.replace('x','.')
                        except:
                            pass
                    for struct in self.similarity_table:
                        try:
                            alt_temp = parse.fetch_residues_from_pdb(struct, [gn])
                            if reference_dict[ref_seg][ref_res]==PDB.Polypeptide.three_to_one(
                                                                    alt_temp[gn_][0].get_parent().get_resname()):
                                orig_res = main_pdb_array[ref_seg][str(ref_res).replace('x','.')]
                                alt_res = alt_temp[gn_]
                                superpose = sp.RotamerSuperpose(orig_res, alt_res)
                                new_atoms = superpose.run()
                                if superpose.backbone_rmsd>0.5:
                                    continue
                                main_pdb_array[ref_seg][str(ref_res).replace('x','.')] = new_atoms
                                template_dict[temp_seg][temp_res] = reference_dict[ref_seg][ref_res]
                                non_cons_res_templates[gn] = struct
                                switched_count+=1
                                no_match = False
                                if 'x' not in ref_res:
                                    num_in_loop = parse.gn_num_extract(ref_res,'|')[1]
                                    seq_num = str(list(Residue.objects.filter(protein_conformation=self.prot_conf,
                                                                              protein_segment__slug=segment))[num_in_loop-1].sequence_number)
                                    self.update_template_source([seq_num],struct,segment,just_rot=True)
                                else:
                                    self.update_template_source([ref_res],struct,segment,just_rot=True)
                                break
                        except:
                            pass
                    if no_match==True:
                        try:
                            if 'free' not in ref_seg:
                                residue = main_pdb_array[ref_seg][str(ref_res).replace('x','.')]
                                main_pdb_array[ref_seg][str(ref_res).replace('x','.')] = residue[0:5]
                                trimmed_residues.append(gn_)
                                trimmed_res_num+=1
                            elif 'free' in ref_seg:
                                trimmed_residues.append(gn_)
                                trimmed_res_num+=1
                        except:
                            logging.warning("Missing atoms in {} at {}".format(self.main_structure,gn))

        self.statistics.add_info('ref_seq_length', ref_length)
        self.statistics.add_info('conserved_num', conserved_count)
        self.statistics.add_info('non_conserved_num', non_cons_count)
        self.statistics.add_info('trimmed_residues_num', trimmed_res_num)
        self.statistics.add_info('non_conserved_switched_num', switched_count)
        self.statistics.add_info('conserved_residues', conserved_residues)
        self.statistics.add_info('non_conserved_residue_templates', non_cons_res_templates)
        self.statistics.add_info('trimmed_residues', trimmed_residues)

        return [main_pdb_array, reference_dict, template_dict, alignment_dict, trimmed_residues]
    
    def write_homology_model_pdb(self, filename, main_pdb_array, trimmed_residues=[]):
        ''' Write PDB file from pdb array to file.
        
            @param filename: str, filename of output file \n
            @param main_pdb_array: OrderedDict(), of atoms of pdb, where keys are generic numbers/residue numbers and
            values are list of atoms. Output of GPCRDBParsingPDB.pdb_array_creator().
        '''
        key = ''
        res_num = 0
        atom_num = 0
        trimmed_resi_nums = OrderedDict()
        with open(filename,'w+') as f:
            for seg_id, segment in main_pdb_array.items():
                trimmed_segment = OrderedDict()
                for key in segment:
                    res_num+=1
                    if key in trimmed_residues:
                        trimmed_segment[key] = res_num
                        if 'x' in segment[key]:
                            f.write("\nTER")
                        if '?' in key:
                            continue
                    if 'x' in segment[key]:
                        f.write("\nTER")
                        continue
                    if '?' in key and '-' in segment[key]:
                        continue
                    if '-term' in seg_id and segment[key]=='-':
                        continue
                    for atom in main_pdb_array[seg_id][key]: 
                        atom_num+=1
                        coord = list(atom.get_coord())
                        coord1 = "%8.3f"% (coord[0])
                        coord2 = "%8.3f"% (coord[1])
                        coord3 = "%8.3f"% (coord[2])
                        if str(atom.get_id())=='CA':
                            if len(key)==4:
                                bfact = "%6.2f"% (float(key))
                            elif '.' not in key:
                                bfact = "%6.2f"% (float(atom.get_bfactor()))
                            else:
                                bfact = " -%4.2f"% (float(key))
                        else:
                            bfact = "%6.2f"% (float(atom.get_bfactor()))
                        occupancy = "%6.2f"% (atom.get_occupancy())
                        template="""
ATOM{atom_num}  {atom}{res} {chain}{res_num}{coord1}{coord2}{coord3}{occupancy}{bfactor}{atom_s}  """
                        context={"atom_num":str(atom_num).rjust(7), "atom":str(atom.get_id()).ljust(4),
                                 "res":atom.get_parent().get_resname(), 
                                 "chain":str(self.main_template_preferred_chain)[0],
                                 "res_num":str(res_num).rjust(4), "coord1":coord1.rjust(12), 
                                 "coord2":coord2.rjust(8), "coord3":coord3.rjust(8), 
                                 "occupancy":str(occupancy).rjust(3),
                                 "bfactor":str(bfact).rjust(4), "atom_s":str(str(atom.get_id())[0]).rjust(12)}
                        f.write(template.format(**context))
                trimmed_resi_nums[seg_id] = trimmed_segment
            f.write("\nTER\nEND")
        return trimmed_resi_nums
                    
    def create_PIR_file(self, ref_temp_alignment, template_file):
        ''' Create PIR file from reference and template alignment (AlignedReferenceAndTemplate).
        
            @param ref_temp_alignment: AlignedReferenceAndTemplate
            @template_file: str, name of template file with path
        '''
        ref_sequence, temp_sequence = '',''
        res_num = 0
        with open(template_file,'r') as f:
            lines = f.readlines()
            for line in lines:
                try:
                    pdb_re = re.search('(ATOM[A-Z\s\d]{13}\S{3}\s\S\s+)(\d+)([A-Z\s\d.-]{49,53})',line)
                    start_num = pdb_re.group(2)
                    break
                except:
                    pass
        for ref_seg, temp_seg in zip(ref_temp_alignment.reference_dict, ref_temp_alignment.template_dict):
            for ref_res, temp_res in zip(ref_temp_alignment.reference_dict[ref_seg], 
                                         ref_temp_alignment.template_dict[temp_seg]):
                res_num+=1
                if ref_temp_alignment.reference_dict[ref_seg][ref_res]=='x':
                    ref_sequence+='-'
                else:
                    ref_sequence+=ref_temp_alignment.reference_dict[ref_seg][ref_res]
                if ref_temp_alignment.template_dict[temp_seg][temp_res]=='x':
                    temp_sequence+='-'
                else:
                    temp_sequence+=ref_temp_alignment.template_dict[temp_seg][temp_res]
        with open("./structure/PIR/"+self.uniprot_id+"_"+self.state+".pir", 'w+') as output_file:
            template="""
>P1;{temp_file}
structure:{temp_file}:{start}:{chain}:{res_num}:{chain}::::
{temp_sequence}*

>P1;{uniprot}
sequence:{uniprot}::::::::
{ref_sequence}*
            """
            context={"temp_file":template_file,
                     "start":start_num,
                     "chain":self.main_template_preferred_chain,
                     "res_num":res_num,
                     "temp_sequence":temp_sequence,
                     "uniprot":self.uniprot_id,
                     "ref_sequence":ref_sequence}
            output_file.write(template.format(**context))
            
    def run_MODELLER(self, pir_file, template, reference, number_of_models, output_file_name, atom_dict=None):
        ''' Build homology model with MODELLER.
        
            @param pir_file: str, file name of PIR file with path \n
            @param template: str, file name of template with path \n
            @param reference: str, Uniprot code of reference sequence \n
            @param number_of_models: int, number of models to be built \n
            @param output_file_name: str, name of output file
        '''
        log.none()
        env = environ(rand_seed=80851) #!!random number generator
        
        if atom_dict==None:
            a = automodel(env, alnfile = pir_file, knowns = template, sequence = reference, 
                          assess_methods=(assess.DOPE))
        else:
            a = HomologyMODELLER(env, alnfile = pir_file, knowns = template, sequence = reference, 
                                 assess_methods=(assess.DOPE), atom_selection=atom_dict)
        
        a.starting_model = 1
        a.ending_model = number_of_models
        a.md_level = refine.slow
        path = "./structure/homology_models/{}_{}".format(self.reference_entry_name,self.state)
        if not os.path.exists(path):
            os.mkdir(path)
        a.make()

        # Get a list of all successfully built models from a.outputs
        ok_models = [x for x in a.outputs if x['failure'] is None]

        # Rank the models by DOPE score
        key = 'DOPE score'
        if sys.version_info[:2] == (2,3):
            # Python 2.3's sort doesn't have a 'key' argument
            ok_models.sort(lambda a,b: cmp(a[key], b[key]))
        else:
            ok_models.sort(key=lambda a: a[key])
        
        # Get top model
        m = ok_models[0]
#        print("Top model: %s (DOPE score %.3f)" % (m['name'], m[key]))        
        
        for file in os.listdir("./"):
            if file==m['name']:
                os.rename("./"+file, "./structure/homology_models/{}_{}/".format(self.reference_entry_name,
                                                                                 self.state)+output_file_name)
            elif file.startswith(self.uniprot_id):
                os.remove("./"+file)


class SilentModeller(object):
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, *args):
        sys.stdout.close()
        sys.stdout = self._stdout

        
class HomologyMODELLER(automodel):
    def __init__(self, env, alnfile, knowns, sequence, assess_methods, atom_selection):
        super(HomologyMODELLER, self).__init__(env, alnfile=alnfile, knowns=knowns, sequence=sequence, 
                                               assess_methods=assess_methods)
        self.atom_dict = atom_selection
        
    def select_atoms(self):
        selection_out = []
        for seg_id, segment in self.atom_dict.items():
            for gn, atom in segment.items():
                selection_out.append(self.residues[str(atom)])
        return selection(selection_out)
        
    def make(self):
        with SilentModeller():
            super(HomologyMODELLER, self).make()


class Loops(object):
    ''' Class to handle loops in GPCR structures.
    '''
    def __init__(self, reference_protein, loop_label, loop_template_structures, main_structure, helix_end_mods, segment_order):
        self.segment_order = segment_order
        self.reference_protein = reference_protein
        self.prot_conf = ProteinConformation.objects.get(protein=self.reference_protein)
        self.loop_label = loop_label
        self.loop_template_structures = loop_template_structures
        self.main_structure = main_structure
        self.helix_end_mods = helix_end_mods
        self.loop_output_structure = None
        self.new_label = None
        self.aligned = False
    
    def fetch_loop_residues(self, main_pdb_array, superpose_modded_loop=False):
        ''' Fetch list of Atom objects of the loop when there is an available template. Returns an OrderedDict().
        '''
        if (self.loop_label=='ECL2' and (self.loop_template_structures==None or 'ECL2_mid' in 
            self.loop_template_structures and self.loop_template_structures['ECL2_mid']==None)):
            return None
        if self.loop_template_structures!=None:
            parse = GPCRDBParsingPDB()            
            seg_list = self.segment_order
            prev_seg = seg_list[seg_list.index(self.loop_label)-1]
            next_seg = seg_list[seg_list.index(self.loop_label)+1]
            orig_before_gns = [i.replace('.','x') for i in list(main_pdb_array[prev_seg].keys())[-4:]]
            orig_after_gns = [j.replace('.','x') for j in list(main_pdb_array[next_seg].keys())[:4]]
            last_before_gn = orig_before_gns[-1]
            first_after_gn = orig_after_gns[0]
            if self.loop_label=='ECL2':
                try:
                    ref_res = Residue.objects.filter(protein_conformation__protein=self.reference_protein,
                                                     protein_segment__slug='ECL2')
                    r_first = list(ref_res)[0].sequence_number
                    r_last = list(ref_res)[-1].sequence_number
                    r_x50 = ref_res.get(generic_number__label='45x50').sequence_number
                except:
                    pass
            if (self.loop_label=='ECL2' and 'ECL2_1' not in self.loop_template_structures) or self.loop_label!='ECL2' or superpose_modded_loop==True:
                for template in self.loop_template_structures:
                    output = OrderedDict()
                    try:
                        if (template==self.main_structure or template=='aligned') and superpose_modded_loop==False:
                            if self.helix_end_mods!=None and (len(self.helix_end_mods['removed'][prev_seg][1])==0 and
                                len(self.helix_end_mods['removed'][next_seg][0])==0):
                                if template=='aligned':
                                    self.aligned = True
                                else:
                                    self.aligned = False
                                try:
                                    loop_res = [r.sequence_number for r in list(Residue.objects.filter(
                                                                                protein_conformation=self.main_structure.protein_conformation,
                                                                                protein_segment__slug=self.loop_label))]
                                    ref_loop_res = list(Residue.objects.filter(protein_conformation=self.prot_conf,
                                                                               protein_segment__slug=self.loop_label))
                                    at_least_one_gn = False
                                    try:
                                        for i in ref_loop_res:
                                            i.generic_number__label
                                            at_least_one_gn = True
                                    except:
                                        pass
                                    if at_least_one_gn==True:
                                        inter_array = parse.fetch_residues_from_pdb(self.main_structure,loop_res)
                                    else:
                                        inter_array = parse.fetch_residues_from_pdb(self.main_structure,loop_res,
                                                                                    just_nums=True)
                                    self.loop_output_structure = self.main_structure
                                    for id_, atoms in inter_array.items():
                                        output[str(id_)] = atoms
                                    return output
                                except:
                                    continue
                            else:
                                print('Warning: need to superpose aligned {}'.format(self.loop_label))
                                return self.fetch_loop_residues(main_pdb_array,superpose_modded_loop=True)
                        else:
                            prot_conf = ProteinConformation.objects.get(protein=self.reference_protein)
                            if self.loop_label=='ICL4' and len(list(Residue.objects.filter(protein_conformation=prot_conf,protein_segment__slug='ICL4')))<3:
                                raise Exception()
                            b_num = Residue.objects.get(protein_conformation=template.protein_conformation,
                                                        generic_number__label=last_before_gn).sequence_number
                            a_num = Residue.objects.get(protein_conformation=template.protein_conformation,
                                                        generic_number__label=first_after_gn).sequence_number
                            before4 = Residue.objects.filter(protein_conformation=template.protein_conformation, 
                                                             sequence_number__in=[b_num,b_num-1,b_num-2,b_num-3])
                            after4 = Residue.objects.filter(protein_conformation=template.protein_conformation, 
                                                             sequence_number__in=[a_num,a_num+1,a_num+2,a_num+3])
                            loop_residues = Residue.objects.filter(protein_conformation=template.protein_conformation,
                                                                   sequence_number__in=list(range(b_num+1,a_num)))
                            if len(loop_residues)!=len(Residue.objects.filter(protein_conformation=self.prot_conf,protein_segment__slug=self.loop_label)):
                                raise Exception()
                            before_gns = [x.sequence_number for x in before4]
                            mid_nums = [x.sequence_number for x in loop_residues]
                            after_gns = [x.sequence_number for x in after4]
                            alt_residues = parse.fetch_residues_from_pdb(template, before_gns+mid_nums+after_gns)
                            orig_residues = parse.fetch_residues_from_pdb(self.main_structure, 
                                                                          orig_before_gns+orig_after_gns)
                            superpose = sp.LoopSuperpose(orig_residues, alt_residues)
                            new_residues = superpose.run()
                            key_list = list(new_residues.keys())[4:-4]
                            for key in key_list:
                                output[key] = new_residues[key]
                            self.loop_output_structure = template
                            return output
                    except:
                        continue
            else:
                output,ECL2_1,ECL2_mid,ECL2_2 = OrderedDict(),OrderedDict(),OrderedDict(),OrderedDict()
                no_first_temp, no_second_temp = True,True
                main_temp_seq = Residue.objects.filter(protein_conformation=self.main_structure.protein_conformation,
                                                       protein_segment__slug=self.loop_label)
                for mid_template in self.loop_template_structures['ECL2_mid']:
                    if mid_template==self.main_structure:
                        ECL2_mid = parse.fetch_residues_from_pdb(self.main_structure,['45x50','45x51','45x52'])
                        x50 = main_temp_seq.get(generic_number__label='45x50').sequence_number
                        break
                o1 = parse.fetch_residues_from_array(main_pdb_array[prev_seg],orig_before_gns)
                o2 = parse.fetch_residues_from_pdb(self.main_structure,['45x50','45x51','45x52'])
                orig_residues1 = parse.add_two_ordereddict(o1,o2)

                if self.loop_template_structures['ECL2_1']==None:
                    no_first_temp=True
                else:
                    for first_temp in self.loop_template_structures['ECL2_1']:
                        if first_temp==self.main_structure:
                            try:
                                ECL2_1 = parse.fetch_residues_from_pdb(self.main_structure,list(range(list(main_temp_seq)[0].sequence_number,x50)))
                                no_first_temp=False
                                break
                            except:
                                continue
                        else:
                            try:
                                b_num = Residue.objects.get(protein_conformation=first_temp.protein_conformation,
                                                            generic_number__label=last_before_gn).sequence_number
                                before4 = Residue.objects.filter(protein_conformation=first_temp.protein_conformation, 
                                                                 sequence_number__in=[b_num,b_num-1,b_num-2,b_num-3])
                                alt_mid1 = Residue.objects.filter(protein_conformation=first_temp.protein_conformation,
                                                                  protein_segment__slug=self.loop_label, 
                                                                  generic_number__label__in=['45x50','45x51','45x52'])
                                alt1_x50 = alt_mid1.get(generic_number__label='45x50').sequence_number
                                loop_res1 = Residue.objects.filter(protein_conformation=first_temp.protein_conformation,
                                                                   sequence_number__in=list(range(b_num, alt1_x50)))
                                before_gns = [x.sequence_number for x in before4]
                                mid_gns1 = [x.sequence_number for x in loop_res1]
                                alt_residues1 = parse.fetch_residues_from_pdb(first_temp,before_gns+mid_gns1+['45x50','45x51','45x52'])
                                superpose = sp.LoopSuperpose(orig_residues1,alt_residues1,ECL2=True,part=1)
                                new_residues = superpose.run()
                                key_list = list(new_residues.keys())[4:-3]
                                for key in key_list:
                                    ECL2_1["1_"+key] = new_residues[key]
                                no_first_temp=False
                                break
                            except:
                                no_first_temp=True

                if no_first_temp==True:
                    for i in range(1,r_x50-r_first+1):
                        ECL2_1['1_'+str(i)]='x'
                    first_temp=None
                o11 = parse.fetch_residues_from_pdb(self.main_structure,['45x50','45x51','45x52'])
                o22 = parse.fetch_residues_from_array(main_pdb_array[next_seg],orig_after_gns)
                orig_residues2 = parse.add_two_ordereddict(o11,o22)
                if self.loop_template_structures['ECL2_2']==None:
                    no_second_temp=True
                else:
                    for second_temp in self.loop_template_structures['ECL2_2']:
                        if second_temp==self.main_structure:
                            try:
                                ECL2_2 = parse.fetch_residues_from_pdb(self.main_structure,list(range(x50+3,list(main_temp_seq)[-1].sequence_number+1)))
                                no_second_temp=False
                                break
                            except:
                                continue
                        else:
                            try:                                
                                a_num = Residue.objects.get(protein_conformation=second_temp.protein_conformation,
                                                            generic_number__label=first_after_gn).sequence_number
                                after4 = Residue.objects.filter(protein_conformation=second_temp.protein_conformation, 
                                                                sequence_number__in=[a_num,a_num+1,a_num+2,a_num+3])
                                alt_mid2 = Residue.objects.filter(protein_conformation=second_temp.protein_conformation,
                                                                  protein_segment__slug=self.loop_label, 
                                                                  generic_number__label__in=['45x50','45x51','45x52'])
                                alt2_x50 = alt_mid2.get(generic_number__label='45x50').sequence_number
                                loop_res2 = Residue.objects.filter(protein_conformation=second_temp.protein_conformation,
                                                                   sequence_number__in=list(range(alt2_x50+3, a_num)))
                                mid_gns2 = [x.sequence_number for x in loop_res2]
                                after_gns = [x.sequence_number for x in after4]
                                alt_residues2 = parse.fetch_residues_from_pdb(second_temp,['45x50','45x51','45x52']+mid_gns2+after_gns)
                                superpose = sp.LoopSuperpose(orig_residues2,alt_residues2,ECL2=True,part=2)
                                new_residues = superpose.run()
                                key_list = list(new_residues.keys())[3:-4]
                                for key in key_list:
                                    ECL2_2["2_"+key] = new_residues[key]
                                no_second_temp=False
                                break
                            except:
                                no_second_temp=True
                if no_second_temp==True:
                    for j in range(1,r_last-r_x50-1):
                        ECL2_2['2_'+str(j)]='x'
                    second_temp=None
                output['ECL2_1'] = ECL2_1
                output['ECL2_mid'] = ECL2_mid
                output['ECL2_2'] = ECL2_2
                self.loop_output_structure = [first_temp,mid_template,second_temp]
                return output
            if len(output.keys())==0:
                return None
        else:
            return None
                    
    def insert_loop_to_arrays(self, loop_output_structure, main_pdb_array, loop_template, reference_dict, 
                              template_dict, alignment_dict):
        ''' Updates the homology model with loop segments. Inserts previously fetched lists of loop Atom objects to 
            the proper arrays, dictionaries.
            
            @param loop_output_structure: Structure object of loop template.
            @param main_pdb_array: nested OrderedDict(), output of GPCRDBParsingPDB().pdb_array_creator().
            @param loop_template: OrderedDict() of loop template with lists of Atom objects as values.
            @param reference_dict: reference dictionary of AlignedReferenceTemplate.
            @param template_dict: template dictionary of AlignedReferenceTemplate.
            @param alignment_dict: alignment dictionary of AlignedReferenceTemplate.
        '''
        shorter_ref, shorter_temp = False, False        
        try:
            for r,t in zip(reference_dict[self.loop_label],template_dict[self.loop_label]):
                if reference_dict[self.loop_label][r]=='-':
                    shorter_ref = True
                elif template_dict[self.loop_label][t]=='-':
                    shorter_temp = True
        except:
            pass
        if loop_template!=None and loop_output_structure!=self.main_structure:
            loop_keys = list(loop_template.keys())[1:-1]
            continuous_loop = False
            self.main_pdb_array = self.discont_loop_insert_to_pdb(main_pdb_array, loop_template, loop_output_structure)           
        elif loop_template!=None and loop_output_structure==self.main_structure or self.aligned==True and (shorter_ref==True or shorter_temp==True):
            loop_keys = list(loop_template.keys())
            continuous_loop = True
            temporary_dict = OrderedDict()
            try:
                if len(loop_keys)<len(template_dict[self.loop_label]):
                    counter=0
                    for i in template_dict[self.loop_label]:
                        if i.replace('x','.') in loop_keys:
                            temporary_dict[i.replace('x','.')] = loop_template[i.replace('x','.')]
                        else:
                            temporary_dict['gap{}'.format(str(counter))] = '-'
                        counter+=1
                    loop_template = temporary_dict
            except:
                pass
            self.main_pdb_array = self.cont_loop_insert_to_pdb(main_pdb_array, loop_template)
        else:
            self.main_pdb_array = main_pdb_array
        
        if loop_template!=None:
            temp_ref_dict, temp_temp_dict, temp_aligned_dict = OrderedDict(),OrderedDict(),OrderedDict()
            if continuous_loop==True:
                if shorter_ref==True and shorter_temp==False:
                    ref_residues = list(reference_dict[self.loop_label].values())
                elif shorter_ref==True and shorter_temp==True:
                    ref_residues = list(reference_dict[self.loop_label].values())
                elif shorter_ref==False and shorter_temp==True:
                    ref_residues = list(reference_dict[self.loop_label].values())
                else:
                    ref_residues = [x.amino_acid for x in Residue.objects.filter(protein_conformation__protein=self.reference_protein,
                                                           protein_segment__slug=self.loop_label)]
            else:
                ref_residues = list(Residue.objects.filter(protein_conformation__protein=self.reference_protein,
                                                           protein_segment__slug=self.loop_label))
            for ref_seg, temp_seg, aligned_seg in zip(reference_dict, template_dict, alignment_dict):
                if ref_seg[0]=='T' and self.segment_order.index(self.loop_label)-self.segment_order.index(ref_seg[:4])==1:
                    temp_ref_dict[ref_seg] = reference_dict[ref_seg]
                    temp_temp_dict[temp_seg] = template_dict[temp_seg]
                    temp_aligned_dict[aligned_seg] = alignment_dict[aligned_seg]
                    input_residues = list(loop_template.keys())
                    ref_loop_seg, temp_loop_seg, aligned_loop_seg = OrderedDict(),OrderedDict(),OrderedDict()
                    if continuous_loop==True:
                        l_res=0
                        for r_res, r_id in zip(ref_residues, input_residues):
                            l_res+=1
                            try:
                                loop_gn = Residue.objects.get(protein_conformation=self.main_structure.protein_conformation, 
                                                              generic_number__label=r_id.replace('.','x')).generic_number.label
                            except:
                                try:
                                    Residue.objects.get(protein_conformation=self.main_structure.protein_conformation, 
                                                        sequence_number=r_id)
        #### Possible bug
                                    loop_gn = self.loop_label+'|'+str(l_res)
                                except:
                                    loop_gn = self.loop_label+'?'+str(l_res)
                            ref_loop_seg[loop_gn] = r_res
                            try:
                                temp_loop_seg[loop_gn] = PDB.Polypeptide.three_to_one(loop_template[r_id][0].get_parent().get_resname())
                            except:
                                temp_loop_seg[loop_gn] = '-'
                            if ref_loop_seg[loop_gn]==temp_loop_seg[loop_gn]:                        
                                aligned_loop_seg[loop_gn] = ref_loop_seg[loop_gn]
                            elif ref_loop_seg[loop_gn]=='-' or temp_loop_seg[loop_gn]=='-':
                                aligned_loop_seg[loop_gn] = '-'    
                            else:
                                aligned_loop_seg[loop_gn] = '.'    
                        self.new_label = self.loop_label+'_cont'
                        temp_ref_dict[self.loop_label+'_cont'] = ref_loop_seg
                        temp_temp_dict[self.loop_label+'_cont'] = temp_loop_seg
                        temp_aligned_dict[self.loop_label+'_cont'] = aligned_loop_seg
                    else:
                        l_res=1
                        ref_loop_seg[self.loop_label+'?'+'1'] = ref_residues[0].amino_acid
                        temp_loop_seg[self.loop_label+'?'+'1'] = 'x'
                        aligned_loop_seg[self.loop_label+'?'+'1'] = 'x'
                        for r_res, r_id in zip(ref_residues[1:-1], input_residues[1:-1]):
                            l_res+=1
                            try:
                                loop_gn = Residue.objects.get(protein_conformation=loop_output_structure.protein_conformation, 
                                                              sequence_number=r_id).generic_number.label
                            except:
                                loop_gn = self.loop_label+'|'+str(l_res)
                            ref_loop_seg[loop_gn] = r_res.amino_acid
                            temp_loop_seg[loop_gn] = PDB.Polypeptide.three_to_one(loop_template[r_id][0].get_parent().get_resname())
                            if ref_loop_seg[loop_gn]==temp_loop_seg[loop_gn]:                        
                                aligned_loop_seg[loop_gn] = ref_loop_seg[loop_gn]
                            else:
                                aligned_loop_seg[loop_gn] = '.'
                        ref_loop_seg[self.loop_label+'?'+str(l_res+1)] = ref_residues[-1].amino_acid
                        temp_loop_seg[self.loop_label+'?'+str(l_res+1)] = 'x'
                        aligned_loop_seg[self.loop_label+'?'+str(l_res+1)] = 'x'
                        self.new_label = self.loop_label+'_dis'
                        temp_ref_dict[self.loop_label+'_dis'] = ref_loop_seg
                        temp_temp_dict[self.loop_label+'_dis'] = temp_loop_seg
                        temp_aligned_dict[self.loop_label+'_dis'] = aligned_loop_seg
                else:
                    temp_ref_dict[ref_seg] = reference_dict[ref_seg]
                    temp_temp_dict[temp_seg] = template_dict[temp_seg]
                    temp_aligned_dict[aligned_seg] = alignment_dict[aligned_seg]
            self.reference_dict = temp_ref_dict
            self.template_dict = temp_temp_dict
            self.alignment_dict = temp_aligned_dict
            try:
                del self.reference_dict[self.loop_label]
                del self.template_dict[self.loop_label]
                del self.alignment_dict[self.loop_label]
            except:
                pass
        else:
            self.reference_dict = reference_dict
            self.template_dict = template_dict
            self.alignment_dict = alignment_dict
            try:
                del self.reference_dict[self.loop_label]
                del self.template_dict[self.loop_label]
                del self.alignment_dict[self.loop_label]
            except:
                pass
        return self
    
    def insert_ECL2_to_arrays(self, loop_output_structure, main_pdb_array, loop_template, reference_dict, 
                              template_dict, alignment_dict):
        temp_array = OrderedDict()
        # first part
        if loop_output_structure[0]!=None:
            if loop_output_structure[0]==self.main_structure:
                temp_array = self.cont_loop_insert_to_pdb(main_pdb_array, loop_template['ECL2_1'], ECL2='')
            else:
                temp_array = self.discont_loop_insert_to_pdb(main_pdb_array, loop_template['ECL2_1'], 
                                                             loop_output_structure, ECL2='')
        else:
            temp_array = self.gap_ECL2(main_pdb_array,loop_template['ECL2_1'],break_chain=True)
        # middle part
        for key, res in loop_template['ECL2_mid'].items():
            temp_array['ECL2'][key] = res
        # second part
        l_res = len(temp_array['ECL2'])
        if loop_output_structure[2]!=None:
            if loop_output_structure[2]==self.main_structure:
                for key, res in loop_template['ECL2_2'].items():
                    l_res+=1
                    if '.' in key:
                        temp_array['ECL2'][key] = res
                    else:
                        temp_array['ECL2'][self.loop_label+'|'+str(l_res)] = res
            else:
                loop_keys = list(loop_template['ECL2_2'].keys())[1:-1]
                temp_array['ECL2'][self.loop_label+'?'+str(l_res+1)] = 'x'
                l_res+=1
                for key in loop_keys:
                    l_res+=1
                    temp_array['ECL2'][self.loop_label+'|'+str(l_res)] = loop_template['ECL2_2'][key]
                temp_array['ECL2'][self.loop_label+'?'+str(l_res+1)] = 'x'
        else:
            for key, res in loop_template['ECL2_2'].items():
                l_res+=1
                temp_array['ECL2'][self.loop_label+'?'+str(l_res)] = '-'
        self.main_pdb_array = temp_array
        temp_ref_dict, temp_temp_dict, temp_aligned_dict = OrderedDict(),OrderedDict(),OrderedDict()
        ref_residues = list(Residue.objects.filter(protein_conformation__protein=self.reference_protein, 
                                                   protein_segment__slug='ECL2'))
        for ref_seg, temp_seg, aligned_seg in zip(reference_dict, template_dict, alignment_dict):
            if ref_seg[0]=='T' and self.segment_order.index(self.loop_label)-self.segment_order.index(ref_seg[:4])==1:
                temp_ref_dict[ref_seg] = reference_dict[ref_seg]
                temp_temp_dict[temp_seg] = template_dict[temp_seg]
                temp_aligned_dict[aligned_seg] = alignment_dict[aligned_seg]
                temp_ref_dict['ECL2'],temp_temp_dict['ECL2'],temp_aligned_dict['ECL2'] = OrderedDict(),OrderedDict(),OrderedDict()
                for ref, key in zip(ref_residues, self.main_pdb_array['ECL2']):
                    temp_ref_dict['ECL2'][key.replace('.','x')] = ref.amino_acid
                    try:
                        temp_temp_dict['ECL2'][key.replace('.','x')] = PDB.Polypeptide.three_to_one(
                                                        self.main_pdb_array['ECL2'][key][0].get_parent().get_resname())
                    except:
                        temp_temp_dict['ECL2'][key.replace('.','x')] = self.main_pdb_array['ECL2'][key.replace('.','x')]
                    if temp_ref_dict['ECL2'][key.replace('.','x')]==temp_temp_dict['ECL2'][key.replace('.','x')]:
                        temp_aligned_dict['ECL2'][key.replace('.','x')] = temp_ref_dict['ECL2'][key.replace('.','x')]
                    elif temp_temp_dict['ECL2'][key.replace('.','x')]=='x':
                        temp_aligned_dict['ECL2'][key.replace('.','x')] = 'x'
                    elif temp_temp_dict['ECL2'][key.replace('.','x')]=='-':
                        temp_aligned_dict['ECL2'][key.replace('.','x')] = '-'
                    else:
                        temp_aligned_dict['ECL2'][key.replace('.','x')] = '.'
            else:
                temp_ref_dict[ref_seg] = reference_dict[ref_seg]
                temp_temp_dict[temp_seg] = template_dict[temp_seg]
                temp_aligned_dict[aligned_seg] = alignment_dict[aligned_seg]
        self.reference_dict = temp_ref_dict
        self.template_dict = temp_temp_dict
        self.alignment_dict = temp_aligned_dict       
        return self       
        
    def gap_ECL2(self, main_pdb_array, loop_template, break_chain=False):
        temp_array, temp_loop = OrderedDict(), OrderedDict()
        for seg_label, gns in main_pdb_array.items():
            if self.segment_order.index(self.loop_label)-self.segment_order.index(seg_label[:4])==1:
                temp_array[seg_label] = gns
                l_res = 0
                for key in loop_template:
                    l_res+=1
                    temp_loop[self.loop_label+'?'+str(l_res)] = '-'
                temp_array[self.loop_label] = temp_loop
            else:
                temp_array[seg_label] = gns
        return temp_array
                
    def cont_loop_insert_to_pdb(self, main_pdb_array, loop_template, ECL2=None):
        temp_array, temp_loop = OrderedDict(), OrderedDict()
        for seg_label, gns in main_pdb_array.items():
            if self.segment_order.index(self.loop_label)-self.segment_order.index(seg_label[:4])==1:
                temp_array[seg_label] = gns
                l_res = 0
                for key in loop_template:
                    l_res+=1
                    if '.' in key:
                        temp_loop[key] = loop_template[key]
                    elif 'gap' in key:
                        temp_loop[self.loop_label+'?'+str(l_res)] = loop_template[key]
                    else:
                        temp_loop[self.loop_label+'|'+str(l_res)] = loop_template[key]   
                if ECL2!=None:
                    temp_array[self.loop_label] = temp_loop
                else:                             
                    temp_array[self.loop_label+'_cont'] = temp_loop
            else:
                temp_array[seg_label] = gns
        return temp_array
        
    def discont_loop_insert_to_pdb(self, main_pdb_array, loop_template, loop_output_structure, ECL2=None):
        temp_array, temp_loop = OrderedDict(), OrderedDict()
        loop_keys = list(loop_template.keys())[1:-1]
        for seg_label, gns in main_pdb_array.items():
            if self.segment_order.index(self.loop_label)-self.segment_order.index(seg_label[:4])==1:
                temp_array[seg_label] = gns
                l_res = 1
                temp_loop[self.loop_label+'?'+'1'] = 'x'
                for key in loop_keys:
                    l_res+=1
                    try:
                        loop_gn = Residue.objects.get(protein_conformation=loop_output_structure.protein_conformation, 
                                                      sequence_number=key).generic_number.label.replace('x','.')
                        temp_loop[loop_gn] = loop_template[key]
                    except:
                        temp_loop[self.loop_label+'|'+str(l_res)] = loop_template[key]
                temp_loop[self.loop_label+'?'+str(l_res+1)] = 'x'
                if ECL2!=None:
                    temp_array[self.loop_label] = temp_loop
                else:                    
                    temp_array[self.loop_label+'_dis'] = temp_loop
            else:
                temp_array[seg_label] = gns
        return temp_array
        
    def insert_gaps_for_loops_to_arrays(self, main_pdb_array, reference_dict, template_dict, alignment_dict):
        ''' When there is no template for a loop region, this function inserts gaps for that region into the main 
            template, fetches the reference residues and inserts these into the arrays. This allows for Modeller to
            freely model these loop regions.
            
            @param main_pdb_array: nested OrderedDict(), output of GPCRDBParsingPDB().pdb_array_creator().
            @param reference_dict: reference dictionary of AlignedReferenceTemplate.
            @param template_dict: template dictionary of AlignedReferenceTemplate.
            @param alignment_dict: alignment dictionary of AlignedReferenceTemplate.
        '''
        residues = Residue.objects.filter(protein_conformation__protein=self.reference_protein, 
                                          protein_segment__slug=self.loop_label)
        temp_pdb_array = OrderedDict()
        for seg_id, seg in main_pdb_array.items():
            if self.segment_order.index(self.loop_label)-self.segment_order.index(seg_id[:4])==1:
                temp_loop = OrderedDict()
                count=0
                temp_pdb_array[seg_id] = seg
                for r in residues:
                    count+=1
                    temp_loop[self.loop_label+'?'+str(count)] = '-'
                temp_pdb_array[self.loop_label+'_free'] = temp_loop
                self.new_label = self.loop_label+'_free'
            else:
                temp_pdb_array[seg_id] = seg
        self.main_pdb_array = temp_pdb_array
        temp_ref_dict, temp_temp_dict, temp_aligned_dict = OrderedDict(), OrderedDict(), OrderedDict()
        for ref_seg, temp_seg, aligned_seg in zip(reference_dict, template_dict, alignment_dict):
            if self.segment_order.index(self.loop_label)-self.segment_order.index(ref_seg[:4])==1:
                temp_ref_loop, temp_temp_loop, temp_aligned_loop = OrderedDict(), OrderedDict(), OrderedDict()
                temp_ref_dict[ref_seg] = reference_dict[ref_seg]
                temp_temp_dict[temp_seg] = template_dict[temp_seg]
                temp_aligned_dict[aligned_seg] = alignment_dict[aligned_seg]
                count=0
                for r in residues:
                    count+=1
                    temp_ref_loop[self.loop_label+'?'+str(count)] = r.amino_acid
                    temp_temp_loop[self.loop_label+'?'+str(count)] = '-'
                    temp_aligned_loop[self.loop_label+'?'+str(count)] = '.'
                temp_ref_dict[self.loop_label+'_free'] = temp_ref_loop
                temp_temp_dict[self.loop_label+'_free'] = temp_temp_loop
                temp_aligned_dict[self.loop_label+'_free'] = temp_aligned_loop
            else:
                temp_ref_dict[ref_seg] = reference_dict[ref_seg]
                temp_temp_dict[temp_seg] = template_dict[temp_seg]
                temp_aligned_dict[aligned_seg] = alignment_dict[aligned_seg]
        self.reference_dict = temp_ref_dict
        self.template_dict = temp_temp_dict
        self.alignment_dict = temp_aligned_dict


class Bulges(object):
    ''' Class to handle bulges in GPCR structures.
    '''
    def __init__(self, gn):
        self.gn = gn
        self.bulge_templates = []
        self.template = None
    
    def find_bulge_template(self, similarity_table, bulge_in_reference):
        ''' Searches for bulge template, returns residues of template (5 residues if the bulge is in the reference, 4
            residues if the bulge is in the template). 
            
            @param gn: str, Generic number of bulge, e.g. 1x411 \n
            @param similarity_table: OrderedDict(), table of structures ordered by preference.
            Output of HomologyModeling().create_similarity_table(). \n
            @param bulge_in_reference: boolean, Set it to True if the bulge is in the reference, set it to False if the
            bulge is in the template.
        '''
        gn = self.gn
        parse = GPCRDBParsingPDB()
        for structure, value in similarity_table.items():
            this_anomaly = ProteinAnomaly.objects.filter(generic_number__label=gn)
            anomaly_list = structure.protein_anomalies.all().prefetch_related()
            if bulge_in_reference==True:
                try:
                    for anomaly in this_anomaly:
                        if anomaly in anomaly_list:
                            gn_list = [parse.gn_indecer(gn,'x',-2),parse.gn_indecer(gn,'x',-1),gn,
                                       parse.gn_indecer(gn,'x',+1),parse.gn_indecer(gn,'x',+2)]
                            alt_bulge = parse.fetch_residues_from_pdb(structure, gn_list)
                            self.template = structure
                            return alt_bulge
                except:
                    pass
            elif bulge_in_reference==False:
                try:
                    suitable_temp = []
                    for anomaly in this_anomaly:
                        if anomaly not in anomaly_list:
                            pass
                        else:
                            suitable_temp.append('no')
                    if 'no' not in suitable_temp:
                        gn_list = [parse.gn_indecer(gn,'x',-2), parse.gn_indecer(gn,'x',-1),
                                   parse.gn_indecer(gn,'x',+1), parse.gn_indecer(gn,'x',+2)]
                        alt_bulge = parse.fetch_residues_from_pdb(structure, gn_list)
                        self.template = structure
                        return alt_bulge
                except:
                    pass
        return None
            
            
class Constrictions(object):
    ''' Class to handle constrictions in GPCRs.
    '''
    def __init__(self, gn):
        self.gn = gn
        self.constriction_templates = []
        self.template = None
    
    def find_constriction_template(self, similarity_table, constriction_in_reference):
        ''' Searches for constriction template, returns residues of template (4 residues if the constriction is in the 
            reference, 5 residues if the constriction is in the template). 
            
            @param gn: str, Generic number of constriction, e.g. 7x44 \n
            @param similarity_table: OrderedDict(), table of structures ordered by preference.
            Output of HomologyModeling().create_similarity_table(). \n
            @param constriction_in_reference: boolean, Set it to True if the constriction is in the reference, set it 
            to False if the constriction is in the template.
        '''
        gn = self.gn
        parse = GPCRDBParsingPDB()
        for structure, value in similarity_table.items():
            this_anomaly = ProteinAnomaly.objects.filter(generic_number__label=gn)
            anomaly_list = structure.protein_anomalies.all().prefetch_related()
            if constriction_in_reference==True:
                try:
                    for anomaly in this_anomaly:
                        if anomaly in anomaly_list:
                            gn_list = [parse.gn_indecer(gn,'x',-2),parse.gn_indecer(gn,'x',-1),
                                       parse.gn_indecer(gn,'x',+1),parse.gn_indecer(gn,'x',+2)]
                            alt_const = parse.fetch_residues_from_pdb(structure, gn_list)
                            self.template = structure
                            return alt_const
                except:
                    pass
            elif constriction_in_reference==False:
                try:
                    suitable_temp = []
                    for anomaly in this_anomaly:
                        if anomaly not in anomaly_list:
                            pass
                        else:
                            suitable_temp.append('no')
                    if 'no' not in suitable_temp:
                        gn_list = [parse.gn_indecer(gn,'x',-2), parse.gn_indecer(gn,'x',-1),gn,
                                   parse.gn_indecer(gn,'x',+1), parse.gn_indecer(gn,'x',+2)]
                        alt_const = parse.fetch_residues_from_pdb(structure, gn_list)
                        self.template = structure
                        return alt_const
                except:
                    pass              
        return None
        
        
class GPCRDBParsingPDB(object):
    ''' Class to manipulate cleaned pdb files of GPCRs.
    '''
    def __init__(self):
        self.segment_coding = OrderedDict([(1,'TM1'),(2,'TM2'),(3,'TM3'),(4,'TM4'),(5,'TM5'),(6,'TM6'),(7,'TM7'),(8,'H8')])
    
    def gn_num_extract(self, gn, delimiter):
        ''' Extract TM number and position for formatting.
        
            @param gn: str, Generic number \n
            @param delimiter: str, character between TM and position (usually 'x')
        '''
        try:
            split = gn.split(delimiter)
            return int(split[0]), int(split[1])
        except:
            try:
                split = gn.split(delimiter)
                return split[0], int(split[1])
            except:
                return '/', '/'
            
    def gn_comparer(self, gn1, gn2, protein_conformation):
        '''
        '''
        res1 = Residue.objects.get(protein_conformation=protein_conformation, generic_number__label=gn1)
        res2 = Residue.objects.get(protein_conformation=protein_conformation, generic_number__label=gn2)
        return res1.sequence_number-res2.sequence_number
            
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

    def fetch_residues_from_pdb(self, structure, generic_numbers, modify_bulges=False, just_nums=False):
        ''' Fetches specific lines from pdb file by generic number (if generic number is
            not available then by residue number). Returns nested OrderedDict()
            with generic numbers as keys in the outer dictionary, and atom names as keys
            in the inner dictionary.
            
            @param structure: Structure, Structure object where residues should be fetched from \n
            @param generic_numbers: list, list of generic numbers to be fetched \n
            @param modify_bulges: boolean, set it to true when used for bulge switching. E.g. you want a 5x461
            residue to be considered a 5x46 residue. 
        '''
        output = OrderedDict()
        atoms_list = []
        for gn in generic_numbers:
            rotamer=None
            if 'x' in str(gn):      
                rotamer = list(Rotamer.objects.filter(structure__protein_conformation=structure.protein_conformation, 
                        residue__generic_number__label=gn, structure__preferred_chain=structure.preferred_chain))
            else:
                rotamer = list(Rotamer.objects.filter(structure__protein_conformation=structure.protein_conformation, 
                        residue__sequence_number=gn, structure__preferred_chain=structure.preferred_chain))
                if just_nums==False:
                    try:
                        gn = Residue.objects.get(protein_conformation=structure.protein_conformation,
                                                 sequence_number=gn).generic_number.label
                    except:
                        pass
            if len(rotamer)>1:
                for i in rotamer:
                    if i.pdbdata.pdb.startswith('COMPND')==False:
                        rotamer = i
                        break
            else:
                rotamer = rotamer[0]
            io = StringIO(rotamer.pdbdata.pdb)
            rota_struct = PDB.PDBParser(QUIET=True).get_structure('structure', io)[0]
            for chain in rota_struct:
                for residue in chain:
                    for atom in residue:
                        atoms_list.append(atom)
                    if modify_bulges==True and len(gn)==5:
                        output[gn.replace('x','.')[:-1]] = atoms_list
                    else:
                        try:
                            output[gn.replace('x','.')] = atoms_list
                        except:
                            output[str(gn)] = atoms_list
                    atoms_list = []
        return output
        
    def fetch_residues_from_array(self, main_pdb_array_segment, list_of_gns):
        array = OrderedDict()
        for i in list_of_gns:
            array[i.replace('x','.')] = main_pdb_array_segment[i.replace('x','.')]
        return array
        
    def add_two_ordereddict(self, dict1, dict2):
        for i,j in dict2.items():
            dict1[i] = j
        return dict1

    def pdb_array_creator(self, structure=None, filename=None):
        ''' Creates an OrderedDict() from the pdb of a Structure object where residue numbers/generic numbers are 
            keys for the residues, and atom names are keys for the Bio.PDB.Residue objects.
            
            @param structure: Structure, Structure object of protein. When using structure, leave filename=None. \n
            @param filename: str, filename of pdb to be parsed. When using filename, leave structure=None).
        '''
        if structure!=None and filename==None:
            io = StringIO(structure.pdb_data.pdb)
        else:
            io = filename
        gn_array = []
        residue_array = []
        pdb_struct = PDB.PDBParser(QUIET=True).get_structure('structure', io)[0]

        assign_gn = as_gn.GenericNumbering(structure=pdb_struct)
        pdb_struct = assign_gn.assign_generic_numbers()
        
        pref_chain = structure.preferred_chain
        for residue in pdb_struct[pref_chain]:
            try:
                if -9.1 < residue['CA'].get_bfactor() < 9.1:
                    gn = str(residue['CA'].get_bfactor())
                    if len(gn.split('.')[1])==1:
                        gn = gn+'0'
                    if gn[0]=='-':
                        gn = gn[1:]+'1'
                    gn_array.append(gn)
                    residue_array.append(residue.get_list())
                else:
                    gn_array.append(str(residue.get_id()[1]))
                    residue_array.append(residue.get_list())
            except:
                logging.warning("Unable to parse {} in {}".format(residue, structure))
        output = OrderedDict()
        for num, label in self.segment_coding.items():
            output[label] = OrderedDict()
        counter=0
        if len(gn_array)!=len(residue_array):
            raise AssertionError()
        for gn, res in zip(gn_array,residue_array):
            if '.' in gn:
                seg_num = int(gn.split('.')[0])
                seg_label = self.segment_coding[seg_num]
                if seg_num==8 and len(output['TM7'])==0:
                    continue
                else:
                    output[seg_label][gn] = res
            else:
                try:
                    found_res = Residue.objects.get(protein_conformation=structure.protein_conformation,
                                                    sequence_number=gn)
                    found_gn = str(found_res.generic_number.label).replace('x','.')
                    if -9.1 < float(found_gn) < 9.1:
                        seg_label = self.segment_coding[int(found_gn.split('.')[0])]
                        output[seg_label][found_gn] = res
                except:
                    pass
            counter+=1
        return output
   
   
class CreateStatistics(object):
    ''' Statistics dictionary for HomologyModeling.
    '''
    def __init__(self, reference):
        self.reference = reference
        self.info_dict = OrderedDict()
    
    def __repr__(self):
        return "<{} \n {} \n>".format(self.reference, self.info_dict)
        
    def items(self):
        ''' Returns the OrderedDict().items().
        '''
        return self.info_dict.items()
    
    def add_info(self, info_name, info):
        ''' Adds new information to the statistics dictionary.
        
            @param info_name: str, info name as dictionary key
            @param info: object, any object as value
        '''
        self.info_dict[info_name] = info

