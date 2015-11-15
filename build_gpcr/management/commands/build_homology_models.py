from django.core.management.base import BaseCommand

from protein.models import Protein, ProteinSegment, ProteinConformation, ProteinAnomaly, ProteinState
from residue.models import Residue
from structure.models import Structure, PdbData, Rotamer, StructureModel, StructureModelLoopTemplates, StructureModelAnomalies, StructureModelResidues
from common.alignment import Alignment, AlignedReferenceTemplate
import structure.structural_superposition as sp
import structure.assign_generic_numbers_gpcr as as_gn
from structure.calculate_RMSD import Validation

import Bio.PDB as PDB
from modeller import *
from modeller.automodel import *
from collections import OrderedDict
import os
import logging
import numpy as np
from io import StringIO
import sys
import multiprocessing
import pprint
import re
from datetime import datetime


startTime = datetime.now()
l = multiprocessing.Lock()

def homology_model_multiprocessing(receptor):
    Homology_model = HomologyModeling(receptor, 'Inactive', ['Inactive'])
    alignment = Homology_model.run_alignment()
    if alignment!=None:
        Homology_model.build_homology_model(alignment)#, switch_bulges=False, switch_constrictions=False, switch_rotamers=False)    
        Homology_model.upload_to_db()
        logger = logging.getLogger('homology_modeling')
        l.acquire()
        logger.info('Model for {} successfully built.'.format(receptor))
        l.release()
        
class Command(BaseCommand):    
    def handle(self, *args, **options):
      
        receptor_list = [ 'gpr15_human']
        if os.path.isfile('./logs/homology_modeling.log'):
            os.remove('./logs/homology_modeling.log')
        logger = logging.getLogger('homology_modeling')
        hdlr = logging.FileHandler('./logs/homology_modeling.log')
        formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
        hdlr.setFormatter(formatter)
        logger.addHandler(hdlr) 
        logger.setLevel(logging.INFO)        
        
        pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
        for i in receptor_list:
            pool.apply_async(homology_model_multiprocessing, [i])           
        pool.close()
        pool.join()

        print('\n###############################')
        print('Total runtime: ',datetime.now() - startTime)
        print('###############################\n')
        

class HomologyModeling(object):
    ''' Class to build homology models for GPCRs. 
    
        @param reference_entry_name: str, protein entry name \n
        @param state: str, endogenous ligand state of reference \n
        @param query_states: list, list of endogenous ligand states to be applied for template search, 
        default: same as reference
    '''
    segment_coding = {1:'TM1',2:'TM2',3:'TM3',4:'TM4',5:'TM5',6:'TM6',7:'TM7'}
    def __init__(self, reference_entry_name, state, query_states):
        self.reference_entry_name = reference_entry_name
        self.state = state
        self.query_states = query_states
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
        self.logger = logging.getLogger('homology_modeling')
        l.acquire()
        self.logger.info('Building model for {} {}.'.format(self.reference_protein, self.state))
        l.release()        
        
    def __repr__(self):
        return "<{}, {}>".format(self.reference_entry_name, self.state)

    def upload_to_db(self):
        # upload StructureModel        
        state=ProteinState.objects.get(name=self.state)
        hommod, created = StructureModel.objects.update_or_create(protein=self.reference_protein, state=state, 
                                                                  main_template=self.main_structure, 
                                                                  pdb=self.format_final_model())
                                                                  
        # upload StructureModelLoopTemplates
        for loop,template in self.statistics.info_dict['loops'].items():
            seg = ProteinSegment.objects.get(slug=loop[:4])
            StructureModelLoopTemplates.objects.update_or_create(homology_model=hommod,template=template,segment=seg)
            
        # upload StructureModelAnomalies
        ref_bulges = self.statistics.info_dict['reference_bulges']
        temp_bulges = self.statistics.info_dict['template_bulges']
        ref_const = self.statistics.info_dict['reference_constrictions']
        temp_const = self.statistics.info_dict['template_constrictions']

        if ref_bulges!=[]:        
            for r_b in ref_bulges:
                print(r_b.keys()[0], r_b.values()[0])
                
        # upload StructureModelResidues
        for gn, res in self.statistics.info_dict['conserved_residues'].items():
            if gn[0] not in ['E','I']:
                res = Residue.objects.get(protein_conformation__protein=self.reference_protein, 
                                          generic_number__label=gn)
                res_temp = Residue.objects.get(protein_conformation=self.main_structure.protein_conformation, 
                                               generic_number__label=gn)
            else:
                res = Residue.objects.filter(protein_conformation__protein=self.reference_protein, 
                                             protein_segment__slug=gn.split('|')[0])[int(gn.split('|')[1])-1]
                res_temp = Residue.objects.filter(protein_conformation=self.main_structure.protein_conformation, 
                                                  protein_segment__slug=gn.split('|')[0])[int(gn.split('|')[1])-1]
            rotamer = Rotamer.objects.filter(structure=self.main_structure, residue=res_temp)
            rotamer = self.right_rotamer_select(rotamer)
            StructureModelResidues.objects.update_or_create(homology_model=hommod, sequence_number=res.sequence_number,
                                                            residue=res, rotamer=rotamer, template=self.main_structure,
                                                            origin='conserved', segment=res.protein_segment)            
        for gn, temp in self.statistics.info_dict['non_conserved_residue_templates'].items():
            res = Residue.objects.get(protein_conformation__protein=self.reference_protein, generic_number__label=gn)
            res_temp = Residue.objects.get(protein_conformation=self.main_structure.protein_conformation,
                                           generic_number__label=gn)
            rotamer = Rotamer.objects.filter(structure=self.main_structure, residue=res_temp)
            rotamer = self.right_rotamer_select(rotamer)
            StructureModelResidues.objects.update_or_create(homology_model=hommod, sequence_number=res.sequence_number, 
                                                            residue=res, rotamer=rotamer, template=temp, 
                                                            origin='switched', segment=res.protein_segment)
        for gn in self.statistics.info_dict['trimmed_residues']:
            if gn[0] not in ['E','I']:
                gn = gn.replace('.','x')
                res = Residue.objects.get(protein_conformation__protein=self.reference_protein, 
                                          generic_number__label=gn)
            else:
                gn = gn.replace('?','|')
                res = Residue.objects.filter(protein_conformation__protein=self.reference_protein, 
                                             protein_segment__slug=gn.split('|')[0])[int(gn.split('|')[1])-1]
            StructureModelResidues.objects.update_or_create(homology_model=hommod, sequence_number=res.sequence_number,
                                                            residue=res, rotamer__isnull=True, template__isnull=True,
                                                            origin='free', segment=res.protein_segment)
                                   
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
        self.starting_res_num = list(Residue.objects.filter(protein_segment=2, 
                                     protein_conformation__protein=self.reference_protein))[0].sequence_number
        resnum = self.starting_res_num
        with open ('./structure/homology_models/{}_{}/modeller_test.pdb'.format(self.uniprot_id, self.state), 'r+') as f:
            pdblines = f.readlines()
            out_list = []
            prev_num = 1
            for line in pdblines:
                try:
                    pdb_re = re.search('(ATOM[A-Z\s\d]+\S{3}\s+)(\d+)([A-Z\s\d.-]+)',line)
                    if int(pdb_re.group(2))>prev_num:
                        resnum+=1
                        prev_num = int(pdb_re.group(2))
                    whitespace = (len(str(resnum))-len(pdb_re.group(2)))*-1
                    if whitespace==0:
                        out_line = pdb_re.group(1)+str(resnum)+pdb_re.group(3)
                    else:
                        out_line = pdb_re.group(1)[:whitespace]+str(resnum)+pdb_re.group(3)
                    out_list.append(out_line)
                except:
                    out_list.append(line)
#            io = StringIO(''.join(out_list))
#            pdb_struct = PDB.PDBParser(PERMISSIVE=True).get_structure('structure', io)[0]
#            assign_gn = as_gn.GenericNumbering(structure=pdb_struct)
#            pdb_struct = assign_gn.assign_generic_numbers()
#            assign_gn.save_gn_to_pdb()
#            outio = PDB.PDBIO()
#            outio.set_structure(pdb_struct)
#            outio.save('./structure/homology_models/{}_{}/modeller_test_ready.pdb'.format(self.uniprot_id, self.state))
#        with open('./structure/homology_models/{}_{}/modeller_test_GPCRDB.pdb'.format(self.uniprot_id, self.state),'r+') as f:
#            pdbdata = f.read()
        return ''.join(out_list)

    def run_alignment(self, core_alignment=True, query_states='default', 
                      segments=['TM1','TM2','TM3','TM4','TM5','TM6','TM7'], order_by='similarity'):
        ''' Creates pairwise alignment between reference and target receptor(s).
            Returns Alignment object.
            
            @param segments: list, list of segments to use, e.g.: ['TM1','IL1','TM2','EL1'] \n
            @param order_by: str, order results by identity, similarity or simscore
        '''
        if query_states=='default':
            query_states=self.query_states
        alignment = AlignedReferenceTemplate(self.reference_protein, segments, query_states, order_by)       
        if core_alignment==True:
            print('Alignment: ',datetime.now() - startTime)
            enhanced_alignment = alignment.enhance_best_alignment()
            print('Enhanced alignment: ',datetime.now() - startTime)
            if enhanced_alignment==None:
                return None
            self.segments = segments
            self.main_structure = alignment.main_template_structure
            self.similarity_table = alignment.similarity_table
            self.similarity_table_all = self.run_alignment(core_alignment=False, 
                                                           query_states=['Inactive','Active']).similarity_table
            self.main_template_preferred_chain = str(self.main_structure.preferred_chain)[0]
            self.statistics.add_info("main_template", self.main_structure)
            self.statistics.add_info("preferred_chain", self.main_template_preferred_chain)
            for loop in ['ICL1','ECL1','ICL2','ECL2','ICL3','ECL3']:
                loop_alignment = AlignedReferenceTemplate(self.reference_protein, [loop], ['Inactive','Active'], 
                                                          order_by='similarity', 
                                                          provide_main_template_structure=self.main_structure,
                                                          provide_similarity_table=self.similarity_table_all)
                self.loop_template_table[loop] = loop_alignment.loop_table
            self.statistics.add_info('similarity_table', self.similarity_table)
            self.statistics.add_info('loops',self.loop_template_table)
            print('Loop alignment: ',datetime.now() - startTime)
        return alignment
        
    def build_homology_model(self, ref_temp_alignment, switch_bulges=True, switch_constrictions=True, loops=True, 
                             switch_rotamers=True):
        ''' Function to identify and switch non conserved residues in the alignment. Optionally,
            it can identify and switch bulge and constriction sites too. 
            
            @param ref_temp_alignment: AlignedReferenceAndTemplate, alignment of reference and main template with 
            alignment string. \n
            @param switch_bulges: boolean, identify and switch bulge sites. Default = True.
            @param switch_constrictions: boolean, identify and switch constriction sites. Default = True.
        '''
        a = ref_temp_alignment
        ref_bulge_list, temp_bulge_list, ref_const_list, temp_const_list = [],[],[],[]
        parse = GPCRDBParsingPDB()
        main_pdb_array = parse.pdb_array_creator(structure=self.main_structure)
        
        print('Create main_pdb_array: ',datetime.now() - startTime)
        # loops
        if loops==True:
            loop_stat = OrderedDict()
            for label, structures in self.loop_template_table.items():
                loop = Loops(self.reference_protein, label, structures, self.main_structure)
                loop_template = loop.fetch_loop_residues()
                loop_insertion = loop.insert_loop_to_arrays(loop.loop_output_structure, main_pdb_array, loop_template, 
                                                            a.reference_dict, a.template_dict, a.alignment_dict)
                main_pdb_array = loop_insertion.main_pdb_array
                a.reference_dict = loop_insertion.reference_dict
                a.template_dict = loop_insertion.template_dict
                a.alignment_dict = loop_insertion.alignment_dict
                if loop.new_label!=None:
                    loop_stat[loop.new_label] = loop.loop_output_structure
                else:
                    loop_stat[label] = loop.loop_output_structure
            self.statistics.add_info('loops', loop_stat)
        print('Integrate loops: ',datetime.now() - startTime)   
        # bulges and constrictions
        if switch_bulges==True or switch_constrictions==True:
            for ref_seg, temp_seg, aligned_seg in zip(a.reference_dict, a.template_dict, a.alignment_dict):
                for ref_res, temp_res, aligned_res in zip(a.reference_dict[ref_seg], a.template_dict[temp_seg], 
                                                          a.alignment_dict[aligned_seg]):
                    gn = ref_res
                    gn_TM = parse.gn_num_extract(gn, 'x')[0]
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
        print('Integrate bulges/constrictions: ',datetime.now() - startTime)
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
        path = "./structure/homology_models/{}_{}/".format(self.uniprot_id,self.state)
        if not os.path.exists(path):
            os.mkdir(path)
        self.write_homology_model_pdb(
                                "./structure/homology_models/{}_{}/pre_switch.pdb".format(self.uniprot_id, self.state), 
                                main_pdb_array, a)        
        print('Check inconsistencies: ',datetime.now() - startTime)
        # inserting loops for free modeling
        for label, template in loop_stat.items():
            if template==None:
                modeling_loops = Loops(self.reference_protein, label, self.similarity_table_all, self.main_structure)
                modeling_loops.insert_gaps_for_loops_to_arrays(main_pdb_array, a.reference_dict, a.template_dict,
                                                               a.alignment_dict)
                main_pdb_array = modeling_loops.main_pdb_array
                a.reference_dict = modeling_loops.reference_dict
                a.template_dict = modeling_loops.template_dict
                a.alignment_dict = modeling_loops.alignment_dict
        print('Free loops: ',datetime.now() - startTime)
        # non-conserved residue switching
        if switch_rotamers==True:
            non_cons_switch = self.run_non_conserved_switcher(main_pdb_array,a.reference_dict,a.template_dict,
                                                              a.alignment_dict)
            main_pdb_array = non_cons_switch[0]
            a.reference_dict = non_cons_switch[1]
            a.template_dict = non_cons_switch[2]
            a.alignment_dict = non_cons_switch[3]
            trimmed_residues = non_cons_switch[4]
        else:
            trimmed_residues=[]
            for seg_id, seg in main_pdb_array.items():
                for key in seg:
                    if a.reference_dict[seg_id][str(key).replace('.','x')]!='-':
                        trimmed_residues.append(key)
        print('Rotamer switching: ',datetime.now() - startTime)
        # write to file
        path = "./structure/homology_models/{}_{}/".format(self.uniprot_id,self.state)
        if not os.path.exists(path):
            os.mkdir(path)
        trimmed_res_nums = self.write_homology_model_pdb(path+self.uniprot_id+"_post.pdb", main_pdb_array, 
                                                         a, trimmed_residues=trimmed_residues)                                                         

        # Model with MODELLER
        self.create_PIR_file(a, path+self.uniprot_id+"_post.pdb")
        self.run_MODELLER("./structure/PIR/"+self.uniprot_id+"_"+self.state+".pir", path+self.uniprot_id+"_post.pdb", 
                          self.uniprot_id, 1, "modeller_test.pdb", atom_dict=trimmed_res_nums)
        
        with open('./structure/homology_models/{}_Inactive/{}.stat.txt'.format(self.uniprot_id, self.uniprot_id), 'w') as stat_file:
            for label, info in self.statistics.items():
                stat_file.write('{} : {}\n'.format(label, info))
                
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
        res_list = []
        for ref_seg, temp_seg, aligned_seg in zip(reference_dict, template_dict, alignment_dict):
            for ref_res, temp_res, aligned_res in zip(reference_dict[ref_seg], template_dict[temp_seg], 
                                                      alignment_dict[aligned_seg]):
                if reference_dict[ref_seg][ref_res]!='-' and reference_dict[ref_seg][ref_res]!='/':
                    ref_length+=1
                if (alignment_dict[aligned_seg][aligned_res]!='.' and
                    alignment_dict[aligned_seg][aligned_res]!='x' and 
                    alignment_dict[aligned_seg][aligned_res]!='-'):
                    conserved_count+=1
                    conserved_residues[ref_res] = alignment_dict[aligned_seg][aligned_res]
                
                gn = ref_res    
                if (alignment_dict[aligned_seg][aligned_res]=='.' and 
                    reference_dict[ref_seg][gn]!=template_dict[temp_seg][gn]) and reference_dict[ref_seg][gn]!='x':
                    non_cons_count+=1
                    gn_ = str(ref_res).replace('x','.')
                    no_match = True
                    for struct in self.similarity_table:
                        try:
                            alt_temp = parse.fetch_residues_from_pdb(struct, [gn])
                            if reference_dict[ref_seg][gn]==PDB.Polypeptide.three_to_one(
                                                                    alt_temp[gn_][0].get_parent().get_resname()):
                                orig_res = main_pdb_array[ref_seg][gn_]
                                alt_res = parse.fetch_residues_from_pdb(struct,[gn])[gn_]
                                superpose = sp.RotamerSuperpose(orig_res, alt_res)
                                new_atoms = superpose.run()
                                if superpose.backbone_rmsd>0.3:
                                    continue
                                main_pdb_array[ref_seg][gn_] = new_atoms
                                template_dict[temp_seg][gn] = reference_dict[ref_seg][gn]
                                switched_count+=1                     
                                non_cons_res_templates[ref_res] = struct
                                res=Residue.objects.get(protein_conformation__protein=self.reference_protein, generic_number__label=ref_res)
                                res_list.append(res.sequence_number)
                                no_match = False
                                break
                        except:
                            pass
                    if no_match==True:
                        try:
                            if 'free' not in ref_seg:
                                residue = main_pdb_array[ref_seg][gn_]
                                main_pdb_array[ref_seg][gn_] = residue[0:5]
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
    
    def write_homology_model_pdb(self, filename, main_pdb_array, ref_temp_alignment, trimmed_residues=[]):
        ''' Write PDB file from pdb array to file.
        
            @param filename: str, filename of output file \n
            @param main_pdb_array: OrderedDict(), of atoms of pdb, where keys are generic numbers/residue numbers and
            values are list of atoms. Output of GPCRDBParsingPDB.pdb_array_creator().
            @param ref_temp_alignment: AlignedReferenceAndTemplate, only writes residues that are in ref_temp_alignment.
        '''
        key = ''
#        self.starting_res_num = list(Residue.objects.filter(protein_segment=2, protein_conformation__protein=self.reference_protein))[0].sequence_number
#        res_num = self.starting_res_num-1
        res_num = 0
        atom_num = 0
        trimmed_resi_nums = OrderedDict()
        with open(filename,'w+') as f:
            for seg_id, segment in main_pdb_array.items():
                trimmed_segment = OrderedDict()
                for key in segment:
                    if str(key).replace('.','x') :#in ref_temp_alignment.reference_dict[seg_id]:
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
structure:{temp_file}:1:{chain}:{res_num}:{chain}::::
{temp_sequence}*

>P1;{uniprot}
sequence:{uniprot}::::::::
{ref_sequence}*
            """
            context={"temp_file":template_file,
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
        path = "./structure/homology_models/{}".format(reference+"_"+self.state)
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
                os.rename("./"+file, "./structure/homology_models/{}_{}/".format(self.uniprot_id,
                                                                                 self.state)+output_file_name)
            elif file.startswith(self.uniprot_id):
                os.remove("./"+file)#, "./structure/homology_models/{}_{}/".format(self.uniprot_id,self.state)+file)


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
                print(self.residues[str(atom)], str(atom))
                selection_out.append(self.residues[str(atom)])
        return selection(selection_out)
        
    def make(self):
        with SilentModeller():
            super(HomologyMODELLER, self).make()


class Loops(object):
    ''' Class to handle loops in GPCR structures.
    '''
    def __init__(self, reference_protein, loop_label, loop_template_structures, main_structure):
        self.segment_order = {'TM1':1, 'ICL1':1.5, 'TM2':2, 'ECL1':2.5, 'TM3':3, 'ICL2':3.5, 'TM4':4, 'ECL2':4.5, 
                              'TM5':5, 'ICL3':5.5, 'TM6':6, 'ECL3':6.5, 'TM7':7}
        self.reference_protein = reference_protein
        self.loop_label = loop_label
        self.loop_template_structures = loop_template_structures
        self.main_structure = main_structure
        self.loop_output_structure = None
        self.new_label = None
    
    def fetch_loop_residues(self):
        ''' Fetch list of Atom objects of the loop when there is an available template. Returns an OrderedDict().
        '''
        if self.loop_template_structures!=None:
            parse = GPCRDBParsingPDB()
            segment = ProteinSegment.objects.get(slug=self.loop_label)
            orig_before_residues = Residue.objects.filter(
                                        protein_conformation=self.main_structure.protein_conformation, 
                                        protein_segment__id=segment.id-1)
            orig_after_residues = Residue.objects.filter(
                                        protein_conformation=self.main_structure.protein_conformation, 
                                        protein_segment__id=segment.id+1)
            orig_before_gns = []
            orig_after_gns = []
            for res in orig_before_residues.reverse():
                if len(orig_before_gns)<4:
                    orig_before_gns.append(res.generic_number.label)
            for res in orig_after_residues:
                if len(orig_after_gns)<4:
                    orig_after_gns.append(res.generic_number.label)
            orig_before_gns = list(reversed(orig_before_gns))
            last_before_gn = orig_before_gns[-1]
            first_after_gn = orig_after_gns[0]
            for template in self.loop_template_structures:
                output = OrderedDict()
                try:
                    if template==self.main_structure:
                        try:
                            loop_res = [r.sequence_number for r in list(Residue.objects.filter(
                                                                        protein_conformation=template.protein_conformation,
                                                                        protein_segment__slug=self.loop_label))]
                            inter_array = parse.fetch_residues_from_pdb(template,loop_res)
                            self.loop_output_structure = self.main_structure
                            for id_, atoms in inter_array.items():
                                output[str(id_)] = atoms
                            return output
                        except:
                            continue
                    else:
                        b_num = Residue.objects.get(protein_conformation=template.protein_conformation,
                                                    generic_number__label=last_before_gn).sequence_number
                        a_num = Residue.objects.get(protein_conformation=template.protein_conformation,
                                                    generic_number__label=first_after_gn).sequence_number
                        before8 = Residue.objects.filter(protein_conformation=template.protein_conformation, 
                                                         sequence_number__in=[b_num,b_num-1,b_num-2,b_num-3])
#                                                                              b_num-4,b_num-5,b_num-6,b_num-7])
                        after8 = Residue.objects.filter(protein_conformation=template.protein_conformation, 
                                                         sequence_number__in=[a_num,a_num+1,a_num+2,a_num+3])
#                                                                              a_num+4,a_num+5,a_num+6,a_num+7])
                        loop_residues = Residue.objects.filter(protein_conformation=template.protein_conformation,
                                                               sequence_number__in=list(range(b_num+1,a_num)))
                        before_gns = [x.sequence_number for x in before8]
                        mid_nums = [x.sequence_number for x in loop_residues]
                        after_gns = [x.sequence_number for x in after8]
                        alt_residues_temp = parse.fetch_residues_from_pdb(template, before_gns+mid_nums+after_gns)
                        alt_residues = OrderedDict()
                        for id_, atoms in alt_residues_temp.items():
                            if '.' not in str(id_):
                                alt_residues[str(id_)] = atoms
                            else:
                                alt_residues[id_] = atoms                            
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
        temp_array = OrderedDict()
        temp_loop = OrderedDict()
        if loop_template!=None and loop_output_structure!=self.main_structure:
            loop_keys = list(loop_template.keys())[1:-1]
            continuous_loop = False
            for seg_label, gns in main_pdb_array.items():
                if self.segment_order[self.loop_label]-self.segment_order[seg_label[:4]]==0.5:
                    temp_array[seg_label] = gns
                    l_res = 1
                    temp_loop[self.loop_label+'?'+'1'] = 'x'
                    for key in loop_keys:
                        l_res+=1
                        temp_loop[self.loop_label+'|'+str(l_res)] = loop_template[key]
                    temp_loop[self.loop_label+'?'+str(l_res+1)] = 'x'
                    temp_array[self.loop_label+'_dis'] = temp_loop
                else:
                    temp_array[seg_label] = gns
            self.main_pdb_array = temp_array
            
        elif loop_template!=None and loop_output_structure==self.main_structure:
            loop_keys = list(loop_template.keys())
            continuous_loop = True
            for seg_label, gns in main_pdb_array.items():
                if self.segment_order[self.loop_label]-self.segment_order[seg_label[:4]]==0.5:
                    temp_array[seg_label] = gns
                    l_res = 0
                    for key in loop_keys:
                        l_res+=1
                        temp_loop[self.loop_label+'|'+str(l_res)] = loop_template[key]
                    temp_array[self.loop_label+'_cont'] = temp_loop
                else:
                    temp_array[seg_label] = gns
            self.main_pdb_array = temp_array          
        else:
            self.main_pdb_array = main_pdb_array
        if loop_template!=None:
            temp_ref_dict, temp_temp_dict, temp_aligned_dict = OrderedDict(),OrderedDict(),OrderedDict()
            ref_residues = list(Residue.objects.filter(protein_conformation__protein=self.reference_protein, 
                                                  protein_segment__slug=self.loop_label))
            for ref_seg, temp_seg, aligned_seg in zip(reference_dict, template_dict, alignment_dict):
                if ref_seg[0]=='T' and self.segment_order[self.loop_label]-self.segment_order[ref_seg[:4]]==0.5:
                    temp_ref_dict[ref_seg] = reference_dict[ref_seg]
                    temp_temp_dict[temp_seg] = template_dict[temp_seg]
                    temp_aligned_dict[aligned_seg] = alignment_dict[aligned_seg]
                    input_residues = list(loop_template.keys())
                    ref_loop_seg, temp_loop_seg, aligned_loop_seg = OrderedDict(),OrderedDict(),OrderedDict()
                    if continuous_loop==True:
                        l_res=0
                        for r_res, r_id in zip(ref_residues, input_residues):
                            l_res+=1
                            ref_loop_seg[self.loop_label+'|'+str(l_res)] = r_res.amino_acid
                            temp_loop_seg[self.loop_label+'|'+str(l_res)] = PDB.Polypeptide.three_to_one(loop_template[r_id][0].get_parent().get_resname())
                            if ref_loop_seg[self.loop_label+'|'+str(l_res)]==temp_loop_seg[self.loop_label+'|'+str(l_res)]:                        
                                aligned_loop_seg[self.loop_label+'|'+str(l_res)] = ref_loop_seg[self.loop_label+'|'+str(l_res)]
                            else:
                                aligned_loop_seg[self.loop_label+'|'+str(l_res)] = '.'    
                        self.new_label = self.loop_label+'_cont'
                        temp_ref_dict[self.loop_label+'_cont'] = ref_loop_seg
                        temp_temp_dict[self.loop_label+'_cont'] = temp_loop_seg
                        temp_aligned_dict[self.loop_label+'_cont'] = aligned_loop_seg
                    else:
                        l_res=1
                        ref_loop_seg[self.loop_label+'?'+'1'] = ref_residues[0].amino_acid
                        temp_loop_seg[self.loop_label+'?'+'1'] = 'x'
                        aligned_loop_seg[self.loop_label+'?'+'1'] = '.'
                        for r_res, r_id in zip(ref_residues[1:-1], input_residues[1:-1]):
                            l_res+=1
                            ref_loop_seg[self.loop_label+'|'+str(l_res)] = r_res.amino_acid
                            temp_loop_seg[self.loop_label+'|'+str(l_res)] = PDB.Polypeptide.three_to_one(loop_template[r_id][0].get_parent().get_resname())
                            if ref_loop_seg[self.loop_label+'|'+str(l_res)]==temp_loop_seg[self.loop_label+'|'+str(l_res)]:                        
                                aligned_loop_seg[self.loop_label+'|'+str(l_res)] = ref_loop_seg[self.loop_label+'|'+str(l_res)]
                            else:
                                aligned_loop_seg[self.loop_label+'|'+str(l_res)] = '.'
                        ref_loop_seg[self.loop_label+'?'+str(l_res+1)] = ref_residues[-1].amino_acid
                        temp_loop_seg[self.loop_label+'?'+str(l_res+1)] = 'x'
                        aligned_loop_seg[self.loop_label+'?'+str(l_res+1)] = '.'
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
        else:
            self.reference_dict = reference_dict
            self.template_dict = template_dict
            self.alignment_dict = alignment_dict
        return self
            
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
            if self.segment_order[self.loop_label]-self.segment_order[seg_id[:4]]==0.5:
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
            if self.segment_order[self.loop_label]-self.segment_order[ref_seg[:4]]==0.5:
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
        self.segment_coding = {1:'TM1',2:'TM2',3:'TM3',4:'TM4',5:'TM5',6:'TM6',7:'TM7'}
    
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

    def fetch_residues_from_pdb(self, structure, generic_numbers, modify_bulges=False):
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
            if len(rotamer)>1:
                for i in rotamer:
                    if i.pdbdata.pdb.startswith('COMPND')==False:
                        rotamer = i
                        break
            else:
                rotamer = rotamer[0]
            io = StringIO(rotamer.pdbdata.pdb)
            rota_struct = PDB.PDBParser().get_structure('structure', io)[0]
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
                            output[gn] = atoms_list
                    atoms_list = []
        return output

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
        residue_array = OrderedDict()
        pdb_struct = PDB.PDBParser(PERMISSIVE=True).get_structure('structure', io)[0]

        assign_gn = as_gn.GenericNumbering(structure=pdb_struct)
        pdb_struct = assign_gn.assign_generic_numbers()
        
        pref_chain = structure.preferred_chain
        for residue in pdb_struct[pref_chain]:
            try:
#                print(residue,residue['CA'].get_bfactor(), residue['N'].get_bfactor())
                if -8.1 < residue['CA'].get_bfactor() < 8.1:
                    gn = str(residue['CA'].get_bfactor())
                    
                    if gn[0]=='-':
                        gn = gn[1:]+'1'
                    elif len(gn)==3:
                        gn = gn+'0'
                    residue_array[gn] = residue.get_list()
                else:
                    residue_array[str(residue.get_id()[1])] = residue.get_list()
            except:
                logging.warning("Unable to parse {} in {}".format(residue, structure))
        
        output = OrderedDict()
        for num, label in self.segment_coding.items():
            output[label] = OrderedDict()
        counter=0
        for gn, res in residue_array.items():
            if '.' in gn:
                seg_label = self.segment_coding[int(gn[0])]
                output[seg_label][gn] = res
            else:
                try:
                    found_res = Residue.objects.get(protein_conformation__protein=structure.protein_conformation.protein.parent,
                                            sequence_number=gn)
                    found_gn = str(found_res.generic_number.label).replace('x','.')
                    if -8.1 < float(found_gn) < 8.1:
                        seg_label = self.segment_coding[int(found_gn[0])]
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

