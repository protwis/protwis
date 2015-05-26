from django.core.management.base import BaseCommand

from protein.models import Protein, ProteinSegment, ProteinConformation
from residue.models import Residue
from structure.models import Structure, PdbData, Rotamer
from common.alignment import Alignment, AlignedReferenceTemplate
import structure.structural_superposition as sp
import structure.assign_generic_numbers_gpcr as as_gn
from build.management.commands.build_structures import Command as rota

import Bio.PDB as PDB
from modeller import *
from modeller.automodel import *
from collections import OrderedDict
import os
import logging
import numpy as np
from io import StringIO
import re
import pprint


class Command(BaseCommand):
    
    def handle(self, *args, **options):
        count=0
        s = [struct.protein_conformation.protein.parent.entry_name for struct in Structure.objects.all()]
        for protein in ProteinConformation.objects.all():
            if protein.protein.entry_name in s and count < 0:
                Homology_model = HomologyModeling(protein.protein.entry_name, 'Inactive', ['Inactive'])
                alignment = Homology_model.run_alignment()
                Homology_model.run_non_conserved_switcher(alignment)
    
                self.stdout.write(Homology_model.statistics, ending='')
                count+=1

        Homology_model = HomologyModeling('gp139_human', 'Inactive', ['Inactive'])
        alignment = Homology_model.run_alignment()
        Homology_model.run_non_conserved_switcher(alignment)        
                    
        val = Validation()
        struct = Structure.objects.get(protein_conformation__protein__name="4ib4")
        print(val.PDB_RMSD(StringIO(struct.pdb_data.pdb),
                "./structure/homology_models/{}_Inactive/pre_switch.pdb".format(Homology_model.uniprot_id),
                assign_gns=[1]))
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
        self.reference_protein = Protein.objects.get(entry_name=self.reference_entry_name)
        self.uniprot_id = self.reference_protein.accession
        self.reference_sequence = self.reference_protein.sequence
        self.statistics.add_info('uniprot_id',self.uniprot_id)
        self.segments = []
        self.similarity_table = OrderedDict()
        self.similarity_table_all = OrderedDict()
        self.main_structure = None
        self.main_template_preferred_chain = ''
        if os.path.isfile('./structure/homology_modeling.log'):
            os.remove('./structure/homology_modeling.log')
        logging.basicConfig(filename='./structure/homology_modeling.log',level=logging.WARNING)
        
    def __repr__(self):
        return "<{}, {}>".format(self.reference_entry_name, self.state)
   
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
        alignment.enhance_best_alignment()
        if core_alignment==True:
            self.segments = segments
            self.main_structure = alignment.main_template_structure
            self.similarity_table = alignment.similarity_table
            self.main_template_preferred_chain = str(self.main_structure.preferred_chain)[0]
            self.statistics.add_info("main_template", self.main_structure)
            self.statistics.add_info("preferred_chain", self.main_template_preferred_chain)
            loops = OrderedDict()
            for loop in ['ICL1','ECL1','ICL2','ECL2','ICL3','ECL3']:
                loop_alignment = AlignedReferenceTemplate(self.reference_protein, [loop], ['Inactive','Active'], 
                                                          order_by='similarity', 
                                                          provide_main_template_structure=self.main_structure)
                loops[loop] = loop_alignment.main_template_structure
        return alignment
        
    def run_non_conserved_switcher(self, ref_temp_alignment, switch_bulges=True, switch_constrictions=True):
        ''' Function to identify and switch non conserved residues in the alignment. Optionally,
            it can identify and switch bulge and constriction sites too. 
            
            @param ref_temp_alignment: AlignedReferenceAndTemplate, alignment of reference and main template with 
            alignment string. \n
            @param switch_bulges: boolean, identify and switch bulge sites. Default = True.
            @param switch_constrictions: boolean, identify and switch constriction sites. Default = True.
        '''
        ref_length = 0
        conserved_count = 0
        non_cons_count = 0
        switched_count = 0
        ref_bulge_list, temp_bulge_list, ref_const_list, temp_const_list = [],[],[],[]
        parse = GPCRDBParsingPDB()
        main_pdb_array = parse.pdb_array_creator(structure=self.main_structure)

        # bulges and constrictions
        if switch_bulges==True or switch_constrictions==True:
            self.similarity_table_all = self.run_alignment(core_alignment=False, 
                                                           query_states=['Inactive','Active']).similarity_table
            for ref_res,temp_res,aligned_res in zip(ref_temp_alignment.reference_dict,ref_temp_alignment.template_dict, 
                                                      ref_temp_alignment.alignment_dict):
                gn = ref_res
                gn_TM = parse.gn_num_extract(gn, 'x')[0]
                gn_num = parse.gn_num_extract(gn, 'x')[1]
                
                if ref_temp_alignment.alignment_dict[aligned_res]=='-':
                    if (ref_temp_alignment.reference_dict[ref_res]=='-' and 
                        ref_temp_alignment.reference_dict[parse.gn_indecer(gn,'x',-1)] not in 
                        ['-','/'] and ref_temp_alignment.reference_dict[parse.gn_indecer(gn,'x',+1)] not in ['-','/']): 
    
                        # bulge in template
                        if len(str(gn_num))==3:
                            if switch_bulges==True:
                                try:
                                    Bulge = Bulges(gn)
                                    bulge_template = Bulge.find_bulge_template(self.similarity_table_all, 
                                                                               bulge_in_reference=False)
                                    bulge_site = parse.fetch_residues_from_pdb(self.main_structure,
                                                                               [parse.gn_indecer(gn,'x',-2),
                                                                                parse.gn_indecer(gn,'x',-1),gn,
                                                                                parse.gn_indecer(gn,'x',+1),
                                                                                parse.gn_indecer(gn,'x',+2)])
                                    superpose = sp.BulgeConstrictionSuperpose(bulge_site, bulge_template)
                                    new_residues = superpose.run()
                                    switch_res = 0
                                    for gen_num, atoms in bulge_template.items():
                                        if switch_res!=0 and switch_res!=3:
                                            gn__ = gen_num.replace('.','x')
                                            main_pdb_array[gen_num] = new_residues[gen_num]
                                            ref_temp_alignment.template_dict[gn__] = PDB.Polypeptide.three_to_one(
                                                                                      atoms[0].get_parent().get_resname())
                                            if ref_temp_alignment.template_dict[gn__]==ref_temp_alignment.reference_dict[gn__]:
                                                ref_temp_alignment.alignment_dict[gn__]=ref_temp_alignment.template_dict[gn__]
                                            else:
                                                ref_temp_alignment.alignment_dict[gn__]='.'
                                        switch_res+=1
                                    del main_pdb_array[gn.replace('x','.')]
                                    del ref_temp_alignment.reference_dict[gn]
                                    del ref_temp_alignment.template_dict[gn]
                                    del ref_temp_alignment.alignment_dict[gn]
                                    temp_bulge_list.append({gn:Bulge.template})
                                except:
                                    temp_bulge_list.append({gn:None})
                                
                        # constriction in reference
                        else:
                            if switch_constrictions==True:
                                try:
                                    Const = Constrictions(gn)
                                    constriction_template = Const.find_constriction_template(self.similarity_table_all,
                                                                                        constriction_in_reference=True)
                                    constriction_site = parse.fetch_residues_from_pdb(self.main_structure,
                                                                                      [parse.gn_indecer(gn,'x',-2),
                                                                                       parse.gn_indecer(gn,'x',-1),gn,
                                                                                       parse.gn_indecer(gn,'x',+1),
                                                                                       parse.gn_indecer(gn,'x',+2)])
                                    superpose = sp.BulgeConstrictionSuperpose(constriction_site, constriction_template)
                                    new_residues = superpose.run()                                  
                                    switch_res = 0
                                    for gen_num, atoms in constriction_template.items():
                                        if switch_res!=0 and switch_res!=3:
                                            gn__ = gen_num.replace('.','x')
                                            main_pdb_array[gen_num] = new_residues[gen_num]
                                            ref_temp_alignment.template_dict[gn__] = PDB.Polypeptide.three_to_one(
                                                                                      atoms[0].get_parent().get_resname())
                                            if ref_temp_alignment.template_dict[gn__]==ref_temp_alignment.reference_dict[gn__]:
                                                ref_temp_alignment.alignment_dict[gn__]=ref_temp_alignment.template_dict[gn__]
                                        switch_res+=1
                                    ref_const_list.append({parse.gn_indecer(gn, 'x', -1)+'-'+parse.gn_indecer(gn, 
                                                                                                'x', +1):Const.template})
                                    del main_pdb_array[gn.replace('x','.')]
                                    del ref_temp_alignment.reference_dict[gn]
                                    del ref_temp_alignment.template_dict[gn]
                                    del ref_temp_alignment.alignment_dict[gn]
                                except:
                                    ref_const_list.append({parse.gn_indecer(gn, 'x', -1)+'-'+parse.gn_indecer(gn, 
                                                                                                'x', +1):None})
                    elif (ref_temp_alignment.template_dict[temp_res]=='-' and 
                          ref_temp_alignment.template_dict[parse.gn_indecer(gn,'x',-1)] not in 
                          ['-','/'] and ref_temp_alignment.template_dict[parse.gn_indecer(gn,'x',+1)] not in ['-','/']): 
                        
                        # bulge in reference
                        if len(str(gn_num))==3:
                            if switch_bulges==True:
                                try:
                                    Bulge = Bulges(gn)
                                    bulge_template = Bulge.find_bulge_template(self.similarity_table_all,
                                                                               bulge_in_reference=True)
                                    bulge_site = parse.fetch_residues_from_pdb(self.main_structure,
                                                                               [parse.gn_indecer(gn,'x',-2),
                                                                                parse.gn_indecer(gn,'x',-1),
                                                                                parse.gn_indecer(gn,'x',+1),
                                                                                parse.gn_indecer(gn,'x',+2)])
                                    superpose = sp.BulgeConstrictionSuperpose(bulge_site, bulge_template)
                                    new_residues = superpose.run()
                                    switch_res = 0
                                    for gen_num, atoms in bulge_template.items():
                                        if switch_res!=0 and switch_res!=4:
                                            gn__ = gen_num.replace('.','x')
                                            main_pdb_array[gen_num] = new_residues[gen_num]
                                            ref_temp_alignment.template_dict[gn__] = PDB.Polypeptide.three_to_one(
                                                                                      atoms[0].get_parent().get_resname())
                                            if ref_temp_alignment.template_dict[gn__]==ref_temp_alignment.reference_dict[gn__]:
                                                ref_temp_alignment.alignment_dict[gn__]=ref_temp_alignment.template_dict[gn__]
                                        switch_res+=1
                                    ref_bulge_list.append({gn:Bulge.template})
                                except:
                                    ref_bulge_list.append({gn:None})
                                
                        # constriction in template
                        else:
                            if switch_constrictions==True:
                                try:
                                    Const = Constrictions(gn)
                                    constriction_template = Const.find_constriction_template(self.similarity_table_all,
                                                                                        constriction_in_reference=False)
                                    constriction_site = parse.fetch_residues_from_pdb(self.main_structure,
                                                                                      [parse.gn_indecer(gn,'x',-2),
                                                                                       parse.gn_indecer(gn,'x',-1),
                                                                                       parse.gn_indecer(gn,'x',+1),
                                                                                       parse.gn_indecer(gn,'x',+2)])
                                    superpose = sp.BulgeConstrictionSuperpose(constriction_site, constriction_template)
                                    new_residues = superpose.run()
                                    switch_res = 0
                                    for gen_num, atoms in constriction_template.items():
                                        if switch_res!=0 and switch_res!=4:
                                            gn__ = gen_num.replace('.','x')
                                            main_pdb_array[gen_num] = new_residues[gen_num]
                                            ref_temp_alignment.template_dict[gn__] = PDB.Polypeptide.three_to_one(
                                                                                      atoms[0].get_parent().get_resname())
                                            if ref_temp_alignment.template_dict[gn__]==ref_temp_alignment.reference_dict[gn__]:
                                                ref_temp_alignment.alignment_dict[gn__]=ref_temp_alignment.template_dict[gn__]
                                        switch_res+=1
                                    temp_const_list.append({parse.gn_indecer(gn, 'x', -1)+'-'+parse.gn_indecer(gn, 
                                                                                                'x', +1):Const.template})
                                except:
                                    temp_const_list.append({parse.gn_indecer(gn, 'x', -1)+'-'+parse.gn_indecer(gn, 
                                                                                                'x', +1):None})
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
                for key, value in main_pdb_array.items():
                    out_pdb_array[key] = value                
                    if str(key)+'1' in bulge_gns:
                        out_pdb_array[str(key)+'1'] = main_pdb_array[str(key)+'1']
                main_pdb_array = out_pdb_array
            
            if temp_const_list!=[]:
                out_pdb_array = OrderedDict()
                const_gns = []
                for const in temp_const_list:
                    if list(const.values())[0]!=None:
                        gn = parse.gn_indecer(list(const.keys())[0].split('-')[0].replace('x','.'), '.', +1)
                        const_gns.append(gn)
                for key, value in main_pdb_array.items():
                    out_pdb_array[key] = value
                    if parse.gn_indecer(key, '.', +1) in const_gns:
                        out_pdb_array[parse.gn_indecer(key, '.', +1)] = main_pdb_array[parse.gn_indecer(key, '.', +1)]
                main_pdb_array = out_pdb_array
        
        # check for inconsitencies with db
        pdb_db_inconsistencies = []
        for gn, res in ref_temp_alignment.template_dict.items():
            try:
                pdb_res = PDB.Polypeptide.three_to_one(main_pdb_array[gn.replace('x','.')][0].get_parent().get_resname())
                if res!=pdb_res:
                    ref_temp_alignment.template_dict[gn] = pdb_res
            except:
                if res!='/':
                    pdb_db_inconsistencies.append({gn:ref_temp_alignment.reference_dict[gn]})
                    ref_temp_alignment.template_dict[gn] = 'x'
        
        self.statistics.add_info('pdb_db_inconsistencies', pdb_db_inconsistencies)
        path = "./structure/homology_models/{}_{}/".format(self.uniprot_id,self.state)
        if not os.path.exists(path):
            os.mkdir(path)
        self.write_homology_model_pdb(
                                "./structure/homology_models/{}_{}/pre_switch.pdb".format(self.uniprot_id, self.state), 
                                main_pdb_array, ref_temp_alignment)        
        
        # non-conserved residues
        non_cons_res_templates = OrderedDict()
        conserved_residues = OrderedDict()
        for ref_res, temp_res, aligned_res in zip(ref_temp_alignment.reference_dict, ref_temp_alignment.template_dict, 
                                                  ref_temp_alignment.alignment_dict):
            if ref_temp_alignment.reference_dict[ref_res]!='-' and ref_temp_alignment.reference_dict[ref_res]!='/':
                ref_length+=1
            if (ref_temp_alignment.alignment_dict[aligned_res]!='.' and 
                ref_temp_alignment.alignment_dict[aligned_res]!='/' and 
                ref_temp_alignment.alignment_dict[aligned_res]!='x' and 
                ref_temp_alignment.alignment_dict[aligned_res]!='-'):
                conserved_count+=1
                conserved_residues[ref_res] = ref_temp_alignment.alignment_dict[aligned_res]
            
            gn = ref_res

            if (ref_temp_alignment.alignment_dict[aligned_res]=='.' and 
                ref_temp_alignment.reference_dict[gn]!=ref_temp_alignment.template_dict[gn]):
                non_cons_count+=1
                residues = Residue.objects.filter(generic_number__label=ref_res)
                proteins_w_this_gn = [res.protein_conformation.protein.parent for res in 
                                        residues if str(res.amino_acid)==ref_temp_alignment.reference_dict[ref_res]]
                proteins_w_this_gn = list(set(proteins_w_this_gn))
                gn_ = ref_res.replace('x','.')
                for struct in self.similarity_table:
                    if struct.protein_conformation.protein.parent in proteins_w_this_gn:
                        try:
                            alt_temp = parse.fetch_residues_from_pdb(struct, [gn])
                            if ref_temp_alignment.reference_dict[gn]==PDB.Polypeptide.three_to_one(
                                                                        alt_temp[gn_][0].get_parent().get_resname()):
                                orig_res = main_pdb_array[gn_]
                                alt_res = parse.fetch_residues_from_pdb(struct,[gn])[gn_]
                                superpose = sp.RotamerSuperpose(orig_res, alt_res)
                                new_atoms = superpose.run()
                                if superpose.backbone_rmsd>0.3:
                                    continue
                                main_pdb_array[gn_] = new_atoms
                                ref_temp_alignment.template_dict[gn] = ref_temp_alignment.reference_dict[gn]
                                switched_count+=1                     
                                non_cons_res_templates[ref_res] = struct
                                break
                        except:
                            pass
                    else:
                        try:
                            residue = main_pdb_array[gn_]
                            main_pdb_array[gn_] = residue[0:4]
                        except:
                            logging.warning("Missing atoms in {} at {}".format(self.main_structure,gn))

        self.statistics.add_info('ref_seq_length', ref_length)
        self.statistics.add_info('conserved_num', conserved_count)
        self.statistics.add_info('non_conserved_num', non_cons_count)
        self.statistics.add_info('non_conserved_switched_num', switched_count)
        self.statistics.add_info('conserved_residues', conserved_residues)
        self.statistics.add_info('non_conserved_residue_templates', non_cons_res_templates)
        self.statistics.add_info('similarity_table', self.similarity_table)
        
        # write to file
        path = "./structure/homology_models/{}_{}/".format(self.uniprot_id,self.state)
        if not os.path.exists(path):
            os.mkdir(path)
        self.write_homology_model_pdb(path+self.uniprot_id+"_post.pdb", main_pdb_array, ref_temp_alignment)
        return ref_temp_alignment       
    
    def write_homology_model_pdb(self, filename, main_pdb_array, ref_temp_alignment):
        ''' Write PDB file from pdb array to file.
        
            @param filename: str, filename of output file \n
            @param main_pdb_array: OrderedDict(), of atoms of pdb, where keys are generic numbers/residue numbers and
            values are list of atoms. Output of GPCRDBParsingPDB.pdb_array_creator().
            @param ref_temp_alignment: AlignedReferenceAndTemplate, only writes residues that are in ref_temp_alignment.
        '''
        res_num=0
        atom_num=0
        with open(filename,'w+') as f:
            for key in main_pdb_array:
                if '.' in str(key) and str(key).replace('.','x') in ref_temp_alignment.reference_dict:
                    res_num+=1
                    segment = int(str(key).split('.')[0])
                    try:
                        segment_pre+=0
                    except:
                        segment_pre = 100
                    if segment>segment_pre:
                        f.write("\nTER")
                    for atom in main_pdb_array[key]:              
                        atom_num+=1
                        coord = list(atom.get_coord())
                        coord1 = "%8.3f"% (coord[0])
                        coord2 = "%8.3f"% (coord[1])
                        coord3 = "%8.3f"% (coord[2])
                        if str(atom.get_id())=='CA':
                            if len(key)==4:
                                bfact = "%6.2f"% (float(key))
                            else:
                                bfact = "-%5.2f"% (float(key))
                        else:
                            bfact = "%6.2f"% (float(atom.get_bfactor()))
                        occupancy = "%6.2f"% (atom.get_occupancy())
                        template="""
ATOM{atom_num}  {atom}{res} {chain}{res_num}{coord1}{coord2}{coord3}{occupancy}{bfactor}{atom_s}  """
                        context={"atom_num":str(atom_num).rjust(7), "atom":str(atom.get_id()).ljust(4),
                                 "res":atom.get_parent().get_resname(), "chain":str(self.main_template_preferred_chain)[0],
                                 "res_num":str(res_num).rjust(4), "coord1":coord1.rjust(12), 
                                 "coord2":coord2.rjust(8), "coord3":coord3.rjust(8), 
                                 "occupancy":str(occupancy).rjust(3),
                                 "bfactor":str(bfact).rjust(4), "atom_s":str(str(atom.get_id())[0]).rjust(12)}
                        f.write(template.format(**context))
                    segment_pre=segment
            f.write("\nTER\nEND")
    
    def create_PIR_file(self, ref_temp_alignment, template_file):
        ''' Create PIR file from reference and template alignment (AlignedReferenceAndTemplate).
        
            @param ref_temp_alignment: AlignedReferenceAndTemplate
            @template_file: str, name of template file with path
        '''
        ref_sequence, temp_sequence = '',''
        for ref_res, temp_res in zip(ref_temp_alignment.reference_dict, ref_temp_alignment.template_dict):
            if ref_temp_alignment.reference_dict[ref_res]=='x':
                ref_sequence+='-'
            else:
                ref_sequence+=ref_temp_alignment.reference_dict[ref_res]
            if ref_temp_alignment.template_dict[temp_res]=='x':
                temp_sequence+='-'
            else:
                temp_sequence+=ref_temp_alignment.template_dict[temp_res]
        with open("./structure/PIR/"+self.uniprot_id+"_"+self.state+".pir", 'w+') as output_file:
            template="""
>P1;{temp_file}
structure:{temp_file}:1:{chain}:210:{chain}::::
{temp_sequence}*

>P1;{uniprot}
sequence:{uniprot}::::::::
{ref_sequence}*
            """
            context={"temp_file":template_file,
                     "chain":self.main_template_preferred_chain,
                     "temp_sequence":temp_sequence,
                     "uniprot":self.uniprot_id,
                     "ref_sequence":ref_sequence}
            output_file.write(template.format(**context))
            
    def run_MODELLER(self, pir_file, template, reference, number_of_models, output_file_name):
        ''' Build homology model with MODELLER.
        
            @param pir_file: str, file name of PIR file with path \n
            @param template: str, file name of template with path \n
            @param reference: str, Uniprot code of reference sequence \n
            @param number_of_models: int, number of models to be built \n
            @param output_file_name: str, name of output file
        '''
        log.verbose()
        env = environ(rand_seed=80851) #!!random number generator
        
        a = automodel(env, alnfile = pir_file, knowns = template, sequence = reference)
        a.starting_model = 1
        a.ending_model = number_of_models
        a.md_level = refine.very_slow
        path = "./structure/homology_models/{}".format(reference+"_"+self.state)
        if not os.path.exists(path):
            os.mkdir(path)
        a.make()
        for file in os.listdir("./"):
            if file.startswith(self.uniprot_id) and file.endswith(".pdb"):
                os.rename("./"+file, "./structure/homology_models/{}_{}/".format(self.uniprot_id,self.state)+output_file_name)
            elif file.startswith(self.uniprot_id):
                os.rename("./"+file, "./structure/homology_models/{}_{}/".format(self.uniprot_id,self.state)+file)
            


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
        if bulge_in_reference==True:
            matches = Residue.objects.filter(generic_number__label=gn)
        elif bulge_in_reference==False:
            excludees = Residue.objects.filter(generic_number__label=gn)
            excludee_proteins = list(OrderedDict.fromkeys([res.protein_conformation.protein.parent.entry_name 
                                        for res in excludees if res.protein_conformation.protein.parent!=None]))
            matches = Residue.objects.filter(generic_number__label=gn[:-1])
        for structure, value in similarity_table.items():  
            protein_object = Protein.objects.get(id=structure.protein_conformation.protein.parent.id)
            try:                            
                for match in matches:
                    if bulge_in_reference==True:
                        if match.protein_conformation.protein.parent==protein_object:
                            self.bulge_templates.append(structure)
                    elif bulge_in_reference==False:
                        if (match.protein_conformation.protein.parent==protein_object and 
                            match.protein_conformation.protein.parent.entry_name not in excludee_proteins):
                            self.bulge_templates.append(structure)
            except:
                pass
        mod_bulge = False
        for temp in self.bulge_templates:
            try:
                if bulge_in_reference==True:
                    alt_bulge = parse.fetch_residues_from_pdb(temp, 
                                                              [parse.gn_indecer(gn,'x',-2),
                                                               parse.gn_indecer(gn,'x',-1),gn,
                                                               parse.gn_indecer(gn,'x',+1),
                                                               parse.gn_indecer(gn,'x',+2)])
                elif bulge_in_reference==False:
                    gn_list = [parse.gn_indecer(gn,'x',-2), parse.gn_indecer(gn,'x',-1),
                               parse.gn_indecer(gn,'x',+1), parse.gn_indecer(gn,'x',+2)]
                    for gn_ in gn_list:
                        if len(Residue.objects.filter(generic_number__label=gn_+'1').filter(
                                                                    protein_conformation=temp.protein_conformation))>0:
                            if int(gn[:-1].split('x')[1])-int(gn_.split('x')[1])==1:
                                gn_list[0] = gn_+'1'
                            elif int(gn[:-1].split('x')[1])-int(gn_.split('x')[1])==-1:
                                gn_list[2] = gn_+'1'
                                gn_list[3] = parse.gn_indecer(gn_,'x',+1)
                            elif int(gn[:-1].split('x')[1])-int(gn_.split('x')[1])==-2:
                                gn_list[3] = gn_+'1'
                            mod_bulge = True
                    if mod_bulge==True:
                        alt_bulge = parse.fetch_residues_from_pdb(temp, gn_list, modify_bulges=True)
                    else:
                        alt_bulge = parse.fetch_residues_from_pdb(temp, gn_list)
                self.template = temp              
                break
            except:
                self.template = None               
        return alt_bulge
            
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
        if constriction_in_reference==True:
            excludees = Residue.objects.filter(generic_number__label=gn)
            excludee_proteins = list(OrderedDict.fromkeys([res.protein_conformation.protein.parent.entry_name 
                                        for res in excludees if res.protein_conformation.protein.parent!=None]))
            matches = Residue.objects.filter(generic_number__label=parse.gn_indecer(gn,'x',-1))
        elif constriction_in_reference==False:
            matches = Residue.objects.filter(generic_number__label=gn)
        for structure, value in similarity_table.items():  
            protein_object = Protein.objects.get(id=structure.protein_conformation.protein.parent.id)
            try:                            
                for match in matches:
                    if constriction_in_reference==True:                        
                        if (match.protein_conformation.protein.parent==protein_object and 
                            match.protein_conformation.protein.parent.entry_name not in excludee_proteins):
                            self.constriction_templates.append(structure)
                    elif constriction_in_reference==False:
                        if match.protein_conformation.protein.parent==protein_object:
                            self.constriction_templates.append(structure)
            except:
                pass
        for temp in self.constriction_templates:
            try:
                if constriction_in_reference==True:
                    alt_bulge = parse.fetch_residues_from_pdb(temp,
                                                              [parse.gn_indecer(gn,'x',-2),
                                                               parse.gn_indecer(gn,'x',-1),
                                                               parse.gn_indecer(gn,'x',+1),
                                                               parse.gn_indecer(gn,'x',+2)])
                elif constriction_in_reference==False:
                    alt_bulge = parse.fetch_residues_from_pdb(temp,
                                                              [parse.gn_indecer(gn,'x',-2),
                                                               parse.gn_indecer(gn,'x',-1),gn,
                                                               parse.gn_indecer(gn,'x',+1),
                                                               parse.gn_indecer(gn,'x',+2)])
                self.template = temp              
                break
            except:
                self.template = None               
        return alt_bulge
        
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
            rotamer = Rotamer.objects.get(structure__protein_conformation=structure.protein_conformation, 
                                                                          residue__generic_number__label=gn)
            io = StringIO(rotamer.pdbdata.pdb)
            rota_struct = PDB.PDBParser().get_structure('structure', io)[0]
            for chain in rota_struct:
                for residue in chain:
                    for atom in residue:
                        atoms_list.append(atom)
                    if modify_bulges==True and len(gn)==5:
                        output[gn.replace('x','.')[:-1]] = atoms_list
                    else:
                        output[gn.replace('x','.')] = atoms_list
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
        
        for chain in pdb_struct:
            for residue in chain:
                try:
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

class Validation():
    ''' Class to validate homology models.
    '''
    def __init__(self):
        pass
    
    def PDB_RMSD(self, pdb_file1, pdb_file2, assign_gns=[1,2]):
        ''' Calculates root-mean-square deviation between the coordinates of two model PDB files. The two files
            must have the same number of atoms.
            
            @param pdb_file1: str, file name of first file with path \n
            @param pdb_file2: str, file name of second file with path
        '''
        array1, array2 = np.array([0,0,0]), np.array([0,0,0])
        parser = PDB.PDBParser()
                
        pdb1 = parser.get_structure('struct1', pdb_file1)[0]
        pdb2 = parser.get_structure('struct2', pdb_file2)[0]
        
        if assign_gns==[1]:
            assign_gn1 = as_gn.GenericNumbering(structure=pdb1)
            pdb1 = assign_gn1.assign_generic_numbers()
        elif assign_gns==[2]:
            assign_gn2 = as_gn.GenericNumbering(structure=pdb2)
            pdb2 = assign_gn2.assign_generic_numbers()
        elif assign_gns==[1,2]:
            assign_gn1 = as_gn.GenericNumbering(structure=pdb1)
            pdb1 = assign_gn1.assign_generic_numbers()
            assign_gn2 = as_gn.GenericNumbering(structure=pdb2)
            pdb2 = assign_gn2.assign_generic_numbers()    
        
        pdb_array1, pdb_array2 = OrderedDict(), OrderedDict()
        
        for chain1, chain2 in zip(pdb1, pdb2):
            for residue1, residue2 in zip(chain1, chain2):
                if -8.1 < residue1['CA'].get_bfactor() < 8.1:
                    pdb_array1[residue1['CA'].get_bfactor()] = residue1
                if -8.1 < residue2['CA'].get_bfactor() < 8.1:
                    pdb_array2[residue2['CA'].get_bfactor()] = residue2

        orig_atomlist, temp_atomlist = [], []               
        for gn1, res1 in pdb_array1.items():
            for gn2, res2 in pdb_array2.items():
                if gn1==gn2 and res1.get_resname()==res2.get_resname():
                    for atom1 in res1:
                        for atom2 in res2:
                            if atom1.get_id()==atom2.get_id():
                                orig_atomlist.append(atom1)
                                temp_atomlist.append(atom2)
                elif gn1==gn2:
                    for atom1 in res1:
                        for atom2 in res2:
                            if atom1.get_id()==atom2.get_id() and atom1.get_id() in ['N','CA','C','O']:
                                orig_atomlist.append(atom1)
                                temp_atomlist.append(atom2)
        superpose = sp.RotamerSuperpose(orig_atomlist, temp_atomlist)
        temp_atomlist = superpose.run()
        for a1, a2 in zip(orig_atomlist, temp_atomlist):
            array1 = np.vstack((array1, list(a1.get_coord())))
            array2 = np.vstack((array2, list(a2.get_coord())))
        rmsd = np.sqrt(((array1[1:]-array2[1:])**2).mean())
        return rmsd
