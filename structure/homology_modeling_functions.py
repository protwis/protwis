from django.conf import settings

from protein.models import Protein, ProteinConformation, ProteinAnomaly, ProteinState, ProteinSegment
from residue.models import Residue
from residue.functions import dgn, ggn
from structure.models import *
from structure.functions import HSExposureCB, PdbStateIdentifier, StructureSeqNumOverwrite, update_template_source, compare_and_update_template_source
from common.alignment import AlignedReferenceTemplate, GProteinAlignment
from common.definitions import *
from common.models import WebLink
from signprot.models import SignprotComplex
import structure.structural_superposition as sp
import structure.assign_generic_numbers_gpcr as as_gn

import Bio.PDB as PDB
from collections import OrderedDict
import os
import logging
import pprint
from io import StringIO, BytesIO
import sys
import re
import math
import yaml
import json
import traceback
import subprocess
from copy import deepcopy
import zipfile
import shutil
from modeller import *
from modeller.automodel import *


class SignprotFunctions(object):
    def __init__(self):
        pass

    def get_subtypes_with_templates(self):
        return SignprotComplex.objects.all().values_list('protein__entry_name', flat=True)

    def get_subfamilies_with_templates(self, receptor_family):
        subfams = []
        for i in SignprotComplex.objects.filter(structure__protein_conformation__protein__family__parent__parent__parent__name=receptor_family):
            if i.protein.family.parent.name not in subfams:
                subfams.append(i.protein.family.parent.name)
        return subfams

    def get_receptor_families_with_templates(self):
        return SignprotComplex.objects.all().values_list('structure__protein_conformation__protein__family__parent__parent__parent__name', flat=True).distinct()

    def get_subfam_subtype_dict(self, subfamilies, receptor_family):
        d = {}
        for s in subfamilies:
            ordered_prots, non_ordered_prots = [], []
            prots = [i.entry_name for i in Protein.objects.filter(family__parent__name=s, species__common_name='Human', accession__isnull=False)]
            for p in prots:
                if p=='gnal_human':
                    continue
                if len(SignprotComplex.objects.filter(protein__entry_name=p, structure__protein_conformation__protein__family__parent__parent__parent__name=receptor_family))>0:
                    ordered_prots.append(p)
                else:
                    non_ordered_prots.append(p)
            d[s] = ordered_prots+non_ordered_prots
        return d

    def get_other_subtypes_in_subfam(self, protein):
        this_prot = Protein.objects.get(entry_name=protein)
        return [i.entry_name for i in Protein.objects.filter(family=this_prot.family.parent, species__common_name='Human', accession__isnull=False).exclude(entry_name=protein)]


class GPCRDBParsingPDB(object):
    ''' Class to manipulate cleaned pdb files of GPCRs.
    '''
    def __init__(self):
        self.segment_coding = OrderedDict([(1,'TM1'),(2,'TM2'),(3,'TM3'),(4,'TM4'),(5,'TM5'),(6,'TM6'),(7,'TM7'),(8,'H8')])

    @staticmethod
    def parse_rotamer_pdb(rotamer):
        atoms_list = []
        io = StringIO(rotamer.pdbdata.pdb)
        rota_struct = PDB.PDBParser(QUIET=True).get_structure('structure', io)[0]
        for chain in rota_struct:
            for residue in chain:
                for atom in residue:
                    atoms_list.append(atom)
        return atoms_list

    
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
        res1 = Residue.objects.get(protein_conformation=protein_conformation, display_generic_number__label=dgn(gn1,protein_conformation))
        res2 = Residue.objects.get(protein_conformation=protein_conformation, display_generic_number__label=dgn(gn2,protein_conformation))
        return res1.sequence_number-res2.sequence_number
            
    def gn_indecer(self, gn, delimiter, direction):
        ''' Get an upstream or downstream generic number from reference generic number.
        
            @param gn: str, Generic number \n
            @param delimiter: str, character between TM and position (usually 'x') \n 
            @param direction: int, n'th position from gn (+ or -)
        '''
        split = self.gn_num_extract(gn, delimiter)
        if len(str(split[1]))==2:
            return str(split[0])+delimiter+str(split[1]+direction)
        elif len(str(split[1]))==3:
            if direction<0:
                direction += 1
            return str(split[0])+delimiter+str(int(str(split[1])[:2])+direction)

    @staticmethod
    def find_before_4_gns(section, reference_dict, distorted_residues, segment_labels):
        minus_offset, minus_gns = [], []
        first_index = list(reference_dict[section[0][0]]).index(section[0][1])
        first_seg = section[0][0]
        ### If not 4 residues in this segment, start from prev segment
        if first_index<4:
            first_seg = segment_labels[segment_labels.index(section[0][0])-1]
            for fi in range(first_index-1, -1, -1):
                minus_offset.append([section[0][0], list(reference_dict[section[0][0]])[fi]])
            if len(reference_dict[first_seg])<4:
                for fi in range(len(reference_dict[first_seg])-1, -1, -1):
                    minus_offset.append([first_seg, list(reference_dict[first_seg])[fi]])
                first_seg = segment_labels[segment_labels.index(first_seg)-1]
            first_index = len(reference_dict[first_seg])
        first_index_orig = first_index
        minus_gns = []

        ### Find before 4
        for j in range(1, 5):
            while list(reference_dict[first_seg])[first_index-j] in distorted_residues[first_seg]:
                minus_offset.append([first_seg,list(reference_dict[first_seg])[first_index-j]])
                first_index-=1
                if first_index_orig>0 and first_index-j<0:
                    first_seg = segment_labels[segment_labels.index(first_seg)-1]
                    first_index = len(reference_dict[first_seg])
                    j = 1
            minus_gns.append(list(reference_dict[first_seg])[first_index-j])

        minus_gns.reverse()
        if len(minus_offset)>0:
            for m_s, m_gn in minus_offset:
                section = [[m_s, m_gn]]+section

        return section, minus_gns, first_seg

    @staticmethod
    def find_after_4_gns(section, reference_dict, distorted_residues, segment_labels):
        ''' Find after 4 residues for superpositioning'''
        plus_offset, plus_gns = [], []
        last_index = list(reference_dict[section[-1][0]]).index(section[-1][1])
        last_seg = section[-1][0]
        last_seg_orig = last_seg

        ### If not 4 residues in this segment, end in next segment
        if len(reference_dict[section[-1][0]])-last_index<5:
            last_seg = segment_labels[segment_labels.index(section[-1][0])+1]
            for li in range(last_index+1, len(reference_dict[section[-1][0]])):
                plus_offset.append([section[-1][0], list(reference_dict[section[-1][0]])[li]])
            if len(reference_dict[last_seg])<4:
                for li in range(0, len(reference_dict[last_seg])+1, 1):
                    plus_offset.append([last_seg, list(reference_dict[last_seg])[li]])
                last_seg = segment_labels[segment_labels.index(last_seg)+1]
            last_index = 0
        if last_seg!=last_seg_orig:
            r = range(0,4)
        else:
            r = range(1,5)
        for k in r:
            while list(reference_dict[last_seg])[last_index+k] in distorted_residues[last_seg]:
                plus_offset.append([last_seg,list(reference_dict[last_seg])[last_index+k]])
                last_index+=1
            plus_gns.append(list(reference_dict[last_seg])[last_index+k])
        
        if len(plus_offset)>0:
            for p_s, p_gn in plus_offset:
                section.append([p_s, p_gn])

        return section, plus_gns, last_seg

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
                        residue__display_generic_number__label=dgn(gn,structure.protein_conformation), 
                        structure__preferred_chain=structure.preferred_chain))
            else:
                residue = Residue.objects.get(protein_conformation=structure.protein_conformation, sequence_number=gn)
                rotamer = list(Rotamer.objects.filter(structure__protein_conformation=structure.protein_conformation, 
                        residue=residue, structure__preferred_chain=structure.preferred_chain))
                if just_nums==False:
                    try:
                        gn = ggn(Residue.objects.get(protein_conformation=structure.protein_conformation,
                                                    sequence_number=gn).display_generic_number.label)
                    except:
                        pass
            if len(rotamer)>1:
                for i in rotamer:
                    if i.pdbdata.pdb.startswith('COMPND')==False:
                        if i.pdbdata.pdb[21] in structure.preferred_chain:
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
        output = OrderedDict()
        for i,j in dict1.items():
            output[i] = j
        for i,j in dict2.items():
            output[i] = j
        return output

    def right_rotamer_select(self, rotamer):
        ''' Filter out compound rotamers.
        '''
        if len(rotamer)>1:
            for i in rotamer:
                if i.pdbdata.pdb.startswith('COMPND')==False:
                    rotamer = i
                    break
        else:
            rotamer=rotamer[0]
        return rotamer

    def pdb_array_creator(self, structure=None, filename=None):
        ''' Creates an OrderedDict() from the pdb of a Structure object where residue numbers/generic numbers are 
            keys for the residues, and atom names are keys for the Bio.PDB.Residue objects.
            
            @param structure: Structure, Structure object of protein. When using structure, leave filename=None. \n
            @param filename: str, filename of pdb to be parsed. When using filename, leave structure=None).
        '''
        # seq_nums_overwrite_cutoff_dict = {'4PHU':2000, '4LDL':1000, '4LDO':1000, '4QKX':1000, '5JQH':1000, '5TZY':2000, '5KW2':2000}
        if structure!=None and filename==None:
            io = StringIO(structure.pdb_data.pdb)
        else:
            io = filename
        gn_array = []
        residue_array = []
        # pdb_struct = PDB.PDBParser(QUIET=True).get_structure(structure.pdb_code.index, io)[0]

        residues = Residue.objects.filter(protein_conformation=structure.protein_conformation)
        gn_list = []
        for i in residues:
            try:
                gn_list.append(ggn(i.display_generic_number.label).replace('x','.'))
            except:
                pass

        ssno = StructureSeqNumOverwrite(structure)
        ssno.seq_num_overwrite('pdb')
        if len(ssno.pdb_wt_table)>0:
            residues = residues.filter(protein_segment__slug__in=['TM1','TM2','TM3','TM4','TM5','TM6','TM7','H8']).order_by('sequence_number')
            output = OrderedDict()
            for r in residues:
                if r.protein_segment.slug==None:
                    continue
                if r.protein_segment.slug not in output:
                    output[r.protein_segment.slug] = OrderedDict()
                rotamer = Rotamer.objects.filter(residue=r)
                rotamer = self.right_rotamer_select(rotamer)
                rota_io = StringIO(rotamer.pdbdata.pdb)
                p = PDB.PDBParser()
                parsed_rota = p.get_structure('rota', rota_io)
                for chain in parsed_rota[0]:
                    for res in chain:
                        atom_list = []
                        for atom in res:
                            # Skip hydrogens
                            if atom.get_id().startswith('H'):
                                continue
                            if atom.get_id()=='N':
                                bw, gn = r.display_generic_number.label.split('x')
                                atom.set_bfactor(bw)
                            elif atom.get_id()=='CA':
                                bw, gn = r.display_generic_number.label.split('x')
                                gn = "{}.{}".format(bw.split('.')[0], gn)
                                if len(gn.split('.')[1])==3:
                                    gn = '-'+gn[:-1]
                                atom.set_bfactor(gn)
                            atom_list.append(atom)
                        output[r.protein_segment.slug][ggn(r.display_generic_number.label).replace('x','.')] = atom_list
            return output
        else:
            assign_gn = as_gn.GenericNumbering(pdb_file=io, pdb_code=structure.pdb_code.index, sequence_parser=True)
            pdb_struct = assign_gn.assign_generic_numbers_with_sequence_parser()
            pref_chain = structure.preferred_chain
            parent_prot_conf = ProteinConformation.objects.get(protein=structure.protein_conformation.protein.parent)
            parent_residues = Residue.objects.filter(protein_conformation=parent_prot_conf)
            last_res = list(parent_residues)[-1].sequence_number
            if len(pref_chain)>1:
                pref_chain = pref_chain[0]
            for residue in pdb_struct[pref_chain]:
                if 'CA' in residue and -9.1 < residue['CA'].get_bfactor() < 9.1:
                    use_resid = False
                    gn = str(residue['CA'].get_bfactor())
                    if len(gn.split('.')[1])==1:
                        gn = gn+'0'
                    if gn[0]=='-':
                        gn = gn[1:]+'1'
                    # Exceptions
                    if structure.pdb_code.index=='3PBL' and residue.get_id()[1]==331:
                        use_resid = True
                    elif structure.pdb_code.index=='6QZH' and residue.get_id()[1]==1434:
                        use_resid = True
                    elif structure.pdb_code.index=='7M3E':
                        use_resid = True
                    #################################################
                    elif gn in gn_list:
                        gn_array.append(gn)
                        residue_array.append(residue.get_list())
                    else:
                        use_resid = True
                    if use_resid:
                        gn_array.append(str(residue.get_id()[1]))
                        residue_array.append(residue.get_list())
            output = OrderedDict()
            for num, label in self.segment_coding.items():
                output[label] = OrderedDict()
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
                        found_res, found_gn = None, None
                        try:
                            found_res = Residue.objects.get(protein_conformation=structure.protein_conformation,
                                                            sequence_number=gn)
                        except:
                            # Exception for res 317 in 5VEX, 5VEW
                            if structure.pdb_code.index in ['5VEX','5VEW'] and gn=='317' and res[0].get_parent().get_resname()=='CYS':
                                found_res = Residue.objects.get(protein_conformation=parent_prot_conf,
                                                                sequence_number=gn)
                            #####################################
                        found_gn = str(ggn(found_res.display_generic_number.label)).replace('x','.')

                        # Exception for res 318 in 5VEX, 5VEW
                        if structure.pdb_code.index in ['5VEX','5VEW'] and gn=='318' and res[0].get_parent().get_resname()=='ILE' and found_gn=='5.47':
                            found_gn = '5.48'
                        #####################################
                        if -9.1 < float(found_gn) < 9.1:
                            if len(res)==1:
                                continue
                            if int(gn)>last_res:
                                continue
                            seg_label = self.segment_coding[int(found_gn.split('.')[0])]
                            output[seg_label][found_gn] = res
                    except:
                        if res[0].get_parent().get_resname()=='YCM' or res[0].get_parent().get_resname()=='CSD':
                            try:
                                found_res = Residue.objects.get(protein_conformation=parent_prot_conf, sequence_number=gn)
                            except:
                                continue
                            if found_res.protein_segment.slug[0] not in ['T','H']:
                                continue
                            try:
                                found_gn = str(ggn(found_res.display_generic_number.label)).replace('x','.')
                            except:
                                found_gn = str(gn)
                            output[found_res.protein_segment.slug][found_gn] = res
        return output

    @staticmethod
    def create_g_alpha_pdb_array(signprot_complex):
        segments = ProteinSegment.objects.filter(proteinfamily='Alpha')
        residues = Residue.objects.filter(protein_conformation__protein__entry_name=signprot_complex.structure.pdb_code.index.lower()+'_a')
        pdb_array = OrderedDict()
        parse = GPCRDBParsingPDB()
        for s in segments:
            if s.slug not in pdb_array:
                pdb_array[s.slug] = OrderedDict()
            for r in residues.filter(protein_segment=s):
                try:
                    rotamers = Rotamer.objects.filter(structure=signprot_complex.structure, residue__display_generic_number__label=r.display_generic_number.label)
                    if len(rotamers)==0:
                        raise Exception()
                    rotamer = parse.right_rotamer_select(rotamers)
                    p = PDB.PDBParser(QUIET=True).get_structure('structure', StringIO(rotamer.pdbdata.pdb))[0]
                    atoms = []
                    for chain in p:
                        for res in chain:
                            for atom in res:
                                atoms.append(atom)
                except:
                    atoms = 'x'
                pdb_array[r.protein_segment.slug][r.display_generic_number.label] = atoms
        return pdb_array


class ImportHomologyModel():
    ''' Imports receptor homology model for complex model building pipeline. The idea is to save time by not rerunning the receptor modeling
        when separate complex models are built with different subunits for the same receptor.
    '''
    def __init__(self, receptor, sign_prot=None, zip_path='./structure/complex_models_zip/'):
        self.receptor = receptor
        self.sign_prot = sign_prot
        self.zip_path = zip_path
        if sign_prot:
            self.segments = [i.slug for i in ProteinSegment.objects.filter(proteinfamily=sign_prot)]
        else:
            self.segments = [i.slug for i in ProteinSegment.objects.filter(proteinfamily='GPCR')]
        self.path_to_pdb = ''

    def find_files(self):
        ''' Locates files of one receptor model that can be used as source.
        '''
        sf = SignprotFunctions()
        other_signprots = sf.get_other_subtypes_in_subfam(self.sign_prot)
        gprots_with_structure = sf.get_subtypes_with_templates()
        if self.sign_prot in gprots_with_structure:
            return None
        if not os.path.exists(self.zip_path):
            os.mkdir(self.zip_path)
        files = os.listdir(self.zip_path)
        for f in files:
            if f.endswith('.zip'):
                modelname = f.split('.')[0]
                if self.receptor not in modelname:
                    continue
                found_other_sf = False
                for o in other_signprots:
                    if o in modelname:
                        found_other_sf = True
                if not found_other_sf:
                    continue
                mod_dir = self.zip_path+modelname
                if not os.path.exists(mod_dir):
                    os.mkdir(mod_dir)
                zip_mod = zipfile.ZipFile(self.zip_path+f, 'r')
                zip_mod.extractall(mod_dir)
                zip_mod.close()

                name_list = modelname.split('_')
                if name_list[3] not in ['Inactive','Active','Intermediate'] and name_list[4]!='refined':
                    self.complex = True
                    self.revise_xtal = False
                    gpcr_class = name_list[0][-1]
                    gpcr_prot = '{}_{}'.format(name_list[1],name_list[2].split('-')[0])
                    sign_prot = '{}_{}'.format(name_list[2].split('-')[1], name_list[3])
                    main_structure = name_list[4]
                    build_date = name_list[5]

                    p = PDB.PDBParser()
                    self.path_to_pdb = os.sep.join([self.zip_path, modelname, modelname+'.pdb'])
                    model = p.get_structure('receptor', os.sep.join([self.zip_path, modelname, modelname+'.pdb']))[0]['R']
                    with open(os.sep.join([self.zip_path, modelname, modelname+'.templates.csv']), 'r') as templates_file:
                        templates = templates_file.readlines()
                    with open(os.sep.join([self.zip_path, modelname, modelname+'.template_similarities.csv']), 'r') as sim_file:
                        similarities = sim_file.readlines()
                    return model, templates, similarities
        return None

    def parse_model(self, modelchain):
        ''' Parses model file.
        '''
        reference_dict, main_pdb_array = OrderedDict(), OrderedDict()
        for i in self.segments:
            reference_dict[i] = OrderedDict()
            main_pdb_array[i] = OrderedDict()
        resis = Residue.objects.filter(protein_conformation__protein__entry_name=self.receptor)
        for res in modelchain:
            seqnum = res.get_id()[1]
            this_res = resis.get(sequence_number=seqnum)
            try:
                if self.sign_prot=='Alpha':
                    reference_dict[this_res.protein_segment.slug][this_res.display_generic_number.label] = this_res.amino_acid
                else:
                    reference_dict[this_res.protein_segment.slug][ggn(this_res.display_generic_number.label)] = this_res.amino_acid
            except AttributeError:
                reference_dict[this_res.protein_segment.slug][str(this_res.sequence_number)] = this_res.amino_acid
            atoms_list = []
            for atom in res:
                if atom.element!='H':
                    atoms_list.append(atom)
            try:
                if self.sign_prot=='Alpha':
                    main_pdb_array[this_res.protein_segment.slug][this_res.display_generic_number.label] = atoms_list
                else:
                    main_pdb_array[this_res.protein_segment.slug][ggn(this_res.display_generic_number.label).replace('x','.')] = atoms_list
            except AttributeError:
                main_pdb_array[this_res.protein_segment.slug][str(this_res.sequence_number)] = atoms_list

        return reference_dict, reference_dict, reference_dict, main_pdb_array

    def parse_template_source(self, source_file):
        ''' Parses template source file.
        '''
        structure_dict = {}
        template_source = OrderedDict()
        for i in self.segments:
            template_source[i] = OrderedDict()
        for line in source_file[1:]:
            split_line = line.split(',')
            if split_line[2]!='None':
                key = split_line[2]
            else:
                key = split_line[1]
            if split_line[3] not in structure_dict:
                if split_line[3]!='None':
                    s1 = Structure.objects.get(pdb_code__index=split_line[3])
                    structure_dict[split_line[3]] = s1
                else:
                    s1 = None
            else:
                s1 = structure_dict[split_line[3]]
            if split_line[4][:-1] not in structure_dict:
                if split_line[4][:-1]!='None':
                    s2 = Structure.objects.get(pdb_code__index=split_line[4][:-1])
                    structure_dict[split_line[4][:-1]] = s2
                else:
                    s2 = None
            else:
                s2 = structure_dict[split_line[3]]
            if split_line[0] in template_source:
                template_source[split_line[0]][key] = [s1, s2]
        return template_source

    def parse_similarities(self, source_file):
        ''' Parses template similarities file.
        '''
        sim_dict = OrderedDict()
        for line in source_file[1:]:
            split_line = line.split(',')
            s = Structure.objects.get(pdb_code__index=split_line[0])
            sim_dict[s] = int(split_line[1])
        return sim_dict

    def find_disulfides(self):
        ''' Finds disulfide bridges.
        '''
        if self.path_to_pdb=='':
            raise AssertionError('Run find_files function first to locate the model source file.')
        disulfide_pairs = []
        with open(self.path_to_pdb, 'r') as f:
            lines = f.readlines()
            c=1
            for line in lines:
                if c>10:
                    break
                if line.startswith('SSBOND'):
                    pdb_re = re.search('SSBOND\s+\d+\s+CYS\sR\s+(\d+)\s+CYS\sR\s+(\d+)', line)
                    num1 = pdb_re.group(1)
                    num2 = pdb_re.group(2)
                    res1, res2 = list(Residue.objects.filter(protein_conformation__protein__entry_name=self.receptor, sequence_number__in=[num1, num2]))
                    if res1.display_generic_number!=None:
                        gn1 = ggn(res1.display_generic_number.label)
                    else:
                        gn1 = str(res1.sequence_number)
                    if res2.display_generic_number!=None:
                        gn2 = ggn(res2.display_generic_number.label)
                    else:
                        gn2 = str(res2.sequence_number)
                    disulfide_pairs.append([gn1, gn2])
                c+=1
        print('FIND disulfides: {}'.format(disulfide_pairs))
        return disulfide_pairs


class DummyStructure():
    def __init__(self, preferred_chain):
        self.preferred_chain = preferred_chain
        self.pdb_code = DummyPDBCode()
        self.protein_conformation = ProteinConformation(protein=Protein(entry_name='AF', parent=Protein(entry_name='AF')))

    def __str__(self):
        return self.pdb_code.index


class DummyPDBCode():
    def __init__(self):
        self.index = "AF"


class Remodeling():
    def __init__(self, model_file, gaps={}, receptor=None, signprot=None, icl3_delete=[]):
        self.model_file = model_file
        self.struct = PDB.PDBParser(QUIET=True).get_structure('structure', self.model_file)[0]
        self.pir_file = './structure/PIR/{}'.format(self.model_file.split('/')[-1].replace('.pdb','.pir'))
        self.gaps = gaps
        self.remark_lines = []
        self.receptor = receptor
        self.signprot = signprot
        self.icl3_delete = icl3_delete

    def find_clashes(self):
        hse = HSExposureCB(self.struct, radius=11, check_chain_breaks=True, check_knots=True, receptor=self.receptor, signprot=self.signprot)
        self.gaps = hse.remodel_resis
        print(hse.chain_breaks)

    def make_pirfile(self):
        #fetch remark lines
        with open(self.model_file, 'r') as f:
            lines = f.readlines()
            for line in lines[:5]:
                if line.startswith('REMARK'):
                    self.remark_lines.append(line)
        text = 'REMARK    {} Remodelled loop regions: '.format(len(self.remark_lines)+1)
        for chain, gaps in self.gaps.items():
            for g in gaps:
                text+=chain+' '+str(g[0])+'-'+str(g[1])+' '
        self.remark_lines.append(text+'\n')
        ref_seq, temp_seq = '', ''
        resnum = 0
        self.chains, self.first_resis = [], []
        for chain in self.struct:
            self.chains.append(chain.get_id())
            first_resi = False
            for resi in chain:
                if not first_resi:
                    first_resi = resi.get_id()[1]
                    self.first_resis.append(resi.get_id()[1])
                ref_seq+=PDB.Polypeptide.three_to_one(resi.get_resname())
                if chain.get_id() in self.gaps:
                    for g in self.gaps[chain.get_id()]:
                        if g[0]<=resi.get_id()[1]<=g[1]:
                            temp_seq+='-'
                        else:
                            temp_seq+=PDB.Polypeptide.three_to_one(resi.get_resname())
                else:
                    temp_seq+=PDB.Polypeptide.three_to_one(resi.get_resname())
                resnum = resi.get_id()[1]
            last_chain_id = chain.get_id()
            ref_seq+='/'
            temp_seq+='/'
        # Remove to be remodeled residues from model
        for chain, gaps in self.gaps.items():
            for g in gaps:
                for i in range(g[0],g[1]+1):
                    if i in self.struct[chain]:
                        self.struct[chain].detach_child(self.struct[chain][i].get_id())

        with open(self.pir_file, 'w+') as output_file:
            template="""
>P1;{model_file}
structure:{model_file}:{start}:{first_chain}:{res_num}:{last_chain}::::
{temp_sequence}*

>P1;{uniprot}
sequence:{uniprot}::::::::
{ref_sequence}*
            """
            context={"model_file":self.model_file,
                     "start":self.first_resis[0],
                     "first_chain":self.chains[0],
                     "last_chain":last_chain_id,
                     "res_num":resnum,
                     "temp_sequence":temp_seq,
                     "uniprot":'uniprot',
                     "ref_sequence":ref_seq}
            output_file.write(template.format(**context))
        io = PDB.PDBIO()
        io.set_structure(self.struct)
        io.save(self.model_file)

    def run(self):
        log.none()
        env = environ()
        print(self.gaps)
        m = LoopRemodel(env,
            alnfile = self.pir_file,
            # inimodel=self.model_file,   # initial model of the target
            knowns=self.model_file,
            sequence='uniprot',                 # code of the target
            assess_methods=assess.DOPE,
            gaps=self.gaps,
            model_chains=self.chains,
            start_resnums=self.first_resis,      # assess loops with DOPE
            icl3_delete=self.icl3_delete)

        m.starting_model= 1           # index of the first loop model
        m.ending_model  = 3           # index of the last loop model
        m.md_level = refine.fast  # loop refinement method

        m.make()

        ok_models = [x for x in m.outputs if x['failure'] is None]
        key = 'DOPE score'
        if sys.version_info[:2] == (2,3):
            ok_models.sort(lambda a,b: cmp(a[key], b[key]))
        else:
            ok_models.sort(key=lambda a: a[key])
        ready_model = ok_models[0]

        for mfile in os.listdir("./"):
            if mfile==ready_model['name']:
                with open("./"+mfile, 'r') as f:
                    lines = f.readlines()
                with open(self.model_file, 'w') as f:
                    for l1 in self.remark_lines:
                        f.write(l1)
                    for l2 in lines:
                        if not l2.startswith('REMARK') and not l2.startswith('EXPDTA'):
                            f.write(l2)
            elif mfile.startswith('uniprot'):
                os.remove("./"+mfile)


class LoopRemodel(automodel):
    def __init__(self, env, alnfile, knowns, sequence, assess_methods, gaps=[], model_chains=[], start_resnums=[], icl3_delete=[]):
        super(LoopRemodel, self).__init__(env, alnfile=alnfile, knowns=knowns, sequence=sequence, 
                                               assess_methods=assess_methods)
        self.alnfile = alnfile
        self.gaps = gaps
        self.model_chains = model_chains
        self.start_resnums = start_resnums
        self.icl3_delete = icl3_delete
        
    def special_patches(self, aln):
        # Rename chains and renumber the residues in each
        self.rename_segments(segment_ids=self.model_chains,
                             renumber_residues=self.start_resnums)

    def select_atoms(self):
        selection_out = []

        for chain, gaps in self.gaps.items():
            for g in gaps:
                if chain=='R' and chain in self.icl3_delete:
                    end_range = g[0]+10
                else:
                    end_range = g[1]+1
                for i in range(g[0],end_range):
                    selection_out.append(self.residues[str(i)+':{}'.format(chain)])
        return selection(selection_out)

    # def make(self):
    #     with SilentModeller():
    #         super(LoopRemodel, self).make()


class LoopRemodel2(loopmodel):

    def __init__(self, env, alnfile, knowns, sequence, assess_methods, gaps=[], model_chains=[], start_resnums=[], icl3_delete=[]):
        super(LoopRemodel2, self).__init__(env, alnfile=alnfile, knowns=knowns, sequence=sequence,
                                               assess_methods=assess_methods)

        self.alnfile = alnfile
        self.gaps = gaps
        self.model_chains = model_chains
        self.start_resnums = start_resnums
        self.icl3_delete = icl3_delete
        
    def special_patches(self, aln):
        # Rename chains and renumber the residues in each
        self.rename_segments(segment_ids=self.model_chains,
                             renumber_residues=self.start_resnums)

    def select_loop_atoms(self):
        selection_list = []
        for chain, gap_resis in self.gaps.items():
            # ICL3
            if chain=='R':
                for loop_section in gap_resis:
                    end_range = loop_section[1]-len(self.icl3_delete["R"])
                    selection_list.append(self.residue_range("{}:R".format(loop_section[0]),"{}:R".format(end_range)))
        if len(selection_list)==1:
            return selection(selection_list[0])
        else:
            ### FIXME
            print("WARNING: Need to fix multiple loop remodeling")


class SilentModeller(object):
    ''' No text to console.
    '''
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, *args):
        sys.stdout.close()
        sys.stdout = self._stdout

