from django.conf import settings

from protein.models import Protein, ProteinConformation, ProteinAnomaly, ProteinState, ProteinSegment
from residue.models import Residue
from residue.functions import dgn, ggn
from structure.models import *
from structure.functions import HSExposureCB, PdbStateIdentifier, HomologyModelingSupportFunctions, update_template_source, compare_and_update_template_source
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
import traceback
import subprocess
from copy import deepcopy

gprotein_segments = ProteinSegment.objects.filter(proteinfamily='Alpha')
gprotein_segment_slugs = [i.slug for i in gprotein_segments]
atom_num_dict = {'E':9, 'S':6, 'Y':12, 'G':4, 'A':5, 'V':7, 'M':8, 'L':8, 'I':8, 'T':7, 'F':11, 'H':10, 'K':9, 
                         'D':8, 'C':6, 'R':11, 'P':7, 'Q':9, 'N':8, 'W':14, '-':0}

class SignprotModeling():
    def __init__(self, main_structure, signprot, template_source, trimmed_residues, alignment, main_pdb_array):
        self.main_structure = main_structure
        self.signprot = signprot
        self.template_source = template_source
        self.trimmed_residues = trimmed_residues
        self.a = alignment
        self.main_pdb_array = main_pdb_array
        self.target_signprot = None

    def run(self):
        parse = GPCRDBParsingPDB()
        self.signprot_complex = SignprotComplex.objects.get(structure=self.main_structure)
        structure_signprot= self.signprot_complex.protein
        if self.signprot!=False:
            self.target_signprot = Protein.objects.get(entry_name=self.signprot)
        else:
            self.target_signprot = self.signprot_complex.protein
        self.signprot_protconf = ProteinConformation.objects.get(protein=self.target_signprot)
        sign_a = GProteinAlignment()
        sign_a.run_alignment(self.target_signprot, structure_signprot)
        io = StringIO(self.main_structure.pdb_data.pdb)
        assign_cgn = as_gn.GenericNumbering(pdb_file=io, pdb_code=self.main_structure.pdb_code.index, sequence_parser=True, signprot=structure_signprot)
        signprot_pdb_array = assign_cgn.assign_cgn_with_sequence_parser(self.signprot_complex.alpha)
        new_array = OrderedDict()

        # Initiate complex part of template source
        source_resis = Residue.objects.filter(protein_conformation__protein=self.target_signprot)
        for res in source_resis:
            if res.protein_segment.slug not in self.template_source:
                self.template_source[res.protein_segment.slug] = OrderedDict()
            if res.protein_segment.category=='loop':
                self.template_source[res.protein_segment.slug][str(res.sequence_number)] = [None, None]
            else:
                self.template_source[res.protein_segment.slug][res.display_generic_number.label] = [self.main_structure, self.main_structure]

        # Superimpose missing regions H1 - hfs2
        alt_complex_struct = None
        segs_for_alt_complex_struct = []
        if self.main_structure.pdb_code.index!='3SN6':
            segs_for_alt_complex_struct = ['H1', 'h1ha', 'HA', 'hahb', 'HB', 'hbhc', 'HC', 'hchd', 'HD', 'hdhe', 'HE', 'hehf', 'HF', 'hfs2']
            alt_complex_struct = Structure.objects.get(pdb_code__index='3SN6')
            io = StringIO(alt_complex_struct.pdb_data.pdb)
            alt_signprot_complex = SignprotComplex.objects.get(structure__pdb_code__index='3SN6')
            alt_assign_cgn = as_gn.GenericNumbering(pdb_file=io, pdb_code='3SN6', sequence_parser=True, signprot=alt_signprot_complex.protein)
            alt_signprot_pdb_array = alt_assign_cgn.assign_cgn_with_sequence_parser(alt_signprot_complex.alpha)
            before_cgns = ['G.HN.50', 'G.HN.51', 'G.HN.52', 'G.HN.53']
            after_cgns =  ['G.H5.03', 'G.H5.04', 'G.H5.05', 'G.H5.06']
            orig_residues1 = parse.fetch_residues_from_array(signprot_pdb_array['HN'], before_cgns)
            orig_residues2 = parse.fetch_residues_from_array(signprot_pdb_array['H5'], after_cgns)
            orig_residues = parse.add_two_ordereddict(orig_residues1, orig_residues2)

            alt_residues1 = parse.fetch_residues_from_array(alt_signprot_pdb_array['HN'], before_cgns)
            alt_residues2 = parse.fetch_residues_from_array(alt_signprot_pdb_array['H5'], after_cgns)
            alt_middle = OrderedDict()
            for s in segs_for_alt_complex_struct:
                alt_middle = parse.add_two_ordereddict(alt_middle, alt_signprot_pdb_array[s])
                self.template_source = update_template_source(self.template_source, list(self.template_source[s].keys()), alt_complex_struct, s)

            alt_residues = parse.add_two_ordereddict(parse.add_two_ordereddict(alt_residues1, alt_middle), alt_residues2)
            del_list = []
            for r, t in alt_middle.items():
                if t=='x':
                    del_list.append(r)
            for r in del_list:
                del alt_residues[r]
            
            superpose = sp.LoopSuperpose(orig_residues, alt_residues)
            new_residues = superpose.run()
            key_list = list(new_residues.keys())[4:-4]
            for key in key_list:
                seg = key.split('.')[1]
                signprot_pdb_array[seg][key] = new_residues[key]

            # alt local loop alignment
            alt_sign_a = GProteinAlignment()
            alt_sign_a.run_alignment(self.target_signprot, alt_signprot_complex.protein, segments=segs_for_alt_complex_struct)
            for alt_seg in segs_for_alt_complex_struct:
                sign_a.reference_dict[alt_seg] = alt_sign_a.reference_dict[alt_seg]
                sign_a.template_dict[alt_seg] = alt_sign_a.template_dict[alt_seg]
                sign_a.alignment_dict[alt_seg] = alt_sign_a.alignment_dict[alt_seg]

            # fix h1ha and hahb and hbhc
            if self.target_signprot.entry_name!='gnas2_human':
                h1ha = Residue.objects.filter(protein_conformation__protein=alt_signprot_complex.protein, protein_segment__slug='h1ha')
                h1ha_dict, hahb_dict = OrderedDict(), OrderedDict()
                for h in h1ha:
                    h1ha_dict[h.generic_number.label] = 'x'
                signprot_pdb_array['h1ha'] = h1ha_dict
                right_order = sorted(list(signprot_pdb_array['hahb'].keys()), key=lambda x: (x))
                for r in right_order:
                    hahb_dict[r] = signprot_pdb_array['hahb'][r]
                signprot_pdb_array['hahb'] = hahb_dict
                
            # Let Modeller model buffer regions
            self.trimmed_residues.append('s1h1_6')
            self.trimmed_residues.append('G.S2.01')
            self.trimmed_residues.append('G.S2.02')
            self.trimmed_residues.append('s4h3_4')
            self.trimmed_residues.append('s4h3_5')

        # New loop alignments for signprot. If length differs between ref and temp, buffer is created in the middle of the loop
        loops = [i.slug for i in ProteinSegment.objects.filter(proteinfamily='Alpha', category='loop')]
        loops_to_model = []
        for r_seg, t_seg, a_seg in zip(sign_a.reference_dict, sign_a.template_dict, sign_a.alignment_dict):
            if r_seg in loops:
                loop_length = len(sign_a.reference_dict[r_seg])
                ref_loop = [i for i in list(sign_a.reference_dict[r_seg].values()) if i not in ['x','-']]
                ref_keys = [i for i in list(sign_a.reference_dict[r_seg].keys()) if i not in ['x','-']]
                ref_loop_residues = Residue.objects.filter(protein_conformation__protein=self.target_signprot, protein_segment__slug=r_seg)
                temp_loop = [i for i in list(sign_a.template_dict[t_seg].values()) if i not in ['x','-']]
                temp_keys = [i for i in list(sign_a.template_dict[t_seg].keys()) if i not in ['x','-']]
                if alt_complex_struct and r_seg in segs_for_alt_complex_struct:
                    temp_loop_residues = Residue.objects.filter(protein_conformation__protein=alt_signprot_complex.protein, protein_segment__slug=r_seg)
                else:
                    temp_loop_residues = Residue.objects.filter(protein_conformation__protein=structure_signprot, protein_segment__slug=r_seg)
                ref_out, temp_out, align_out = OrderedDict(), OrderedDict(), OrderedDict()
                # ref is longer
                if len(ref_loop)>len(temp_loop):
                    mid_temp = math.ceil(len(temp_loop)/2)
                    j = 0
                    for i in range(0, loop_length):
                        key = r_seg+'_'+str(i+1)
                        if i+1<=mid_temp:
                            temp_out[key] = temp_loop[i]
                            self.template_source = compare_and_update_template_source(self.template_source, r_seg, signprot_pdb_array, i, ref_loop_residues[i].display_generic_number.label, 
                                                                                      ref_loop_residues[i].sequence_number, segs_for_alt_complex_struct, alt_complex_struct, self.main_structure)
                        elif mid_temp<i+1<=loop_length-mid_temp+1:
                            if i+1==loop_length-mid_temp+1 and len(temp_loop)%2==0:
                                temp_out[key] = temp_loop[mid_temp+j]
                                self.template_source = compare_and_update_template_source(self.template_source, r_seg, signprot_pdb_array, mid_temp+j, ref_loop_residues[i].display_generic_number.label, 
                                                                                          ref_loop_residues[i].sequence_number, segs_for_alt_complex_struct, alt_complex_struct, self.main_structure)
                                j+=1
                            else:
                                temp_out[key.replace('_','?')] = '-'
                                self.template_source = compare_and_update_template_source(self.template_source, r_seg, signprot_pdb_array, mid_temp+j, ref_loop_residues[i].display_generic_number.label, 
                                                                                          ref_loop_residues[i].sequence_number, segs_for_alt_complex_struct, alt_complex_struct, self.main_structure)
                        else:
                            temp_out[key] = temp_loop[mid_temp+j]
                            self.template_source = compare_and_update_template_source(self.template_source, r_seg, signprot_pdb_array, mid_temp+j, ref_loop_residues[i].display_generic_number.label, 
                                                                                      ref_loop_residues[i].sequence_number, segs_for_alt_complex_struct, alt_complex_struct, self.main_structure)
                            j+=1
                    for i, j in enumerate(list(sign_a.reference_dict[r_seg].values())):
                        key = r_seg+'_'+str(i+1)
                        try:
                            temp_out[key]
                            ref_out[key] = j
                        except:
                            ref_out[key.replace('_','?')] = j
                        i+=1
                # temp is longer
                elif len(ref_loop)<len(temp_loop):
                    mid_ref = math.ceil(len(ref_loop)/2)
                    j = 0
                    for i in range(0, loop_length):
                        key = r_seg+'_'+str(i+1)
                        if i+1<=mid_ref:
                            ref_out[key] = ref_loop[i]
                            self.template_source = compare_and_update_template_source(self.template_source, r_seg, signprot_pdb_array, i, temp_loop_residues[i].display_generic_number.label, 
                                                                                      ref_loop_residues[i].sequence_number, segs_for_alt_complex_struct, alt_complex_struct, self.main_structure)
                        elif mid_ref<i+1<=loop_length-mid_ref+1:
                            if i+1==loop_length-mid_ref+1 and len(ref_loop)%2==0:
                                ref_out[key] = ref_loop[mid_ref+j]
                                self.template_source = compare_and_update_template_source(self.template_source, r_seg, signprot_pdb_array, mid_ref+j, temp_loop_residues[i].display_generic_number.label, 
                                                                                          ref_loop_residues[mid_ref+j].sequence_number, segs_for_alt_complex_struct, alt_complex_struct, self.main_structure)
                                j+=1
                            else:
                                ref_out[key.replace('_','?')] = '-'
                                self.template_source = compare_and_update_template_source(self.template_source, r_seg, signprot_pdb_array, mid_ref+j, temp_loop_residues[i].display_generic_number.label, 
                                                                                          ref_loop_residues[mid_ref+j].sequence_number, segs_for_alt_complex_struct, alt_complex_struct, self.main_structure)
                        else:
                            ref_out[key] = ref_loop[mid_ref+j]
                            self.template_source = compare_and_update_template_source(self.template_source, r_seg, signprot_pdb_array, mid_ref+j, temp_loop_residues[i].display_generic_number.label, 
                                                                                      ref_loop_residues[mid_ref+j].sequence_number, segs_for_alt_complex_struct, alt_complex_struct, self.main_structure)
                            j+=1
                    for i, j in enumerate(list(sign_a.template_dict[t_seg].values())):
                        key = r_seg+'_'+str(i+1)
                        try:
                            ref_out[key]
                            temp_out[key] = j
                        except:
                            temp_out[key.replace('_','?')] = j
                        i+=1
                    loops_to_model.append(r_seg)
                # ref and temp length equal
                else:
                    c = 1
                    for i, j in zip(list(sign_a.reference_dict[r_seg].values()), list(sign_a.template_dict[t_seg].values())):
                        ref_out[r_seg+'_'+str(c)] = i
                        temp_out[r_seg+'_'+str(c)] = j
                        self.template_source = compare_and_update_template_source(self.template_source, r_seg, signprot_pdb_array, c-1, temp_loop_residues[c-1].display_generic_number.label, 
                                                                                  ref_loop_residues[c-1].sequence_number, segs_for_alt_complex_struct, alt_complex_struct, self.main_structure)
                        c+=1
                c = 1

                # update alignment dict
                for i, j in zip(list(ref_out.values()), list(temp_out.values())):
                    key = r_seg+'_'+str(c)
                    if i=='-' or j=='-':
                        align_out[key.replace('_','?')] = '-'
                    elif i!=j:
                        align_out[key] = '.'
                    elif i==j:
                        align_out[key] = i
                    c+=1
                # update pdb array
                new_pdb_array = OrderedDict()
                atoms_list = list(signprot_pdb_array[t_seg].values())
                j = 0
                for t_c, t in temp_out.items():
                    jplus1 = False
                    if t!='-':
                        for i in range(j, len(atoms_list)):
                            if atoms_list[j]!='-':
                                new_pdb_array[t_c] = atoms_list[j]
                                jplus1 = True
                                break
                        if jplus1:
                            j+=1
                    else:
                        new_pdb_array[t_c] = 'x'
                        # j+=1

                # pprint.pprint(new_pdb_array)
                # for i,j in new_pdb_array.items():
                #     try:
                #         print(i, PDB.Polypeptide.three_to_one(j[0].get_parent().get_resname()))
                #     except:
                #         print(i, j)

                # update dictionary keys with '?' if no backbone template
                ref_out_final, temp_out_final, align_out_final, new_pdb_array_final = OrderedDict(), OrderedDict(), OrderedDict(), OrderedDict()
                # self.template_source[r_seg] = OrderedDict()
                for i,j in new_pdb_array.items():
                    if '?' not in i and j=='x':
                        ref_out_final[i.replace('_','?').replace('.','?')] = ref_out[i]
                        temp_out_final[i.replace('_','?').replace('.','?')] = temp_out[i]
                        align_out_final[i.replace('_','?').replace('.','?')] = align_out[i]
                        new_pdb_array_final[i.replace('_','?').replace('.','?')] = new_pdb_array[i]
                    else:
                        ref_out_final[i] = ref_out[i]
                        temp_out_final[i] = temp_out[i]
                        align_out_final[i] = align_out[i]
                        new_pdb_array_final[i] = new_pdb_array[i]
                sign_a.reference_dict[r_seg] = ref_out_final
                sign_a.template_dict[t_seg] = temp_out_final
                sign_a.alignment_dict[a_seg] = align_out_final
                signprot_pdb_array[r_seg] = new_pdb_array_final
                
                align_loop = list(sign_a.alignment_dict[a_seg].values())

        self.a.reference_dict = deepcopy(self.a.reference_dict)
        self.a.template_dict = deepcopy(self.a.template_dict)
        self.a.alignment_dict = deepcopy(self.a.alignment_dict)

        for seg, values in sign_a.reference_dict.items():
            new_array[seg] = OrderedDict()
            # self.template_source[seg] = OrderedDict()
            final_values = deepcopy(values)
            for key, res in values.items():
                try:
                    if signprot_pdb_array[seg][key]=='x':
                        new_array[seg][key] = 'x'
                        self.template_source = update_template_source(self.template_source, [key], None, seg)
                    else:
                        new_array[seg][key] = signprot_pdb_array[seg][key]
                except:
                    if res!='-':
                        new_array[seg][key] = '-'
                        self.template_source = update_template_source(self.template_source, [key], None, seg)
            self.a.reference_dict[seg] = final_values
        for seg, values in sign_a.template_dict.items():
            for key, res in values.items():
                if new_array[seg][key]=='x':
                    sign_a.template_dict[seg][key] = 'x'
                else:
                    if new_array[seg][key]=='-':
                        sign_a.template_dict[seg][key] = '-'
                    else:
                        pdb_res = PDB.Polypeptide.three_to_one(new_array[seg][key][0].get_parent().get_resname())
                        if pdb_res!=sign_a.template_dict[seg][key]:
                            sign_a.template_dict[seg][key] = pdb_res
            self.a.template_dict[seg] = sign_a.template_dict[seg]

        for seg, values in sign_a.alignment_dict.items():
            for key, res in values.items():
                if new_array[seg][key]=='x':
                    values[key] = 'x'
            self.a.alignment_dict[seg] = values
        signprot_pdb_array = new_array

        for seg, values in signprot_pdb_array.items():
            self.main_pdb_array[seg] = values

        delete_HN_begin = []
        for i in self.a.reference_dict['HN']:
            if i=='G.HN.30':
                break
            delete_HN_begin.append(i)

        for d in delete_HN_begin:
            del self.a.reference_dict['HN'][d]
            try:
                del self.a.template_dict['HN'][d]
            except:
                pass
            try:
                del self.a.alignment_dict['HN'][d]
            except:
                pass
            del self.main_pdb_array['HN'][d]
            try:
                del self.template_source['HN'][d]
            except:
                pass

        # add residues to model to self.trimmed_residues
        gprot_segments = [i.slug for i in ProteinSegment.objects.filter(proteinfamily='Alpha')]
        for i,j in self.a.reference_dict.items():
            if i in gprot_segments:
                for k,l in j.items():
                    if '?' in k or self.main_pdb_array[i][k] in ['-','x']:
                        self.trimmed_residues.append(k)
                    if i in loops_to_model:
                        self.trimmed_residues.append(k)

        # custom mods
        long_HG_prots = Protein.objects.filter(family__name='Gs')
        if structure_signprot in long_HG_prots and self.target_signprot not in long_HG_prots:
            self.trimmed_residues.append('G.HG.08')
            self.trimmed_residues.append('G.HG.09')
            self.trimmed_residues.append('G.HG.12')
            self.trimmed_residues.append('G.HG.13')
            self.trimmed_residues.append('G.HG.14')
            self.trimmed_residues.append('G.HG.16')
            self.trimmed_residues.append('G.HG.17')
        if structure_signprot!=self.target_signprot or alt_signprot_complex.protein not in [None, self.target_signprot]:
            # hbhc
            hbhc_keys = list(self.a.reference_dict['hbhc'].keys())
            self.trimmed_residues.append(hbhc_keys[2])
            self.trimmed_residues.append(hbhc_keys[3])
            self.trimmed_residues.append(hbhc_keys[-3])
            self.trimmed_residues.append(hbhc_keys[-2])
            # H1
            self.trimmed_residues.append('G.H1.07')
            self.trimmed_residues.append('G.H1.08')
        if 'hgh4' in loops_to_model:
            self.trimmed_residues.append('G.H4.01')
            self.trimmed_residues.append('G.H4.02')
            self.trimmed_residues.append('G.H4.03')

        # Add mismatching residues to trimmed residues for modeling
        for seg, val in self.a.alignment_dict.items():
            if seg in gprotein_segment_slugs:
                for key, res in val.items():
                    if res=='.':
                        self.trimmed_residues.append(key)
        # Add residues with missing atom coordinates to trimmed residues for modeling
        for seg, val in self.main_pdb_array.items():
            if seg in gprotein_segment_slugs:
                for key, atoms in val.items():
                    if atoms not in ['-','x']:
                        if atom_num_dict[PDB.Polypeptide.three_to_one(atoms[0].get_parent().get_resname())]>len(atoms):
                            self.trimmed_residues.append(key)

        # Add Beta and Gamma chains
        p = PDB.PDBParser(QUIET=True).get_structure('structure', StringIO(self.main_structure.pdb_data.pdb))[0]
        beta = p[self.signprot_complex.beta_chain]
        gamma = p[self.signprot_complex.gamma_chain]
        self.a.reference_dict['Beta'] = OrderedDict()
        self.a.template_dict['Beta'] = OrderedDict()
        self.a.alignment_dict['Beta'] = OrderedDict()
        self.main_pdb_array['Beta'] = OrderedDict()
        self.template_source['Beta'] = OrderedDict()
        self.a.reference_dict['Gamma'] = OrderedDict()
        self.a.template_dict['Gamma'] = OrderedDict()
        self.a.alignment_dict['Gamma'] = OrderedDict()
        self.main_pdb_array['Gamma'] = OrderedDict()
        self.template_source['Gamma'] = OrderedDict()
        for b_res in beta:
            key = str(b_res.get_id()[1])
            self.a.reference_dict['Beta'][key] = PDB.Polypeptide.three_to_one(b_res.get_resname())
            self.a.template_dict['Beta'][key] = PDB.Polypeptide.three_to_one(b_res.get_resname())
            self.a.alignment_dict['Beta'][key] = PDB.Polypeptide.three_to_one(b_res.get_resname())
            atoms = [atom for atom in b_res]
            self.main_pdb_array['Beta'][key] = atoms
            self.template_source['Beta'][key] = [self.main_structure, self.main_structure]
        for g_res in gamma:
            key = str(g_res.get_id()[1])
            self.a.reference_dict['Gamma'][key] = PDB.Polypeptide.three_to_one(g_res.get_resname())
            self.a.template_dict['Gamma'][key] = PDB.Polypeptide.three_to_one(g_res.get_resname())
            self.a.alignment_dict['Gamma'][key] = PDB.Polypeptide.three_to_one(g_res.get_resname())
            atoms = [atom for atom in g_res]
            self.main_pdb_array['Gamma'][key] = atoms
            self.template_source['Gamma'][key] = [self.main_structure, self.main_structure]

        # raise AssertionError
        # for i,j,k,l in zip(sign_a.reference_dict, sign_a.template_dict, sign_a.alignment_dict, signprot_pdb_array):
        #     pprint.pprint(self.template_source[i])
        #     for v,b,n,m in zip(sign_a.reference_dict[i], sign_a.template_dict[j], sign_a.alignment_dict[k], signprot_pdb_array[l]):
        #         print(v, b, n, m, sign_a.reference_dict[i][v], sign_a.template_dict[j][b], sign_a.alignment_dict[k][n], signprot_pdb_array[l][m])


class SignprotFunctions(object):
    def __init__(self):
        pass

    def get_subfamilies_with_templates(self):
        subfams = []
        for i in SignprotComplex.objects.all():
            if i.protein.family.name not in subfams:
                subfams.append(i.protein.family.name)
        return subfams

    def get_subfam_subtype_dict(self, subfamilies):
        d = {}
        for s in subfamilies:
            d[s] = [i.entry_name for i in Protein.objects.filter(family__name=s, species__common_name='Human')]
        return d


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
                rotamer = list(Rotamer.objects.filter(structure__protein_conformation=structure.protein_conformation, 
                        residue__sequence_number=gn, structure__preferred_chain=structure.preferred_chain))
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

    def pdb_array_creator(self, structure=None, filename=None):
        ''' Creates an OrderedDict() from the pdb of a Structure object where residue numbers/generic numbers are 
            keys for the residues, and atom names are keys for the Bio.PDB.Residue objects.
            
            @param structure: Structure, Structure object of protein. When using structure, leave filename=None. \n
            @param filename: str, filename of pdb to be parsed. When using filename, leave structure=None).
        '''
        seq_nums_overwrite_cutoff_dict = {'4PHU':2000, '4LDL':1000, '4LDO':1000, '4QKX':1000, '5JQH':1000, '5TZY':2000, '5KW2':2000}
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

        assign_gn = as_gn.GenericNumbering(pdb_file=io, pdb_code=structure.pdb_code.index, sequence_parser=True)
        pdb_struct = assign_gn.assign_generic_numbers_with_sequence_parser()
        pref_chain = structure.preferred_chain
        parent_prot_conf = ProteinConformation.objects.get(protein=structure.protein_conformation.protein.parent)
        parent_residues = Residue.objects.filter(protein_conformation=parent_prot_conf)
        last_res = list(parent_residues)[-1].sequence_number
        if len(pref_chain)>1:
            pref_chain = pref_chain[0]
        for residue in pdb_struct[pref_chain]:
            try:
                if -9.1 < residue['CA'].get_bfactor() < 9.1:
                    gn = str(residue['CA'].get_bfactor())
                    if len(gn.split('.')[1])==1:
                        gn = gn+'0'
                    if gn[0]=='-':
                        gn = gn[1:]+'1'
                    # Exception for 3PBL 331, gn get's assigned wrong
                    if structure.pdb_code.index=='3PBL' and residue.get_id()[1]==331:
                        raise Exception()
                    #################################################
                    if gn in gn_list:
                        if int(residue.get_id()[1])>1000:
                            if structure.pdb_code.index in seq_nums_overwrite_cutoff_dict and int(residue.get_id()[1])>=seq_nums_overwrite_cutoff_dict[structure.pdb_code.index]:
                                gn_array.append(gn)
                                residue_array.append(residue.get_list())
                            else:
                                raise Exception()
                        else:
                            gn_array.append(gn)
                            residue_array.append(residue.get_list())
                    else:
                        raise Exception()
                else:
                    raise Exception()
            except:
                if structure!=None and structure.pdb_code.index in seq_nums_overwrite_cutoff_dict:
                    if int(residue.get_id()[1])>seq_nums_overwrite_cutoff_dict[structure.pdb_code.index]:
                        gn_array.append(str(int(str(residue.get_id()[1])[1:])))
                    else:
                        gn_array.append(str(residue.get_id()[1]))
                else:
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
                        found_res = Residue.objects.get(protein_conformation=parent_prot_conf, sequence_number=gn)
                        if found_res.protein_segment.slug[0] not in ['T','H']:
                            continue
                        try:
                            found_gn = str(ggn(found_res.display_generic_number.label)).replace('x','.')
                        except:
                            found_gn = str(gn)
                        output[found_res.protein_segment.slug][found_gn] = res

        return output

