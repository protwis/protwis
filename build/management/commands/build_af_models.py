from django.core.management.base import BaseCommand, CommandError
from django.core.exceptions import MultipleObjectsReturned
from django.conf import settings
from django.db import connection
from django.utils.text import slugify
from django.db import IntegrityError

from build.management.commands.base_build import Command as BaseBuild
from build.management.commands.build_ligand_functions import get_or_create_ligand, match_id_via_unichem
from protein.models import (Protein, ProteinConformation, ProteinState, ProteinAnomaly, ProteinAnomalyType,
    ProteinSegment)
from residue.models import ResidueGenericNumber, ResidueNumberingScheme, Residue, ResidueGenericNumberEquivalent
from common.models import WebLink, WebResource, Publication
from common.tools import test_model_updates
from structure.models import Structure, StructureType, StructureStabilizingAgent,PdbData, Rotamer, Fragment
from construct.functions import *

from contactnetwork.models import *
import contactnetwork.interaction as ci
from contactnetwork.cube import compute_interactions

from Bio.PDB import PDBParser,PPBuilder
from Bio import pairwise2

from structure.assign_generic_numbers_gpcr import GenericNumbering
from structure.functions import StructureBuildCheck, ParseAFModelsCSV
from ligand.models import Ligand, LigandType, LigandRole, LigandPeptideStructure
from interaction.models import *
from interaction.views import runcalculation_2022, regexaa, check_residue, extract_fragment_rotamer
from residue.functions import dgn

import django.apps
import logging
import os
import re
import yaml
import time
import gc
import pandas as pd
from collections import OrderedDict
from datetime import datetime, date
import json
from urllib.request import urlopen
from Bio.PDB import parse_pdb_header
from Bio.PDB.Selection import *

# import traceback

_mapping_tag = yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG
def dict_constructor(loader, node):
    return OrderedDict(loader.construct_pairs(node))

def represent_ordereddict(dumper, data):
    value = []

    for item_key, item_value in data.items():
        node_key = dumper.represent_data(item_key)
        node_value = dumper.represent_data(item_value)

        value.append((node_key, node_value))

    return yaml.nodes.MappingNode(u'tag:yaml.org,2002:map', value)

yaml.add_representer(OrderedDict, represent_ordereddict)
yaml.add_constructor(_mapping_tag, dict_constructor)

class Command(BaseBuild):
    help = 'Reads source data and creates pdb structure records'

    def add_arguments(self, parser):
        parser.add_argument('-p', '--proc',
            type=int,
            action='store',
            dest='proc',
            default=1,
            help='Number of processes to run')
        parser.add_argument('-s', '--structure',
            dest='structure',
            help='Structure to import (PDB ID)',
            nargs='+')
        parser.add_argument('-u', '--purge',
            action='store_true',
            dest='purge',
            default=False,
            help='Purge existing records')
        parser.add_argument('--skip_cn',
            action='store_true',
            default=False,
            help='Skip building contact network for test build')
        parser.add_argument('-i', '--incremental',
            action='store_true',
            dest='incremental',
            default=False,
            help='Incremental update to structures for small live update')

    tracker = {}
    all_models = django.apps.apps.get_models()[6:]
    test_model_updates(all_models, tracker, initialize=True)

    # source file directory
    pdb_data_dir = os.sep.join([settings.DATA_DIR, 'structure_data', 'pdbs'])

    ### USE below to fix seg ends
    xtal_seg_end_file = os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation', 'mod_xtal_segends.yaml'])
    with open(xtal_seg_end_file, 'r') as f:
        xtal_seg_ends = yaml.load(f, Loader=yaml.Loader)

    xtal_anomalies_file = os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation', 'all_anomalities.yaml'])
    with open(xtal_anomalies_file, 'r') as f2:
        xtal_anomalies = yaml.load(f2, Loader=yaml.Loader)

    xtal_representatives = os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation', 'xtal_representatives.yaml'])
    with open(xtal_representatives, 'r') as f3:
        xtal_representatives = yaml.load(f3, Loader=yaml.Loader)

    non_xtal_seg_end_file = os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation', 'non_xtal_segends.yaml'])
    with open(non_xtal_seg_end_file, 'r') as f:
        non_xtal_seg_ends = yaml.load(f, Loader=yaml.Loader)

    s = ProteinSegment.objects.all()
    segments = {}
    for segment in s:
        segments[segment.slug] = segment

    parsed_pdb = None

    construct_errors, rotamer_errors, contactnetwork_errors, interaction_errors = [],[],[],[]

    with open(os.sep.join([settings.DATA_DIR, 'residue_data', 'unnatural_amino_acids.yaml']), 'r') as f_yaml:
        raw_uaa = yaml.safe_load(f_yaml)
        unnatural_amino_acids = {}
        for i, j in raw_uaa.items():
            unnatural_amino_acids[i] = j

    def handle(self, *args, **options):
        # delete any existing structure data
        if options['purge']:
            try:
                self.purge_structures()
                self.tracker = {}
                test_model_updates(self.all_models, self.tracker, initialize=True)
            except Exception as msg:
                print(msg)
                self.logger.error(msg)

        if options['skip_cn']:
            self.run_contactnetwork=False
        else:
            self.run_contactnetwork=True

        self.parsed_structures = ParseAFModelsCSV()
        self.parsed_structures.parse_ligands()

        if options['structure']:
            self.parsed_structures.complexes = [i for i in self.parsed_structures.complexes if i in options['structure'] or i.lower() in options['structure']]

        if options['incremental']:
            self.incremental_mode = True
        else:
            self.incremental_mode = False

        try:
            self.logger.info('CREATING STRUCTURES')
            # run the function twice (once for representative structures, once for non-representative)
            iterations = 1
            for i in range(1,iterations+1):
                self.prepare_input(options['proc'], self.parsed_structures.complexes, i)
            test_model_updates(self.all_models, self.tracker, check=True)
            self.logger.info('COMPLETED CREATING STRUCTURES')
        except Exception as msg:
            print(msg)
            self.logger.error(msg)

    @staticmethod
    def queryset_iterator(qs, batchsize = 5000, gc_collect = True):
        # See https://www.guguweb.com/2020/03/27/optimize-django-memory-usage/
        iterator = qs.values_list('pk', flat=True).order_by('pk').distinct().iterator()
        eof = False
        while not eof:
            primary_key_buffer = []
            try:
                while len(primary_key_buffer) < batchsize:
                    primary_key_buffer.append(iterator.__next__())
            except StopIteration:
                eof = True
            for obj in qs.filter(pk__in=primary_key_buffer).order_by('pk').iterator():
                yield obj
            if gc_collect:
                gc.collect()

    #Adapted from Gaspar's Original code from build_structures
    def create_rotamers(self, struct, pdb_path, d):
        wt_lookup = {} #used to match WT seq_number to WT residue record
        pdbseq = {} #used to keep track of pdbseq residue positions vs index in seq
        ref_positions = {} #WT postions in alignment
        mapped_seq = {} # index in contruct, tuple of AA and WT [position,AA]
        debug = False
        preferred_chain = struct.preferred_chain
        wt_pdb_lookup = []

        if len(preferred_chain.split(','))>1: #if A,B
            preferred_chain = preferred_chain.split(',')[0]


        AA = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D',
             'CYS':'C', 'GLN':'Q', 'GLU':'E', 'GLY':'G',
             'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K',
             'MET':'M', 'PHE':'F', 'PRO':'P', 'SER':'S',
             'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V',
             'YCM':'C', 'CSD':'C', 'TYS':'Y', 'SEP':'S', 'TPO':'T'} #non-standard AAs

        atom_num_dict = {'E':9, 'S':6, 'Y':12, 'G':4, 'A':5, 'V':7, 'M':8, 'L':8, 'I':8, 'T':7, 'F':11, 'H':10, 'K':9,
                         'D':8, 'C':6, 'R':11, 'P':7, 'Q':9, 'N':8, 'W':14}


        entry_name = d['construct_crystal']['uniprot']

        s = self.parsed_pdb

        chain = s[preferred_chain] #select only one chain (avoid n-mer receptors)

        ppb=PPBuilder()
        seq = ''
        i = 1

        # Checking Na+ atom in xtal
        parent_prot_conf = ProteinConformation.objects.get(protein=struct.protein_conformation.protein)
        try:
            wt_2x50 = Residue.objects.get(protein_conformation=parent_prot_conf, display_generic_number__label=dgn('2x50',parent_prot_conf))
            wt_3x39 = Residue.objects.get(protein_conformation=parent_prot_conf, display_generic_number__label=dgn('3x39',parent_prot_conf))
            if wt_2x50.amino_acid=='D' and wt_3x39.amino_acid=='S':
                if chain[wt_2x50.sequence_number].get_resname()=='ASP' and chain[wt_3x39.sequence_number].get_resname()=='SER':
                    v_2x50 = chain[wt_2x50.sequence_number]['OD1'].get_vector()
                    v_3x39 = chain[wt_3x39.sequence_number]['OG'].get_vector()
                    all_resis = uniqueify(chain)
                    for r in all_resis:
                        id_ = r.get_id()
                        if id_[0]=='H_ NA':
                            v_na = r['NA'].get_vector()
                            d_2x50 = (v_na-v_2x50).norm()
                            d_3x39 = (v_na-v_3x39).norm()
                            if d_2x50<3 and d_3x39<3:
                                struct.sodium = True
                                struct.save()
                                break
        except:
            pass

        check_1000 = 0
        prev_id = 0
        bigjump = False
        all_pdb_residues_in_chain = 0
        for pp in ppb.build_peptides(chain, aa_only=False): #remove >1000 pos (fusion protein / gprotein)
            for i,res in enumerate(pp,1 ):
                all_pdb_residues_in_chain += 1
                residue_id = res.get_full_id()

        for pp in ppb.build_peptides(chain, aa_only=False): #remove >1000 pos (fusion protein / gprotein)
            for i,res in enumerate(pp,1 ):
                id = res.id
                residue_id = res.get_full_id()
                prev_id = id[1]

        i = 1
        for pp in ppb.build_peptides(chain, aa_only=False):
            seq += str(pp.get_sequence()) #get seq from fasta (only chain A)
            for residue in pp:
                residue_id = residue.get_full_id()
                chain = residue_id[2]
                if chain not in pdbseq:
                    pdbseq[chain] = {}
                pos = residue_id[3][1]

                if residue.resname != "NH2": # skip amidation of peptide
                    pdbseq[chain][pos] = [i, AA[residue.resname]]
                    i += 1

        parent_seq_protein = str(struct.protein_conformation.protein.sequence)
        # print(structure.protein_conformation.protein.parent.entry_name)
        rs = Residue.objects.filter(protein_conformation__protein=struct.protein_conformation.protein).prefetch_related('display_generic_number','generic_number','protein_segment')
        #retrieve segment ends from wild type
        if struct.protein_conformation.protein.entry_name in self.non_xtal_seg_ends:
            seg_ends = self.non_xtal_seg_ends[struct.protein_conformation.protein.entry_name]
            # k = next(iter(seg_ends))
            # seg_ends.pop(k)
        else:
            seg_ends = []
            # print('No SEG ENDS info for {}'.format(structure.pdb_code.index))
            print('No SEG ENDS info for {}'.format(struct.protein_conformation.protein.entry_name))
        parent_seq = ""
        for i,r in enumerate(rs, 1): #required to match WT position to a record (for duplication of GN values)
            # wt_lookup[r.sequence_number] = r
            wt_lookup[i] = r
            parent_seq += r.amino_acid

        if len(wt_lookup)==0:
            print("No residues for", struct.protein_conformation.protein.entry_name)
            return None
        #align WT with structure seq -- make gaps penalties big, so to avoid too much overfitting
        pw2 = pairwise2.align.localms(parent_seq, seq, 3, -4, -5, -2)

        gaps = 0
        unmapped_ref = {}
        ref_seq, temp_seq = str(pw2[0][0]), str(pw2[0][1])

        for i, r in enumerate(ref_seq, 1): #loop over alignment to create lookups (track pos)
            if r == "-":
                gaps += 1
            if r != "-":
                ref_positions[i] = [i-gaps,r]
            elif r == "-":
                ref_positions[i] = [None,'-']

            if temp_seq[i-1]=='-':
                unmapped_ref[i-gaps] = '-'

        gaps = 0
        for i, r in enumerate(temp_seq, 1): #make second lookup
            # print(i,r,ref_seq[i-1]) #print alignment for sanity check
            if r == "-":
                gaps += 1
            if r != "-":
                mapped_seq[i-gaps] = [r,ref_positions[i]]

        pdb = struct.pdb_data.pdb
        protein_conformation=struct.protein_conformation
        temp = ''
        check = 0
        errors = 0
        mismatch_seq = 0
        match_seq = 0
        not_matched = 0
        matched_by_pos = 0
        aa_mismatch = 0
        generic_change = 0
        generic_change_from = []

        pdblines_temp = pdb.splitlines()
        pdblines = []
        for line in pdblines_temp: #Get rid of all odd records
            if line.startswith('ATOM') or (line[17:20] in ['YCM','CSD','TYS','SEP'] and line.startswith('HETATM')):
                pdblines.append(line)
            # Only build residues for the first model
            if line.startswith('MODEL        2'):
                break
        pdblines.append('') #add a line to not "run out"
        rotamer_bulk = []
        rotamer_data_bulk = []
        residues_bulk = []

        for i,line in enumerate(pdblines):
            # print(line)
            if line.startswith('ATOM') or (line[17:20] in ['YCM','CSD','TYS','SEP'] and line.startswith('HETATM')):
                # if line[17:20] in ['YCM','CSD','TYS','SEP']: # sanity check for non-standard helix residues
                #     print(line)
                chain = line[21]
                if preferred_chain and chain!=preferred_chain: #If perferred is defined and is not the same as the current line, then skip
                    pass
                else:
                    nextline = pdblines[i+1]
                    residue_number = line[22:26].strip()
                    if (check==0 or nextline[22:26].strip()==check) and nextline.startswith('TER')==False and (nextline.startswith('ATOM')==True or nextline.startswith('HETATM')==True): #If this is either the begining or the same as previous line add to current rotamer
                        temp += line + "\n"
                        #print('same res',pdb.splitlines()[i+1])
                    else: #if this is a new residue
                        #print(pdb.splitlines()[i+1][22:26].strip(),check)
                        temp += line + "\n"
                        #(int(check.strip())<2000 or structure.pdb_code.index=="4PHU") and
                        residue = Residue()
                        residue.sequence_number = int(check.strip())
                        residue.amino_acid = AA[residue_name.upper()]
                        residue.protein_conformation = protein_conformation
                        try:
                            seq_num_pos = pdbseq[chain][residue.sequence_number][0]
                        except:
                            # print('failed residue',pdb_path,residue.sequence_number)
                            temp = "" #start new line for rotamer
                            check = pdblines[i+1][22:26].strip()
                            continue
                        # print('hi',seq_num_pos,residue.sequence_number)
                        if seq_num_pos in mapped_seq:
                            # print(int(check.strip()),seq_num_pos) #sanity check
                            if mapped_seq[seq_num_pos][1][0]==None:
                                # print('no match found') #sanity check
                                # print(residue.sequence_number,residue.amino_acid) #sanity check
                                residue.display_generic_number = None
                                residue.generic_number = None
                                residue.protein_segment = None
                                not_matched +=1
                            else:
                                wt_r = wt_lookup[mapped_seq[seq_num_pos][1][0]]
                                # print(seq_num_pos, mapped_seq[seq_num_pos], wt_r, residue.sequence_number, residue.amino_acid) #sanity check for mapped resis
                                if residue.sequence_number!=wt_r.sequence_number and residue.amino_acid!=wt_r.amino_acid and residue.sequence_number in wt_lookup: #if pos numbers not work -- see if the pos number might be in WT and unmapped
                                    if wt_lookup[residue.sequence_number].amino_acid==residue.amino_acid:
                                        if residue.sequence_number in unmapped_ref: #WT was not mapped, so could be it
                                            # print(residue.sequence_number,residue.amino_acid) #sanity check
                                            # print('wrongly matched, better match on pos+aa',residue.sequence_number,residue.amino_acid,wt_r.sequence_number,wt_r.amino_acid)
                                            wt_r = wt_lookup[residue.sequence_number]
                                            matched_by_pos +=1
                                            match_seq += 1
                                        else:
                                            mismatch_seq += 1
                                            ### REPLACE seq number with WT to fix odd PDB annotation. FIXME kinda dangerous, but best way to ensure consistent GN numbering
                                            residue.sequence_number = wt_r.sequence_number
                                            #print('could have been matched, but already aligned to another position',residue.sequence_number,residue.amino_acid,wt_r.sequence_number,wt_r.amino_acid)
                                    else:
                                        # print('WT pos not same AA, mismatch',residue.sequence_number,residue.amino_acid,wt_r.sequence_number,wt_r.amino_acid)
                                        wt_pdb_lookup.append(OrderedDict([('WT_POS',wt_r.sequence_number), ('PDB_POS',residue.sequence_number), ('AA','.')]))
                                        mismatch_seq += 1
                                        aa_mismatch += 1
                                elif residue.sequence_number!=wt_r.sequence_number:
                                    # print('WT pos not same pos, mismatch',residue.sequence_number,residue.amino_acid,wt_r.sequence_number,wt_r.amino_acid)
                                    wt_pdb_lookup.append(OrderedDict([('WT_POS',wt_r.sequence_number), ('PDB_POS',residue.sequence_number), ('AA',wt_r.amino_acid)]))
                                    if residue.sequence_number in unmapped_ref:
                                        #print('residue.sequence_number',residue.sequence_number,'not mapped though')
                                        if residue.amino_acid == wt_lookup[residue.sequence_number].amino_acid:
                                            #print('they are same amino acid!')
                                            wt_r = wt_lookup[residue.sequence_number]
                                            mismatch_seq -= 1
                                    mismatch_seq += 1
                                    ### REPLACE seq number with WT to fix odd PDB annotation. FIXME kinda dangerous, but best way to ensure consistent GN numbering
                                    ### 2019.01.18 DISABLED underneat, to be sure that sequence number can be found in DB correctly.
                                    # residue.sequence_number = wt_r.sequence_number

                                if residue.amino_acid!=wt_r.amino_acid:
                                    if debug: print('aa mismatch',residue.sequence_number,residue.amino_acid,wt_r.sequence_number,wt_r.amino_acid)
                                    aa_mismatch += 1

                                else:
                                    match_seq += 1
                                if wt_r.generic_number is not None:
                                    residue.display_generic_number = wt_r.display_generic_number
                                    residue.generic_number = wt_r.generic_number
                                else:
                                    residue.display_generic_number = None
                                    residue.generic_number = None
                                    # print('no GN')

                                residue.protein_segment = wt_r.protein_segment
                                residue.missing_gn = False
                                # print(residue, residue.protein_segment, residue.display_generic_number, wt_r, wt_r.protein_segment, wt_r.display_generic_number)
                                if len(seg_ends):
                                    if residue.protein_segment.slug=='TM1':
                                        if seg_ends['1b']!='-' and seg_ends['1e']!='-':
                                            if residue.sequence_number<seg_ends['1b']:
                                                residue.protein_segment = self.segments['N-term']
                                                residue.display_generic_number = None
                                                residue.generic_number = None
                                            elif residue.sequence_number>seg_ends['1e']:
                                                residue.protein_segment = self.segments['ICL1']
                                                residue.display_generic_number = None
                                                residue.generic_number = None
                                    elif residue.protein_segment.slug=='ICL1':
                                        if seg_ends['i1b']!='-' and seg_ends['i1e']!='-':
                                            if residue.sequence_number<seg_ends['i1b'] and residue.sequence_number<=seg_ends['1e']:
                                                residue.protein_segment = self.segments['TM1']
                                            # elif residue.sequence_number>seg_ends['i1e']:
                                            #     residue.protein_segment = self.segments['TM2']
                                            elif (residue.sequence_number>=seg_ends['i1b'] and residue.sequence_number<=seg_ends['i1e']) and residue.generic_number is None:
                                                if debug: print("Missing GN in loop!",residue.sequence_number)
                                                residue.missing_gn = True
                                            elif (residue.sequence_number<seg_ends['i1b'] and residue.sequence_number>seg_ends['i1e']):
                                                if debug: print("Remove former GN for ",residue.generic_number)
                                                residue.display_generic_number = None
                                                residue.generic_number = None
                                        else:
                                            if residue.generic_number is not None:
                                                if debug: print("Remove former GN for ",residue.generic_number)
                                                residue.display_generic_number = None
                                                residue.generic_number = None
                                            if residue.sequence_number<=seg_ends['1e']:
                                                residue.protein_segment = self.segments['TM1']
                                            elif residue.sequence_number>=seg_ends['2b']:
                                                residue.protein_segment = self.segments['TM2']
                                    elif residue.protein_segment.slug=='TM2':
                                        if seg_ends['2b']!='-' and seg_ends['2e']!='-':
                                            if residue.sequence_number<seg_ends['2b']:
                                                residue.protein_segment = self.segments['ICL1']
                                                residue.display_generic_number = None
                                                residue.generic_number = None
                                            elif residue.sequence_number>seg_ends['2e']:
                                                residue.protein_segment = self.segments['ECL1']
                                                residue.display_generic_number = None
                                                residue.generic_number = None
                                    elif residue.protein_segment.slug=='ECL1':
                                        if seg_ends['e1b']!='-' and seg_ends['e1e']!='-' and seg_ends['2e']!='-':
                                            if residue.sequence_number<seg_ends['e1b'] and residue.sequence_number<=seg_ends['2e']:
                                                residue.protein_segment = self.segments['TM2']
                                            # elif residue.sequence_number>seg_ends['e1e']:
                                            #     residue.protein_segment = self.segments['TM3']
                                            elif (residue.sequence_number>=seg_ends['e1b'] and residue.sequence_number<=seg_ends['e1e']) and residue.generic_number is None:
                                                if debug: print("Missing GN in loop!",residue.sequence_number)
                                                residue.missing_gn = True
                                            elif (residue.sequence_number<seg_ends['e1b'] and residue.sequence_number>seg_ends['e1e']):
                                                if debug: print("Remove former GN for ",residue.generic_number)
                                                residue.display_generic_number = None
                                                residue.generic_number = None
                                        else:
                                            if residue.generic_number is not None:
                                                if debug: print("Remove former GN for ",residue.generic_number)
                                                residue.display_generic_number = None
                                                residue.generic_number = None
                                            if residue.sequence_number<=seg_ends['2e']:
                                                residue.protein_segment = self.segments['TM2']
                                            elif residue.sequence_number>=seg_ends['3b']:
                                                residue.protein_segment = self.segments['TM3']
                                    elif residue.protein_segment.slug=='TM3':
                                        if seg_ends['3b']!='-' and seg_ends['3e']!='-':
                                            if residue.sequence_number<seg_ends['3b']:
                                                residue.protein_segment = self.segments['ECL1']
                                                residue.display_generic_number = None
                                                residue.generic_number = None
                                            elif residue.sequence_number>seg_ends['3e']:
                                                residue.protein_segment = self.segments['ICL2']
                                                residue.display_generic_number = None
                                                residue.generic_number = None
                                    elif residue.protein_segment.slug=='ICL2':
                                        if seg_ends['i2b']!='-' and seg_ends['i2e']!='-':
                                            # if residue.sequence_number<seg_ends['i2b']:
                                            #     residue.protein_segment = self.segments['TM3']
                                            if residue.sequence_number>seg_ends['i2e'] and residue.sequence_number>=seg_ends['4b']:
                                                residue.protein_segment = self.segments['TM4']
                                            elif (residue.sequence_number>=seg_ends['i2b'] and residue.sequence_number<=seg_ends['i2e']) and residue.generic_number is None:
                                                if debug: print("Missing GN in loop!",residue.sequence_number)
                                                residue.missing_gn = True
                                            elif (residue.sequence_number<seg_ends['i2b'] and residue.sequence_number>seg_ends['i2e']):
                                                if debug: print("Remove former GN for ",residue.generic_number)
                                                residue.display_generic_number = None
                                                residue.generic_number = None
                                        else:
                                            if residue.generic_number is not None:
                                                if debug: print("Remove former GN for ",residue.generic_number)
                                                residue.display_generic_number = None
                                                residue.generic_number = None
                                            if residue.sequence_number<=seg_ends['3e']:
                                                residue.protein_segment = self.segments['TM3']
                                            elif residue.sequence_number>=seg_ends['4b']:
                                                residue.protein_segment = self.segments['TM4']
                                    elif residue.protein_segment.slug=='TM4':
                                        if seg_ends['4b']!='-' and seg_ends['4e']!='-':
                                            if residue.sequence_number<seg_ends['4b']:
                                                residue.protein_segment = self.segments['ICL2']
                                                residue.display_generic_number = None
                                                residue.generic_number = None
                                            elif residue.sequence_number>seg_ends['4e']:
                                                residue.protein_segment = self.segments['ECL2']
                                                residue.display_generic_number = None
                                                residue.generic_number = None
                                    elif residue.protein_segment.slug=='ECL2':
                                        if seg_ends['e2b']!='-' and seg_ends['e2e']!='-' and seg_ends['4e']!='-':
                                            if residue.sequence_number<seg_ends['e2b'] and residue.sequence_number<=seg_ends['4e']:
                                                residue.protein_segment = self.segments['TM4']
                                            elif residue.sequence_number>seg_ends['e2e'] and residue.sequence_number>=seg_ends['5b']:
                                                residue.protein_segment = self.segments['TM5']
                                        else:
                                            if residue.sequence_number<=seg_ends['4e']:
                                                residue.protein_segment = self.segments['TM4']
                                            elif residue.sequence_number>=seg_ends['5b']:
                                                residue.protein_segment = self.segments['TM5']
                                    elif residue.protein_segment.slug=='TM5':
                                        if seg_ends['5b']!='-' and seg_ends['5e']!='-':
                                            if residue.sequence_number<seg_ends['5b']:
                                                residue.protein_segment = self.segments['ECL2']
                                                residue.display_generic_number = None
                                                residue.generic_number = None
                                            elif residue.sequence_number>seg_ends['5e']:
                                                residue.protein_segment = self.segments['ICL3']
                                                residue.display_generic_number = None
                                                residue.generic_number = None
                                    elif residue.protein_segment.slug=='ICL3':
                                        if residue.sequence_number<=seg_ends['5e']:
                                            residue.protein_segment = self.segments['TM5']
                                        elif residue.sequence_number>=seg_ends['6b']:
                                            residue.protein_segment = self.segments['TM6']
                                    elif residue.protein_segment.slug=='TM6':
                                        if seg_ends['6b']!='-' and seg_ends['6e']!='-':
                                            if residue.sequence_number<seg_ends['6b']:
                                                residue.protein_segment = self.segments['ICL3']
                                                residue.display_generic_number = None
                                                residue.generic_number = None
                                            elif residue.sequence_number>seg_ends['6e']:
                                                residue.protein_segment = self.segments['ECL3']
                                                residue.display_generic_number = None
                                                residue.generic_number = None
                                    elif residue.protein_segment.slug=='TM7':
                                        if seg_ends['7b']!='-' and seg_ends['7e']!='-':
                                            if residue.sequence_number<seg_ends['7b']:
                                                residue.protein_segment = self.segments['ECL3']
                                                residue.display_generic_number = None
                                                residue.generic_number = None
                                            elif residue.sequence_number>seg_ends['7e']:
                                                if seg_ends['8b']!='-':
                                                    if residue.sequence_number<seg_ends['8b']:
                                                        residue.protein_segment = self.segments['ICL4']
                                                        residue.display_generic_number = None
                                                        residue.generic_number = None
                                                else:
                                                    residue.protein_segment = self.segments['C-term']
                                                    residue.display_generic_number = None
                                                    residue.generic_number = None
                                    elif residue.protein_segment.slug=='H8':
                                        if seg_ends['8b']!='-' and seg_ends['8e']!='-':
                                            if residue.sequence_number<seg_ends['8b']:
                                                residue.protein_segment = self.segments['ICL4']
                                                residue.display_generic_number = None
                                                residue.generic_number = None
                                            elif residue.sequence_number>seg_ends['8e']:
                                                residue.protein_segment = self.segments['C-term']
                                                residue.display_generic_number = None
                                                residue.generic_number = None
                                        else:
                                            residue.protein_segment = self.segments['C-term']
                                            residue.display_generic_number = None
                                            residue.generic_number = None
                                    elif residue.protein_segment.slug=='C-term':
                                        if seg_ends['8e']!='-':
                                            if residue.sequence_number<=seg_ends['8e']:
                                                residue.protein_segment = self.segments['H8']

                                    residue.sequence_number = int(check.strip())
                                    #HAVE TO RESET SO IT FITS FOR INTERACTION SCRIPT

                                    if residue.protein_segment != wt_r.protein_segment:
                                        residue.display_generic_number = None
                                        residue.generic_number = None
                                        generic_change += 1
                                        # print(residue.protein_segment.slug[0:2])
                                        if residue.protein_segment.slug[0:2]=="TM" or 1==1:
                                            if debug: print(struct.protein_conformation.protein.entry_name,residue.amino_acid,"XTAL POS",residue.sequence_number, "WT POS",wt_r.sequence_number,"XTAL:",residue.protein_segment,"WT:",wt_r.protein_segment)
                                            pass

                            # print('aa ',residue.sequence_number,residue.amino_acid,residue.display_generic_number, residue.protein_segment)
                            if residue.protein_segment==None:
                                pass
                            else:
                                #print('inserted',residue.sequence_number) #sanity check
                                # residue.save()
                                residues_bulk.append(residue)
                                rotamer_data, created = PdbData.objects.get_or_create(pdb=temp)
                                #rotamer_data_bulk.append(PdbData(pdb=temp))
                                missing_atoms = False
                                if rotamer_data.pdb.startswith('COMPND'):
                                    lines = len(rotamer_data.pdb.split('\n'))-2
                                else:
                                    lines = len(rotamer_data.pdb.split('\n'))
                                if lines<atom_num_dict[residue.amino_acid]:
                                    missing_atoms = True
                                rotamer_data_bulk.append([rotamer_data, missing_atoms])

                        temp = "" #start new line for rotamer
                        check = pdblines[i+1][22:26].strip()

                    check = pdblines[i+1][22:26].strip()
                chain = line[21]
                residue_name = line[17:20].title() #use title to get GLY to Gly so it matches

        bulked_res = Residue.objects.filter(protein_conformation__protein=struct.protein_conformation.protein)
        bulked_rot = rotamer_data_bulk #Rotamer.objects.filter(structure__protein_conformation__protein=struct.protein_conformation.protein)
        rotamer_bulk = []
        for i,res in enumerate(bulked_res):
            rotamer_bulk.append(Rotamer(residue=res, structure=struct, pdbdata=bulked_rot[i][0],
                                    missing_atoms=bulked_rot[i][1]))

        if not os.path.exists(os.sep.join([settings.DATA_DIR, 'structure_data', 'wt_pdb_lookup'])):
            os.makedirs(os.sep.join([settings.DATA_DIR, 'structure_data', 'wt_pdb_lookup']))
        wt_pdb_lookup_folder = os.sep.join([settings.DATA_DIR, 'structure_data', 'wt_pdb_lookup', str(struct) + '.json'])
        if len(wt_pdb_lookup)>0:
            with open(wt_pdb_lookup_folder, 'w') as f2:
                json.dump(wt_pdb_lookup, f2)
        elif os.path.exists(os.sep.join([settings.DATA_DIR, 'structure_data', 'wt_pdb_lookup', str(struct) + '.json'])):
            os.remove(os.sep.join([settings.DATA_DIR, 'structure_data', 'wt_pdb_lookup', str(struct) + '.json']))
        return None

    def build_contact_network(self, location, protein, ligand):
        try:
            # interacting_pairs, distances  = compute_interactions(pdb_code, save_to_db=True)
            compute_interactions(location, protein=protein, lig=ligand, do_peptide_ligand=True, save_to_db=True, file_input=True)
        except:
            print('Error with computing interactions (%s)' % (location))
            self.logger.error('Error with computing interactions (%s)' % (location))
            return

    def purge_structures(self):
        Structure.objects.filter(structure_type__slug='af-peptide').delete()
        ResidueFragmentInteraction.objects.filter(structure_ligand_pair__structure__structure_type__slug='af-peptide').delete()
        # ResidueFragmentInteractionType.objects.all().delete()
        StructureLigandInteraction.objects.filter(structure__structure_type__slug='af-peptide').delete()
        #Remove previous Rotamers/Residues to prepare repopulate
        Fragment.objects.filter(structure__structure_type__slug='af-peptide').delete()
        Rotamer.objects.filter(structure__structure_type__slug='af-peptide').delete()
        # PdbData.objects.all().delete()

    @staticmethod
    def parsecalculation(sd, data, molecule, ignore_ligand_preset=False):
        module_dir = '/tmp/interactions/'
        pdb_id = sd['pdb']
        pdb_name = sd['location'].split('/')[-1]
        complex_name = sd['location'].split('/')[-1].split('-r')[0]
        gpcrdb_id = molecule['gpcrdb id']
        pdb_location = module_dir + 'pdbs/' + complex_name + '/' + pdb_name
        web_resource = WebResource.objects.get(slug='pdb')
        web_link, _ = WebLink.objects.get_or_create(index=pdb_id)
        structure = Structure.objects.filter(pdb_code=web_link)
        if structure.exists():
            structure = Structure.objects.get(pdb_code=web_link)

            if structure.pdb_data is None:
                if os.path.isfile(pdb_location):
                    pdbdata, created = PdbData.objects.get_or_create(pdb=open(pdb_location, 'r').read())  # does this close the file?
                else:
                    print('quitting due to no pdb in filesystem')
                    quit()
                structure.pdb_data = pdbdata
                structure.save()

            protein = structure.protein_conformation
            lig_key = list(data.keys())[0]
            #/tmp/interactions/pdbs/oprd_mouse/oprd_mouse-1643-rank0.pdb
            prot_pep = sd['location'].split('/')[-1].split('-r')[0]
            # /tmp/interactions/results/ranked_0/interaction
            f = module_dir + "results/" + prot_pep + "/interaction/" + prot_pep + "_" + lig_key + ".pdb"
            print(f)

            if os.path.isfile(f):
                pdbdata, created = PdbData.objects.get_or_create(pdb=open(f, 'r').read())  # does this close the file?
            else:
                print('quitting due to no pdb for fragment in filesystem', f)
                quit()

            struct_lig_interactions = StructureLigandInteraction.objects.filter(pdb_reference=lig_key, ligand_id=gpcrdb_id, structure=structure, annotated=True) #, pdb_file=None
            if struct_lig_interactions.exists():  # if the annotated exists
                try:
                    struct_lig_interactions = struct_lig_interactions.get()
                    struct_lig_interactions.pdb_file = pdbdata
                    ligand = struct_lig_interactions.ligand
                except Exception as msg:
                    print('error with duplication structureligand',lig_key,msg)
                    quit() #not sure about this quit
            elif StructureLigandInteraction.objects.filter(pdb_reference=lig_key, structure=structure).exists():
                try:
                    struct_lig_interactions = StructureLigandInteraction.objects.filter(pdb_reference=lig_key, structure=structure).get()
                    struct_lig_interactions.pdb_file = pdbdata
                except StructureLigandInteraction.DoesNotExist: #already there
                    struct_lig_interactions = StructureLigandInteraction.objects.filter(pdb_reference=lig_key, structure=structure, pdb_file=pdbdata).get()
                ligand = struct_lig_interactions.ligand
            else:  # create ligand and pair
                print(pdb_id, "Skipping interactions with ", pdb_id)
                pass

            struct_lig_interactions.save()

            ResidueFragmentInteraction.objects.filter(structure_ligand_pair=struct_lig_interactions).delete()

            for interaction in data[lig_key]['interactions']:
                aa = interaction[0]
                if aa[-1] != structure.preferred_chain:
                    continue
                aa, pos, _ = regexaa(aa)
                residue = check_residue(protein, pos, aa)
                f = interaction[1]
                fragment, rotamer = extract_fragment_rotamer(f, residue, structure, ligand)
                if fragment is not None:
                    interaction_type, created = ResidueFragmentInteractionType.objects.get_or_create(
                                                slug=interaction[2],
                                                name=interaction[3],
                                                type=interaction[4], direction=interaction[5])
                    fragment_interaction, created = ResidueFragmentInteraction.objects.get_or_create(
                                                    structure_ligand_pair=struct_lig_interactions,
                                                    interaction_type=interaction_type,
                                                    fragment=fragment, rotamer=rotamer)
        else:
            print('Something went wrong and we passed')
            pass

    def main_func(self, positions, iteration, count, lock):
        # setting up processes
        complexes = self.parsed_structures.complexes
        complexes = list(set(complexes)) #removing duplicates
        complexes = complexes[730:]
        while count.value < len(complexes):
            print('******************************************')
            cmpx = complexes[count.value]
            print('ITERATION NUMBER {} of {}'.format(count.value, len(complexes)-1))
            print('******************************************')
            print('PARSING DATA FOR COMPLEX: {}'.format(cmpx))
            count.value +=1

            sd = self.parsed_structures.structures[cmpx]
            # print('Building structure {0} with peptide {1}'.format(sd['protein'], sd['peptide_id']))
            self.logger.info('Building structure {0} with peptide {1}'.format(sd['protein'], sd['peptide_id']))

            representative = False
            # does the construct exist?
            try:
                con = Protein.objects.get(entry_name=sd['protein'].lower())
            except Protein.DoesNotExist:
                print('BIG ERROR Construct {} does not exists, skipping!'.format(sd['name'].lower()))
                self.logger.error('Construct {} does not exists, skipping!'.format(sd['name'].lower()))
                continue

            # get the PDB file and save to DB
            sd['pdb'] = 'AFM_' + cmpx.replace('-','_').upper()

            # create a structure record
            # check if there is a ligand
            try:
                struct = Structure.objects.get(protein_conformation__protein=con, pdb_code__index=sd['pdb'], structure_type__slug='af-peptide')
            except Structure.DoesNotExist:
                struct = Structure()

            # protein state
            if 'state' not in sd:
                self.logger.warning('State not defined, using default state {}'.format(settings.DEFAULT_PROTEIN_STATE))
                state = settings.DEFAULT_STATE.title()
            else:
                state = sd['state']
            state_slug = slugify(state)

            try:
                ps, created = ProteinState.objects.get_or_create(slug=state_slug, defaults={'name': state})
                if created:
                    self.logger.info('Created protein state {}'.format(ps.name))
            except IntegrityError:
                ps = ProteinState.objects.get(slug=state_slug)

            struct.representative = representative
            struct.state = ps
            struct.author_state = ps

            # protein conformation
            try:
                struct.protein_conformation = ProteinConformation.objects.get(protein=con)
            except ProteinConformation.DoesNotExist:
                self.logger.error('Protein conformation for construct {} does not exists'.format(con))
                continue

            if struct.protein_conformation.state is not state:
                ProteinConformation.objects.filter(protein=con).update(state=ps)


            pdb_path = sd['location']
            if not os.path.isfile(pdb_path):
                print('Generated model file for protein {} is not available, skipping.'.format(sd['protein']))
                self.logger.info('Generated model file for protein {} is not available, skipping.'.format(sd['protein']))
                continue
            else:
                try:
                    with open(pdb_path, 'r') as pdb_file:
                        pdbdata_raw = pdb_file.read()
                except FileNotFoundError:
                    print('File {} does not exist. Skipping'.format(pdb_path))
                    continue

            pdbdata, created = PdbData.objects.get_or_create(pdb=pdbdata_raw)
            struct.pdb_data = pdbdata

            self.parsed_pdb = PDBParser(PERMISSIVE=True, QUIET=True).get_structure('ref', pdb_path)[0]

            if 'pdb' in sd:
                web_resource = WebResource.objects.get(slug='pdb')
                struct.pdb_code, created = WebLink.objects.get_or_create(index=sd['pdb'], web_resource=web_resource)
            else:
                self.logger.error('PDB code not specified for structure {}, skipping!'.format(sd['pdb']))
                continue

            # structure type
            sd['structure_method'] = sd['model']
            if 'structure_method' in sd and sd['structure_method']:
                structure_type = sd['structure_method'].capitalize()
                structure_type_slug = slugify(sd['structure_method'])
                try:
                    st, created = StructureType.objects.get_or_create(slug=structure_type_slug, defaults={'name': structure_type})
                    if created:
                        self.logger.info('Created structure type {}'.format(st))
                except IntegrityError:
                    st = StructureType.objects.get(slug=structure_type_slug)
                struct.structure_type = st
            else:
                self.logger.warning('No structure type specified in PDB file {}'.format(sd['pdb']))

            try:
                if 'ligand' in sd and sd['ligand'] and sd['ligand']!='None':
                    ligand = sd['ligand'][0]
                    l = Ligand.objects.get(id=ligand['gpcrdb id'])
                else:
                    l = None
                    print('No ligand inserted in the complex {}.'.format(sd['protein'].lower()))
            except:
                print('BIG ERROR Construct {} does not exists, skipping!'.format(sd['name'].lower()))

#################################################

            # insert into plain text fields
            if 'preferred_chain' in sd:
                struct.preferred_chain = sd['preferred_chain']
            else:
                self.logger.warning('Preferred chain not specified for structure {}'.format(sd['pdb']))
            if 'resolution' in sd:
                struct.resolution = float(sd['resolution'])
            else:
                struct.resolution = None
                self.logger.warning('Resolution not specified for structure {}. Setting as null'.format(sd['pdb']))
            if 'publication_date' in sd:
                struct.publication_date = sd['publication_date']
            else:
                sd['publication_date'] = build_date = date.today()
                struct.publication_date = sd['publication_date']
                self.logger.warning('Publication date not specified for structure {}. Defaulting at today'.format(sd['pdb']))

            struct.annotated = True
            struct.refined = False
            struct.stats_text = None

            #Resolution and
            # save structure before adding M2M relations
            struct.save()

            # ligands
            peptide_chain = ""
            pdb_reference = 'pep'
            if 'ligand' in sd and sd['ligand'] and sd['ligand']!='None':
                if isinstance(sd['ligand'], list):
                    ligands = sd['ligand']
                else:
                    ligands = [sd['ligand']]
                for ligand in ligands:
                    l = False
                    peptide_chain = ""
                    if ligand['chain']!='':
                        peptide_chain = ligand['chain']
                        # ligand['name'] = 'pep'
                    if ligand['name'] and ligand['name'] != 'None': # some inserted as none.
                        ligand['type'] = ligand['type'].lower()
                        # use annoted ligand type or default type
                        if ligand['type']:
                            lt, created = LigandType.objects.get_or_create(slug=slugify(ligand['type']),
                                defaults={'name': ligand['type']})
                        else:
                            lt, created = LigandType.objects.get_or_create(
                                slug=slugify(default_ligand_type), defaults={'name': default_ligand_type})

                        if 'title' in ligand and ligand['title']:
                            ligand_title = ligand['title']
                        else:
                            ligand_title = ligand['name']

                        l = Ligand.objects.get(id=ligand['gpcrdb id'])

                        # Create LigandPeptideStructure object to store chain ID for peptide ligands - supposed to b TEMP
                        if ligand['type'] in ['peptide','protein']:
                            lps, created = LigandPeptideStructure.objects.get_or_create(ligand=l, structure=struct, chain=peptide_chain)
                    else:
                        continue

                    # structure-ligand interaction
                    if l and ligand['role']:
                        role_slug = slugify(ligand['role'])
                        try:
                            lr, created = LigandRole.objects.get_or_create(slug=role_slug,
                            defaults={'name': ligand['role']})
                            if created:
                                self.logger.info('Created ligand role {}'.format(ligand['role']))
                        except IntegrityError:
                            lr = LigandRole.objects.get(slug=role_slug)
                        i, created = StructureLigandInteraction.objects.get_or_create(ligand=l, ligand_role=lr, structure=struct, pdb_file=pdbdata, annotated=True, defaults={'pdb_reference': 'pep'})
                        if i.pdb_reference != pdb_reference:
                            i.pdb_reference = pdb_reference
                            i.save()

####################################################
            # protein anomalies
            scheme = struct.protein_conformation.protein.residue_numbering_scheme
            try:
                anomaly_entry = self.xtal_anomalies[sd['protein']]
                segment_codes = {'1':'TM1','12':'ICL1','2':'TM2','23':'ECL1','3':'TM3','34':'ICL2','4':'TM4','5':'TM5','6':'TM6','7':'TM7'}
                all_bulges, all_constrictions = OrderedDict(), OrderedDict()
                for key, val in anomaly_entry.items():
                    if 'x' not in val:
                        continue
                    if key[0] not in segment_codes:
                        continue
                    segment = segment_codes[key.split('x')[0]]
                    if len(key.split('x')[1])==3:
                        try:
                            all_bulges[segment] = all_bulges[segment]+[val]
                        except:
                            all_bulges[segment] = [val]
                    else:
                        try:
                            all_constrictions[segment] = all_constrictions[segment]+[val]
                        except:
                            all_constrictions[segment] = [val]

                if len(all_bulges)>0:
                    pa_slug = 'bulge'
                    try:
                        pab, created = ProteinAnomalyType.objects.get_or_create(slug=pa_slug, defaults={
                            'name': 'Bulge'})
                        if created:
                            self.logger.info('Created protein anomaly type {}'.format(pab))
                    except IntegrityError:
                        pab = ProteinAnomalyType.objects.get(slug=pa_slug)

                    for segment, bulges in all_bulges.items():
                        for bulge in bulges:
                            try:
                                gn, created = ResidueGenericNumber.objects.get_or_create(label=bulge,
                                    scheme=scheme, defaults={'protein_segment': ProteinSegment.objects.get(
                                    slug=segment)})
                                if created:
                                    self.logger.info('Created generic number {}'.format(gn))
                            except IntegrityError:
                                gn =  ResidueGenericNumber.objects.get(label=bulge, scheme=scheme)

                            try:
                                pa, created = ProteinAnomaly.objects.get_or_create(anomaly_type=pab,
                                    generic_number=gn)
                                if created:
                                    self.logger.info('Created protein anomaly {}'.format(pa))
                            except IntegrityError:
                                pa, created = ProteinAnomaly.objects.get(anomaly_type=pab, generic_number=gn)

                            struct.protein_anomalies.add(pa)

                if len(all_constrictions)>0:
                    pa_slug = 'constriction'
                    try:
                        pac, created = ProteinAnomalyType.objects.get_or_create(slug=pa_slug, defaults={
                            'name': 'Constriction'})
                        if created:
                            self.logger.info('Created protein anomaly type {}'.format(pac))
                    except IntegrityError:
                        pac = ProteinAnomalyType.objects.get(slug=pa_slug)

                    for segment, constrictions in all_constrictions.items():
                        for constriction in constrictions:
                            try:
                                gn, created = ResidueGenericNumber.objects.get_or_create(label=constriction,
                                    scheme=scheme, defaults={'protein_segment': ProteinSegment.objects.get(
                                    slug=segment)})
                                if created:
                                    self.logger.info('Created generic number {}'.format(gn))
                            except IntegrityError:
                                gn = ResidueGenericNumber.objects.get(label=constriction, scheme=scheme)

                            try:
                                pa, created = ProteinAnomaly.objects.get_or_create(anomaly_type=pac,
                                    generic_number=gn)
                                if created:
                                    self.logger.info('Created protein anomaly {}'.format(pa))
                            except IntegrityError:
                                pa, created = ProteinAnomaly.objects.get(anomaly_type=pac, generic_number=gn)

                            struct.protein_anomalies.add(pa)
            except KeyError:
                pass

            # save structure
            struct.save()
            #Delete previous interaction data to prevent errors.
            ResidueFragmentInteraction.objects.filter(structure_ligand_pair__structure=struct, structure_ligand_pair__ligand=l).delete()
            #Remove previous Rotamers/Residues to prepare repopulate
            Fragment.objects.filter(structure=struct, ligand=l).delete()
            Rotamer.objects.filter(structure=struct, pdbdata=pdbdata).delete()
            # Residue.objects.filter(protein_conformation=struct.protein_conformation).delete()

####################################################
            d = {}

            try:
                current = time.time()
                #protein = Protein.objects.filter(entry_name=s.protein_conformation).get()
                d = fetch_pdb_info(sd['location'], con, model=True)
                #delete before adding new
                #Construct.objects.filter(name=d['construct_crystal']['pdb_name']).delete()
                # add_construct(d)
                end = time.time()
                diff = round(end - current,1)
                print('construction calculations done for {}. {} seconds.'.format(struct.protein_conformation.protein.entry_name, diff))
            except Exception as msg:
                print(msg)
                print('ERROR WITH CONSTRUCT FETCH {}'.format(sd['pdb']))
                self.logger.error('ERROR WITH CONSTRUCT FETCH for {}'.format(sd['pdb']))
                self.construct_errors.append(struct)

            try:
                current = time.time()
                self.create_rotamers(struct, pdb_path, d)
                # residue_errors = sbc.check_rotamers(s.pdb_code.index)
                # if len(residue_errors)>0:
                #     raise Exception('Error with rotamer check: {}'.format(residue_errors))
                end = time.time()
                diff = round(end - current,1)
                print('Create resides/rotamers done for {}. {} seconds.'.format(struct.protein_conformation.protein.entry_name, diff))
            except Exception as msg:
                print(msg)
                print('ERROR WITH ROTAMERS {}'.format(sd['pdb']))
                self.logger.error('Error with rotamers for {}'.format(sd['pdb']))
                self.rotamer_errors.append(struct)

            try:
                struct.protein_conformation.generate_sites()
            except:
                pass

            try:
                current = time.time()
                # compute_interactions(sd['location'], protein=struct, lig=l, do_peptide_ligand=True, save_to_db=True, file_input=True)
                self.build_contact_network(sd['location'], protein=struct, ligand=l)
                end = time.time()
                diff = round(end - current,1)
                print('Create contactnetwork done for {}.'.format(struct.protein_conformation.protein.entry_name))
            except Exception as msg:
                print(msg)
                print('ERROR WITH CONTACTNETWORK {}'.format(sd['pdb']))
                self.logger.error('Error with contactnetwork for {}'.format(sd['pdb']))

            for ligand in ligands:
                try:
                    current = time.time()
                    peptide_chain = ""
                    if ligand['chain']!='':
                        peptide_chain = ligand['chain']
                    data_results = runcalculation_2022(sd['location'].strip('.pdb'), peptide_chain, file_input=True)
                    self.parsecalculation(sd, data_results, ligand)
                    end = time.time()
                    diff = round(end - current,1)
                    print('Interaction calculations done for {}. {} seconds.'.format(struct.protein_conformation.protein.entry_name, diff))
                except Exception as msg:
                    print(msg)
                    # print(traceback.format_exc())
                    print('ERROR WITH INTERACTIONS {}'.format(sd['pdb']))
                    self.logger.error('Error parsing interactions output for {}'.format(sd['pdb']))
                    self.interaction_errors.append(struct)
