from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.utils.text import slugify

from django.db.models import  F, Q

from common.definitions import GPCR_CLASS_SLUG_PREFIX
from protein.models import Protein, ProteinFamily
from residue.models import Residue
from ligand.models import Ligand, LigandRole, AssayExperiment, Endogenous_GTP

from signprot.views import CouplingBrowser



import os
import logging
import hashlib
import base64
import sqlite3
import csv
import copy
import re
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


TABLE_NAME = 'ligand' #cannot be "U0"
db_file_path = os.path.join(settings.DATA_DIR,"structure_data","ligand_sequence_hash.sqlite3")
csv_cleaned_seq_filename = "cleaned_seqs.csv"
af_scripts_csv_filename = "af_scripts_csv.csv"
hash_names_csv_filename = "hash_names.csv"

max_buffer_size = 1000

HUMAN_G_PROT_BETA = 'gbb1_human' # UniprotKB protein entry name of the human G protein Beta
HUMAN_G_PROT_GAMMA = 'gbg2_human' # UniprotKB protein entry name of the human G protein Gamma

SEGMENT_SPACER = 7
CTER_CUTOFF = 10

GPCR_CLASSES = ['A','B1','F','T2']

# Use an int to force an error
HUMAN_G_PROT_ALFA_DICT = {
('-', '-'): 'gnai1_human', #  {'cml2_human', 'c5ar2_human', 't2r14_human', 'lgr6_human', 'casr_human', 'c3ar_human', 'ta2r4_human', 'ackr3_human', 'ta2r1_human', 'o51e2_human', 'gpr15_human'}
('G12/13', '-'): 'gna12_human',
('G12/13', 'G11'): 10, 
('G12/13', 'G12'): 'gna12_human',
('G12/13', 'G13'): 'gna13_human',
('GPa1 family', '-'): None, # Yeast only
('Gi/o', '-'): 'gnai1_human', #g37l1_human
('Gi/o', 'G11'): 10,
('Gi/o', 'G15'): 10,
('Gi/o', 'Gi1'): 'gnai1_human',
('Gi/o', 'Gi2'): 'gnai2_human',
('Gi/o', 'GoA'): 'gnao_human',
('Gi/o', 'GoB'): 'gnao_human',
('Gi/o', 'Gz'): 'gnaz_human',
('Gq/11', '-'): 'gnaq_human', # fzd6_human: gnai1_human (Gi/o at same rank). mtlr_human: gna12_human (G12/G13 at same rank)
('Gq/11', 'G11'): 'gna11_human',
('Gq/11', 'G12'): 10, 
('Gq/11', 'G14'): 'gna14_human',
('Gq/11', 'G15'): 'gna15_human',
('Gq/11', 'Gq'): 'gnaq_human',
('Gq/11', 'Gz'): 10, 
('Gs', '-'): 'gnas2_human',
('Gs', 'G13'): 20, #LPAR6 wrong. According to source is cannonical sequence of GNA13_HUMAN. And G12 might be stronger.
('Gs', 'Gi1'): 10, 
('Gs', 'GsL'): 'gnas2_human',
('Gs', 'GsS'): 'gnas2_human', #FIX because it is difficult to use non-canonical isoforms.
}

GPCR_SEGMENTS_SLUGS = ['N-term','C-term','ICL3','H8']

MAX_UNKNOWN_AA_FRACTION = 0.5 # if equal to MAX_UNKNOWN_AA_FRACTION is also discarded
MIN_PEPTIDE_LENGTH = 5

valid_aa_pattern_standard_only = re.compile(r'^[ACDEFGHIKLMNPQRSTVWY]+$')
valid_aa_pattern_non_standard_included = re.compile(r'^[ACDEFGHIKLMNPQRSTVWYBXZJUO]+$')

AF_SCRIPTS_CSV_FIELDNAMES = ['gpcr','sequence','c-term','n-term','icl3','clas',
                             'modification','new_seq','starting','plus_numbering','mod_seq','x_start','x_end']

CLASS_NAME_PREFIX = 'Class '


def get_GPCR_ProteinFamilies_from_shorter_name(names_list):

    # 1) Get GPCR class ProteinFamilies slugs
    prefix = CLASS_NAME_PREFIX
    slug_prefix = GPCR_CLASS_SLUG_PREFIX
    name = names_list[0]
    myQs = Q(name__startswith=prefix+name)
    for name in names_list[1:]:
        myQ = Q(name__startswith=prefix+name)
        myQs = myQs | myQ
    q = ProteinFamily.objects.filter(slug__startswith=slug_prefix).filter(myQs).values_list('slug',flat=True)
    slugs = list(q)
    
    # 2) Get GPCR ProteinFamilies from slugs
    slug = slugs[0]
    myQs = Q(slug__startswith=slug)
    for slug in slugs:
        myQ = Q(slug__startswith=slug)
        myQs = myQs | myQ
    q2 = ProteinFamily.objects.filter(slug__startswith=slug_prefix).filter(myQs).values_list('id',flat=True)
    return q2

def get_stimulatory_peptide_like_ligand_AssayExperiment_obj():
    stimulatory_ligand_value_types = ['EC50','pEC50']
    q = AssayExperiment.objects.all().select_related('ligand')
    q = q.filter(value_type__in=stimulatory_ligand_value_types)
    q = q.filter(ligand__ligand_type__slug__in=['peptide','protein'])
    q = q.exclude(ligand__sequence=None).exclude(ligand__sequence='')

    #ChEMBL
    qc = AssayExperiment.objects.exclude(document_chembl_id=None).exclude(assay_description__icontains='antagonist')
    qc = qc.filter(value_type='AC50',assay_description__contains='agonist').select_related('ligand')
    qc = qc.filter(ligand__ligand_type__slug__in=['peptide','protein'])
    qc = qc.exclude(ligand__sequence=None).exclude(ligand__sequence='')
    
    qlf = q|qc

    return qlf

def get_inhibitory_peptide_like_ligand_AssayExperiment_obj():
    inhibitory_ligand_value_types = ['pA2','pKB','pKB','IC50','pIC50']
    q = AssayExperiment.objects.all().select_related('ligand')
    q = q.filter(value_type__in=inhibitory_ligand_value_types)
    q = q.filter(ligand__ligand_type__slug__in=['peptide','protein'])
    q = q.exclude(ligand__sequence=None).exclude(ligand__sequence='')


    #ChEMBL
    qc1 = AssayExperiment.objects.exclude(document_chembl_id=None)
    qc1 = qc1.filter(value_type='AC50',assay_description__icontains='antagonist').select_related('ligand')
    qc1 = qc1.filter(ligand__ligand_type__slug__in=['peptide','protein'])
    qc1 = qc1.exclude(ligand__sequence=None).exclude(ligand__sequence='')

    qc2 = AssayExperiment.objects.exclude(document_chembl_id=None)
    qc2= qc2.filter(value_type='AC50',assay_description__icontains='inhibit').select_related('ligand')
    qc2 = qc2.filter(ligand__ligand_type__slug__in=['peptide','protein'])
    qc2 = qc2.exclude(ligand__sequence=None).exclude(ligand__sequence='')
    
    qf = qc1|qc2

    qlf = q|qf

    return qlf

def get_human_gpcr_AssayExperiment_obj():
    qp = AssayExperiment.objects.annotate(species=F('protein__species__common_name'))
    qp = qp.annotate(protein_family_id=F('protein__family__id'))
    qp = qp.filter(protein_family_id__in=get_GPCR_ProteinFamilies_from_shorter_name(GPCR_CLASSES))
    qp = qp.filter(species__iexact='Human')
    return qp

def get_endogenous_ligand_receptor_id_pairs_set():
    q = Endogenous_GTP.objects.all().values('ligand_id','receptor_id')
    return set([(r['ligand_id'],r['receptor_id'],) for r in q])

def get_ligand_receptor_id_pairs(state,remove_endogenous=False):
    if state == 'active':
        ligand_AssayExperiment_obj = get_stimulatory_peptide_like_ligand_AssayExperiment_obj()
    elif state == 'inactive':
        ligand_AssayExperiment_obj = get_inhibitory_peptide_like_ligand_AssayExperiment_obj()

    human_gpcr_AssayExperiment_obj = get_human_gpcr_AssayExperiment_obj()

    q = ligand_AssayExperiment_obj & human_gpcr_AssayExperiment_obj

    q = q.annotate(entry_name=F('protein__entry_name'),ligand_name=F('ligand__name'))

    fields = ['ligand_id','protein_id','entry_name','ligand_name']
    if remove_endogenous:
        q = [r for r in q.values(*fields) if (r['ligand_id'],r['protein_id'],) not in get_endogenous_ligand_receptor_id_pairs_set()]
    else:
        q = q.values(*fields)
    return q

def get_human_primary_transducer_g_protein_alfa(options=None):
    view_instance = CouplingBrowser()
    tab_fields, header = view_instance.tab_fields(view_instance.subunit_filter, view_instance.families)

    primary_transducers_dict = {}
    for v in tab_fields.values():
        entry_name = v['protein']['entryname']
        sets = v['sets']
        if 'GproteinDb' in sets:
            GproteinDb = sets['GproteinDb']
        elif len(list(sets.keys())) == 1:
            GproteinDb = sets[list(sets.keys())[0]]
        else:
            raise ValueError('Invalid row in GproteinDb Couplings data table view context data.'+\
                             'Multiple biosensors and no merge data row:\n'+str(v))  
        prim_fam = GproteinDb['prim_fam']
        prim_subtype = GproteinDb['prim_subtype']

        #harcoded FIXES:
        if entry_name == 'lpar6_human' and (prim_fam, prim_subtype,) == ('Gs', 'G13',):
            prim_fam = 'G12/13'
        elif entry_name == 'cckar_human' and (prim_fam, prim_subtype,) == ('Gq/11', 'G12',):
            prim_fam = 'G12/13'
        elif entry_name == 'c5ar1_human' and (prim_fam, prim_subtype,) == ('Gi/o', 'G15',):
            prim_fam = 'Gq/11'
        elif entry_name == 'ogr1_human' and (prim_fam, prim_subtype,) == ('G12/13', 'G11'):
            prim_fam = 'Gq/11'

        if prim_fam is None:
            prim_fam = '-'
        
        if prim_fam != '-':
            percent_of_primary_family = GproteinDb[slugify(prim_fam)]['percent_of_primary_family']
        else:
            percent_of_primary_family = '-'

        if prim_subtype is None:
            prim_subtype = '-'

        if prim_subtype != '-':
            percent_of_primary_subtype = GproteinDb[prim_subtype]['percent_of_primary_subtype']
        else:
            percent_of_primary_subtype = '-'

        if options is not None:
            if options['primary_transducers'] or options['primary_transducers_dict']:
                if prim_fam is None:
                    prim_fam = '-'

        if entry_name in primary_transducers_dict: # e.g. calcr_human is duplicated
            old_percent_of_primary_family, old_percent_of_primary_subtype = primary_transducers_dict[entry_name]['percentages']

            if old_percent_of_primary_family == '-':
                old_percent_of_primary_family = 0
            if old_percent_of_primary_subtype == '-':
                old_percent_of_primary_subtype = 0

            new_percent_of_primary_family = percent_of_primary_family
            new_percent_of_primary_subtype = percent_of_primary_subtype

            if new_percent_of_primary_family == '-':
                new_percent_of_primary_family = 0
            if new_percent_of_primary_subtype == '-':
                new_percent_of_primary_subtype = 0

            # 1. If there is only one primary transducer family that has no subtype data, a default
            #    subtype is chosen even if there is subtype data for other families.
            #
            # 2. If there are two primary transducer families (average percentage is the same),
            #    the family with the higher subtype percentage is chosen. The lack of subtype percentage is considered as 0%.
            #
            # 3. If there are two primary transducer subtypes, one of them is chosen randomly.
            #

            if old_percent_of_primary_family != new_percent_of_primary_family:
                print('WARNING: duplicated entry',entry_name)
            if old_percent_of_primary_family == new_percent_of_primary_family and \
               old_percent_of_primary_subtype == new_percent_of_primary_subtype:
                print('WARNING: entry with two or more primary subtypes',entry_name)
            if new_percent_of_primary_family < old_percent_of_primary_family or \
               new_percent_of_primary_family == old_percent_of_primary_family and \
               old_percent_of_primary_subtype >= new_percent_of_primary_subtype:
                continue
        else:
            primary_transducers_dict[entry_name] = {}
        primary_transducers_dict[entry_name]['family_and_subtype'] = (prim_fam,prim_subtype,)
        primary_transducers_dict[entry_name]['percentages'] = (percent_of_primary_family,percent_of_primary_subtype,)

    return primary_transducers_dict

def get_GPCR_segments(entry_names):
    
    # get default conformation. Assumes only a single conformation per protein
    resq = Residue.objects.annotate(protein_entry_name=F('protein_conformation__protein__entry_name'))
    resq = resq.annotate(protein_segment_slug=F('protein_segment__slug'))
    resq = resq.filter(protein_segment_slug__in=GPCR_SEGMENTS_SLUGS)
    resq = resq.filter(protein_entry_name__in=entry_names).order_by('protein_entry_name','protein_segment_slug','sequence_number')
    resq = resq.values('protein_entry_name','protein_segment_slug','sequence_number','amino_acid')

    # concatenate residue records to obtain segment sequences
    protein_segments_dict = {}
    prev_entry_name = None
    for res in resq:
        entry_name = res['protein_entry_name']
        seq_num = res['sequence_number']
        seg_slug = res['protein_segment_slug']
        aa = res['amino_acid']
        if entry_name not in protein_segments_dict:
            seg_slug_dict  = protein_segments_dict[entry_name] = {}
            prev_seq_num = 0
            prev_seg_slug = None
            seq = ''
            not_first = False

        if prev_seg_slug != seg_slug:
            if seg_slug not in seg_slug_dict and not_first:
                seg_slug_dict[seg_slug] = {'sequence':seq,
                                                     'starting_seq_num':start_seq_num,'ending_seq_num':seq_num}
                not_first = True
            start_seq_num = seq_num 
            prev_seq_num = 0
            seq = ''

        if prev_seq_num+1 < seq_num and prev_seg_slug == seg_slug and prev_entry_name == entry_name:
            for i in range(prev_seq_num,seq_num-1):
                seq += 'X'
                print('WARNING: missing residues for '+entry_name+' '+seg_slug+'.')
        if prev_seq_num >= seq_num:
            raise ValueError("Unexpected residue number %i in %s %s." %(seq_num,entry_name,seg_slug))
        seq += aa
        not_first = True
        prev_seq_num = seq_num
        prev_seg_slug = seg_slug
        prev_entry_name = entry_name

        seg_slug_dict[seg_slug] = {'sequence':seq,
                                                'starting_seq_num':start_seq_num,'ending_seq_num':seq_num}
    return protein_segments_dict 
        
        

def is_valid_protein_sequence_fasta_line(sequence,allow_non_standard_amino_acids=False):
    new_sequence = sequence.strip()

    if allow_non_standard_amino_acids:
        valid_aa_pattern = valid_aa_pattern_non_standard_included
    else:
        valid_aa_pattern = valid_aa_pattern_standard_only
        
    return bool(valid_aa_pattern.fullmatch(sequence))

def create_AF_Bio_SeqRecord(seq,chain,name):                

    record = SeqRecord(
                Seq(seq),
                id="Chain "+chain+" "+name,
                description="",
            )
    return record

def write_fasta_file_from_SeqRecord(fasta_path,records,overwrite=False):
    if os.path.exists(fasta_path) and not overwrite:
        raise ValueError("'{}' already exists.".format(fasta_path))
    with open(fasta_path, "w") as fastafile:
        SeqIO.write(records, fastafile, "fasta")
           
def modification_checker(n_ter, c_ter,modification):

    """

    """

    try:
        if modification in ['c', 'b2_trunc', 'b2_st','b2_st_cter', 'b2_trunc_cter', 'c_cter','cter']:
            gpcr_n_len = 0
        elif modification in ['unmod', 'icl3','icl3_cter']:
            gpcr_n_len = len(n_ter)
        else:
            gpcr_n_len = 0
    except:
        gpcr_n_len = 0

    try:
        gpcr_c_len = len(c_ter)
    except:
        gpcr_c_len = 0

    return gpcr_n_len, gpcr_c_len


class Command(BaseCommand):
    help = 'Generate AlphaFold 2 input files. '

    def add_arguments(self, parser):
        super(Command, self).add_arguments(parser=parser)
        parser.add_argument('--output', default=False, action='store', help='Output file.')
        parser.add_argument('--overwrite', default=False, action='store_true', help='Overwrite output directory.')
        parser.add_argument('--verbose', default=False, action='store_true', help='Print progress in stdout.')
        parser.add_argument('--primary-transducers', default=False, action='store_true', help='Print primary transducers.')
        parser.add_argument('--primary-transducers-dict', default=False, action='store_true', help='Print primary transducers dict.')



    logger = logging.getLogger(__name__)


    def handle(self, *args, **options):    
        output_root_path = options['output']
        overwrite = options['overwrite']

        if options['verbose']: print('Retrieving primary transducers...')
        primary_transducers = get_human_primary_transducer_g_protein_alfa(options)

        if options['primary_transducers'] or options['primary_transducers_dict']:
            transducers_set = set()
            for v in primary_transducers.values():
                transducers_set.add(v['family_and_subtype'])
            transducers_list = sorted(sorted(list(transducers_set),key=lambda x: x[1]),key=lambda x: x[0])
            
            if options['primary_transducers']:
                for t in transducers_list:
                    print(t)
            if options['primary_transducers_dict']:
                str_to_append = ": '',"
                print('{')
                for i in range(0,len(transducers_list)-1):
                    print(transducers_list[i].__str__()+str_to_append)
                print(transducers_list[-1].__str__()+str_to_append)
                print('}')
            return

        output_folder = options["output"]
        con = sqlite3.connect(db_file_path)
        
        sql_query = 'SELECT "table_name_0"."sequence_hash", "table_name_0"."sequence_hash_col"'+ \
        'FROM "%s" AS "table_name_0"' % (TABLE_NAME)  + \
        'WHERE "table_name_0"."gpcrdb_pk" = ?'


        ligand_receptor_id_pairs_dict = {k: get_ligand_receptor_id_pairs(k,remove_endogenous=True) for k in ['active','inactive']}

        g_prot_b = HUMAN_G_PROT_BETA
        g_prot_g = HUMAN_G_PROT_GAMMA

        protein_entry_names = set((g_prot_b,g_prot_g,))
        ligand_ids = set()
        
        entry_names_with_issues = {}
        missing_transducer_data_entry_names = set()
        ligand_hash_base_dict = {} # "hash base" is the sequence hash without collision ID
        ligand_id_2_hash_base_dict = {}
        # hash_entry_names_dict = {}
        
        # Open and prepare for later CSV file for storing cleaned ligands
        cleaned_seq_csvfile_path = os.path.join(output_folder,csv_cleaned_seq_filename)
        if os.path.exists(cleaned_seq_csvfile_path) and not overwrite:
            raise ValueError("'{}' already exists.".format(cleaned_seq_csvfile_path))
        with open(cleaned_seq_csvfile_path, 'w', newline='') as cleaned_seq_csvfile:
            cleaned_seq_csvfile_fieldnames = ['ligand_id','ligand_name','old_hash','old_hash_col','old_sequence','cleaned_seq_hash',
                                             'cleaned_seq_hash_col', 'cleaned_sequence']
            cleaned_seq_csvwriter = csv.DictWriter(cleaned_seq_csvfile, fieldnames=cleaned_seq_csvfile_fieldnames)
            cleaned_seq_csvwriter.writeheader()

            
            for k,v in ligand_receptor_id_pairs_dict.items():
                p_to_remove = []
                hashes_dict_set = {}
                for p_index,p in enumerate(v):
                    entry_name = p['entry_name'] #here there are only UniprotKB entry names for the receptors and not for the ligands.
                    ligand_id = p['ligand_id']

                    # Store sequence hashes in ligand_receptor_id_pairs_dict

                    cur = con.cursor()
                    cur.execute(sql_query,(ligand_id,))
                    con.commit()
                    row = cur.fetchone()
                    my_hash_base, my_hash_col = row
                    p['hash_base'] = my_hash_base # "hash base" is the sequence hash without collision ID
                    p['hash_col'] = my_hash_col
                    my_hash = my_hash_base+my_hash_col
                    p['hash'] = my_hash

                    # Store sequence hashes and ligand names in dictionaries for later easy access
                    if my_hash_base not in ligand_hash_base_dict:
                        ligand_hash_base_dict[my_hash_base] =  {}
                    if ligand_id not in ligand_hash_base_dict[my_hash_base]:
                        ligand_hash_base_dict[my_hash_base] = {}
                    ligand_hash_base_dict[my_hash_base][ligand_id] = {'hash_col': my_hash_col}
                    ligand_id_2_hash_base_dict[ligand_id] = {'hash_base':my_hash_base,'hash_col': my_hash_col,
                                                             'ligand_name':p['ligand_name']}
                    
                    # if my_hash not in hash_entry_names_dict:
                    #     hash_entry_names_dict[my_hash] = set()
                    # hash_entry_names_dict[my_hash].add(entry_name)
                    
                    
                    #remove duplicated ligand sequences for the same receptor
                    if entry_name not in hashes_dict_set:
                        hashes_dict_set[entry_name] = set()

                    if my_hash not in hashes_dict_set[entry_name]:
                        hashes_dict_set[entry_name].add(my_hash)
                    else:
                        p_to_remove.append(p_index) 
                        continue
                        pass
                    cur.close()

                    # Associate transducers to the ligand-GPCR pairs and assign folder names
                    if k == 'active':
                        if entry_name in primary_transducers:
                            family_and_subtype_pair = primary_transducers[entry_name]['family_and_subtype']
                            g_prot_a = HUMAN_G_PROT_ALFA_DICT[family_and_subtype_pair]
                            
                        else:
                            g_prot_a = None
                            print('WARNING: NO transducer data for', entry_name)
                            
                        if g_prot_a is None:
                            g_prot_a = HUMAN_G_PROT_ALFA_DICT[('-','-',)]
                            missing_transducer_data_entry_names.add(entry_name)


                        if isinstance(g_prot_a, int):
                            raise ValueError("Invalid family and subtype pair: "+family_and_subtype_pair.__str__())    
                        protein_entry_names.add(g_prot_a)

                        foldername = entry_name+'-'+'hashedseq'+'['+my_hash+']'+'-'+g_prot_a+'_'+g_prot_b+'_'+g_prot_g

                        p['g_prot_a'] = g_prot_a
                    elif k == 'inactive':
                        foldername = entry_name+'-'+'hashedseq'+'['+my_hash+']'
                    ligand_ids.add(ligand_id)
                    protein_entry_names.add(entry_name)
                    p['folder_name'] = foldername
                    
                for i in sorted(p_to_remove, reverse=True):
                    del ligand_receptor_id_pairs_dict[k][i]
            del hashes_dict_set


            print('active count:',len(ligand_receptor_id_pairs_dict['active']))
            print('inactive count:',len(ligand_receptor_id_pairs_dict['inactive']))

            #Get ligand sequences
            q = Ligand.objects.filter(id__in=list(ligand_ids))
            q = list(q.values('id','sequence'))
            ligand_sequences =  {r['id']:r['sequence'] for r in q }
    

            #Ligand clean up, sequences fixes and sequence hashes generation
            ligands_to_discard = set()
            n_ter_acylation_re = re.compile(r'^Ac[-]')
            c_ter_re = re.compile(r'[-]OH$')
            c_ter_amination_re = re.compile(r'[-]NH2$')
            abu_aa_re = re.compile(r'\(Abu\)')
            cleaned_ligands_to_update = {}
            
            for id,seq in ligand_sequences.items():
                if not is_valid_protein_sequence_fasta_line(seq,allow_non_standard_amino_acids=False):
                    new_seq = seq
                    if not is_valid_protein_sequence_fasta_line(seq,allow_non_standard_amino_acids=True):
                        new_seq = n_ter_acylation_re.sub('',seq)
                        new_seq = c_ter_re.sub('',new_seq)
                        new_seq = c_ter_amination_re.sub('',new_seq)
                        new_seq = abu_aa_re.sub('A',new_seq) # replace α-Aminobutyric acid  by Ala

                        if seq == "Ac-FKPLAAaR-OH":
                            new_seq = "FKPLAAAR"
                        elif seq == "Ac-FKPLA(Abu)aR-OH":
                            new_seq = "FKPLAAAR"
                        elif seq == "L-Phe-Phe-Phe":
                            new_seq = "FFF"
                        elif seq == "NWTPNAALYLFGPQa":
                            new_seq = "NWTPNAALYLFGPQ"
                        elif seq == "H-EGTFISDYSIAMDKIK(C16-diacid)QQDFVNWLLAQKGKKNDWKHN-OH":
                            new_seq = "EGTFISDYSIAMDKIKQQDFVNWLLAQKGKKNDWKHN"
                        elif seq == "xYIQNCXLX":
                            new_seq = "CYIQNCPLP"
                        elif seq == "hydrocinnamoyl-XPAWR":   
                            new_seq = "FAPGWR"                    
                        elif not is_valid_protein_sequence_fasta_line(new_seq,allow_non_standard_amino_acids=True):
                            print("WARNING: Ligand ID: %i Invalid sequence: '%s'. Original sequence: '%s' ." % (id,new_seq,seq))
                            continue

                    if is_valid_protein_sequence_fasta_line(new_seq,allow_non_standard_amino_acids=True):
                        if len(new_seq) < MIN_PEPTIDE_LENGTH:
                            ligands_to_discard.add(id)
                            continue
                        if float(new_seq.count('X')) / float(len(new_seq)) >= MAX_UNKNOWN_AA_FRACTION:
                            ligands_to_discard.add(id)
                            continue

                        if seq == "RAAPYGVRLSGREZIRAZIFTSGGSRW":
                            new_seq = seq.replace('Z','M')
                        elif seq == 'Ac-LEGREKVRAQIUUEGXSTWSURKK-NH2':
                            new_seq = 'LEGREKVRAQIAAEGMSTWSARKK'


                        u_count = new_seq.count('U')
                        if new_seq.count('U') > 0:
                            print("WARNING: ligand %i %s with %i Us.  Original sequence: '%s'" % (id,new_seq,u_count,seq))
                        new_seq = new_seq.replace('U','C') # replace selenocystein by cystein
                        new_seq = new_seq.replace('U','C') # replace selenocystein by cystein
                        z_count = new_seq.count('Z')
                        b_count = new_seq.count('B')
                        j_count = new_seq.count('J')
                        o_count = new_seq.count('O')
                        if new_seq.count('Z') > 0 or new_seq.count('B') > 0 or new_seq.count('J') > 0:
                            print("WARNING: ligand %i %s with %i Zs, %i Bs and %i Js. Original sequence: %s" % (id,new_seq,z_count,b_count,j_count, seq))
                        if new_seq.count('O') > 0:
                            print("WARNING: ligand %i %s with %i Os.  Original sequence: '%s'" % (id,new_seq,o_count,seq))

                        new_seq = new_seq.replace('Z','Q') # replace ambiguous E or Q by Q
                        new_seq = new_seq.replace('B','N') # replace ambiguous D or N by N
                        new_seq = new_seq.replace('J','I') # replace ambiguous L or I by I
                        new_seq = new_seq.replace('X','A') # replace ambiguous by A
                    else:
                        raise ValueError("Ligand ID: %i Invalid sequence: '%s' . Original sequence: '%s'" % (id, new_seq, seq))

                    if not is_valid_protein_sequence_fasta_line(new_seq,allow_non_standard_amino_acids=False):
                        raise ValueError("Ligand ID: %i Invalid sequence: '%s' . Original sequence: '%s'" % (id, new_seq, seq))

                    my_new_hash = base64.b32encode(hashlib.md5(new_seq.encode()).digest()).decode().strip('=')
                    cleaned_ligands_to_update[id] = {'new_seq':new_seq,'new_hash':my_new_hash}


            if options['test_ligands']:
                test_ids = [1918,381,906,1373]
                len_test_ids = len(test_ids)
                for i,id in enumerate(test_ids):
                    mixed_ligand_id = test_ids[len_test_ids - i - 1]
                    if i < (len_test_ids / 2):
                        this_hash = ligand_id_2_hash_base_dict[mixed_ligand_id]['hash_base']
                    cleaned_ligands_to_update[id] = {'new_seq':ligand_sequences[mixed_ligand_id],'new_hash':this_hash}
                    


            print('active count:',len(ligand_receptor_id_pairs_dict['active']))
            print('inactive count:',len(ligand_receptor_id_pairs_dict['inactive']))
                
            # filter ligands that do not meet criteria from ligand_receptor_id_pairs_dict
            for k,v in ligand_receptor_id_pairs_dict.items():
                ligand_receptor_id_pairs_dict[k] = [p for p in v if p['ligand_id'] not in ligands_to_discard]
            print('ligands_to_discard',ligands_to_discard)
            print('active count:',len(ligand_receptor_id_pairs_dict['active']))
            print('inactive count:',len(ligand_receptor_id_pairs_dict['inactive']))

            # Add fixed sequences to ligand_receptor_id_pairs_dict
            new_col_count_dict = {}
            data_to_update_dict = {}
            hash_brand_new_cleaned_sequence_dict = {}
            duplicated_seq_to_before_cleaning_ligand_ids = {}
            after_cleaning_hash_2_ligand_ids = {}
            
            for cleaned_ligand_id in cleaned_ligands_to_update:
                ligand_name = ligand_id_2_hash_base_dict[cleaned_ligand_id]['ligand_name']
                my_new_hash = cleaned_ligands_to_update[cleaned_ligand_id]['new_hash']
                cleaned_sequence = cleaned_ligands_to_update[cleaned_ligand_id]['new_seq']
                cleaned_seq_csv_row = {'ligand_id':cleaned_ligand_id,'cleaned_sequence':cleaned_sequence,
                                       'ligand_name':ligand_name}
                cleaned_seq_csv_row['old_hash'] = ligand_id_2_hash_base_dict[cleaned_ligand_id]['hash_base']
                cleaned_seq_csv_row['old_hash_col'] = ligand_id_2_hash_base_dict[cleaned_ligand_id]['hash_col']
                cleaned_seq_csv_row['old_sequence'] = ligand_sequences[cleaned_ligand_id]
                if my_new_hash in ligand_hash_base_dict:
                    ligand_hash_base_dict_current_hash = ligand_hash_base_dict[my_new_hash]
                    cleaned_seq_csv_row['old_hash'] = my_new_hash
                    for ligand_id in ligand_hash_base_dict_current_hash:
                        my_hash_col = ligand_hash_base_dict_current_hash[ligand_id]['hash_col']
                        sequence = ligand_sequences[ligand_id]
                        if cleaned_sequence == sequence:
                            if cleaned_ligand_id != ligand_id:
                                # duplicate between before-cleaning set and after-cleaning set of sequences  
                                if cleaned_ligand_id not in duplicated_seq_to_before_cleaning_ligand_ids:
                                    duplicated_seq_to_before_cleaning_ligand_ids[cleaned_ligand_id] = []
                                duplicated_seq_to_before_cleaning_ligand_ids[cleaned_ligand_id].append(ligand_id)

                                #anotate duplicate in cleaned sequences CSV file
                                cleaned_seq_csv_row['cleaned_seq_hash'] = my_new_hash
                                cleaned_seq_csv_row['cleaned_seq_hash_col'] = my_hash_col
                                cleaned_seq_csvwriter.writerow(cleaned_seq_csv_row)
                        else:
                            # collision with before-cleaning set of sequences
                            if my_new_hash not in new_col_count_dict:
                                col_count = new_col_count_dict[my_new_hash] = 0
                            else:
                                new_col_count_dict[my_new_hash] += 1
                                col_count = new_col_count_dict[my_new_hash]
                            cleaned_seq_new_hash = my_new_hash+'zco'
                            cleaned_seq_new_hash_col = hex(col_count)[2:][::-1]
                            cleaned_seq_csv_row['cleaned_seq_hash'] = cleaned_seq_new_hash
                            cleaned_seq_csv_row['cleaned_seq_hash_col'] = cleaned_seq_new_hash_col
                            cleaned_seq_csvwriter.writerow(cleaned_seq_csv_row)
                            data_to_update_dict[cleaned_ligand_id] = {'cleaned_seq':cleaned_sequence,
                                                                    'cleaned_seq_hash':cleaned_seq_new_hash,
                                                                    'cleaned_seq_hash_col': cleaned_seq_new_hash_col}
                else:
                    # Check for collisions among after-cleaning set of sequences
                    # and set up new special hashes for cleaned-fixed sequences not in the set of before-cleaning
                    # sequences.

                    is_duplicate = False
                    cleaned_seq_new_hash = my_new_hash+'zcn'
                    if my_new_hash not in hash_brand_new_cleaned_sequence_dict:
                        # first sequence hash instance
                        hash_col_count = 0
                        hash_brand_new_cleaned_sequence_dict[my_new_hash] = {'seqs':{cleaned_sequence:hash_col_count},
                                                                             'hash_col_count':hash_col_count}
                    else:
                        hash_brand_new_cleaned_sequence_dict_current = hash_brand_new_cleaned_sequence_dict[my_new_hash]
                        if cleaned_sequence not in hash_brand_new_cleaned_sequence_dict_current['seqs']:
                            #collision
                            hash_brand_new_cleaned_sequence_dict_current['hash_col_count'] += 1
                            hash_col_count  = hash_brand_new_cleaned_sequence_dict_current['hash_col_count']
                            hash_brand_new_cleaned_sequence_dict_current['seqs'][cleaned_sequence] = hash_col_count
                        else:
                            #duplicate
                            hash_col_count  = hash_brand_new_cleaned_sequence_dict_current['hash_col_count']
                            is_duplicate = True
                            
                    cleaned_seq_new_hash_col = hex(hash_col_count)[2:][::-1]
                    cleaned_seq_csv_row['cleaned_seq_hash'] = cleaned_seq_new_hash
                    cleaned_seq_csv_row['cleaned_seq_hash_col'] = cleaned_seq_new_hash_col
                    cleaned_seq_csvwriter.writerow(cleaned_seq_csv_row)

                    data_to_update_dict[cleaned_ligand_id] = {'cleaned_seq':cleaned_sequence,
                                                'cleaned_seq_hash':cleaned_seq_new_hash,
                                                'cleaned_seq_hash_col': cleaned_seq_new_hash_col}
                    cleaned_seq_new_hash_hash_col = cleaned_seq_new_hash+cleaned_seq_new_hash_col
                    if not is_duplicate:
                        after_cleaning_hash_2_ligand_ids[cleaned_seq_new_hash_hash_col] = set()
                    after_cleaning_hash_2_ligand_ids[cleaned_seq_new_hash_hash_col].add(cleaned_ligand_id)

        # update pair data with cleaned seqs data
        for k,v in ligand_receptor_id_pairs_dict.items():
            for p in v:
                ligand_id = p['ligand_id']
                if ligand_id in cleaned_ligands_to_update and ligand_id not in duplicated_seq_to_before_cleaning_ligand_ids:
                    data_to_update_dict_current = data_to_update_dict[ligand_id]
                    for k2,v2 in data_to_update_dict_current.items():
                        p[k2] = v2
                    str_2_replace = 'hashedseq'+'['+p['hash']+']'
                    my_hash = p['cleaned_seq_hash']+p['cleaned_seq_hash_col']
                    str_2_replace_4 = 'hashedseq'+'['+my_hash+']'
                    p['folder_name'] = p['folder_name'].replace(str_2_replace,str_2_replace_4)

        print('active count:',len(ligand_receptor_id_pairs_dict['active']))
        print('inactive count:',len(ligand_receptor_id_pairs_dict['inactive']))

        # removing from ligand_receptor_id_pairs_dict duplicated GPCR-ligand pairs


        # Get sets of ligand IDs with duplicated sequences respect other IDs
        # Ligand ID sets of cleaned duplicate sequences respect of sequences before cleaning 

        temp_list = []
        for ids in duplicated_seq_to_before_cleaning_ligand_ids.values():
            temp_list += ids
        duplicated_seq_to_before_cleaning_non_cleaned_ligand_ids = set(temp_list)
        del temp_list

        duplicated_seq_to_after_cleaning_ligand_ids = [
            ids for ids in after_cleaning_hash_2_ligand_ids.values() if len(ids) > 1
        ]

        # Ligand ID sets of cleaned duplicate sequences respect of sequences after cleaning 
        temp_list = []
        for ids in duplicated_seq_to_after_cleaning_ligand_ids:
            temp_list += ids
        duplicated_seq_to_after_cleaning_ligand_ids_set = set(temp_list)
        del temp_list




        # Get receptor entry names for the duplicated cleaned ligands

        for k,pairs in ligand_receptor_id_pairs_dict.items():

            dup_seq_to_b_clean_lig_ids_entry_names = {'cleaned':{},'non-cleaned':{}}
            dup_seq_to_a_clean_lig_entry_names_2_lig_id = {}
            for p in pairs:
                ligand_id = p['ligand_id']
                entry_name = p['entry_name']
                if ligand_id in duplicated_seq_to_before_cleaning_ligand_ids:
                    clean = 'cleaned'
                elif ligand_id in duplicated_seq_to_before_cleaning_non_cleaned_ligand_ids:
                    clean = 'non-cleaned'
                elif ligand_id in duplicated_seq_to_after_cleaning_ligand_ids_set:
                    if entry_name not in dup_seq_to_a_clean_lig_entry_names_2_lig_id:
                        dup_seq_to_a_clean_lig_entry_names_2_lig_id[entry_name] = set()
                    dup_seq_to_a_clean_lig_entry_names_2_lig_id[entry_name].add(ligand_id)
                    continue
                else:
                    continue
                if ligand_id not in dup_seq_to_b_clean_lig_ids_entry_names[clean]:
                    dup_seq_to_b_clean_lig_ids_entry_names[clean][ligand_id] = set()
                dup_seq_to_b_clean_lig_ids_entry_names[clean][ligand_id].add(entry_name)


            # Remove the duplicated cleaned ligand - receptor pairs
            # PARTIALLY TESTED
            ligand_receptor_id_pairs_to_remove = set()
            for cleaned_ligand_id,cleaned_ligand_entry_names in dup_seq_to_b_clean_lig_ids_entry_names['cleaned'].items():
                for ligand_id in duplicated_seq_to_before_cleaning_ligand_ids[cleaned_ligand_id]:
                    entry_names = dup_seq_to_b_clean_lig_ids_entry_names['non-cleaned'][ligand_id]
                    dup_entry_names = cleaned_ligand_entry_names.intersection(entry_names)
                for entry_name in list(dup_entry_names):
                    ligand_receptor_id_pairs_to_remove.add((cleaned_ligand_id,entry_name,))
            # TESTED
            for entry_name, ligand_ids in dup_seq_to_a_clean_lig_entry_names_2_lig_id.items():
                for dup_seq_ligand_ids in duplicated_seq_to_after_cleaning_ligand_ids:
                    dup_ligand_ids = ligand_ids.intersection(dup_seq_ligand_ids)
                if len(dup_ligand_ids) > 1:
                    ligand_ids_to_remove = sorted(list(dup_ligand_ids))[1:]
                    for ligand_id in ligand_ids_to_remove:
                        ligand_receptor_id_pairs_to_remove.add((ligand_id,entry_name,))

            ligand_receptor_id_pairs_dict[k] = [
                p 
                for p in pairs if (p['ligand_id'],p['entry_name'],) not in ligand_receptor_id_pairs_to_remove
            ]       
            print('ligand_receptor_id_pairs_to_remove:',ligand_receptor_id_pairs_to_remove)
        print('active count:',len(ligand_receptor_id_pairs_dict['active']))
        print('inactive count:',len(ligand_receptor_id_pairs_dict['inactive']))

        if len(entry_names_with_issues.keys()) > 0: print("WARNING: Receptors with transducer data issues: ",entry_names_with_issues)
        if len(missing_transducer_data_entry_names) > 0:
            missing_transducer_data_complexes_count_by_entry_name = {}
            for p in ligand_receptor_id_pairs_dict['active']:
                entry_name = p['entry_name']
                if entry_name in missing_transducer_data_complexes_count_by_entry_name:
                    missing_transducer_data_complexes_count_by_entry_name[entry_name] += 1
                elif entry_name in missing_transducer_data_entry_names:
                    missing_transducer_data_complexes_count_by_entry_name[entry_name] = 1
                
            print("WARNING: %i receptors with missing transducer data: " % (sum([v for v in missing_transducer_data_complexes_count_by_entry_name.values()])),
                  missing_transducer_data_complexes_count_by_entry_name)
            
        #Get protein sequences
        q = Protein.objects.filter(entry_name__in=list(protein_entry_names))
        q = list(q.values('entry_name','sequence'))
        orig_sequences = {r['entry_name']:r['sequence'] for r in q }
        if len(q) < len(protein_entry_names):
            missing_seqs = []
            for p in protein_entry_names:
                if p not in orig_sequences:
                    missing_seqs.append()
            raise ValueError("Missing sequences for: %s" %(','.join(missing_seqs)))
        del q
       
        #Get segments
        gpcr_entry_names = set()
        for k,v in ligand_receptor_id_pairs_dict.items():
            for p in v:
                gpcr_entry_names.add(p['entry_name'])

        gpcr_segments = get_GPCR_segments(gpcr_entry_names)

        # GPCR classes
        q = Protein.objects.filter(entry_name__in=list(gpcr_entry_names))
        gpcr_class_dict = {r.entry_name:r.get_protein_class_from_slug(very_short=True) for r in q }
        del q

        #Process sequences and retrieve CSV fields
        AF_SCRIPTS_CSV_FIELDNAMES = ['gpcr','sequence','c-term','n-term','icl3','clas',
                             'modification','new_seq','starting','plus_numbering','mod_seq','x_start','x_end']
        processed_sequences = {}

        for entry_name in list(gpcr_entry_names):
            new_seq = mod_seq = seq = orig_sequences[entry_name]
            gpcr_class = gpcr_class_dict[entry_name]
            modification = 'unmod'
            plus_numbering = 0
            starting = None
            x_end = len(seq)
            x_start = 0
            
            residues = {i:a for i,a in enumerate(seq,start=1)}
            modifications = []
            gpcr_segments_current = gpcr_segments[entry_name]


            if 'N-term' in gpcr_segments_current:
                nter_seg = gpcr_segments_current['N-term']
                nter_seq = nter_seg['sequence']
                x_start = nter_seg['ending_seq_num']
            else:
                nter_seg = None
                nter_seq = ''
                x_start = 0
            if 'C-term' in gpcr_segments_current:
                cter_seg = gpcr_segments_current['C-term']
                cter_seq = cter_seg['sequence']
            else:
                cter_seg = None
                cter_seq = ''
            if 'ICL3' in gpcr_segments_current:
                icl3_seg = gpcr_segments_current['ICL3']
                icl3_seq = icl3_seg['sequence']
            else:
                icl3_seg = None
                icl3_seq = ''

            if gpcr_class == 'A':
                if icl3_seg is not None:
                    # Process ICL3
                    icl3_len = icl3_seg['ending_seq_num'] - icl3_seg['starting_seq_num'] + 1
                    if icl3_len > SEGMENT_SPACER*2:
                        modifications.append('icl3')
                        plus_numbering = icl3_len - SEGMENT_SPACER*2 #The length of the fragment of sequence removed
                        starting = icl3_seg['starting_seq_num'] - 1 + SEGMENT_SPACER #TO CHECK, it might need a -1.
                        residues = {i:r for i,r in residues.items() if not (i > starting and i < plus_numbering + starting)}
                        mod_seq = new_seq = ''.join([x[1] for x in sorted([(i,r,) for i,r in residues.items()],key = lambda x : x[0])])

            if cter_seg is not None:
                cter_start = cter_seg['starting_seq_num']
                cter_end = cter_seg['ending_seq_num']
                cter_len = cter_end - cter_start + 1
                if CTER_CUTOFF is not None:
                    if cter_len > CTER_CUTOFF:
                        modifications.append('cter')
                        last_res = cter_start - 1 + CTER_CUTOFF
                        residues = {i:r for i,r in residues.items() if not (i > last_res)}
                        mod_seq = ''.join([x[1] for x in sorted([(i,r,) for i,r in residues.items()],key = lambda x : x[0])])
                        x_end = last_res

            if len(modifications) > 0:
                modification = '_'.join(modifications)

            # from pae_mean.py
            n_len, c_len = modification_checker(nter_seq, cter_seq,modification)
            x_start = n_len
            x_end = len(mod_seq) - c_len

            processed_sequences[entry_name] = {'gpcr':entry_name,'sequence':seq,'c-term':cter_seq,'n-term':nter_seq,
                                               'icl3':icl3_seq,
                                               'new_seq':new_seq,'starting':starting,
                                               'plus_numbering': plus_numbering,'mod_seq':mod_seq,
                                               'modification':modification,'clas':gpcr_class,'x_start':x_start,'x_end':x_end}

        #CSV
        af_scripts_csv_csvfile_path = os.path.join(output_folder,af_scripts_csv_filename)
        if os.path.exists(af_scripts_csv_csvfile_path) and not overwrite:
            raise ValueError("'{}' already exists.".format(af_scripts_csv_csvfile_path))
        with open(af_scripts_csv_csvfile_path, 'w', newline='') as af_scripts_csv_csvfile:
            af_scripts_csv_fieldnames = AF_SCRIPTS_CSV_FIELDNAMES
            af_scripts_csv_csvwriter = csv.DictWriter(af_scripts_csv_csvfile, fieldnames=af_scripts_csv_fieldnames)
            af_scripts_csv_csvwriter.writeheader()
            for entry_name in sorted(processed_sequences.keys()):
                af_scripts_csv_csvwriter.writerow(processed_sequences[entry_name])
        self.logger.info('Alphafold inputs generated.')

        hash_names_csvfile_path = os.path.join(output_folder,hash_names_csv_filename)
        if os.path.exists(hash_names_csvfile_path) and not overwrite:
            raise ValueError("'{}' already exists.".format(hash_names_csvfile_path))
        with open(hash_names_csvfile_path, 'w', newline='') as hash_names_csvfile:
            hash_names_csv_fieldnames = ['hash','ligand_name']
            hash_names_csv_csvwriter = csv.DictWriter(hash_names_csvfile, fieldnames=hash_names_csv_fieldnames)
            hash_names_csv_csvwriter.writeheader()

            for k,v in ligand_receptor_id_pairs_dict.items():
                for p in v:
                    #Create output folders
                    
                    foldername = p['folder_name']
                    print(foldername)
                    complex_path = os.path.join(output_root_path,foldername)
                    os.makedirs(complex_path,exist_ok=overwrite)

                    #get ligand data
                    ligand_id = p['ligand_id']
                    ligand_name = ligand_id_2_hash_base_dict[ligand_id]['ligand_name']
                    if 'cleaned_seq_hash' in p:
                        my_hash = p['cleaned_seq_hash']+p['cleaned_seq_hash_col']
                        ligand_seq = p['cleaned_seq']
                    else:
                        my_hash = p['hash']
                        ligand_seq = ligand_sequences[ligand_id]

                    hash_names_csv_csvwriter.writerow({'hash':my_hash,'ligand_name':ligand_name})
                    records = []
                    # receptor
                    record = create_AF_Bio_SeqRecord(processed_sequences[entry_name]['sequence'],'A',entry_name)
                    records.append(record)
                    fasta_path = os.path.join(complex_path,entry_name+".fasta")
                    write_fasta_file_from_SeqRecord(fasta_path,record,overwrite=overwrite)
                    
                    # G protein alpha
                    record = create_AF_Bio_SeqRecord(orig_sequences[g_prot_a],'B',g_prot_a)
                    records.append(record)
                    fasta_path = os.path.join(complex_path,g_prot_a+".fasta")
                    write_fasta_file_from_SeqRecord(fasta_path,record,overwrite=overwrite)
                    if k == 'active':
                        # G protein beta
                        record = create_AF_Bio_SeqRecord(orig_sequences[g_prot_b],'C',g_prot_b)
                        records.append(record)
                        fasta_path = os.path.join(complex_path,g_prot_b+".fasta")
                        write_fasta_file_from_SeqRecord(fasta_path,record,overwrite=overwrite)
                        # G protein gamma
                        record = create_AF_Bio_SeqRecord(orig_sequences[g_prot_g],'D',g_prot_g)
                        records.append(record)
                        fasta_path = os.path.join(complex_path,g_prot_g+".fasta")
                        write_fasta_file_from_SeqRecord(fasta_path,record,overwrite=overwrite)

                    #ligand
                    record = create_AF_Bio_SeqRecord(ligand_seq,'E',my_hash)
                    records.append(record)
                    fasta_path = os.path.join(complex_path,my_hash+".fasta")
                    write_fasta_file_from_SeqRecord(fasta_path,record,overwrite=overwrite)

                    namefile = os.path.join(complex_path,my_hash+".name.txt")
                    if os.path.exists(namefile) and not overwrite:
                        raise ValueError("'{}' already exists.".format(namefile))
                    with open(namefile , "w") as f:
                        print(ligand_name,file=f)

                    fasta_path = os.path.join(complex_path,foldername+".fasta")
                    write_fasta_file_from_SeqRecord(fasta_path,records,overwrite=overwrite)

 