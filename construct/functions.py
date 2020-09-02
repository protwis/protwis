from django.conf import settings
from django.utils.text import slugify
from django.db import IntegrityError, connection
from protein.models import Protein, ProteinConformation
from residue.models import Residue
from structure.models import Structure
from construct.models import *

from ligand.models import Ligand, LigandType, LigandRole
from ligand.functions import get_or_make_ligand

from common.tools import fetch_from_web_api
from urllib.parse import quote
from string import Template
from urllib.request import urlopen
from operator import itemgetter
from itertools import groupby
import time
import json
import datetime
from collections import OrderedDict
import pickle
import logging
import os

AA_three = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
# to override some faulty PDB DBREF entries
uniprot_convert_table = {'Q548Y0_HUMAN':'OX2R_HUMAN'}

# def look_for_value(d,k):
#     ### look for a value in dict if found, give back, otherwise None

def fetch_pdb_info(pdbname,protein,new_xtal=False, ignore_gasper_annotation=False):
    # ignore_gaspar_annotation skips PDB_RANGE edits that mark missing residues as deleted, which messes up constructs.

    if not protein:
        if pdbname=='6ORV':
            with open(os.sep.join([settings.DATA_DIR, 'structure_data', 'pdbs', '6ORV.pdb']), 'r') as pdb6orv:
                pdbdata_raw = pdb6orv.read()
        else:
            url = 'https://www.rcsb.org/pdb/files/%s.pdb' % pdbname
            pdbdata_raw = urlopen(url).read().decode('utf-8')
        # figure out what protein this is
        for line in pdbdata_raw.split('\n'):
            if line.startswith('DBREF'):
                line = line.split()
                if len(line)>7:
                    if line[7] in uniprot_convert_table:
                        uniprot = uniprot_convert_table[line[7]]
                    else:
                        uniprot = line[7]
                    try:
                        p = Protein.objects.get(entry_name=uniprot.lower())
                        slug = p.family.slug
                        print(slug,slug.startswith("00"))
                        if slug.startswith("00"):
                            print('PROTEIN GPCR')
                            protein = p
                    except:
                        pass

    SIFT_exceptions = {'5GLI':[395,403], '5GLH':[395,401]}
    logger = logging.getLogger('build')
    #d = {}
    d = OrderedDict()
    d['construct_crystal'] = {}
    d['construct_crystal']['pdb'] = pdbname
    d['construct_crystal']['pdb_name'] = 'auto_'+pdbname
    try:
        d['construct_crystal']['uniprot'] = protein.parent.entry_name
        d['protein'] = protein.parent.name
        d['wt_seq'] = protein.parent.sequence
    except:
        # if above fails, it's likely due to there being no child, and we are looking at top level
        d['construct_crystal']['uniprot'] = protein.entry_name
        d['protein'] = protein.name
        d['wt_seq'] = protein.sequence


    d['contact_info'] = {}
    d['contact_info']['name_cont'] = 'gpcrdb'
    d['contact_info']['pi_email'] = 'info@gpcrdb.org'
    d['contact_info']['pi_name'] = 'gpcrdb'
    d['contact_info']['url'] = 'gpcrdb.org'
    d['contact_info']['date'] = time.strftime('%m/%d/%Y')
    d['contact_info']['address'] = ''

    d['pdb'] = pdbname
    d['links'] = []
    d['xml_not_observed'] = []
    d['xml_segments'] = []

    pos_in_wt = list(range(1,len(d['wt_seq'])+1))

    # GET PDB FILE TO GET INITIAL VALUES - remove known WT that do not exist
    pdbdata_raw = None
    try:
        structure = Structure.objects.filter(pdb_code__index=d['construct_crystal']['pdb'].upper()).get()
        if structure.pdb_data.pdb:
            pdbdata_raw = structure.pdb_data.pdb
    except:
        pass
    if not pdbdata_raw:
        pdb_data_dir = os.sep.join([settings.DATA_DIR, 'structure_data', 'pdbs'])
        pdb_path = os.sep.join([pdb_data_dir, pdbname + '.pdb'])
        if not os.path.isfile(pdb_path):
            url = 'https://www.rcsb.org/pdb/files/%s.pdb' % pdbname
            pdbdata_raw = urlopen(url).read().decode('utf-8')
            with open(pdb_path, 'w') as f:
                f.write(pdbdata_raw)
        else:
            with open(pdb_path, 'r') as pdb_file:
                pdbdata_raw = pdb_file.read()

    if new_xtal==False:
        try:
            structure = Structure.objects.filter(pdb_code__index=d['construct_crystal']['pdb'].upper()).get()

            if 1==1: #update pdbs
                structure.pdb_data.pdb = pdbdata_raw
                structure.pdb_data.save()

            pdb_file = structure.pdb_data.pdb
        except:
            pdb_file = pdbdata_raw
        # print(pdb_file)
    else:
        pdb_file = pdbdata_raw
    pdb_range = []
    uniprot_code = ''

    # do uniprot_code check to see if they only label with PDB code
    for line in pdb_file.split('\n'):
        if line.startswith('DBREF'):
            line = line.split()
            if len(line)>7:
                if line[7] in uniprot_convert_table:
                    uniprot = uniprot_convert_table[line[7]]
                else:
                    uniprot = line[7]
                if uniprot == d['construct_crystal']['uniprot'].upper():
                    uniprot_code = line[6]

    for line in pdb_file.split('\n'):
        if line.startswith('DBREF'):
            # print(line)
            line = line.split()
            if len(line)<8:
                continue
            if line[7] in uniprot_convert_table:
                uniprot = uniprot_convert_table[line[7]]
            else:
                uniprot = line[7]
            start = line[8]
            end = line[9]
            # print(line,uniprot,d['construct_crystal']['uniprot'].upper()) #show DBREF
            if uniprot == d['construct_crystal']['uniprot'].upper() or (uniprot==pdbname.upper() and uniprot_code==''):
                uniprot_code = line[6]
                # print(line)
                pdb_range += range(int(start),int(end)+1)
        elif line.startswith('SEQADV'):
            line = line.split()
            # if it is relevant to correct uniprot
            if line[7]=='DELETION':
                if line[4]==uniprot_code:
                        #remove from pdb_range
                        # print(line)
                        # pdb_range += [int(line[6])]
                        if int(line[6]) in pdb_range:
                            pdb_range.remove(int(line[6]))
            elif len(line)>10 and line[10]=='MUTATION':
                # print("mutation",line)
                pass
            else:
                # print('uknonwn',line)
                pass

    #replace if we could create via pdb_file
    dbref_found = False
    # print(pdb_range)
    # print(list(set(pdb_range)))
    # print("pos_in_wt",pos_in_wt)

    # 5T1A: The final model included 295 residues (37–225 and 241–319) of the 360 residues of CCR2
    # The sequence of human CCR2 isoform B (Uniprot ID P41597-2) was engineered for crystallization by truncation of C-terminal residues 329–360
    # ISOFORM: 314-374: SLFHIALGCR...EASLQDKEGA → RYLSVFFRKH...TGEQEVSAGL
    # DBREF  5T1A A    2   233  UNP    P41597   CCR2_HUMAN       2    233
    # DBREF: DBREF  5T1A A  234   328  UNP    P41597   CCR2_HUMAN     234    328
    if pdbname.upper()=='5T1A':
        pdb_range = list(range(2,226))+list(range(241,329))

    # 4Z9G
    # GLN   103 EXPRESSION TAG
    # DBREF  4Z9G A  103   220  UNP    P34998   CRFR1_HUMAN    103    220
    # DBREF  4Z9G A  223   372  UNP    P34998   CRFR1_HUMAN    252    401 (223-372 in right isoform)
    if pdbname.upper()=='4Z9G':
        pdb_range = list(range(104,221))+list(range(223,402))

    # 5GLI 304 - 310 arent in pdb?
    # ['DBREF', '5GLI', 'A', '63', '416', 'PDB', '5GLI', '5GLI', '63', '416'] <-- unparsed?
    #  The C terminus was truncated after Ser407, and three cysteine residues were mutated to alanine (C396A, C400A and C405A) to avoid heterogeneous palmitoylation
    # T4 lysozyme containing the C54T and C97A mutations52 was introduced into intracellular loop 3, between Lys3035.68 and Leu3116.23 (ETBR-Y5-T4L)
    # a tobacco etch virus (TEV) protease recognition sequence was introduced between Gly57 and Leu66,
    # WT starts at 30 till 57, then TEV, then WT-66, since TEV cleaves, from 66 is in XTAL
    if pdbname.upper()=='5GLI' or pdbname.upper()=='5GLH':
        pdb_range = list(range(66,304))+list(range(311,408))

    # https://www.nature.com/nature/journal/v546/n7657/fig_tab/nature22378_SF1.html
    if pdbname.upper()=='5VEW' or pdbname.upper()=='5VEX':
        pdb_range = list(range(128,205))+list(range(214,258))+list(range(261,432))

    # The refined structure contains receptor residues 19–330 with the segment of residues 230–242 within ICL3 replaced by rubredoxin
    # ['DBREF', '5VBL', 'B', '7', '229', 'UNP', 'P35414', 'APJ_HUMAN', '7', '229']
    # ['DBREF', '5VBL', 'B', '243', '330', 'UNP', 'P35414', 'APJ_HUMAN', '243', '330']
    if pdbname.upper()=='5VBL':
        pdb_range = list(range(19,230))+list(range(243,331))

    # Fix 327 being labelled as there, since it is part of the expression tag
    if pdbname.upper()=='4Z36':
        pdb_range = list(range(2,233))+list(range(249,327))

    # http://www.pnas.org/content/suppl/2014/01/22/1317903111.DCSupplemental/pnas.201317903SI.pdf
    # Amino acids V280-I295 were deleted in the constructs ΔIC3A (3ZEV + 4BV0) and E273-T290 in ΔIC3B (4BUO).
    if pdbname.upper()=='3ZEV':
        pdb_range = list(range(50,280))+list(range(296,391))

    # REMARK 999 THE AUTHORS STATE THAT THE SEQUENCE "NV" IS THE ORIGINAL SEQUENCE
    # REMARK 999 FROM ENDOTHELIUM, AS OPPOSED TO KSL WHICH IS THE GENOMIC SEQUENCE.
    # REMARK 999 THE UNIPROT RECORD WAS CHANGED TO THE GENOMIC SEQUENCE "KSL" AFTER
    # REMARK 999 THE CLONE WAS CREATED AND THE CONSTRUCT WAS NOT CHANGED BECAUSE IT
    # REMARK 999 WAS PERFORMING WELL.
    # SEQADV 3V2W     A       UNP  P21453    LYS   250 SEE REMARK 999
    # SEQADV 3V2W ASN A  251  UNP  P21453    SER   251 SEE REMARK 999
    # SEQADV 3V2W VAL A  252  UNP  P21453    LEU   252 SEE REMARK 999
    if pdbname.upper()=='3V2W' or pdbname.upper()=='3V2Y':
        pdb_range = list(range(2,232))+list(range(244,250))+list(range(251,327))

    # COMPND   6 OTHER_DETAILS: RESIDUES 3-32 AT THE N-TERMINUS AND RESIDUES 244-271
    # COMPND   7  OF THE THIRD INTRACELLULAR LOOP WERE DELETED FROM THE CONSTRUCT.
    # COMPND   8  THE CONSTRUCT WAS TRUNCATED AFTER RESIDUE 367 AND A HEXAHIS TAG
    # COMPND   9  ADDED.
    # if pdbname.upper()=='5A8E':
    #     pdb_range = list(range(33,244))+list(range(272,367))

    # US28 wastruncated by 10 amino acids (1-10) at the N-terminus and 44 amino acids (311-354) at the C-terminus (US28∆N∆C) (Fig. S5B);
    # http://science.sciencemag.org/content/sci/suppl/2015/03/04/347.6226.1113.DC1/Burg.SM.pdf
    if pdbname.upper()=='4XT1' or pdbname.upper()=='4XT3':
        pdb_range = list(range(11,311))

    # Issue with isoform
    if pdbname.upper()=='6NIY':
        pdb_range = list(range(1,475))

    # # Maybe this messes with construct build, temp fix
    # if pdbname.upper()=='6K41':
    #     pdb_range = list(range(11,79))+list(range(85,150))+list(range(164,167))+list(range(170,204))+list(range(361,394))+list(range(404,442))
    # if pdbname.upper()=='6K42':
    #     pdb_range = list(range(11,152))+list(range(162,205))+list(range(358,399))+list(range(401,446))
    # if pdbname.upper()=='6KUY':
    #     pdb_range = list(range(33,173))+list(range(184,228))+list(range(365,444))
    # if pdbname.upper()=='6KUX':
    #     pdb_range = list(range(29,174))+list(range(183,228))+list(range(365,444))

    if not ignore_gasper_annotation:
        ## THIS BLOCK IS DONE BY GASPAR -- IT ALSO REMOVES MISSING RESIDUES, NOT TO BE USED TO WRITE CONSTRUCTS
        # Misannotated DBREF in PDB file
        if pdbname.upper()=='6QZH':
            pdb_range = list(range(52,160))+list(range(165,252))+list(range(262,289))+list(range(297,342))
        if pdbname.upper()=='3SN6':
            pdb_range = list(range(30,366))
        elif pdbname.upper()=='5ZKP':
            pdb_range = list(range(6,124))+list(range(138,217))+list(range(224,316))
        elif pdbname.upper()=='5L7D':
            pdb_range = list(range(58,429))+list(range(446,552))
        elif pdbname.upper()=='5L7I':
            pdb_range = list(range(58,429))+list(range(446,553))
        elif pdbname.upper()=='5WIV':
            pdb_range = list(range(32,177))+list(range(182,228))+list(range(383,465))
        elif pdbname.upper()=='5WIU':
            pdb_range = list(range(34,177))+list(range(182,228))+list(range(383,463))
        elif pdbname.upper()=='5YC8' or pdbname.upper()=='5ZKC':
            pdb_range = list(range(16,215))+list(range(380,459))
        elif pdbname.upper()=='5ZK8' or pdbname.upper()=='5ZK3':
            pdb_range = list(range(18,215))+list(range(383,459))
        elif pdbname.upper()=='5V54':
            pdb_range = list(range(38,192))+list(range(198,240))+list(range(305,389))
        elif pdbname.upper()=='6D32':
            pdb_range = list(range(36,402))+list(range(416,526))
        elif pdbname.upper() in ['6H7N','6H7L','6H7M']:
            pdb_range = list(range(40,359))
        elif pdbname.upper() in ['6H7J','6H7O']:
            pdb_range = list(range(40,358))
        elif pdbname.upper() in ['2YDO']:
            pdb_range = list(range(6,214))+list(range(224,325))
        elif pdbname.upper() in ['3QAK']:
            pdb_range = list(range(3,149))+list(range(158,209))+list(range(222,309))
        elif pdbname.upper() in ['4BVN']:
            pdb_range = list(range(36, 241))+list(range(275,359))
        elif pdbname.upper() in ['4NTJ']:
            pdb_range = list(range(16, 88))+list(range(92,133))+list(range(136,163))+list(range(179,224))+list(range(231,313))
        elif pdbname.upper() in ['4O9R']:
            pdb_range = list(range(192,345))+list(range(356,434))+list(range(441,497))+list(range(505,552))
        elif pdbname.upper() in ['4UHR']:
            pdb_range = list(range(6, 155))+list(range(158,263))+list(range(264,320))
        elif pdbname.upper() in ['4YAY']:
            pdb_range = list(range(12,173))+list(range(177,186))+list(range(190,225))+list(range(235,318))
        elif pdbname.upper() in ['4ZUD']:
            pdb_range = list(range(12,134))+list(range(141,186))+list(range(189,223))+list(range(235,305))
        elif pdbname.upper() in ['5A8E']:
            pdb_range = list(range(36,241))+list(range(275,355))
        elif pdbname.upper() in ['5GLH']:
            pdb_range = list(range(88,130))+list(range(135,207))+list(range(217,304))+list(range(311,402))
        elif pdbname.upper() in ['5VEW','5VEX']:
            pdb_range = list(range(136,204))+list(range(218,258))+list(range(261,373))+list(range(380,423))
        elif pdbname.upper()=='5W0P':
            pdb_range = list(range(1,325))
        elif pdbname.upper()=='6HLO':
            pdb_range = list(range(28,227))+list(range(238,279))+list(range(282,325))
        elif pdbname.upper()=='5UZ7':
            pdb_range = list(range(136,206))+list(range(213,332))+list(range(337,361))+list(range(366,419))
        elif pdbname.upper()=='3OE8':
            pdb_range = list(range(35,230))+list(range(232,306))
        elif pdbname.upper()=='3OE9':
            pdb_range = list(range(35,229))+list(range(236,272))+list(range(274,304))
        elif pdbname.upper()=='3V2W':
            pdb_range = list(range(47,149))+list(range(156,232))+list(range(245,326))
        elif pdbname.upper()=='3V2Y':
            pdb_range = list(range(16,149))+list(range(156,232))+list(range(245,331))
        elif pdbname.upper()=='6G79':
            pdb_range = list(range(45,188))+list(range(197,241))+list(range(305,339))+list(range(345,386))

    # Uncertain about exact cut -- pdb/article do not compliment eachother.
    if pdbname.upper()=='4XEE' or pdbname.upper()=='4XES':
        pdb_range = list(range(43,269))+list(range(297,385))
    if pdb_range:
        dbref_found = True
        for pos in list(range(1,len(d['wt_seq'])+1)):
            # print('check for ',pos)
            if pos in list(set(pdb_range)):
                # print(pos,'there')
                pos_in_wt.remove(int(pos))

        # for pos in list(set(pdb_range)):
        #     pos_in_wt.remove(int(pos))
    else:
        print('NO DB REF!!')
        dbref_found = False
    # print("pos_in_wt",pos_in_wt)

    # To prevent the otherwise overrides from faulty SIFTS
    chain_over_ride = None
    if pdbname=='3SN6':
        chain_over_ride = 'D'
    elif pdbname=='5UZ7':
        chain_over_ride = 'E'
    # elif pdbname=='5VBL': # fixed by pdb_range
    #     chain_over_ride = 'B'
    # print(pdb_range)
    #http://files.gpcrdb.org/uniprot_mapping.txt
    ## get uniprot to name mapping
    uniprot_mapping = cache.get('gpcrdb_uniprot_mapping')
    if not uniprot_mapping:
        url = 'http://files.gpcrdb.org/uniprot_mapping.txt'
        req = urlopen(url)
        uniprot_mapping = req.read().decode('UTF-8')
        rows = ( line.split(' ') for line in uniprot_mapping.split('\n') )
        uniprot_mapping = { row[0]:row[1:] for row in rows }
        cache.set('gpcrdb_uniprot_mapping',uniprot_mapping,60*60*24)

    #errors, fix it.
    uniprot_mapping['P08483'] = ['acm3_rat']
    uniprot_mapping['P42866'] = ['oprm_mouse']

    variants_mapping = {}
    cache_dir = ['sifts', 'xml']
    url = 'http://www.uniprot.org/uniprot/$index.xml'
    insert_info = fetch_from_web_api(url, d['construct_crystal']['uniprot'], cache_dir, xml = True)
    for elm in insert_info.findall('.//{http://uniprot.org/uniprot}feature'):
        if elm.attrib['type']=="sequence variant":
            if 'description' in elm.attrib:
                desc = elm.attrib['description']
            else:
                desc = ''
            if 'id' in elm.attrib:
                var_id = elm.attrib['id']
            else :
                var_id = None
           #  print(desc,var_id)
            try:
                ori = elm.find('{http://uniprot.org/uniprot}original').text
                var = elm.find('{http://uniprot.org/uniprot}variation').text
                pos = elm.find('{http://uniprot.org/uniprot}location')[0].attrib['position']
                if pos not in variants_mapping:
                    variants_mapping[pos] = {}
                if var not in variants_mapping[pos]:
                    variants_mapping[pos][var] = []
                variants_mapping[pos][var].append([desc,var_id])
            except:
                pass
    for elm in insert_info.findall('.//{http://uniprot.org/uniprot}sequence'):
        uniprot_seq = elm.text
        if uniprot_seq:
            import re
            uniprot_seq = re.sub('[\s+]', '', uniprot_seq)
            # print(uniprot_seq)
    # print(variants_mapping)
    # if len(uniprot_seq)!=len(d['wt_seq']): print("gpcrdb seq",len(d['wt_seq']),'uniport len',len(uniprot_seq))

    # Parsing

    #ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/1xyz.xml.gz
    # Alternative : url = 'https://www.rcsb.org/pdb/files/$index.sifts.xml'
    url = 'ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/$index.xml.gz'
    sifts = fetch_from_web_api(url, pdbname.lower(), cache_dir, xml = True)
    d['links'].append(Template(url).substitute(index=quote(str(pdbname.lower()), safe='')))
    d['mutations'] = []
    d['auxiliary'] = OrderedDict()
    d['construct_sequences'] = OrderedDict()
    receptor_seq_ids = []
    other_seq_ids = {}
    receptor_chain = ''
    uniprot_to_name = {}
    mutations_check = []
    if sifts: #success
        insert_position = 'N-term'
        insert_start = 0
        msg_1 = 0
        msg_2 = 0
        sifts_https = False
        for elem in sifts:
            if "https" in elem:
                sifts_https = True
        if sifts_https == True:
            sfits_https = "https"
        else:
            sfits_https = "http"

        for elem in sifts.findall('.//{'+sfits_https+'://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd}segment'):
            receptor = False
            chain = elem.attrib['segId'].split('_')[1]
            for res in elem[0]: #first element is residuelist
                if receptor_chain!='':
                    break #break if found
                for node in res:
                    if node.tag == '{'+sfits_https+'://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd}crossRefDb':
                        source = node.attrib['dbSource']
                        if source=='UniProt':
                            u_id = node.attrib['dbAccessionId']
                            if u_id in uniprot_mapping:
                                receptor_chain = chain
                                break

        prev_elem_name = ""
        pdb_resid_total = []
        pdb_resid_total_accounted = []
        # print(d['wt_seq'])
        for elem in sifts.findall('.//{'+sfits_https+'://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd}segment'):
            if 'segId' not in elem.attrib:
                continue #not receptor
            if elem.attrib['segId']==prev_elem_name:
                # pass
                # print("skip ",elem.attrib['segId'])
                continue
            prev_elem_name = elem.attrib['segId']
            seg_uniprot_ids = []
            max_pos = 0
            min_pos = 99999
            pos_list = []
            uniprot_pos = None
            pos = None
            receptor = False
            u_id_source = 'N/A'
            chain = elem.attrib['segId'].split('_')[1]
            seg_resid_list = []
            elem_seq = ""
            prev_raw_u_id = ""
            raw_u_id = ""
            prev_pos = 0
            prev_receptor = False
            seg_had_receptor = False

            if (chain=="A" or chain=="B") and pdbname.lower()=="4k5y":
                continue
            # print(chain,'chain')
            for res in elem[0]: #first element is residuelist
                u_id = 'N/A'
                pdb_aa = ''
                uniprot_pos = None
                pos = None
                for node in res:
                    if raw_u_id!=prev_raw_u_id:
                        # print("New u_id",raw_u_id,u_id)
                        prev_raw_u_id = raw_u_id
                    if node.tag == '{'+sfits_https+'://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd}crossRefDb':
                        source = node.attrib['dbSource']
                        if source=='UniProt':
                            u_id = node.attrib['dbAccessionId']
                            raw_u_id = node.attrib['dbAccessionId']
                            u_id_source = 'UniProt'
                            if u_id in uniprot_mapping:
                                u_id = uniprot_mapping[u_id][0]
                                receptor = True ## this is receptor element
                                seg_had_receptor = True
                                if receptor_chain=='' or receptor_chain==chain:
                                    receptor_chain = chain
                                elif msg_1==0:
                                    msg_1 = 1
                                    # print('\t', pdbname.lower(),'receptor in many chains?!',chain,receptor_chain)
                                    logger.warning('{} has receptor in many chains {} {}'.format(pdbname.lower(),chain,receptor_chain))
                                insert_position = 'Within Receptor'
                            else:
                                receptor = False
                                if raw_u_id in uniprot_to_name:
                                    u_id = uniprot_to_name[raw_u_id]

                                else:
                                    url = 'http://www.uniprot.org/uniprot/$index.xml'
                                    insert_info = fetch_from_web_api(url, raw_u_id, cache_dir, xml = True)
                                    found_u_id = None
                                    for elm in insert_info.findall('.//{http://uniprot.org/uniprot}feature'):
                                        # GRAB NON RECEPTOR VARIANTS
                                        try:
                                            if elm.attrib['type']=="sequence variant":
                                                desc = elm.attrib['description']
                                                ori = elm.find('{http://uniprot.org/uniprot}original').text
                                                var = elm.find('{http://uniprot.org/uniprot}variation').text
                                                pos = elm.find('{http://uniprot.org/uniprot}location')[0].attrib['position']
                                                # print(raw_u_id,desc,ori,var,pos)
                                        except:
                                            pass

                                    for elm in insert_info.findall('.//{http://uniprot.org/uniprot}recommendedName'):
                                        new_u_id = elm.find('{http://uniprot.org/uniprot}fullName').text
                                        uniprot_to_name[raw_u_id] = new_u_id
                                        u_id = new_u_id
                                        found_u_id = True
                                        break #no need to continue looking
                                    if not found_u_id:
                                        for elm in insert_info.findall('.//{http://uniprot.org/uniprot}submittedName'):
                                            new_u_id = elm.find('{http://uniprot.org/uniprot}fullName').text
                                            uniprot_to_name[raw_u_id] = new_u_id
                                            u_id = new_u_id
                                            found_u_id = True
                                            break #no need to continue looking

                            if u_id not in seg_uniprot_ids:
                                seg_uniprot_ids.append(u_id)
                            uniprot_pos = int(node.attrib['dbResNum'])
                            uniprot_aa = node.attrib['dbResName']
                            if u_id=='crfr1_human':
                                # fix for crfr1_human isoform blah
                                if uniprot_pos>145:
                                    # print('crfr1_human',uniprot_pos,pos,pdb_aa)
                                    uniprot_pos = pos
                            if pdbname == '5UZ7':
                                #Special fix for 5UZ7 due to faulty annotation, there is an offset of 34 at the end of the isoforms
                                if pos:
                                    # print(uniprot_pos,d['wt_seq'][uniprot_pos-1-34],uniprot_aa,pos)
                                    if uniprot_pos>452 and uniprot_aa==d['wt_seq'][uniprot_pos-1-34]:
                                        #these are unmapped at this point, make sure the AA are the same, in case it gets fixed later
                                        uniprot_pos = uniprot_pos-34
                                        pos = uniprot_pos
                                        # Assume we are talking about the non-observed residues
                                        if uniprot_pos not in d['xml_not_observed']:
                                            d['xml_not_observed'].append(uniprot_pos)
                                    else:
                                        uniprot_pos = int(pos)
                                else:
                                    receptor = False
                            if pdbname == '6NIY':
                                #Special fix for 5UZ7 due to faulty annotation, there is an offset of 34 at the end of the isoforms
                                if pos and receptor:
                                    # print(uniprot_pos,d['wt_seq'][uniprot_pos-1-34],uniprot_aa,pos,d['wt_seq'][uniprot_pos-1-18], receptor)
                                    if uniprot_aa==d['wt_seq'][uniprot_pos-1-18]:
                                        uniprot_pos = uniprot_pos-18
                                        pos = uniprot_pos
                                        # print('match 1')
                                        # print('found',pos, uniprot_aa)
                                        # # Assume we are talking about the non-observed residues
                                        # if uniprot_pos not in d['xml_not_observed']:
                                        #     d['xml_not_observed'].append(uniprot_pos)
                                    elif uniprot_pos>192 and uniprot_aa==d['wt_seq'][uniprot_pos-1-34]:
                                        #these are unmapped at this point, make sure the AA are the same, in case it gets fixed later
                                        # print('match 2')
                                        uniprot_pos = uniprot_pos-34
                                        pos = uniprot_pos
                                        # print('found',pos, uniprot_aa)
                                        # Assume we are talking about the non-observed residues
                                        # if uniprot_pos not in d['xml_not_observed']:
                                        #     d['xml_not_observed'].append(uniprot_pos)
                                    else:
                                        #print('no match')
                                        uniprot_pos = int(pos)
                                else:
                                    receptor = False
                            # if pdbname in ['6KUX','6KUY']:
                            #     # Special fix for shift in annotation
                            #     if pos:
                            #         uniprot_pos = uniprot_pos-15

                            # if receptor:
                            #     print(receptor, uniprot_pos, pos,uniprot_aa, u_id,raw_u_id,chain,node.attrib['dbResNum'],d['wt_seq'][uniprot_pos-1])
                        elif source=='PDB' and node.attrib['dbResNum'].lstrip('-').isdigit(): #use instead of isinstance(node.attrib['dbResNum'], int):
                            pos = int(node.attrib['dbResNum'])
                            # print("PDB pos",pos)
                            try:
                                pdb_aa = AA_three[node.attrib['dbResName'].upper()]
                            except:
                                pdb_aa = "X"
                            if receptor:
                                receptor_seq_ids.append(pos)
                            elem_seq += pdb_aa
                            # print(pos,pdb_aa)
                            seg_resid_list.append(pos)
                            if pos not in pdb_resid_total:
                                pdb_resid_total.append(pos)
                            if pos>max_pos: max_pos = pos
                            if pos<min_pos: min_pos = pos
                        elif source=='PDB':
                            #if above fails, it's probably cos its missing residue, still should capture pdb_aa
                            if node.attrib['dbResName'].upper() in AA_three:
                                pdb_aa = AA_three[node.attrib['dbResName'].upper()]
                            else:
                                pdb_aa ="X"

                    elif pdb_aa and node.tag == '{'+sfits_https+'://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd}residueDetail':
                        # print(node.text)
                        if u_id!='N/A' and u_id not in d['construct_sequences']:
                            d['construct_sequences'][u_id] = OrderedDict()
                            d['construct_sequences'][u_id]['seq'] = ''
                            d['construct_sequences'][u_id]['pos'] = []
                            d['construct_sequences'][u_id]['where'] = insert_position
                        if u_id!='N/A' and uniprot_pos and receptor_chain==chain:
                            if uniprot_pos not in d['construct_sequences'][u_id]['pos']:
                                # d['construct_sequences'][u_id][uniprot_pos] = [uniprot_aa ,pdb_aa]
                                d['construct_sequences'][u_id]['seq'] += pdb_aa
                                d['construct_sequences'][u_id]['pos'].append(uniprot_pos)
                            if node.text == 'Engineered mutation':
                                # print(node.attrib['property'],node.text)
                                # print(u_id)
                                # print(uniprot_pos)
                                # print(pdb_aa)
                                # print(node.text)
                                if 'mutations' not in d['construct_sequences'][u_id]:
                                    d['construct_sequences'][u_id]['mutations'] = OrderedDict()
                                d['construct_sequences'][u_id]['mutations'][uniprot_pos] = [uniprot_aa,pdb_aa]

                        if node.text=='Not_Observed' and receptor:
                            # print('not observed!',elem.attrib['segId'],uniprot_pos,pos,pdb_aa)
                            if uniprot_pos:
                                if uniprot_pos not in d['xml_not_observed']:
                                    d['xml_not_observed'].append(uniprot_pos)
                                if not pos:
                                    #if no pos but uniprot looks like pos_in_wt
                                    if uniprot_pos in pdb_range:
                                        pos = uniprot_pos
                            elif pos: #in rare cases a uniprot_pos is not captured, but it is safe to assume that pos then must be correct since we are in receptor
                                if pos not in d['xml_not_observed']:
                                    d['xml_not_observed'].append(pos)

                                    if receptor and pos in pos_in_wt:
                                        #make sure to not get this pos 'deleted'
                                        uniprot_pos = pos
                                        pos_in_wt.remove(pos)
                                        insert_start =  str(pos+1)

                        elif node.attrib['property']=='Annotation' and u_id=='N/A':
                            u_id = node.text
                            if u_id not in seg_uniprot_ids:
                                seg_uniprot_ids.append(u_id)
                        elif receptor and node.attrib['property']=='Annotation' and node.text == 'Engineered mutation': ## only in receptor
                            if (uniprot_pos not in pos_in_wt and dbref_found) or not dbref_found:
                                if {'mut':pdb_aa,'wt':uniprot_aa,'pos':uniprot_pos,'type':node.text} not in d['mutations']: #prevent duplicates
                                    d['mutations'].append({'mut':pdb_aa,'wt':uniprot_aa,'pos':uniprot_pos,'type':node.text})
                                    mutations_check.append(uniprot_pos)
                                    # print({'mut':pdb_aa,'wt':uniprot_aa,'pos':uniprot_pos,'type':node.text},chain,u_id,max_pos)
                        elif receptor and node.attrib['property']=='Annotation' and node.text == 'Conflict': ## only in receptor
                            # if (uniprot_pos not in pos_in_wt and dbref_found) or not dbref_found:
                            #     if {'mut':pdb_aa,'wt':uniprot_aa,'pos':uniprot_pos,'type':node.text} not in d['mutations']: #prevent duplicates
                            #         d['mutations'].append({'mut':pdb_aa,'wt':uniprot_aa,'pos':uniprot_pos,'type':node.text})
                            #         mutations_check.append(uniprot_pos)
                            #         print({'mut':pdb_aa,'wt':uniprot_aa,'pos':uniprot_pos,'type':node.text},chain,u_id,max_pos)
                            # print('conflict!','mut',pdb_aa,'wt',uniprot_aa,'pos',uniprot_pos)
                            pass


                # print(chain,uniprot_pos,pos,receptor)

                if u_id!='N/A':
                    if pos not in pdb_resid_total_accounted:
                        pdb_resid_total_accounted.append(pos)

                if not chain_over_ride:
                    if not receptor and pos:
                        if pos in pos_in_wt and uniprot_seq[pos-1]==pdb_aa:
                            # if uniprot receptor aa hasnt been found already and aa matches
                            if int(pos)<prev_pos:
                                # print('suddenly lower pos number, maybe receptor again?',pos)
                                receptor = True
                                uniprot_pos = pos
                            elif prev_receptor and int(pos-1)==prev_pos:
                                # print('last was receptor and this number is only one higher?',pos)
                                receptor = True
                                uniprot_pos = pos
                        else:
                            pass
                            # print(pos,'not receptor!',pdb_aa)
                    elif not pos:
                        receptor = False

                    if pos:
                        prev_pos = int(pos)
                    if receptor:
                            if not uniprot_pos:
                                uniprot_pos = pos
                            # print(chain,pos,uniprot_pos,uniprot_aa)
                            wt_aa = d['wt_seq'][uniprot_pos-1]
                            prev_receptor = True
                            # if pos==250 or uniprot_pos==250:
                            #     print(pos,uniprot_pos,pdb_aa,d['wt_seq'][uniprot_pos-1],d['wt_seq'][pos-1])
                    # if receptor and uniprot_pos==None :
                    #     if pos<len(d['wt_seq']):
                    #         print("receptor but no uniprot pos?",pos,pdb_aa,u_id)
                    #         print("WT AA ",d['wt_seq'][pos-1])
                    #         wt_aa = d['wt_seq'][pos-1]
                    #         if uniprot_pos in pos_in_wt:
                    #             pos_in_wt.remove(pos)
                    #             insert_start =  str(pos+1)
                    #         pos_list.append(pos)
                    #         if wt_aa!=pdb_aa:
                    #             # mutation!
                    #             if {'mut':pdb_aa,'wt':wt_aa,'pos':pos,'type':'custom_maybe_wrong'} not in d['mutations']: #prevent duplicates
                    #                 d['mutations'].append({'mut':pdb_aa,'wt':wt_aa,'pos':pos,'type':'custom_maybe_wrong'})

                    if not receptor and dbref_found:
                        # if not receptor and we have dbref_found
                        if pos in pdb_range:
                            # if the pdb pos is in the pdb_range then assume its sifts error
                            receptor = True
                            uniprot_pos = pos
                            # print('fixed',pos)

                if uniprot_pos and u_id=='crfr1_human':
                    # fix for crfr1_human isoform blah
                    if uniprot_pos>145:
                        # print('crfr1_human',uniprot_pos,pos,pdb_aa)
                        uniprot_pos = pos

                # print(chain,uniprot_pos,pos,receptor)
                # if receptor and uniprot_pos and pos:
                #     print(receptor, uniprot_pos, pos,uniprot_aa, pdb_aa,d['wt_seq'][uniprot_pos-1],d['wt_seq'][pos-1])
                if uniprot_pos:
                    if pos:
                        pos = int(pos)
                    uniprot_pos = int(uniprot_pos)
                    if receptor and uniprot_pos not in pos_in_wt and raw_u_id==uniprot_code: #FIXME, consider undoing raw_u_id since it doesnt capture anything
                            # print('hi',pos,uniprot_pos)
                        # if pos and pos>1000 and uniprot_pos<1000 and uniprot_pos!=pos-1000 and uniprot_pos!=pos-2000:
                        #     # print('PDB residue number:',pos, 'Receptor Uniprot pos:',uniprot_pos)
                        #     pass
                        #     # not wt likely
                        # else:
                            if uniprot_pos<len(d['wt_seq']):
                                wt_aa = d['wt_seq'][uniprot_pos-1]
                                # print(pdb_aa,wt_aa)
                                if wt_aa!=pdb_aa and pdb_aa:
                                    # mutation!
                                    # print("MUTATION",u_id, uniprot_pos,pos ,uniprot_aa,"|",pdb_aa,"|",wt_aa)
                                    if uniprot_pos not in mutations_check: #prevent duplicates
                                        mut_type = "non_annotated_mutation"
                                        if str(uniprot_pos) in variants_mapping:
                                            mut_type = 'SNP location (Not this AA)'
                                            if pdb_aa in variants_mapping[str(uniprot_pos)]:
                                                mut_type = 'SNP location: '+variants_mapping[str(uniprot_pos)][pdb_aa][0][0] + variants_mapping[str(uniprot_pos)][pdb_aa][0][1]

                                        # if  uniprot_seq[uniprot_pos-1]==pdb_aa:
                                        #     mut_type = "ISOFOR MISMATCH"
                                        if pos and pos>1000 and 1==2:
                                            ## this is probably not a real mutation but an annotation error, ignore
                                            pass
                                        else:
                                            # print("MUT NOT ANNOTATED WT AA ",wt_aa,pdb_aa,uniprot_pos,pos,mut_type,u_id,'uni',uniprot_seq[uniprot_pos-1])
                                            d['mutations'].append({'mut':pdb_aa,'wt':wt_aa,'pos':uniprot_pos,'type':mut_type})
                                            mutations_check.append(uniprot_pos)
                                elif not pdb_aa:
                                    #no pdb_aa seen to be missing
                                    print(uniprot_pos,' not found')
                                    if uniprot_pos not in d['xml_not_observed']:
                                        d['xml_not_observed'].append(uniprot_pos)
                            if uniprot_pos in pdb_range:
                                pdb_range.remove(uniprot_pos)
                            insert_start =  str(uniprot_pos+1)
                    elif receptor:
                        # print('wierd error with position already deleted',uniprot_pos)
                        pass
                    else:
                        if u_id not in other_seq_ids:
                            other_seq_ids[u_id] = []
                        if pos not in other_seq_ids[u_id]:
                            other_seq_ids[u_id].append(pos)
                elif pos:
                    #if this segment doesnt have a uniprot equivilant (linkers/tags)
                    if u_id not in d['construct_sequences']:
                            d['construct_sequences'][u_id] = OrderedDict()
                            d['construct_sequences'][u_id]['seq'] = ''
                            d['construct_sequences'][u_id]['pos'] = []
                            d['construct_sequences'][u_id]['where'] = insert_position
                    if pos not in d['construct_sequences'][u_id]['pos']:
                        # d['construct_sequences'][u_id][uniprot_pos] = [uniprot_aa ,pdb_aa]
                        d['construct_sequences'][u_id]['seq'] += pdb_aa
                        d['construct_sequences'][u_id]['pos'].append(pos)



            ranges = []
            for k, g in groupby(enumerate(pos_list), lambda x:x[0]-x[1]):
                group = list(map(itemgetter(1), g))
                ranges.append((group[0], group[-1]))
                if (group[0], group[-1]) not in ranges:
                    ranges.append((group[0], group[-1]))

            if pdbname in SIFT_exceptions:
                try:
                    if ranges[0][1]==SIFT_exceptions[pdbname][0]:
                        ranges=[(ranges[0][0],SIFT_exceptions[pdbname][1])]
                        actually_present = list(range(max_pos+1,ranges[0][1]+1))
                        seg_resid_list+=actually_present
                        max_pos=SIFT_exceptions[pdbname][1]
                except:
                    pass
            if len(seg_uniprot_ids)>0 and seg_uniprot_ids[0]=='Not_Observed' and 'actually_present' in locals():
                if min_pos==SIFT_exceptions[pdbname][0]+1:
                    min_pos=SIFT_exceptions[pdbname][1]+1
                seg_resid_list = [i for i in seg_resid_list if i not in actually_present]
                pos_in_wt = [i for i in pos_in_wt if i not in actually_present]

            mutations = None

            if receptor==False and u_id_source=='UniProt':
                if seg_uniprot_ids[0] in d['construct_sequences']:
                    if 'mutations' in d['construct_sequences'][seg_uniprot_ids[0]]:
                        mutations = d['construct_sequences'][seg_uniprot_ids[0]]['mutations']

                if insert_info!=False:
                    for elm in insert_info.findall('.//{http://uniprot.org/uniprot}recommendedName'):
                        seg_uniprot_ids[0] = elm.find('{http://uniprot.org/uniprot}fullName').text
            # Custom annotation fixes
            if pdbname in ['6D32','6D35'] and min_pos==416 and max_pos==525:
                seg_uniprot_ids = ['Uncharacterized protein']
            # elif pdbname in ['4LDE'] and min_pos==1029 and max_pos==1342:
            #     seg_uniprot_ids = ['adrb2_human']
            if pdbname=='7BZ2' and chain=='E':
                seg_uniprot_ids = ['adrb2_human']
                receptor = [{'start': 30, 'end': 340, 'origin': 'user'}]
            if pdbname=='6IQL' and chain in ['A','B'] and min_pos==304:
                seg_uniprot_ids = ['drd4_mouse']

            # print([elem.attrib['segId'],seg_uniprot_ids,min_pos,max_pos,ranges,insert_position,seg_resid_list,mutations,seg_had_receptor])
            d['xml_segments'].append([elem.attrib['segId'],seg_uniprot_ids,min_pos,max_pos,ranges,insert_position,seg_resid_list,mutations,seg_had_receptor])

            # print("end of segment",elem.attrib['segId'],seg_uniprot_ids,max_pos)
            if [elem.attrib['segId'],seg_uniprot_ids,min_pos,max_pos,ranges,insert_position,seg_resid_list,mutations,seg_had_receptor] not in d['xml_segments']:
                d['xml_segments'].append([elem.attrib['segId'],seg_uniprot_ids,min_pos,max_pos,ranges,insert_position,seg_resid_list,mutations,seg_had_receptor])

            if receptor == False and receptor_chain==chain: #not receptor, but is in same chain
                if len(seg_uniprot_ids):
                    subtype =seg_uniprot_ids[0]
                else:
                    subtype ='N/A'
                    continue #do not add segments without information
                if subtype == 'Not_Observed':
                    continue #ignore "aux" that are 'not observed'
                if subtype == 'Engineered mutation':
                    continue #ignore "aux" that are 'not observed'
                if subtype == 'S-arrestin':
                    continue #  S-arrestin is not part of the chain
                # print(subtype)
                seq = ''
                if len(seg_uniprot_ids)==1:
                    seq = elem_seq
                d['auxiliary']['aux'+str(len(d['auxiliary']))] = {'type':'auto','subtype':subtype,'presence':'YES','position':insert_position, 'start':insert_start, 'sequence':seq}
            elif receptor == False:
                # print('\t',pdbname.lower(),'Protein in PDB, not part of receptor chain',seg_uniprot_ids,'chain',chain)
                logger.warning('{} Protein in structure, but not part of receptor chain {} {}'.format(pdbname.lower(),seg_uniprot_ids,chain))

        # print(sorted(pdb_resid_total))
        # print(sorted(pdb_resid_total_accounted))
        non_accounted = sorted(list(set(pdb_resid_total) - set(pdb_resid_total_accounted)))
        d['non_accounted'] = []
        for k, g in groupby(enumerate(non_accounted), lambda x:x[0]-x[1]):
            group = list(map(itemgetter(1), g))
            d['non_accounted'].append((group[0], group[-1]))
        # if len(d['non_accounted']):
        #     print("non_accounted",d['non_accounted'])
        d['deletions'] = []
        for k, g in groupby(enumerate(pos_in_wt), lambda x:x[0]-x[1]):
            group = list(map(itemgetter(1), g))
            d['deletions'].append({'start':group[0], 'end':group[-1], 'origin':'user'})
        d['not_observed'] = []
        if len(d['xml_not_observed']):
            # print(d['xml_not_observed'])
            for k, g in groupby(enumerate(sorted(d['xml_not_observed'])), lambda x:x[0]-x[1]):
                group = list(map(itemgetter(1), g))
                d['not_observed'].append((group[0], group[-1]))

        for i,v in d['construct_sequences'].items():
            d['construct_sequences'][i]['ranges'] = []
            for k, g in groupby(enumerate(v['pos']), lambda x:x[0]-x[1]):
                group = list(map(itemgetter(1), g))
                d['construct_sequences'][i]['ranges'].append((group[0], group[-1]))

            d['construct_sequences'][i]['pos'] = ''
    else:
        pass
        # print('failed sifts')


    for aux,v in d['auxiliary'].items():
        if int(v['start'])>max(receptor_seq_ids):
            v['position'] = "C-term"
        # print(max(receptor_seq_ids))

    #https://www.ebi.ac.uk/pdbe/api/pdb/entry/experiment/2RH1
    ## experiment data
    cache_dir = ['pdbe', 'experiment']
    url = 'https://www.ebi.ac.uk/pdbe/api/pdb/entry/experiment/$index'
    pdbe = fetch_from_web_api(url, pdbname, cache_dir)
    d['links'].append(Template(url).substitute(index=quote(str(pdbname), safe='')))
    if pdbe: #success
        r = pdbe[pdbname.lower()][0]
        d['resolution'] = r.get('resolution')
        d['crystal_growth'] = r.get('crystal_growth')
        d['r_factor'] = r.get('r_factor')
        d['experimental_method'] = r.get('experimental_method')
    else:
        pass
        # print('failed pdbe')

    # #https://www.ebi.ac.uk/pdbe/api/pdb/entry/modified_AA_or_NA/2RH1
    # ## modified AA (empty on 2RH1)
    # cache_dir = ['pdbe', 'modified_AA_or_NA']
    # url = 'https://www.ebi.ac.uk/pdbe/api/pdb/entry/modified_AA_or_NA/$index'
    # pdbe_mod = fetch_from_web_api(url, pdbname, cache_dir)
    # d['links'].append(Template(url).substitute(index=quote(str(pdbname), safe='')))

    # if pdbe_mod: #success
    #     print(pdbe_mod)
    # else:
    #     d['modifications3'] = 'None'
    #     print('failed pdbe_mod')

    #http://www.rcsb.org/pdb/explore/jmol.do?structureId=4LDO&json=true
    ## modifications for their jmol -- "hacky" way to get it
    try:
        cache_dir = ['rcsb', 'jmol_modifications']
        url = 'https://www.rcsb.org/pdb/explore/jmol.do?structureId=$index&json=true'
        rcsb_mod = fetch_from_web_api(url, pdbname, cache_dir)
        d['links'].append(Template(url).substitute(index=quote(str(pdbname), safe='')))
        # print(Template(url).substitute(index=quote(str(pdbname), safe='')))
    except:
        print('rscb failed for ',pdbname)
        rcsb_mod = None
    if rcsb_mod: #success
        d['modifications'] = []
        d['modifications2'] = rcsb_mod
        # print(receptor_seq_ids)
        for mod in rcsb_mod['protmod']['domains']:
            t = mod['range'].split(',')
            if t[0].split(':')[1]!=receptor_chain:
                # print('modification not in receptor chain, not interested')
                continue
            if len(t)>1:
                position_type = 'pair'
                position_info = [t[0].split(':')[0],t[1].split(':')[0]]
            elif len(t)==1:
                position_type = 'single'
                position_info = [t[0].split(':')[0],0]
            else:
                print('error',t)
                continue
            # print(mod['id'],pair,mod['description'])
            if mod['id']=='crosslink2': mod['id']="Disulfide bond" #replace non-descript crosslink2
            d['modifications'].append({'position':[position_type,position_info],'type':mod['id'],'remark':mod['description']})
            #{{v.id}} {{v.range}} {{v.description}} {{v.pdbCcId}} <br><br>

    else:
        d['modifications2'] = 'None'
        # print('failed pdbe_mod')

    #https://www.ebi.ac.uk/pdbe/api/pdb/entry/ligand_monomers/2RH1
    cache_dir = ['pdbe', 'ligands']
    url = 'https://www.ebi.ac.uk/pdbe/api/pdb/entry/ligand_monomers/$index'
    pdbe_ligands = fetch_from_web_api(url, pdbname, cache_dir)
    d['links'].append(Template(url).substitute(index=quote(str(pdbname), safe='')))
    # print(Template(url).substitute(index=quote(str(pdbname), safe='')))
    if pdbe_ligands: #success
        d['ligands'] = {}
        for name,pdb in pdbe_ligands.items():
            for ligand in pdb:
                if ligand['chem_comp_id'] not in d['ligands']:
                    d['ligands'][ligand['chem_comp_id']] = {'comp_name':ligand['chem_comp_name'], 'number_of_entries':1}
                else:
                    d['ligands'][ligand['chem_comp_id']]['number_of_entries'] += 1
        # print(d['ligands'])

    else:
        d['ligands'] = 'None'
        # print('failed pdbe_ligands')


    ## NOT NEED - FETCH MUT FROM XML

    # #https://www.ebi.ac.uk/pdbe/api/pdb/entry/mutated_AA_or_NA/2RH1
    # ## mutated AA
    # ### got conflicts, engerineered mutation and expression tag examples
    # cache_dir = ['pdbe', 'mutated_AA_or_NA']
    # url = 'https://www.ebi.ac.uk/pdbe/api/pdb/entry/mutated_AA_or_NA/$index'
    # pdbe_mut = fetch_from_web_api(url, pdbname, cache_dir)
    # d['links'].append(Template(url).substitute(index=quote(str(pdbname), safe='')))

    # if pdbe_mut: #success
    #     r = pdbe_mut[pdbname.lower()]
    #     d['mutations_pdbe'] = []
    #     for mut in r:
    #         mut_from = mut['mutation_details']['from']
    #         mut_to = mut['mutation_details']['to']
    #         mut_type = mut['mutation_details']['type']
    #         construct_seq_number = mut['residue_number']
    #         wt_seq_number = mut['author_residue_number']
    #         t = {'wt':mut_from,'mut':mut_to,'type':mut_type,'c_seq_nr':construct_seq_number,'pos':wt_seq_number}
    #         d['mutations_pdbe'].append(t)
    # else:
    #     print('failed pdbe_mut')


    #https://www.rcsb.org/pdb/rest/das/pdb_uniprot_mapping/alignment?query=2RH1
    ## uniprot mappings
    ### seems to be IDs of stuff then use:
    # http://www.uniprot.org/uniprot/P00720.xml
    cache_dir = ['rcsb', 'pdb_uniprot_mapping']
    url = 'https://www.rcsb.org/pdb/rest/das/pdb_uniprot_mapping/alignment?query=$index'
    uniprot_map = fetch_from_web_api(url, pdbname, cache_dir, xml = True)
    d['links'].append(Template(url).substitute(index=quote(str(pdbname), safe='')))

    if uniprot_map: #success
        inserts = {}
        inserts_fixed = {}
        for block in uniprot_map[0]:
            if block.tag[-5:]!='block':
                continue #only interested in the blocks...
            i = 0
            for segment in block:
                if i==0:
                    construct_range = [segment.attrib['start'],segment.attrib['end']]
                else:
                    insert_range = [segment.attrib['start'],segment.attrib['end']]
                    insert_id = segment.attrib['intObjectId']
                prev_block = segment
                i += 1
            i = inserts.setdefault(insert_id, [])
            i.append({'c':construct_range,'i':insert_range})
        for insert,blocks in inserts.items():

            if insert in uniprot_mapping:
                insert = uniprot_mapping[insert][0]

            inserts_fixed[insert] = {}
            cache_dir = ['uniprot', 'id']
            url = 'http://www.uniprot.org/uniprot/$index.xml'
            insert_info = fetch_from_web_api(url, insert, cache_dir, xml = True)
            d['links'].append(Template(url).substitute(index=quote(str(insert), safe='')))

            if insert_info:
                for elm in insert_info.findall('.//{http://uniprot.org/uniprot}recommendedName'):
                    inserts_fixed[insert]['alt_name'] = elm.find('{http://uniprot.org/uniprot}fullName').text
            else:
                inserts_fixed[insert]['alt_name'] = insert
            # print(insert_info.findall('.//.'))

            blocks_num = len(blocks)
            prev_block = None
            temp = []
            for i, b in enumerate(blocks): #for each block, to glue them together
                if i==0:
                    start = [b['i'][0],b['c'][0]]
                    end = [b['i'][1],b['c'][1]]
                # print(i,b)
                if i<blocks_num-1: #if not last
                    # print('cur',b,'next',blocks[i+1])
                    if int(b['i'][1])==int(blocks[i+1]['i'][0])-1 and int(b['c'][1])==int(blocks[i+1]['c'][0])-1:
                        #if insert is a contination #if construct continues
                        end = [blocks[i+1]['i'][1],blocks[i+1]['c'][1]]
                    else:
                        #gap
                        temp.append({'i_start':start[0],'i_end':end[0],'c_start':start[1],'c_end':end[1]})
                        # temp.append([start,end])
                        start = [blocks[i+1]['i'][0],blocks[i+1]['c'][0]]
                        end = [blocks[i+1]['i'][1],blocks[i+1]['c'][1]]
            temp.append({'i_start':start[0],'i_end':end[0],'c_start':start[1],'c_end':end[1]})
            i = inserts_fixed[insert].setdefault('positions', [])
            i.append(temp)

        d['inserts'] = inserts_fixed


    else:
        pass
        # print('failed uniprot_map')
    return d


def add_construct(d):

    #delete if already name there
    Construct.objects.filter(name = d['construct_crystal']['pdb_name']).delete()

    protein = Protein.objects.filter(entry_name=d['construct_crystal']['uniprot']).get()
    structure = Structure.objects.filter(pdb_code__index=d['construct_crystal']['pdb'].upper()).get()
    protein_conformation = structure.protein_conformation

    construct = Construct()
    construct.protein = protein
    construct.name = d['construct_crystal']['pdb_name'].strip()
    construct.json = json.dumps(d, indent=4, separators=(',', ': '))
    construct.structure = structure

    #CrystalInfo
    crystal = CrystalInfo()
    crystal.resolution = structure.resolution
    crystal.pdb_data = structure.pdb_data
    crystal.pdb_code = structure.pdb_code.index
    crystal.save()

    construct.crystal = crystal

    #Contact INFO
    if 'contact_info' in d:
        construct.contributor, created = ContributorInfo.objects.get_or_create(name = d['contact_info']['name_cont'],
                                                       pi_email = d['contact_info']['pi_email'],
                                                       pi_name = d['contact_info']['pi_name'],
                                                       urls = d['contact_info']['url'],
                                                       date = datetime.datetime.strptime(d['contact_info']['date'], '%m/%d/%Y').strftime('%Y-%m-%d'),
                                                       address = d['contact_info']['address'])

    construct.save()
    #MUTATIONS
    for mutation in d['mutations']:

        if 'type' not in mutation:
            mutation['type'] = ''

        if 'remark' not in mutation:
            mutation['remark'] = ''

        res_wt = Residue.objects.get(protein_conformation__protein=protein_conformation.protein.parent, sequence_number=mutation['pos'])
        # if res_wt.amino_acid != mutation['wt']:
        #     print('aa dont match',construct,mutation['pos'],"annotated wt:", mutation['wt'], "DB wt:",res_wt.amino_acid, "Annotated Mut",mutation['mut'])

        mutation_type, created = ConstructMutationType.objects.get_or_create(slug=slugify(mutation['type']),name=mutation['type'], effect=None)

        #construct=construct, TODO: create a unique one for each mutation per construct to avoid unambiguity
        mut = ConstructMutation.objects.create(construct=construct, sequence_number=mutation['pos'],wild_type_amino_acid=mutation['wt'],mutated_amino_acid=mutation['mut'],remark=mutation['remark'], residue=res_wt)
        mut.effects.add(mutation_type)
        #construct.mutations.add(mut)

    #print(d['raw_data'])
    #make sure to order auxiliary correct
    ip_lookup = {}
    if 'raw_data' in d:
        for name,aux in d['auxiliary'].items():
            id = name.replace('aux','')
            aux['sort'] = 0
            try:
                aux['sort'] = int(d['raw_data']["sort_pos_"+id])
            except:
                pass
        d['auxiliary'] = OrderedDict(sorted(d['auxiliary'].items(), key=lambda x: (x[1]['sort'], x[0])))
        temp = OrderedDict()
        for i, (name,aux) in enumerate(d['auxiliary'].items()):
            old_id = name.replace('aux','')
            temp['aux'+str(i)] = aux
            ip_lookup[old_id] = str(i)
        d['auxiliary'] = temp


    #DELETIONS
    insert_deletions = {}
    for deletion in d['deletions']:
        # if a 'deletion' is a single type and of non-user origin, assume its an insert and the pos is not actually deleted (3odu)
        dele = False
        if 'start' in deletion:
            dele, created = ConstructDeletion.objects.get_or_create(construct=construct, start=deletion['start'],end=deletion['end'])
        else:
            if deletion['origin']=='user':
                dele, created = ConstructDeletion.objects.get_or_create(construct=construct, start=deletion['pos'],end=deletion['pos'])
        # if dele:
        #     construct.deletions.add(dele)
        if deletion['origin']!='user':
            id = deletion['origin'].split('_')[1]
            if id in ip_lookup:
                id = ip_lookup[id]
            insert_deletions[id] = deletion



    #INSERTIONS (AUX)
    for name,aux in d['auxiliary'].items():
        id = name.replace('aux','')
        aux_type, created = ConstructInsertionType.objects.get_or_create(name=aux['type'],subtype=aux['subtype'])
        insert = ConstructInsertion.objects.create(construct=construct, insert_type=aux_type,presence=aux['presence'],position=aux['position']+"_"+id)

        if insert.presence == 'YES' and insert.position.startswith('Within Receptor'):
            #need to fetch range
            if 'start' in aux:
                insert.start = aux['start']
                insert.end = aux['start']
            else:
                if insert_deletions[id]['type'] == 'range':
                    insert.start = insert_deletions[id]['start']
                    insert.end = insert_deletions[id]['end']
                else:
                    insert.start = insert_deletions[id]['pos']
                    insert.end = insert_deletions[id]['pos']
            insert.save()

        # construct.insertions.add(insert)

    #MODIFICATIONS
    for modification in d['modifications']:
        mod, created = ConstructModification.objects.get_or_create(construct=construct, modification=modification['type'],position_type=modification['position'][0],
                                                   pos_start=modification['position'][1][0],
                                                   pos_end=modification['position'][1][1],remark=modification['remark'] )
        # construct.modifications.add(mod)


    #EXPRESSION
    if 'expression' in d:
        if 'expr_method' in d['expression']:

            if 'expr_remark' not in d['expression']:
                d['expression']['expr_remark'] = ''

            if d['expression']['host_cell'] == 'other [See next field]':
                d['expression']['host_cell'] = d['expression']['other_host_cell']

            if d['expression']['host_cell_type'] == 'other [See next field]':
                d['expression']['host_cell_type'] = d['expression']['other_host']

            if d['expression']['expr_method'] == 'Other [In case of E.Coli or Yeast recombinant expression]':
                d['expression']['expr_method'] = d['expression']['expr_other']

            construct.expression,created = ExpressionSystem.objects.get_or_create(expression_method=d['expression']['expr_method'],
                                                            host_cell_type=d['expression']['host_cell_type'],
                                                            host_cell=d['expression']['host_cell'],
                                                            remarks=d['expression']['expr_remark'])



    #solubilization
    if 'solubilization' in d:
        if 'deterg_type' in d['solubilization']:
            c_list = ChemicalList()
            list_name,created  = ChemicalListName.objects.get_or_create(name='Solubilization')
            c_list.name = list_name
            c_list.save()
            for item,value in d['solubilization'].items():
                if item.startswith(('deterg_type')):
                    d_id = ''
                    if item != 'deterg_type': #if it has deterg_type_2 index
                        d_id = "_" + item.split('_')[2]

                    if value == 'other [See next field]':
                        value = d['raw_data']['other_deterg_type'+ d_id]

                    ct, created = ChemicalType.objects.get_or_create(name='detergent')
                    chem, created = Chemical.objects.get_or_create(name=value, chemical_type=ct)
                    if 'deterg_concentr' + d_id in d['solubilization']:
                        cc, created = ChemicalConc.objects.get_or_create(concentration=d['solubilization']['deterg_concentr' + d_id], concentration_unit=d['solubilization']['deterg_concentr_unit' + d_id], chemical=chem)
                    else: #if no concentr is dictionary, then it was inputted before caputring concentration for additinal chemicals
                        cc, created = ChemicalConc.objects.get_or_create(concentration='', concentration_unit='',chemical=chem)
                    c_list.chemicals.add(cc)

            ct, created = ChemicalType.objects.get_or_create(name='additive')
            chem, created = Chemical.objects.get_or_create(name=d['solubilization']['solub_additive'], chemical_type=ct)
            cc, created = ChemicalConc.objects.get_or_create(concentration=d['solubilization']['additive_concentr'], concentration_unit=d['solubilization']['addit_concentr_unit'], chemical=chem)
            c_list.chemicals.add(cc)

            solubilization = Solubilization.objects.create(chemical_list = c_list)

            construct.solubilization = solubilization
            construct.save()

            #Purification
            purification = Purification.objects.create()
            for puri,step in d['solubilization'].items():
                if not puri.startswith(('chem_enz_treatment','sol_remark')):
                    continue
                else:
                    if step == 'other [See next field]':
                        continue #there will be sol_remark instead
                    if step == 'None':
                        continue #dont put in none step
                    s,created = PurificationStep.objects.get_or_create(name=step)
                    purification.steps.add(s)
            construct.purification = purification
    construct.save()

    #CRYSTALLIZATION
    if 'crystallization' in d:
        if 'crystal_type' in d['crystallization']:
            c = Crystallization()


            if d['crystallization']['crystal_method'] == 'other [See next field]':
                d['crystallization']['crystal_method'] = d['raw_data']['other_method']

            if d['crystallization']['crystal_type'] == 'other [See next field]':
                d['crystallization']['crystal_type'] = d['raw_data']['other_crystal_type']

            sub_name = "" if 'lcp_lipid' not in d['crystallization'] else d['crystallization']['lcp_lipid']
            c_type, created = CrystallizationTypes.objects.get_or_create(name=d['crystallization']['crystal_type'], sub_name=sub_name)
            c_method, created = CrystallizationMethods.objects.get_or_create(name=d['crystallization']['crystal_method'])

            c.crystal_type = c_type
            c.crystal_method = c_method
            if 'crystal_remark' in d['crystallization']:
                c.remarks = d['crystallization']['crystal_remark']
            c.temp = d['crystallization']['temperature']

            if d['crystallization']['ph']=='single_ph':
                c.ph_start = d['crystallization']['ph_single']
                c.ph_end = d['crystallization']['ph_single']
            else:
                c.ph_start = d['crystallization']['ph_range_one']
                c.ph_end = d['crystallization']['ph_range_two']


            c.protein_conc = d['crystallization']['protein_concentr']
            c.protein_conc_unit = d['crystallization']['protein_conc_unit']

            c.save()

            #MAKE LISTS
            c_list = ChemicalList()
            list_name,created  = ChemicalListName.objects.get_or_create(name='Additional')
            c_list.name = list_name
            c_list.save()
            if 'chemical_components' in d['crystallization']:
                for chemical in d['crystallization']['chemical_components']:
                    if 'type' not in chemical: #to fix legacy json files
                        chemical['type'] = 'unknown'
                    ct, created = ChemicalType.objects.get_or_create(name=chemical['type'])
                    chem, created = Chemical.objects.get_or_create(name=chemical['component'], chemical_type=ct)
                    cc, created = ChemicalConc.objects.get_or_create(concentration=chemical['value'], concentration_unit=chemical['unit'], chemical=chem)
                    c_list.chemicals.add(cc)
                c.chemical_lists.add(c_list)

            if d['crystallization']['crystal_type']=='lipidic cubic phase': #make list of LCP stuff
                c_list = ChemicalList()
                # c_list.name = d['crystallization']['lcp_lipid']
                list_name,created  = ChemicalListName.objects.get_or_create(name='LCP')
                c_list.name = list_name
                c_list.save()
                ct, created = ChemicalType.objects.get_or_create(name='LCP Lipid additive')
                chem, created = Chemical.objects.get_or_create(name=d['crystallization']['lcp_add'], chemical_type=ct)
                cc, created = ChemicalConc.objects.get_or_create(concentration=d['crystallization']['lcp_conc'], concentration_unit=d['crystallization']['lcp_conc_unit'], chemical=chem)
                c_list.chemicals.add(cc)
                c.chemical_lists.add(c_list)

            #DETERGENT
            if 'detergent' in d['crystallization']:
                c_list = ChemicalList()
                list_name,created  = ChemicalListName.objects.get_or_create(name='Detergent')
                c_list.name = list_name
                c_list.save()
                ct, created = ChemicalType.objects.get_or_create(name='detergent')
                chem, created = Chemical.objects.get_or_create(name=d['crystallization']['detergent'], chemical_type=ct)
                cc, created = ChemicalConc.objects.get_or_create(concentration=d['crystallization']['deterg_conc'], concentration_unit=d['crystallization']['deterg_conc_unit'], chemical=chem)
                c_list.chemicals.add(cc)
                c.chemical_lists.add(c_list)

            #LIPID
            if 'lipid' in d['crystallization']:
                c_list = ChemicalList()
                list_name,created  = ChemicalListName.objects.get_or_create(name='Lipid')
                c_list.name = list_name
                c_list.save()
                ct, created = ChemicalType.objects.get_or_create(name='lipid')
                chem, created = Chemical.objects.get_or_create(name=d['crystallization']['lipid'], chemical_type=ct)
                cc, created = ChemicalConc.objects.get_or_create(concentration=d['crystallization']['lipid_concentr'], concentration_unit=d['crystallization']['lipid_concentr_unit'], chemical=chem)
                c_list.chemicals.add(cc)
                c.chemical_lists.add(c_list)



            #Use ligand function to get ligand if it exists or otherwise create. Lots of checks for inchi/smiles/name
            ligand = get_or_make_ligand(d['construct_crystal']['ligand_id'],d['construct_crystal']['ligand_id_type'],d['construct_crystal']['ligand_name'])
            if 'ligand_activity' not in d['construct_crystal']:
                d['construct_crystal']['ligand_activity'] = 'unknown'
            if ligand and 'ligand_activity' in d['construct_crystal']:
                role_slug = slugify(d['construct_crystal']['ligand_activity'])
                try:
                    lr, created = LigandRole.objects.get_or_create(slug=role_slug,
                    defaults={'name': d['construct_crystal']['ligand_activity']})
                except IntegrityError:
                    lr = LigandRole.objects.get(slug=role_slug)
            if ligand:
                ligand_c = CrystallizationLigandConc()
                ligand_c.construct_crystallization = c
                ligand_c.ligand = ligand
                if lr:
                    ligand_c.ligand_role = lr
                if 'ligand_conc' in d['construct_crystal']:
                    ligand_c.ligand_conc = d['construct_crystal']['ligand_conc']
                if 'ligand_conc_unit' in d['construct_crystal']:
                    ligand_c.ligand_conc_unit = d['construct_crystal']['ligand_conc_unit']
                ligand_c.save()

                c.ligands.add(ligand_c)

            construct.crystallization = c

    construct.save()

def convert_ordered_to_disordered_annotation(d):
    if 'raw_data' not in d:
        d['raw_data'] = {}
    d['raw_data']['pdb'] = d['construct_crystal']['pdb']
    d['raw_data']['uniprot'] = d['construct_crystal']['uniprot']
    if 'pdb_name' not in d['construct_crystal']:
        d['raw_data']['pdb_name'] = d['construct_crystal']['uniprot']+"_construct"
    else:
        d['raw_data']['pdb_name'] = d['construct_crystal']['pdb_name']


    #contributor
    for k,v in d['contact_info'].items():
        d['raw_data'][k] = v

    #growth information
    if 'crystal_growth' in d:
        d['raw_data']['crystal_type'] = d['crystal_growth'][0]['grow_method']

    i = 2 #starts with two for some reason, ask Anna
    insert_starts = {}
    for aux,v in d['auxiliary'].items():
        # print(aux)
        d['raw_data']['protein_type_'+str(i)] = v['type']
        d['raw_data']['position_'+str(i)] = v['position']
        d['raw_data']['presence_'+str(i)] = v['presence']
        d['raw_data']['fusion_prot_'+str(i)] = "Please Select"
        d['raw_data']['aux'+str(i)+'_subtype'] = v['subtype']
        if 'start' in v:
            insert_starts[v['start']] = str(i)
            #         "position": "N-term",
            # "presence": "NO",
            # "type": "signal",
            # "subtype": "hemagglutinin"
        i+=1

    # print(insert_starts)

    i = 2
    for mut in d['mutations']:
        d['raw_data']['mut_aa_'+str(i)] = mut['mut']
        d['raw_data']['wt_aa_'+str(i)] = mut['wt']
        d['raw_data']['aa_no_'+str(i)] = mut['pos']
        d['raw_data']['mut_type_'+str(i)] = ''
        i+=1

    i = 2
    for dele in d['deletions']:
        if str(dele['start']) not in insert_starts.keys():
            d['raw_data']['delet_start_'+str(i)] = dele['start']
            d['raw_data']['delet_end_'+str(i)] = dele['end']
            d['raw_data']['delet_origin_'+str(i)] = dele['origin']
            d['raw_data']['deletion_'+str(i)] = "del_range"
            i+=1
        else:
            i_id = insert_starts[str(dele['start'])]
            d['raw_data']['insert_pos_type_'+str(i_id)] = "ins_range"
            d['raw_data']['ins_start_'+str(i_id)] = dele['start']
            d['raw_data']['ins_end_'+str(i_id)] = dele['end']

    i = 2
    for mod in d['modifications']:
        d['raw_data']['aamod_'+str(i)] = mod['type']
        if mod['position'][0]=='single':
            d['raw_data']['aamod_single_'+str(i)] = mod['position'][1][0]
        else:
            d['raw_data']['aamod_pair_one_'+str(i)] = mod['position'][1][0]
            d['raw_data']['aamod_pair_two_'+str(i)] = mod['position'][1][1]
        d['raw_data']['aamod_position_'+str(i)] = mod['position'][0]
        d['raw_data']['mod_remark_'+str(i)] = mod['remark']
        i+=1

    return d
