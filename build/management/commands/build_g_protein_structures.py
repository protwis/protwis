from django.core.management.base import BaseCommand, CommandError
from build.management.commands.base_build import Command as BaseBuild
from django.conf import settings
from django.db import connection
from django.db import IntegrityError

from protein.models import (Protein, ProteinGProtein,ProteinGProteinPair, ProteinConformation, ProteinState, ProteinFamily, ProteinAlias,
        ProteinSequenceType, Species, Gene, ProteinSource, ProteinSegment)
from residue.models import (ResidueNumberingScheme, ResidueGenericNumber, Residue, ResidueGenericNumberEquivalent)
from signprot.models import SignprotComplex, SignprotStructure, SignprotStructureExtraProteins
from common.models import WebResource, WebLink, Publication
from structure.models import StructureType, StructureStabilizingAgent, Structure

import re
from Bio import pairwise2
from collections import OrderedDict
import logging
import shlex, subprocess
from io import StringIO
from Bio.PDB import PDBParser,PPBuilder
from Bio import pairwise2
import pprint
import json
import yaml
from urllib.request import urlopen, Request

import traceback
import sys, os
import datetime


AA = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D',
     'CYS':'C', 'GLN':'Q', 'GLU':'E', 'GLY':'G',
     'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K',
     'MET':'M', 'PHE':'F', 'PRO':'P', 'SER':'S',
     'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V',
     'YCM':'C', 'CSD':'C', 'TYS':'Y', 'SEP':'S'} #non-standard AAs


class Command(BaseBuild):

    local_uniprot_dir = os.sep.join([settings.DATA_DIR, 'g_protein_data', 'uniprot'])
    local_uniprot_beta_dir = os.sep.join([settings.DATA_DIR, 'g_protein_data', 'uniprot_beta'])
    local_uniprot_gamma_dir = os.sep.join([settings.DATA_DIR, 'g_protein_data', 'uniprot_gamma'])
    with open(os.sep.join([settings.DATA_DIR, 'g_protein_data', 'g_protein_display_names.yaml']), 'r') as y:
        display_name_lookup = yaml.load(y, Loader=yaml.FullLoader)

    def add_arguments(self, parser):
        parser.add_argument('--purge', default=False, action='store_true', help='Purge G protein structures from database')
        parser.add_argument('--only_signprot_structures', default=False, action='store_true', help='Only build SignprotStructure objects')

    def handle(self, *args, **options):
        self.options = options
        if self.options['purge']:
            Residue.objects.filter(protein_conformation__protein__entry_name__endswith='_a', protein_conformation__protein__family__parent__parent__name='Alpha').delete()
            ProteinConformation.objects.filter(protein__entry_name__endswith='_a', protein__family__parent__parent__name='Alpha').delete()
            Protein.objects.filter(entry_name__endswith='_a', family__parent__parent__name='Alpha').delete()

        if not options['only_signprot_structures']:
            # Building protein and protconf objects for g protein structure in complex
            scs = SignprotComplex.objects.all()
            for sc in scs:
                self.logger.info('Protein, ProteinConformation and Residue build for alpha subunit of {} is building'.format(sc))
                try:
                    # Alpha subunit
                    try:
                        alpha_protein = Protein.objects.get(entry_name=sc.structure.pdb_code.index.lower()+'_a')
                    except:
                        alpha_protein = Protein()
                        alpha_protein.entry_name = sc.structure.pdb_code.index.lower()+'_a'
                        alpha_protein.accession = None
                        alpha_protein.name = sc.structure.pdb_code.index.lower()+'_a'
                        alpha_protein.sequence = sc.protein.sequence
                        alpha_protein.family = sc.protein.family
                        alpha_protein.parent = sc.protein
                        alpha_protein.residue_numbering_scheme = sc.protein.residue_numbering_scheme
                        alpha_protein.sequence_type = ProteinSequenceType.objects.get(slug='mod')
                        alpha_protein.source = ProteinSource.objects.get(name='OTHER')
                        alpha_protein.species = sc.protein.species
                        alpha_protein.save()
                    try:
                        alpha_protconf = ProteinConformation.objects.get(protein__entry_name=sc.structure.pdb_code.index.lower()+'_a')
                    except:
                        alpha_protconf = ProteinConformation()
                        alpha_protconf.protein = alpha_protein
                        alpha_protconf.state = ProteinState.objects.get(slug='active')
                        alpha_protconf.save()
                    pdbp = PDBParser(PERMISSIVE=True, QUIET=True)
                    s = pdbp.get_structure('struct', StringIO(sc.structure.pdb_data.pdb))
                    chain = s[0][sc.alpha]
                    nums = []
                    for res in chain:
                        try:
                            res['CA']
                            nums.append(res.get_id()[1])
                        except:
                            pass

                    resis = Residue.objects.filter(protein_conformation__protein=sc.protein)
                    num_i = 0
                    temp_seq2 = ''
                    pdb_num_dict = OrderedDict()
                    # Create first alignment based on sequence numbers
                    for n in nums:
                        if sc.structure.pdb_code.index=='6OIJ' and n<30:
                            nr = n+6
                        else:
                            nr = n
                        pdb_num_dict[n] = [chain[n], resis.get(sequence_number=nr)]
                    # Find mismatches
                    mismatches = []
                    for n, res in pdb_num_dict.items():
                        if AA[res[0].get_resname()]!=res[1].amino_acid:
                            mismatches.append(res)

                    pdb_lines = sc.structure.pdb_data.pdb.split('\n')
                    seqadv = []
                    for l in pdb_lines:
                        if l.startswith('SEQADV'):
                            seqadv.append(l)
                    mutations, shifted_mutations = OrderedDict(), OrderedDict()
                    # Search for annotated engineered mutations in pdb SEQADV
                    for s in seqadv:
                        line_search = re.search('SEQADV\s{1}[A-Z\s\d]{4}\s{1}([A-Z]{3})\s{1}([A-Z]{1})\s+(\d+)[\s\S\d]{5}([\s\S\d]{12})([A-Z]{3})\s+(\d+)(\s\S+)',s)
                        if line_search!=None:
                            if line_search.group(2)==sc.alpha:
                                if line_search.group(4).strip()==sc.protein.accession:
                                    if line_search.group(3)==line_search.group(6):
                                        mutations[int(line_search.group(3))] = [line_search.group(1), line_search.group(5)]
                                    else:
                                        shifted_mutations[int(line_search.group(3))] = [line_search.group(1), line_search.group(5), int(line_search.group(6))]
                                else:
                                    # Exception for 6G79
                                    if line_search.group(3)!=line_search.group(6) and 'CONFLICT' in line_search.group(7):
                                        mutations[int(line_search.group(3))] = [line_search.group(1), line_search.group(5)]
                                    # Exception for 5G53
                                    if line_search.group(4).strip()!=sc.protein.accession:
                                        mutations[int(line_search.group(3))] = [line_search.group(1), line_search.group(5)]
                    remaining_mismatches = []

                    # Check and clear mismatches that are registered in pdb SEQADV as engineered mutation
                    for m in mismatches:
                        num = m[0].get_id()[1]
                        if num in mutations:
                            if m[0].get_resname()!=mutations[num][0] and m[1].amino_acid!=AA[mutations[num][1]]:
                                remaining_mismatches.append(m)
                        elif num in shifted_mutations:
                            remaining_mismatches.append(m)
                        else:
                            remaining_mismatches.append(m)

                    ### sanity check
                    # print(sc)
                    # print(mutations)
                    # print(shifted_mutations)
                    # print(mismatches)
                    # print(remaining_mismatches)
                    # pprint.pprint(pdb_num_dict)

                    # Mismatches remained possibly to seqnumber shift, making pairwise alignment to try and fix alignment
                    if len(remaining_mismatches)>0 and sc.structure.pdb_code.index not in ['6OIJ','6OY9', '6OYA', '6LPB']:
                        ppb = PPBuilder()
                        seq = ''
                        for pp in ppb.build_peptides(chain, aa_only=False):
                            seq += str(pp.get_sequence())
                        pw2 = pairwise2.align.localms(sc.protein.sequence, seq, 2, -1, -.5, -.1)
                        ref_seq, temp_seq = str(pw2[0][0]), str(pw2[0][1])
                        wt_pdb_dict = OrderedDict()
                        pdb_wt_dict = OrderedDict()
                        j, k = 0, 0
                        for i, ref, temp in zip(range(0,len(ref_seq)), ref_seq, temp_seq):
                            if ref!='-' and temp!='-':
                                wt_pdb_dict[resis[j]] = pdb_num_dict[nums[k]]
                                pdb_wt_dict[pdb_num_dict[nums[k]][0]] = resis[j]
                                j+=1
                                k+=1
                            elif ref=='-':
                                wt_pdb_dict[i] = pdb_num_dict[nums[k]]
                                pdb_wt_dict[pdb_num_dict[nums[k]][0]] = i
                                k+=1
                            elif temp=='-':
                                wt_pdb_dict[resis[j]] = i
                                pdb_wt_dict[i] = resis[j]
                                j+=1
                        for i, r in enumerate(remaining_mismatches):
                            # Adjust for shifted residue when residue is a match
                            if r[0].get_id()[1]-remaining_mismatches[i-1][0].get_id()[1]>1:
                                pdb_num_dict[r[0].get_id()[1]-1][1] = pdb_wt_dict[chain[r[0].get_id()[1]-1]]
                            # Adjust for shifted residue when residue is mutated and it's logged in SEQADV
                            if r[0].get_id()[1] in shifted_mutations:
                                pdb_num_dict[r[0].get_id()[1]][1] = resis.get(sequence_number=shifted_mutations[r[0].get_id()[1]][2])
                            # Adjust for shift
                            else:
                                pdb_num_dict[r[0].get_id()[1]][1] = pdb_wt_dict[r[0]]

                    bulked_residues = []
                    for key, val in pdb_num_dict.items():
                        # print(key, val) # sanity check
                        if not isinstance(val[1], int):
                         res_obj = Residue()
                         res_obj.sequence_number = val[0].get_id()[1]
                         res_obj.amino_acid = AA[val[0].get_resname()]
                         res_obj.display_generic_number = val[1].display_generic_number
                         res_obj.generic_number = val[1].generic_number
                         res_obj.protein_conformation = alpha_protconf
                         res_obj.protein_segment = val[1].protein_segment
                         bulked_residues.append(res_obj)
                        else:
                         self.logger.info('Skipped {} as no annotation was present, while building for alpha subunit of {}'.format(val[1], sc))

                    Residue.objects.bulk_create(bulked_residues)
                    self.logger.info('Protein, ProteinConformation and Residue build for alpha subunit of {} is finished'.format(sc))
                except Exception as msg:
                    #print('Protein, ProteinConformation and Residue build for alpha subunit of {} has failed'.format(sc))
                    #print(msg)
                    #print(traceback.format_exc())
                    #exit(0)
                    self.logger.info('Protein, ProteinConformation and Residue build for alpha subunit of {} has failed'.format(sc))


        # Build SignprotStructure objects from non-complex signprots
        g_prot_alphas = Protein.objects.filter(family__slug__startswith='100_001', accession__isnull=False)#.filter(entry_name='gnai1_human')
        complex_structures = SignprotComplex.objects.all().values_list('structure__pdb_code__index', flat=True)
        for a in g_prot_alphas:
            pdb_list = self.get_pdb_ids(a.accession)
            for pdb in pdb_list:
                if pdb not in complex_structures:
                    try:
                        data = self.fetch_gprot_data(pdb, a)
                        if data:
                            self.build_g_prot_struct(a, pdb, data)
                    except Exception as msg:
                        self.logger.error('SignprotStructure of {} {} failed\n{}: {}'.format(a.entry_name, pdb, type(msg), msg))

    def get_pdb_ids(self, uniprot_id):
        pdb_list = []
        data = { "query": { "type": "terminal", "service": "text", "parameters":{ "attribute":"rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession", "operator":"in", "value":[ uniprot_id ] } }, "request_options": { "pager": {"start": 0,"rows": 99999 }}, "return_type": "entry" }
        url = 'http://search.rcsb.org/rcsbsearch/v1/query'
        req = Request(url)
        req.add_header('Content-Type', 'application/json; charset=utf-8')
        jsondata = json.dumps(data)
        jsondataasbytes = jsondata.encode('utf-8')
        req.add_header('Content-Length', len(jsondataasbytes))
        response = urlopen(req, jsondataasbytes)
        rr = response.read()
        if len(rr)==0:
            return []
        out = json.loads(rr)
        for i in out['result_set']:
            pdb_list.append(i['identifier'])
        response.close()
        return pdb_list

    def fetch_gprot_data(self, pdb, alpha_protein):
        data = {}
        beta_uniprots = os.listdir(self.local_uniprot_beta_dir)
        gamma_uniprots = os.listdir(self.local_uniprot_gamma_dir)

        response = urlopen('https://data.rcsb.org/rest/v1/core/entry/{}'.format(pdb))
        json_data = json.loads(response.read())

        data['method'] = json_data['exptl'][0]['method']
        if data['method'].startswith("THEORETICAL") or data['method']=='SOLUTION NMR':
            return None
        if 'pubmed_id' in json_data['rcsb_entry_container_identifiers']:
            data['pubmedId'] = json_data['rcsb_entry_container_identifiers']['pubmed_id']
        else:
            data['pubmedId'] = None
        
        # Format server time stamp to match release date shown on entry pages
        date = datetime.date.fromisoformat(json_data['rcsb_accession_info']['initial_release_date'][:10])
        date += datetime.timedelta(days=1)
        data['release_date'] = datetime.date.isoformat(date)
        data['resolution'] = json_data['rcsb_entry_info']['resolution_combined'][0]
        entities_num = len(json_data['rcsb_entry_container_identifiers']['polymer_entity_ids'])
        data['alpha'] = alpha_protein.accession
        data['alpha_chain'] = None
        data['alpha_coverage'] = None
        data['beta'] = None
        data['beta_chain'] = None
        data['gamma'] = None
        data['gamma_chain'] = None
        data['other'] = []
        for i in range(1,entities_num+1):
            response = urlopen('https://data.rcsb.org/rest/v1/core/polymer_entity/{}/{}'.format(pdb, i))
            json_data = json.loads(response.read())
            if 'uniprot_ids' in json_data['rcsb_polymer_entity_container_identifiers']:
                for j, u_id in enumerate(json_data['rcsb_polymer_entity_container_identifiers']['uniprot_ids']):
                    if u_id+'.txt' in beta_uniprots:
                        data['beta'] = u_id
                        data['beta_chain'] = json_data['rcsb_polymer_entity_container_identifiers']['auth_asym_ids'][j][0]
                    elif u_id+'.txt' in gamma_uniprots:
                        data['gamma'] = u_id
                        data['gamma_chain'] = json_data['rcsb_polymer_entity_container_identifiers']['auth_asym_ids'][j][0]
                    elif u_id==alpha_protein.accession:
                        data['alpha'] = u_id
                        data['alpha_coverage'] = json_data['entity_poly']['rcsb_sample_sequence_length']
                        # pprint.pprint(json_data)
                        try:
                            data['alpha_chain'] = json_data['rcsb_polymer_entity_container_identifiers']['auth_asym_ids'][j][0]
                        except IndexError as e:
                            data['alpha_chain'] = json_data['rcsb_polymer_entity_container_identifiers']['auth_asym_ids'][j-1][0]
                    else:
                        if json_data['rcsb_polymer_entity']['pdbx_description'] not in data['other']:
                            data['other'].append(json_data['rcsb_polymer_entity']['pdbx_description'])
            else:
                if json_data['rcsb_polymer_entity']['pdbx_description'] not in data['other']:
                    data['other'].append(json_data['rcsb_polymer_entity']['pdbx_description'])
        return data   

    def build_g_prot_struct(self, alpha_prot, pdb, data):
        ss = SignprotStructure()
        pdb_code, p_c = WebLink.objects.get_or_create(index=pdb, web_resource=WebResource.objects.get(slug='pdb'))
        pub_date = data['release_date']
        # Structure type
        if 'x-ray' in data['method'].lower():
            structure_type_slug = 'x-ray-diffraction'
        elif 'electron' in data['method'].lower():
            structure_type_slug = 'electron-microscopy'
        else:
            structure_type_slug = '-'.join(data['method'].lower().split(' '))
        try:
            structure_type = StructureType.objects.get(slug=structure_type_slug)
        except StructureType.DoesNotExist as e:
            structure_type, c = StructureType.objects.get_or_create(slug=structure_type_slug, name=data['method'])
            self.logger.info('Created StructureType:'+str(structure_type))
        # Publication
        if data['pubmedId']:
            try:
                pub = Publication.objects.get(web_link__index=data['pubmedId'])
            except Publication.DoesNotExist as e:
                pub = Publication()
                wl, created = WebLink.objects.get_or_create(index=data['pubmedId'], web_resource=WebResource.objects.get(slug='pubmed'))
                pub.web_link = wl
                pub.update_from_pubmed_data(index=data['pubmedId'])
                pub.save()
                self.logger.info('Created Publication:'+str(pub))
        else:
            pub = None
        ss.pdb_code = pdb_code
        ss.structure_type = structure_type
        ss.resolution = data['resolution']
        ss.publication_date = pub_date
        ss.publication = pub
        ss.protein = alpha_prot
        ss.save()
        # Stabilizing agent
        for o in data['other']:
            if len(o)>75:
                continue
            if o=='REGULATOR OF G-PROTEIN SIGNALING 14':
                o = 'Regulator of G-protein signaling 14'
            elif o=='Nanobody 35':
                o = 'Nanobody-35'
            elif o=='ADENYLATE CYCLASE, TYPE V':
                o = 'Adenylate cyclase, type V'
            elif o=='1-phosphatidylinositol-4,5-bisphosphate phosphodiesterase beta-3':
                o = '1-phosphatidylinositol 4,5-bisphosphate phosphodiesterase beta-3'
            stabagent, sa_created = StructureStabilizingAgent.objects.get_or_create(slug=o.replace(' ','-').replace(' ','-'), name=o)
            ss.stabilizing_agents.add(stabagent)
        ss.save()
        # Extra proteins
        # Alpha - ### A bit redundant, consider changing this in the future
        if data['alpha']:
            alpha_sep = SignprotStructureExtraProteins()
            alpha_sep.wt_protein = alpha_prot
            alpha_prot.structure = ss
            alpha_sep.protein_conformation = ProteinConformation.objects.get(protein=alpha_prot)
            alpha_sep.display_name = self.display_name_lookup[alpha_prot.family.name]
            alpha_sep.note = None
            alpha_sep.chain = data['alpha_chain']
            alpha_sep.category = 'G alpha'
            alpha_sep.wt_coverage = round(data['alpha_coverage']/len(alpha_prot.sequence)*100)
            alpha_sep.save()
            ss.extra_proteins.add(alpha_sep)
        # Beta
        if data['beta']:
            beta_prot = Protein.objects.get(accession=data['beta'])
            beta_sep = SignprotStructureExtraProteins()
            beta_sep.wt_protein = beta_prot
            beta_sep.structure = ss
            beta_sep.protein_conformation = ProteinConformation.objects.get(protein=beta_prot)
            beta_sep.display_name = self.display_name_lookup[beta_prot.name]
            beta_sep.note = None
            beta_sep.chain = data['beta_chain']
            beta_sep.category = 'G beta'
            beta_sep.wt_coverage = None
            beta_sep.save()
            ss.extra_proteins.add(beta_sep)
        # Gamma
        if data['gamma']:
            gamma_prot = Protein.objects.get(accession=data['gamma'])
            gamma_sep = SignprotStructureExtraProteins()
            gamma_sep.wt_protein = gamma_prot
            gamma_sep.structure = ss
            gamma_sep.protein_conformation = ProteinConformation.objects.get(protein=gamma_prot)
            gamma_sep.display_name = self.display_name_lookup[gamma_prot.name]
            gamma_sep.note = None
            gamma_sep.chain = data['gamma_chain']
            gamma_sep.category = 'G gamma'
            gamma_sep.wt_coverage = None
            gamma_sep.save()
            ss.extra_proteins.add(gamma_sep)
        ss.save()
        self.logger.info('Created SignprotStructure: {}'.format(ss.pdb_code))

