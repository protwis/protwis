# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 15:50:57 2016

@author: Gaspar Pandy
"""
from build.management.commands.base_build import Command as BaseBuild
from django.db.models import Q

from protein.models import Protein, ProteinConformation, ProteinSequenceType, ProteinSource, ProteinState
from residue.models import Residue
from structure.models import Structure, PdbData, StructureType
from structure.sequence_parser import SequenceParser
from structure.functions import PdbChainSelector, PdbStateIdentifier, ParseStructureCSV, get_pdb_ids
from structure.management.commands.structure_yaml_editor import StructureYaml
from construct.functions import *
from common.models import WebResource, WebLink, Publication
from tools.management.commands.blast_recent_PDB import Command as BlastRecentPDB

import Bio.PDB as PDB
from datetime import datetime
import urllib
import re
import json
import os
import xmltodict
import yaml
import shlex
import subprocess
import pprint
import csv


structs_with_missing_x50 = ['5EM9', '5AER', '3BEF', '3LU9', '3HKI', '3HKJ', '1NRR', '1NRQ', '1NRP', '1NRO', '1NRN', '3QDZ', '2ZPK', '1YTV', '4JQI', '6NI2', '5YD3',
                            '5YD5', '5YD4', '5YY4', '6KVA', '6KVF', '4NUU', '4NUV', '6K3F', '1XWD', '4AY9', '4MQW', '2XWT', '3G04', '4KT1', '4QXE', '4QXF', '4UFR',
                            '4UFS', '4BSU', '4BSR', '4BSS', '4BST', '5II0', '6PFO', '6PGQ', '3AQF', '4RWF', '4RWG', '5V6Y', '6D1U', '6ZHO', '6ZIS', '6UMG', '6V2E',
                            '2XDG', '2QKH', '3C59', '3C5T', '3IOL', '5E94', '4ZGM', '5OTW', '2A83', '3CZF', '4ERS', '4LF3', '3L2J', '3H3G', '3C4M', '5EMB', '4Z8J',
                            '3N94', '3B3I', '3B6S', '3DTX', '3HCV', '5IB1', '5IB3', '5IB2', '5IB5', '5IB4', '5DEG', '5DEF', '1OF2', '1OGT', '2X57', '5AFB', '5CMN',
                            '6VHH', '2BO2', '2BOU', '2BOX', '4DLO', '5K5T', '5K5S', '5FBK', '5FBH', '4PAS', '4MR7', '4MR8', '4MQE', '4MQF', '4MS3', '4MS4', '4MRM',
                            '4MS1', '4MR9', '4F11', '4F12', '6M8R', '6OCP', '3KS9', '4XAQ', '4XAS', '5CNI', '5CNJ', '5KZN', '5KZQ', '3SM9', '4XAR', '5CNK', '5CNM',
                            '6B7H', '3LMK', '6N4Y', '6N4X', '6N50', '3MQ4', '5C5C', '6E5V', '6BSZ', '6BT5', '6C0B', '5BPB', '5BQC', '5UWG', '5BPQ', '5BQE', '6NE1',
                            '5CL1', '5CM4', '5URZ', '5URY', '6O39', '5URV', '4Z33', '6NE2', '6NE4', '6O3B', '6O3A', '5WBS', '5T44', '5UN6', '5UN5', '6NDZ', '5KZV',
                            '5KZY', '5KZZ', '7JQD', '2PUX', '2PV9', '6EXJ', '7JNZ', '7NW3', '4F8K', '1RY1', '2J28', '4UE5', '6O9I', '6O9H', '7ALO', '6SKA', '4DLQ',
                            '5OVP', '6SKE', '5FTT', '5FTU', '4YEB', '4RMK', '4RML', '6JBU', '6IDX', '5KVM', '6V55', '7NJZ', '3N96', '3N93', '3N95', '7D86', '7DO4',
                            '4OAJ', '4LI1', '4LI2', '7R84', '7R85', '7R86', '1EWK', '1EWT', '1EWV', '1ISR', '1ISS', '2E4U', '2E4V', '2E4W', '2E4X', '2E4Y', '1DDV', 
                            '2E4Z', '5X2M', '5X2N', '5X2O', '5X2P', '5X2Q', '7N95', '7N97', '7N9S', '1IJY', '4F0A', '6AHY', '6TFB', '6TFM', '4C79', '4C7A', '7DF9', 
                            '7DFA', '7DFB', '7DFC', '7P8X', '7P93']


class Command(BaseBuild):

    def add_arguments(self, parser):
        super(Command, self).add_arguments(parser=parser)
        parser.add_argument('--classified', help="Use PDB's 'G protein-coupled receptors' classification", default=False,
                            action='store_true')
        parser.add_argument('-r', help="Query specific receptor(s) with UniProt entry names", default=False, type=str, nargs='+')

    def handle(self, *args, **options):
        if options['verbosity'] in [0,1,2,3]:
            self.verbose = True
        else:
            self.verbose = False
        if options['classified']:
            q = QueryPDBClassifiedGPCR()
            q.new_xtals(self.verbose)
        else:
            if options['r']:
                self.uniprots = self.fetch_accession_from_entryname(options['r'])
            else:
                self.uniprots = self.get_all_GPCR_uniprots()
            self.pdbs = ParseStructureCSV().pdb_ids
            self.prepare_input(options['proc'], self.uniprots)

    def main_func(self, positions, iteration, count, lock):
        if not positions[1]:
            uniprot_list = self.uniprots[positions[0]:]
        else:
            uniprot_list = self.uniprots[positions[0]:positions[1]]
        q = QueryPDB(self.uniprots, self.pdbs)

        brp = BlastRecentPDB()
        blast_pdbs = brp.run()
        self.blast_uniprot_dict = {}
        for b in blast_pdbs:
            blast_uniprots = q.pdb_request_by_pdb(b, 'polymer_entity')
            if not blast_uniprots:
                q.csv_list.append(b)
                continue
            for bu in blast_uniprots:
                if bu not in self.blast_uniprot_dict:
                    self.blast_uniprot_dict[bu] = [b]
                else:
                    self.blast_uniprot_dict[bu].append(b)

        consider_list, error_list = [], []
        print('{} number of receptors to check'.format(len(uniprot_list)))

        # uniprot_list = ['P28223']

        for uni in uniprot_list:
            # print(uni)
            q.new_xtals(uni, self.blast_uniprot_dict)
            for i in q.consider_list:
                if i not in consider_list and i not in structs_with_missing_x50:
                    consider_list.append(i)
            for i in q.error_list:
                if i not in error_list:
                    error_list.append(i)
        if self.verbose:
            print('Missing from db: ', q.db_list)
            print('Missing from csv: ', q.csv_list)
            print('Structures with missing x50s: {} structures {}'.format(len(consider_list), consider_list))
            print('Structures with an error: {} structures {}'.format(len(error_list), error_list))

    def fetch_accession_from_entryname(self, listof_entrynames):
        return [i.accession for i in Protein.objects.filter(entry_name__in=listof_entrynames)]

    def get_all_GPCR_uniprots(self):
        try:
            uniprots = [i.accession for i in Protein.objects.filter(accession__isnull=False).filter(family__slug__istartswith='00')]
            if len(uniprots)<100:
                raise Exception()
        except:
            uniprots = [i.split('.')[0] for i in os.listdir('/protwis/data/protwis/gpcr/protein_data/uniprot/')]
        return uniprots

    # Deprecated
    def get_all_yamls(self):
        yamls = [i.split('.')[0] for i in os.listdir('/protwis/data/protwis/gpcr/structure_data/structures/')]
        return yamls


class QueryPDB():
    ''' Queries PDB using GPCRdb protein and structure entries. If those are not available, it uses the structure and uniprot data folders.
    '''
    def __init__(self, uniprots, pdbs):
        self.exceptions = []
        self.uniprots = uniprots
        self.pdbs = pdbs
        self.db_list, self.csv_list = [], []
        self.missing_x50_list = ['4KNG','3N7P','3N7R','3N7S','4HJ0','6DKJ','5OTT','5OTU','5OTV','5OTX','6GB1']
        self.missing_x50_exceptions = ['6TPG','6TPJ']
        self.consider_list, self.error_list = [], []

    def new_xtals(self, uniprot, blast_uniprot_dict):
        ''' List GPCR crystal structures missing from GPCRdb and the csv files. Adds missing structures to DB.
        '''
        # structs = self.pdb_request_by_uniprot(uniprot)
        structs = get_pdb_ids(uniprot)
        try:
            protein = Protein.objects.get(accession=uniprot)
        except:
            protein = None
        try:
            x50s = Residue.objects.filter(protein_conformation__protein=protein, generic_number__label__in=['1x50','2x50','3x50','4x50','5x50','6x50','7x50'])
        except:
            x50s = None
        if uniprot in blast_uniprot_dict:
            for i in blast_uniprot_dict[uniprot]:
                if i not in structs:
                    if structs==['null']:
                        structs = [i]
                    else:
                        structs.append(i)
        if structs!=['null']:
            for s in structs:
                # print(s)
                missing_from_db, missing_csv = False, False
                try:
                    st_obj = Structure.objects.get(pdb_code__index=s)
                except:
                    if s not in self.exceptions:
                        check = self.pdb_request_by_pdb(s, 'entry')
                        if check:
                            self.db_list.append(s)
                            missing_from_db = True

                if s not in self.pdbs and s not in self.exceptions:
                    if s not in self.db_list:
                        check = self.pdb_request_by_pdb(s, 'entry')
                    else:
                        check = True
                    if check:
                        self.csv_list.append(s)
                        missing_csv = True
                if not missing_from_db:
                    continue
                try:
                    pdb_data_dict = fetch_pdb_info(s, protein, new_xtal=True)
                    # pprint.pprint(pdb_data_dict)
                    exp_method = pdb_data_dict['experimental_method']
                    if exp_method=='Electron Microscopy':
                        st_type = StructureType.objects.get(slug='electron-microscopy')
                    elif exp_method=='X-ray diffraction':
                        st_type = StructureType.objects.get(slug='x-ray-diffraction')
                    elif exp_method=='Electron crystallography':
                        st_type = StructureType.objects.get(slug='electron-crystallography')

                    if 'deletions' in pdb_data_dict:
                        for d in pdb_data_dict['deletions']:
                            presentx50s = []
                            for x in x50s:
                                if not d['start']<x.sequence_number<d['end']:
                                    presentx50s.append(x)
                            # Filter out ones without all 7 x50 positions present in the xtal
                            if len(presentx50s)!=7:
                                if s not in self.missing_x50_list:
                                    self.consider_list.append(s)
                                if s not in self.missing_x50_exceptions:
                                    try:
                                        del self.db_list[self.db_list.index(s)]
                                        missing_from_db = False
                                        del self.csv_list[self.csv_list.index(s)]
                                    except:
                                        pass
                    if 'not_observed' in pdb_data_dict:
                        for no in pdb_data_dict['not_observed']:
                            presentx50s = []
                            for x in x50s:
                                if not no[0]<x.sequence_number<no[1]:
                                    presentx50s.append(x)
                            if len(presentx50s)!=7:
                                if s not in self.missing_x50_list:
                                    self.consider_list.append(s)
                                if s not in self.missing_x50_exceptions:
                                    try:
                                        del self.db_list[self.db_list.index(s)]
                                        missing_from_db = False
                                        del self.csv_list[self.csv_list.index(s)]
                                    except:
                                        pass
                    else:
                        print('Warning: no deletions in pdb info, check {}'.format(s))
                        continue

                    if missing_from_db:
                        pref_chain = ''
                        resolution = pdb_data_dict['resolution']
                        pdb_code, created = WebLink.objects.get_or_create(index=s, web_resource=WebResource.objects.get(slug='pdb'))
                        pdbl = PDB.PDBList()
                        pdbl.retrieve_pdb_file(s, pdir='./', file_format="pdb")
                        with open('./pdb{}.ent'.format(s).lower(),'r') as f:
                            lines = f.readlines()
                        pdb_file = ''
                        publication_date, pubmed, doi = '','',''
                        state = ProteinState.objects.get(slug='inactive')
                        new_prot, created = Protein.objects.get_or_create(entry_name=s.lower(), accession=None, name=s.lower(), sequence=pdb_data_dict['wt_seq'], family=protein.family,
                                                                          parent=protein, residue_numbering_scheme=protein.residue_numbering_scheme,
                                                                          sequence_type=ProteinSequenceType.objects.get(slug='mod'), source=ProteinSource.objects.get(name='OTHER'),
                                                                          species=protein.species)
                        new_prot_conf, created = ProteinConformation.objects.get_or_create(protein=new_prot, state=state)
                        for line in lines:
                            if line.startswith('REVDAT   1'):
                                publication_date = line[13:22]
                            if line.startswith('JRNL        PMID'):
                                pubmed = line[19:].strip()
                            if line.startswith('JRNL        DOI'):
                                doi = line[19:].strip()
                            pdb_file+=line
                        pdb_data, created = PdbData.objects.get_or_create(pdb=pdb_file)
                        d = datetime.strptime(publication_date,'%d-%b-%y')
                        publication_date = d.strftime('%Y-%m-%d')
                        try:
                            if doi!='':
                                publication = Publication.get_or_create_from_doi(doi)
                            elif pubmed!='':
                                publication = Publication.get_or_create_from_pubmed(pubmed)
                        except:
                            pass
                        pcs = PdbChainSelector(s, protein)
                        pcs.run_dssp()
                        preferred_chain = pcs.select_chain()

                        # Run state identification

                        # Create yaml files
                        # with open(os.sep.join([settings.DATA_DIR, 'structure_data','constructs', '{}.yaml'.format(pdb_code.index)]), 'w') as construct_file:
                        #     yaml.dump({'name': pdb_code.index.lower(), 'protein': protein.entry_name}, construct_file, indent=4)
                        # with open(os.sep.join([settings.DATA_DIR, 'structure_data','structures','{}.yaml'.format(pdb_code.index)]), 'w') as structure_file:
                        #     struct_yaml_dict = {'construct': pdb_code.index.lower(), 'pdb': pdb_code.index, 'preferred_chain': preferred_chain, 'auxiliary_protein': '',
                        #                         'ligand': {'name': 'None', 'pubchemId': 'None', 'title': 'None', 'role': '.nan', 'type': 'None'}, 'signaling_protein': 'None', 'state': 'Inactive'}
                        #     auxiliary_proteins, ligands = [], []
                        #     if pdb_data_dict['ligands']!='None':
                        #         for key, values in pdb_data_dict['ligands'].items():
                        #             if key in ['SO4','NA','CLR','OLA','OLB','OLC','TAR','NAG','EPE','BU1','ACM','GOL','PEG','PO4','TLA','BOG','CIT','PLM','BMA','MAN','MLI','PGE','SIN','PGO','MES','ZN','NO3','NI','MG','PG4']:
                        #                 continue
                        #             else:
                        #                 ligands.append({'name': key, 'pubchemId': 'None', 'title': pdb_data_dict['ligands'][key]['comp_name'], 'role': '.nan', 'type': 'None'})
                        #         sy = StructureYaml(s+'.yaml')
                        #         bril, by = sy.check_aux_protein('BRIL')
                        #         t4, ty = sy.check_aux_protein('T4-Lysozyme')
                        #         if bril:
                        #             auxiliary_proteins.append('BRIL')
                        #         if t4:
                        #             auxiliary_proteins.append('T4-Lysozyme')
                        #         for key, values in pdb_data_dict['auxiliary'].items():
                        #             if pdb_data_dict['auxiliary'][key]['subtype'] in ['Expression tag', 'Linker']:
                        #                 continue
                        #             else:
                        #                 if pdb_data_dict['auxiliary'][key]['subtype']=='Soluble cytochrome b562':
                        #                     aux_p = 'BRIL'
                        #                 elif pdb_data_dict['auxiliary'][key]['subtype'] in ['Endolysin','T4-Lysozyme']:
                        #                     aux_p = 'T4-Lysozyme'
                        #                 else:
                        #                     aux_p = pdb_data_dict['auxiliary'][key]['subtype']
                        #                 if aux_p not in auxiliary_proteins:
                        #                     auxiliary_proteins.append(aux_p)
                        #         for key, values in pdb_data_dict['construct_sequences'].items():
                        #             if key!=protein.entry_name and key not in struct_yaml_dict['auxiliary_protein']:
                        #                 if 'arrestin' in key:
                        #                     struct_yaml_dict['signaling_protein'] = key
                        #         if len(auxiliary_proteins)>1:
                        #             struct_yaml_dict['auxiliary_protein'] = ', '.join(auxiliary_proteins)
                        #         if len(ligands)>1:
                        #             struct_yaml_dict['ligand'] = ligands
                        #     yaml.dump(struct_yaml_dict, structure_file, indent=4, default_flow_style=False)

                        # # Build residue table for structure
                        # build_structure_command = shlex.split('/env/bin/python3 manage.py build_structures -f {}.yaml'.format(pdb_code.index))
                        # subprocess.call(build_structure_command)

                        # # Check state
                        # struct = Structure.objects.get(pdb_code__index=pdb_code.index)
                        # pi = PdbStateIdentifier(struct)
                        # pi.run()
                        # if pi.state!=None:
                        #     Structure.objects.filter(pdb_code__index=pdb_code.index).update(state=pi.state)
                        #     print(pi.state, pi.activation_value)
                        #     with open('../../data/protwis/gpcr/structure_data/structures/{}.yaml'.format(pdb_code.index), 'r') as yf:
                        #         struct_yaml = yaml.load(yf, Loader=yaml.FullLoader)
                        #     struct_yaml['state'] = pi.state.name
                        #     try:
                        #         struct_yaml['distance'] = round(float(pi.activation_value), 2)
                        #     except:
                        #         struct_yaml['distance'] = None
                        #     with open('../../data/protwis/gpcr/structure_data/structures/{}.yaml'.format(pdb_code.index), 'w') as struct_yaml_file:
                        #         yaml.dump(struct_yaml, struct_yaml_file, indent=4, default_flow_style=False)

                        # # Check sodium pocket
                        # new_prot_conf.sodium_pocket()

                        print('{} added to db (preferred_chain chain: {})'.format(s, preferred_chain))
                except Exception as msg:
                    print(s, msg)
                    self.error_list.append(s)
                    del self.db_list[self.db_list.index(s)]
                    missing_from_db = False
                    del self.csv_list[self.csv_list.index(s)]



    def pdb_request_by_uniprot(self, uniprot_id):
        url = 'https://www.rcsb.org/pdb/rest/search'

        queryText = """
<orgPdbQuery>
    <queryType>org.pdb.query.simple.UpAccessionIdQuery</queryType>
    <description>Simple query for a list of UniprotKB Accession IDs: {}</description>
    <accessionIdList>{}</accessionIdList>
</orgPdbQuery>

        """.format(uniprot_id, uniprot_id)
        req = urllib.request.Request(url, data=bytes(queryText, 'utf-8'))
        f = urllib.request.urlopen(req)
        result = f.read()
        structures = [i.split(':')[0] for i in result.decode('utf-8').split('\n')[:-1]]
        return structures

    def pdb_request_by_pdb(self, pdb, request_type):
        data = {}
        response = urlopen('https://data.rcsb.org/rest/v1/core/entry/{}'.format(pdb))
        json_data = json.loads(response.read())
        response.close()
        if request_type=='entry':
            data['method'] = json_data['exptl'][0]['method']
            if data['method'].startswith("THEORETICAL") or data['method'] in ['SOLUTION NMR','SOLID-STATE NMR']:
                return False
            if 'pubmed_id' in json_data['rcsb_entry_container_identifiers']:
                data['pubmedId'] = json_data['rcsb_entry_container_identifiers']['pubmed_id']
            else:
                data['pubmedId'] = None
            return True
        elif request_type=='polymer_entity':
            entity_list = json_data['rcsb_entry_container_identifiers']['entity_ids']
            uniprot_ids = []
            for i in entity_list:
                try:
                    response2 = urlopen('https://data.rcsb.org/rest/v1/core/polymer_entity/{}/{}'.format(pdb, i))
                except urllib.error.HTTPError:
                    continue
                json_data2 = json.loads(response2.read())
                response2.close()
                try:
                    uniprot_ids+=json_data2['rcsb_polymer_entity_container_identifiers']['uniprot_ids']
                except KeyError:
                    continue
            if len(uniprot_ids)>0:
                return uniprot_ids
            else:
                return False

    def pdb_request_by_pdb_deprecated(self, pdb_code):
        response = urllib.request.urlopen('https://www.rcsb.org/pdb/rest/describePDB?structureId={}'.format(pdb_code.lower()))
        response_mol = urllib.request.urlopen('https://www.rcsb.org/pdb/rest/describeMol?structureId={}'.format(pdb_code.lower()))
        str_des = str(response.read())
        dic = xmltodict.parse(response_mol.read())
        if 'NMR' in str_des or 'extracellular' in str_des:
            return 0
        if pdb_code in ['4QXE','1XWD','4QXF','4MQW','6B7H','6BSZ','6BT5','5OTW','3G04','3KS9','4XAQ','5II0','6N4X']:
            return 0
        polymer = dic['molDescription']['structureId']['polymer']
        description = ''
        if type(polymer)==type([]):
            for mol in polymer:
                if 'receptor' in mol['polymerDescription']['@description'] or 'Rhodopsin' in mol['polymerDescription']['@description'] or 'Smoothened' in mol['polymerDescription']['@description']:
                    description = mol['polymerDescription']['@description']
                if description=='' or int(mol['@length'])<100:
                    pass
                else:
                    try:
                        if polymer['macroMolecule'][0]['accession']['@id'] in self.uniprots:
                            return 1
                        else:
                            raise Exception()
                    except:
                        try:
                            for m in mol['macroMolecule']:
                                try:
                                    if mol['macroMolecule']['accession']['@id'] in self.uniprots:
                                        return 1
                                except:
                                    try:
                                        if m['accession']['@id'] in self.uniprots:
                                            return 1
                                    except:
                                        pass
                            raise Exception()
                        except:
                            pass
            return 0
        else:
            if 'receptor' in polymer['polymerDescription']['@description'] or 'Rhodopsin' in polymer['polymerDescription']['@description'] or 'Smoothened' in polymer['polymerDescription']['@description'] or 'Frizzled' in polymer['polymerDescription']['@description']:
                if int(polymer['@length'])<100:
                    return 0
                if type(polymer['macroMolecule'])==type([]):
                    for mM in polymer['macroMolecule']:
                        if mM['accession']['@id'] in self.uniprots:
                            return 1
                else:
                    if polymer['macroMolecule']['accession']['@id'] in self.uniprots:
                        return 1
                    else:
                        return 0
            else:
                return 0


class QueryPDBClassifiedGPCR():
    ''' Queries PDB using GPCRdb protein and structure entries using the 'G protein-coupled receptors' classification on PDB. Tree node number (<n>248</n>)
        need to be updated when new xtals come in.
    '''
    def __init__(self):
        self.num_struct = None
        self.new_structures = None
        self.new_uniques = None

    def list_xtals(self, verbose=True):
        ''' Lists structures and matching receptors from PDB that are not on GPCRdb yet. '''

        url = 'https://www.rcsb.org/pdb/rest/search'

        queryText = """
<orgPdbQuery>
    <version>head</version>
    <queryType>org.pdb.query.simple.TreeQuery</queryType>
    <description>TransmembraneTree Search for G Protein-Coupled Receptors (GPCRs)</description>
    <t>19</t>
    <n>248</n>
    <nodeDesc>G Protein-Coupled Receptors (GPCRs)</nodeDesc>
</orgPdbQuery>

        """
        req = urllib.request.Request(url, data=bytes(queryText, 'utf-8'))
        f = urllib.request.urlopen(req)
        result = f.read()
        structures = result.decode('utf-8').split('\n')[:-1]
        if len(structures)<159:
            raise AssertionError('Less than 159 structures, change the pdb query.')
        if verbose:
            print('Number of GPCR structures on PDB:',len(structures))
        new_struct = []
        new_unique = []
        for i in structures:
            response = urllib.request.urlopen('https://www.rcsb.org/pdb/rest/describeMol?structureId={}'.format(i.lower()))
            response_des = urllib.request.urlopen('https://www.rcsb.org/pdb/rest/describePDB?structureId={}'.format(i.lower()))
            str_text = str(response.read())
            str_des = str(response_des.read())
            if 'NMR' in str_des:
                continue
            if 'extracellular' in str_des:
                continue
            if i=='1EDN':
                continue
            uniprots = re.findall('accession id="([A-Z0-9]+)"', str_text)
            try:
                s = Protein.objects.get(entry_name=i.lower())
                continue
            except:
                new_struct.append((i, uniprots))
            miss_count = 0
            for j in uniprots:
                try:
                    p = Protein.objects.get(accession=j)
                    try:
                        parent = Protein.objects.filter(parent=p)
                        if len(parent)==0:
                            raise Exception()
                    except:
                        miss_count+=1

                except:
                    pass
            if miss_count==1:
                new_unique.append((i,p))
        if verbose:
            print('\nStructures not on GPCRdb: ',len(new_struct),'\n',new_struct)
            print('\nNew unique structures: ',len(new_unique),'\n',new_unique)
        self.num_struct = len(structures)
        self.new_structures = new_struct
        self.new_uniques = new_unique


def yamls_to_csv():
    s_dir = os.sep.join([settings.DATA_DIR, 'structure_data', 'structures'])
    c_dir = os.sep.join([settings.DATA_DIR, 'structure_data', 'constructs'])
    g_dir = os.sep.join([settings.DATA_DIR, 'g_protein_data'])
    yamls = os.listdir(s_dir)
    constructs = os.listdir(s_dir)
    d = OrderedDict()
    for i in yamls:
        pdb = i.split('.')[0]
        try:
            s_obj = Structure.objects.get(pdb_code__index=pdb)
        except Structure.DoesNotExist:
            continue
        with open(os.sep.join([s_dir, i]), 'r') as f1:
            s_y = yaml.load(f1, Loader=yaml.FullLoader)
            d[pdb] = s_y
        with open(os.sep.join([c_dir, i]), 'r') as f2:
            c_y = yaml.load(f2, Loader=yaml.FullLoader)
            d[pdb]['protein'] = c_y['protein']
        d[pdb]['obj'] = s_obj
        if type(d[pdb]['ligand'])!=type([]):
            d[pdb]['ligand'] = [d[pdb]['ligand']]
    # Order by pub date
    ordered_structs = Structure.objects.filter(pdb_code__index__in=d.keys()).order_by('publication_date', 'pdb_code__index__in').values_list('pdb_code__index', flat=True)
    temp_d = OrderedDict()
    for s in ordered_structs:
        temp_d[s] = d[s]
    d = temp_d
    # structures.csv
    with open(os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation', 'structures.csv']), 'w', newline='') as s_csv:
        struct_w = csv.writer(s_csv, delimiter=',', quotechar="'", quoting=csv.QUOTE_MINIMAL)
        struct_w.writerow(['PDB', 'Receptor_UniProt', 'Method', 'Resolution', 'State', 'ChainID', 'Note', 'Date'])
        for pdb, vals in d.items():
            if vals['obj'].structure_type.name.startswith('X-ray'):
                method = 'X-ray'
            elif vals['obj'].structure_type.name=='Electron microscopy':
                method = 'cryo-EM'
            else:
                method = vals['obj'].structure_type.name
            struct_w.writerow([pdb, vals['protein'], method, vals['obj'].resolution, vals['state'], vals['preferred_chain'], '', vals['obj'].publication_date])
    # ligands.csv
    with open(os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation', 'ligands.csv']), 'w', newline='') as l_csv:
        lig_w = csv.writer(l_csv, delimiter='\t', quotechar="'", quoting=csv.QUOTE_MINIMAL)
        lig_w.writerow(['PDB', 'ChainID', 'Name', 'PubChemID', 'Role', 'Title', 'Type'])
        for pdb, vals in d.items():
            lig = d[pdb]['ligand']
            if isinstance(lig, list):
                for l in lig:
                    if 'chain' in l:
                        chain = l['chain']
                    else:
                        chain = ''
                    if 'title' not in l:
                        title = l['name']
                    else:
                        title = l['title']
                    lig_w.writerow([pdb, chain, l['name'], l['pubchemId'], l['role'], title, l['type']])
    fusion_prots = OrderedDict()
    ramp, grk = OrderedDict(), OrderedDict()
    # nanobodies.csv fusion_proteins.csv ramp.csv grk.csv
    with open(os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation', 'nanobodies.csv']), 'w', newline='') as n_csv:
        nb_w = csv.writer(n_csv, delimiter=',', quotechar="'", quoting=csv.QUOTE_MINIMAL)
        nb_w.writerow(['PDB', 'Name'])
        for pdb, vals in d.items():
            auxs = d[pdb]['auxiliary_protein'].split(',')
            for aux in auxs:
                if aux=='' or aux=='None' or aux.startswith('GABA'):
                    continue
                if aux[0]==' ':
                    aux = aux[1:]
                if aux.startswith('Antibody') or aux.startswith('scFv') or 'Fab' in aux or 'Nanobody' in aux or aux.startswith('Camelid') or aux.startswith('IgG') or aux in ['Sb51','DN13','Anti-RON nanobody','T-cell surface glycoprotein CD4']:
                    if aux.startswith('Nanobody '):
                        aux = aux.split(' ')[0]+'-'+aux.split(' ')[1]
                    elif aux.startswith('Nanobody') and '-' not in aux and aux!='Nanobody':
                        aux = 'Nanobody-'+aux[8:]
                    nb_w.writerow([pdb, aux])
                elif aux in ['BRIL', 'T4-Lysozyme', 'Flavodoxin', 'Rubredoxin', 'GlgA glycogen synthase', 'Glycogen synthase', 'Thioredoxin 1', 'TrxA', 'Sialidase NanA']:
                    fusion_prots[pdb] = aux
                elif aux.startswith('RAMP'):
                    ramp[pdb] = aux
                elif aux.startswith('GRK'):
                    grk[pdb] = aux
                else:
                    print('====',aux)
    with open(os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation', 'fusion_proteins.csv']), 'w', newline='') as f_csv:
        f_w = csv.writer(f_csv, delimiter=',', quotechar="'", quoting=csv.QUOTE_MINIMAL)
        f_w.writerow(['PDB', 'Name'])
        for pdb, name in fusion_prots.items():
            f_w.writerow([pdb, name])
    with open(os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation', 'ramp.csv']), 'w', newline='') as r_csv:
        r_w = csv.writer(r_csv, delimiter=',', quotechar="'", quoting=csv.QUOTE_MINIMAL)
        r_w.writerow(['PDB', 'Name'])
        for pdb, name in ramp.items():
            r_w.writerow([pdb, name])
    with open(os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation', 'grk.csv']), 'w', newline='') as g_csv:
        g_w = csv.writer(g_csv, delimiter=',', quotechar="'", quoting=csv.QUOTE_MINIMAL)
        g_w.writerow(['PDB', 'Name'])
        for pdb, name in grk.items():
            g_w.writerow([pdb, name])
    # g proteins
    with open(os.sep.join([g_dir, 'complex_model_templates.yaml']), 'r') as f2:
        gprots = yaml.load(f2, Loader=yaml.FullLoader)
    with open(os.sep.join([settings.DATA_DIR, 'structure_data', 'extra_protein_notes.yaml']), 'r') as f3:
        extra = yaml.load(f3, Loader=yaml.FullLoader)
    arrestin = OrderedDict()
    # # Order g proteins and extra by pub date
    # temp_gprots, temp_extra = OrderedDict(), OrderedDict()
    # ordered_g_structs = Structure.objects.filter(pdb_code__index__in=gprots.keys()).order_by('publication_date').values_list('pdb_code__index', flat=True)
    # ordered_extra = Structure.objects.filter(pdb_code__index__in=extra.keys()).order_by('publication_date').values_list('pdb_code__index', flat=True)
    # for s in ordered_g_structs:
    #     temp_gprots[s] =  gprots[s]
    # for s in ordered_extra:
    #     temp_extra[s] = extra[s]
    # gprots, extra = temp_gprots, temp_extra
    with open(os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation', 'g_proteins.csv']), 'w', newline='') as gp_csv:
        gp_w = csv.writer(gp_csv, delimiter=',', quotechar="'", quoting=csv.QUOTE_MINIMAL)
        gp_w.writerow(['PDB', 'Alpha_UniProt', 'Alpha_ChainID', 'Beta_UniProt', 'Beta_ChainID', 'Gamma_UniProt', 'Gamma_ChainID', 'Note'])
        gprot_pdbs = []
        for alpha, vals_l in gprots.items():
            for vals in vals_l:
                note = ''
                if vals['pdb'] in extra:
                    note = extra[vals['pdb']]['note']
                if vals['beta']['protein']=='None':
                    beta_prot, beta_chain = '', ''
                else:
                    beta_prot, beta_chain = vals['beta']['protein'], vals['beta']['chain']
                if vals['gamma']['protein']=='None':
                    gamma_prot, gamma_chain = '', ''
                else:
                    gamma_prot, gamma_chain = vals['gamma']['protein'], vals['gamma']['chain']
                gprot_pdbs.append(vals['pdb'])
                gp_w.writerow([vals['pdb'], alpha, vals['alpha'], beta_prot, beta_chain, gamma_prot, gamma_chain, note])
        prot_code = {'GNAS2':'gnas2_human','GNAT1':'gnat1_bovin','GNAT3':'gnat3_bovin'}
        for pdb, vals in extra.items():
            if 'category' in vals:
                if vals['category']=='G alpha':
                    gp_w.writerow([pdb, vals['prot'], vals['chain'], '', '', '', '', vals['note']])
                elif vals['category']=='Arrestin':
                    arrestin[pdb] = vals
    with open(os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation', 'arrestins.csv']), 'w', newline='') as ar_csv:
        ar_w = csv.writer(ar_csv, delimiter=',', quotechar="'", quoting=csv.QUOTE_MINIMAL)
        ar_w.writerow(['PDB', 'UniProt', 'ChainID', 'Note'])
        for pdb, vals in arrestin.items():
            ar_w.writerow([pdb, vals['prot'], vals['chain'], vals['note']])

