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
from structure.functions import PdbChainSelector, PdbStateIdentifier
from structure.management.commands.structure_yaml_editor import StructureYaml
from construct.functions import *
from common.models import WebResource, WebLink, Publication

import Bio.PDB as PDB
from datetime import datetime
import urllib
import re
import os
import xmltodict
import yaml
import shlex
import subprocess


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
            self.yamls = self.get_all_yamls()
            self.prepare_input(options['proc'], self.uniprots)

    def main_func(self, positions, iteration, count, lock):
        if not positions[1]:
            uniprot_list = self.uniprots[positions[0]:]
        else:
            uniprot_list = self.uniprots[positions[0]:positions[1]]
        q = QueryPDB(self.uniprots, self.yamls)
        consider_list = []
        for uni in uniprot_list:
            q.new_xtals(uni)
            for i in q.consider_list:
                if i not in consider_list:
                    consider_list.append(i)
        if self.verbose:
            print('Missing from db: ', q.db_list)
            print('Missing yamls: ', q.yaml_list)
            print('Structures with missing x50s: {} structures {}'.format(len(consider_list), consider_list))

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

    def get_all_yamls(self):
        yamls = [i.split('.')[0] for i in os.listdir('/protwis/data/protwis/gpcr/structure_data/structures/')]
        return yamls


class QueryPDB():
    ''' Queries PDB using GPCRdb protein and structure entries. If those are not available, it uses the structure and uniprot data folders.
    '''
    def __init__(self, uniprots, yamls):
        self.exceptions = []
        self.uniprots = uniprots
        self.yamls = yamls
        self.db_list, self.yaml_list = [], []
        self.missing_x50_list = ['4KNG','3N7P','3N7R','3N7S','4HJ0','6DKJ','5OTT','5OTU','5OTV','5OTX','6GB1']
        self.missing_x50_exceptions = ['6TPG','6TPJ']
        self.consider_list = []

    def new_xtals(self, uniprot):
        ''' List GPCR crystal structures missing from GPCRdb and the yaml files. Adds missing structures to DB.
        '''
        structs = self.pdb_request_by_uniprot(uniprot)
        try:
            protein = Protein.objects.get(accession=uniprot)
        except:
            protein = None
        try:
            x50s = Residue.objects.filter(protein_conformation__protein=protein,generic_number__label__in=['1x50','2x50','3x50','4x50','5x50','6x50','7x50'])
        except:
            x50s = None
        if structs!=['null']:
            for s in structs:
                missing_from_db, missing_yaml = False, False
                try:
                    st_obj = Structure.objects.get(pdb_code__index=s)
                except:
                    if s not in self.exceptions:
                        check = self.pdb_request_by_pdb(s)
                        if check==1:
                            self.db_list.append(s)
                            missing_from_db = True

                if s not in self.yamls and s not in self.exceptions:
                    if s not in self.db_list:
                        check = self.pdb_request_by_pdb(s)
                    else:
                        check = 1
                    if check==1:
                        self.yaml_list.append(s)
                        missing_yaml = True
                if not missing_from_db:
                    continue
                try:
                    pdb_data_dict = fetch_pdb_info(s, protein, new_xtal=True)
                    exp_method = pdb_data_dict['experimental_method']
                    if exp_method=='Electron Microscopy':
                        st_type = StructureType.objects.get(slug='electron-microscopy')
                    elif exp_method=='X-ray diffraction':
                        st_type = StructureType.objects.get(slug='x-ray-diffraction')
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
                                        del self.yaml_list[self.yaml_list.index(s)]
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
                                try:
                                    publication = Publication.objects.get(web_link__index=doi)
                                except Publication.DoesNotExist as e:
                                    p = Publication()
                                    try:
                                        p.web_link = WebLink.objects.get(index=doi, web_resource__slug='doi')
                                    except WebLink.DoesNotExist:
                                        wl = WebLink.objects.create(index=doi,
                                            web_resource = WebResource.objects.get(slug='doi'))
                                        p.web_link = wl
                                    p.update_from_doi(doi=doi)
                                    p.save()
                                    publication = p
                            elif pubmed!='':
                                try:
                                    publication = Publication.objects.get(web_link__index=pubmed)
                                except Publication.DoesNotExist as e:
                                    p = Publication()
                                    try:
                                        p.web_link = WebLink.objects.get(index=pubmed,
                                            web_resource__slug='pubmed')
                                    except WebLink.DoesNotExist:
                                        wl = WebLink.objects.create(index=pubmed,
                                            web_resource = WebResource.objects.get(slug='pubmed'))
                                        p.web_link = wl
                                    p.update_from_pubmed_data(index=pubmed)
                                    p.save()
                                    publication = p
                        except:
                            pass
                        pcs = PdbChainSelector(s, protein)
                        pcs.run_dssp()
                        preferred_chain = pcs.select_chain()

                        # Run state identification

                        # Create yaml files
                        with open(os.sep.join([settings.DATA_DIR, 'structure_data','constructs', '{}.yaml'.format(pdb_code.index)]), 'w') as construct_file:
                            yaml.dump({'name': pdb_code.index.lower(), 'protein': protein.entry_name}, construct_file, indent=4)
                        with open(os.sep.join([settings.DATA_DIR, 'structure_data','structures','{}.yaml'.format(pdb_code.index)]), 'w') as structure_file:
                            struct_yaml_dict = {'construct': pdb_code.index.lower(), 'pdb': pdb_code.index, 'preferred_chain': preferred_chain, 'auxiliary_protein': '',
                                                'ligand': {'name': 'None', 'pubchemId': 'None', 'title': 'None', 'role': '.nan', 'type': 'None'}, 'signaling_protein': 'None', 'state': 'Inactive'}
                            auxiliary_proteins, ligands = [], []
                            if pdb_data_dict['ligands']!='None':
                                for key, values in pdb_data_dict['ligands'].items():
                                    if key in ['SO4','NA','CLR','OLA','OLB','OLC','TAR','NAG','EPE','BU1','ACM','GOL','PEG','PO4','TLA','BOG','CIT','PLM','BMA','MAN','MLI','PGE','SIN','PGO','MES','ZN','NO3','NI','MG','PG4']:
                                        continue
                                    else:
                                        ligands.append({'name': key, 'pubchemId': 'None', 'title': pdb_data_dict['ligands'][key]['comp_name'], 'role': '.nan', 'type': 'None'})
                                sy = StructureYaml(s+'.yaml')
                                bril, by = sy.check_aux_protein('BRIL')
                                t4, ty = sy.check_aux_protein('T4-Lysozyme')
                                if bril:
                                    auxiliary_proteins.append('BRIL')
                                if t4:
                                    auxiliary_proteins.append('T4-Lysozyme')
                                for key, values in pdb_data_dict['auxiliary'].items():
                                    if pdb_data_dict['auxiliary'][key]['subtype'] in ['Expression tag', 'Linker']:
                                        continue
                                    else:
                                        if pdb_data_dict['auxiliary'][key]['subtype']=='Soluble cytochrome b562':
                                            aux_p = 'BRIL'
                                        elif pdb_data_dict['auxiliary'][key]['subtype'] in ['Endolysin','T4-Lysozyme']:
                                            aux_p = 'T4-Lysozyme'
                                        else:
                                            aux_p = pdb_data_dict['auxiliary'][key]['subtype']
                                        if aux_p not in auxiliary_proteins:
                                            auxiliary_proteins.append(aux_p)
                                for key, values in pdb_data_dict['construct_sequences'].items():
                                    if key!=protein.entry_name and key not in struct_yaml_dict['auxiliary_protein']:
                                        if 'arrestin' in key:
                                            struct_yaml_dict['signaling_protein'] = key
                                if len(auxiliary_proteins)>1:
                                    struct_yaml_dict['auxiliary_protein'] = ', '.join(auxiliary_proteins)
                                if len(ligands)>1:
                                    struct_yaml_dict['ligand'] = ligands
                            yaml.dump(struct_yaml_dict, structure_file, indent=4, default_flow_style=False)

                        # Build residue table for structure
                        build_structure_command = shlex.split('/env/bin/python3 manage.py build_structures -f {}.yaml'.format(pdb_code.index))
                        subprocess.call(build_structure_command)

                        # Check state
                        struct = Structure.objects.get(pdb_code__index=pdb_code.index)
                        pi = PdbStateIdentifier(struct)
                        pi.run()
                        if pi.state!=None:
                            Structure.objects.filter(pdb_code__index=pdb_code.index).update(state=pi.state)
                            print(pi.state, pi.activation_value)
                            with open('../../data/protwis/gpcr/structure_data/structures/{}.yaml'.format(pdb_code.index), 'r') as yf:
                                struct_yaml = yaml.load(yf, Loader=yaml.FullLoader)
                            struct_yaml['state'] = pi.state.name
                            try:
                                struct_yaml['distance'] = round(float(pi.activation_value), 2)
                            except:
                                struct_yaml['distance'] = None
                            with open('../../data/protwis/gpcr/structure_data/structures/{}.yaml'.format(pdb_code.index), 'w') as struct_yaml_file:
                                yaml.dump(struct_yaml, struct_yaml_file, indent=4, default_flow_style=False)

                        # Check sodium pocket
                        new_prot_conf.sodium_pocket()

                        print('{} added to db (preferred_chain chain: {})'.format(s, preferred_chain))
                except Exception as msg:
                    print(s, msg)
        


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

    def pdb_request_by_pdb(self, pdb_code):
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
