# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 15:50:57 2016

@author: Gaspar Pandy
"""
from django.core.management.base import BaseCommand

from protein.models import Protein, ProteinConformation, ProteinSequenceType, ProteinSource, ProteinState
from residue.models import Residue
from structure.models import Structure, PdbData, StructureType
from structure.sequence_parser import SequenceParser
from structure.functions import PdbChainSelector
from construct.functions import *
from common.models import WebResource, WebLink, Publication

import Bio.PDB as PDB
from datetime import datetime
import urllib
import re
import os
import xmltodict


class Command(BaseCommand):
        
    def add_arguments(self, parser):
        super(Command, self).add_arguments(parser=parser)
        parser.add_argument('--verbose', help='Print output', default=False, 
                            action='store_true')
        parser.add_argument('--classified', help="Use PDB's 'G protein-coupled receptors' classification", default=False, 
                            action='store_true')

    def handle(self, *args, **options):
        if options['verbose']:
            self.verbose = True
        else:
            self.verbose = False
        if options['classified']:
            q = QueryPDBClassifiedGPCR()
            q.new_xtals(self.verbose)
        else:
            q = QueryPDB()
            q.new_xtals(self.verbose)


class QueryPDB():
    ''' Queries PDB using GPCRdb protein and structure entries. If those are not available, it uses the structure and uniprot data folders.
    '''
    def __init__(self):
        self.exceptions = []
        try:
            self.uniprots = [i.accession for i in Protein.objects.filter(accession__isnull=False)]
            if len(self.uniprots)<100:
                raise Exception()
        except:
            self.uniprots = [i.split('.')[0] for i in os.listdir('/protwis/data/protwis/gpcr/protein_data/uniprot/')]
        self.yamls = [i.split('.')[0] for i in os.listdir('/protwis/data/protwis/gpcr/structure_data/structures/')]

    def new_xtals(self, verbose=False):
        ''' List GPCR crystal structures missing from GPCRdb and the yaml files.
        '''
        db_list, yaml_list = [], []
        # sts = [i.protein_conformation.protein.parent.accession for i in Structure.objects.all()]
        # unis = []
        # for st in sts:
        #     if st not in unis:
        #         unis.append(st)
        # self.uniprots = unis
        # self.uniprots = ['P02699']
        for u in self.uniprots:
            structs = self.pdb_request_by_uniprot(u)
            try:
                protein = Protein.objects.get(accession=u)
            except:
                protein = None
            try:
                x50s = Residue.objects.filter(protein_conformation__protein=protein,generic_number__label__in=['1x50','2x50','3x50','4x50','5x50','6x50','7x50'])
            except:
                x50s = None
            if structs!=['null']:
                for s in structs:
                    missing_from_db = False
                    try:
                        st_obj = Structure.objects.get(pdb_code__index=s)
                    except:
                        if s not in self.exceptions:
                            check = self.pdb_request_by_pdb(s)
                            if check==1:
                                db_list.append(s)
                                missing_from_db = True
                    if s not in self.yamls and s not in self.exceptions:
                        if s not in db_list:
                            check = self.pdb_request_by_pdb(s)
                        else:
                            check = 1
                        if check==1:
                            yaml_list.append(s)
                    if not missing_from_db:
                        continue
                    try:
                        pdb_data_dict = fetch_pdb_info(s, protein)
                        for d in pdb_data_dict['deletions']:
                            presentx50s = []
                            for x in x50s:
                                if not d['start']<x.sequence_number<d['end']:
                                    presentx50s.append(x)                                    
                            # Filter out ones without all 7 x50 positions present in the xtal
                            if len(presentx50s)!=7:
                                try:
                                    del db_list[db_list.index(s)]
                                    missing_from_db = False
                                    del yaml_list[yaml_list.index(s)]
                                except:
                                    pass
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
                            new_prot_conf, created = ProteinConformation.objects.get_or_create(protein=new_prot, state=state, template_structure=None)
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
                            os.remove('./pdb{}.ent'.format(s).lower())

                            # Create new structure object
                            Structure.objects.get_or_create(preferred_chain=preferred_chain, resolution=resolution, publication_date=publication_date, representative='f', pdb_code=pdb_code,
                                                            pdb_data=pdb_data, protein_conformation=new_prot_conf, publication=publication, state=state, 
                                                            structure_type=StructureType.objects.get(slug='x-ray-diffraction'))
                            print('{} added to db (preferred_chain chain: {})'.format(s, preferred_chain))
                    
                    except Exception as msg:
                        print(msg)
        if verbose:
            print('Missing from db: ', db_list)
            print('Missing yamls: ', yaml_list)

        return db_list, yaml_list

    def pdb_request_by_uniprot(self, uniprot_id):
        url = 'http://www.rcsb.org/pdb/rest/search'

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
        response = urllib.request.urlopen('http://www.rcsb.org/pdb/rest/describePDB?structureId={}'.format(pdb_code.lower()))
        response_mol = urllib.request.urlopen('http://www.rcsb.org/pdb/rest/describeMol?structureId={}'.format(pdb_code.lower()))
        str_des = str(response.read())
        dic = xmltodict.parse(response_mol.read())
        if 'NMR' in str_des or 'extracellular' in str_des:
            return 0
        if pdb_code=='AAAA':
            return 0
        polymer = dic['molDescription']['structureId']['polymer']
        if type(polymer)==type([]):
            for mol in polymer:
                if 'receptor' in mol['polymerDescription']['@description'] or 'Rhodopsin' in mol['polymerDescription']['@description']:
                    if int(mol['@length'])<100:
                        return 0
                    else:
                        try:
                            if polymer['macroMolecule'][0]['accession']['@id'] in self.uniprots:
                                return 1
                            else:
                                raise Exception()
                        except:
                            return 0
        else:
            if 'receptor' in polymer['polymerDescription']['@description'] or 'Rhodopsin' in polymer['polymerDescription']['@description']:
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
    
        url = 'http://www.rcsb.org/pdb/rest/search'
    
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
            response = urllib.request.urlopen('http://www.rcsb.org/pdb/rest/describeMol?structureId={}'.format(i.lower()))
            response_des = urllib.request.urlopen('http://www.rcsb.org/pdb/rest/describePDB?structureId={}'.format(i.lower()))
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
