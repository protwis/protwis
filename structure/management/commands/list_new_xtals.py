# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 15:50:57 2016

@author: Gaspar Pandy
"""
from django.core.management.base import BaseCommand

from protein.models import Protein
from structure.models import Structure

import urllib
import re
import os
import xmltodict, pprint


class Command(BaseCommand):
        
    def add_arguments(self, parser):
        super(Command, self).add_arguments(parser=parser)
        parser.add_argument('--verbose', help='Print output', default=True, 
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
            q.list_xtals(self.verbose)
        else:
            q = QueryPDB()
            q.list_xtals(self.verbose)


class QueryPDB():
    ''' Queries PDB using GPCRdb protein and structure entries. If those are not available, it uses the structure and uniprot data folders.
    '''
    def __init__(self):
        self.exceptions = []
        try:
            self.uniprots = Protein.objects.filter(accession__isnull=False)
            if len(self.uniprots)<100:
                raise Exception()
        except:
            self.uniprots = [i.split('.')[0] for i in os.listdir('/protwis/data/protwis/gpcr/protein_data/uniprot/')]
        self.yamls = [i.split('.')[0] for i in os.listdir('/protwis/data/protwis/gpcr/structure_data/structures/')]

    def list_xtals(self, verbose=True):
        ''' List GPCR crystal structures missing from GPCRdb and the yaml files.
        '''      
        db_list, yaml_list = [], []
        for u in self.uniprots:
            structs = self.pdb_request_by_uniprot(u.accession)
            if structs!=['null']:
                for s in structs:
                    try:
                        st_obj = Structure.objects.get(pdb_code__index=s)
                    except:
                        if s not in self.exceptions:
                            check = self.pdb_request_by_pdb(s)
                            if check==1:
                                db_list.append(s)
                    if s not in self.yamls and s not in self.exceptions:
                        if s not in db_list:
                            check = self.pdb_request_by_pdb(s)
                        else:
                            check = 1
                        if check==1:
                            yaml_list.append(s)
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
                if 'receptor' in mol['polymerDescription']['@description']:
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
            if 'receptor' in polymer['polymerDescription']['@description'] and int(polymer['@length'])>100:
                try:
                    if polymer['macroMolecule'][0]['accession']['@id'] in self.uniprots:
                        return 1
                except:
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
