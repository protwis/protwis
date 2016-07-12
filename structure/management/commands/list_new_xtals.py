# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 15:50:57 2016

@author: Gaspar Pandy
"""
from django.core.management.base import BaseCommand

from protein.models import Protein

import urllib
import re


class Command(BaseCommand):
        
    def handle(self, *args, **options):
        q = QueryPDB()
        q.list_xtals()


class QueryPDB():    
    def __init__(self):
        self.num_struct = None
        self.new_structures = None
        self.new_uniques = None
    
    def list_xtals(self):
        ''' Lists structures and matching receptors from PDB that are not on GPCRdb yet. '''
    
        url = 'http://www.rcsb.org/pdb/rest/search'
    
        queryText = """
<orgPdbQuery>
    <version>head</version>
    <queryType>org.pdb.query.simple.TreeQuery</queryType>
    <description>TransmembraneTree Search for G Protein-Coupled Receptors (GPCRs)</description>
    <t>19</t>
    <n>236</n>
    <nodeDesc>G Protein-Coupled Receptors (GPCRs)</nodeDesc>
</orgPdbQuery>
    
        """   
        req = urllib.request.Request(url, data=bytes(queryText, 'utf-8'))
        f = urllib.request.urlopen(req)
        result = f.read()
        structures = result.decode('utf-8').split('\n')[:-1]
        if len(structures)<159:
            raise AssertionError('Less than 159 structures, change the pdb query.')
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
            uniprots = re.findall('accession id="([A-Z0-9]+)"', str_text)
            try:
                s = Protein.objects.get(entry_name=i.lower())
                continue
            except:
                new_struct.append(i)
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
        print('\nStructures not on GPCRdb: ',len(new_struct),'\n',new_struct)
        print('\nNew unique structures: ',len(new_unique),'\n',new_unique)
        self.num_struct = len(structures)
        self.new_structures = new_struct
        self.new_uniques = new_unique