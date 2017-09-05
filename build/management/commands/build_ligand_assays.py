from django.core.management.base import BaseCommand, CommandError
from build.management.commands.base_build import Command as BaseBuild
from django.utils.text import slugify
from django.conf import settings
from ligand.functions import get_or_make_ligand
from common.models import WebResource, WebLink
from protein.models import Protein
from ligand.models import Ligand, LigandProperities, LigandRole, LigandType, ChemblAssay, AssayExperiment
from ligand.models import LigandVendorLink, LigandVendors
from ligand.functions import get_or_make_ligand
from common.tools import fetch_from_web_api
from collections import defaultdict
import requests
from optparse import make_option
import logging
import shlex
import csv
import os
from collections import OrderedDict

class Command(BaseBuild):
    help = 'Reads source data and creates links to other databases'

    def add_arguments(self, parser):
        parser.add_argument('-p', '--proc',
            type=int,
            action='store',
            dest='proc',
            default=1,
            help='Number of processes to run')
        parser.add_argument('--filename', action='append', dest='filename',
            help='Filename to import. Can be used multiple times')

    logger = logging.getLogger(__name__)

    # source file directory
    links_data_dir = os.sep.join([settings.DATA_DIR, 'ligand_data', 'assay_data'])
    dictionary_file = os.sep.join([settings.DATA_DIR, 'ligand_data', 'assay_data', 'dictionary.txt'])
    defaults = {
        'name': 'ChEMBL_compound_ids',
        'url': 'https://www.ebi.ac.uk/chembl/compound/inspect/$index'
        }
    wr, created = WebResource.objects.get_or_create(slug='chembl_ligand', defaults = defaults)
    wr_pubchem = WebResource.objects.get(slug='pubchem')
        
    
    
    def handle(self, *args, **options):
        if options['filename']:
            filenames = options['filename']
        else:
            filenames = os.listdir(self.links_data_dir)

        #chembl_mol_ids,assay_ids = self.extract_chembl_mol_ids(filenames)
        data,header_dict = self.extract_chembl_mol_ids(filenames) #reading the entry file from the data filder.
        self.chembl_cid_dict = self.read_dict_file(self.dictionary_file)
        self.chembl_cid_dict = OrderedDict(self.chembl_cid_dict)
        #found_ids_dict = chembl_cid_dict #testing time. else commented.

        print("size of current dict",len(self.chembl_cid_dict))

        ##### find the chembl_molecular_ids ############
        self.chembl_mol_ids = set()
        for i in data:
            self.chembl_mol_ids.add(i[header_dict['molecule_chembl_id']])
            #print (i[header_dict['pchembl_value']])
        ######## end ##############
        print("amount of ids in dataset",len(self.chembl_mol_ids))
            
        ##notfound_ids_set = self.read_notfound(filenames) # testing
        # chembl_cid_dict, notfound_ids_set = self.find_cid(chembl_mol_ids, chembl_cid_dict)
        # chembl_cid_dict, notfound_ids_set = self.call_pubchem_service(chembl_cid_dict, notfound_ids_set)
        # self.write_dictionary_file(chembl_cid_dict)
        ###print (len(chembl_mol_ids))

        self.prepare_input(options['proc'], self.chembl_mol_ids)

        # self.create_assays(data,header_dict)
        # self.create_assaysExperiment(data,header_dict)
        
        

    
    def create_assaysExperiment(self, data, header):
        for record in data:

            target = record[header['target_chembl_id']]
            assay = record[header['assay_chembl_id']]
            ligand =record[header['molecule_chembl_id']]
            p = Protein.objects.get(web_links__index = target, web_links__web_resource__slug = 'chembl')
            ls = Ligand.objects.filter(properities__web_links__index=ligand, properities__web_links__web_resource__slug = 'chembl_ligand', canonical=True)
            for l in ls:
                if len(ls)>1:
                    print('issue with canonical! give to munk',l,l.pk,ligand)
                    break
            assay_experiments = AssayExperiment.objects.filter( protein=p, ligand=l, assay__assay_id=assay)
            
            if assay_experiments.exists():
                assay_experiment = assay_experiments.get()
            else:
                assay_experiment = AssayExperiment()
                assay_experiment.assay = ChemblAssay.objects.get(web_links__index = assay, web_links__web_resource__slug = 'chembl_assays' )
                assay_experiment.ligand = l
                assay_experiment.protein = p
            
            
            assay_experiment.assay_type = record[header['assay_type']]
            assay_experiment.pchembl_value = record[header['pchembl_value']]
            assay_experiment.assay_description = record[header['assay_description']]
            assay_experiment.published_value = record[header['published_value']]
            assay_experiment.published_relation = record[header['published_relation']]
            assay_experiment.published_type = record[header['published_type']]
            assay_experiment.published_units = record[header['published_units']]
            
            assay_experiment.standard_value = record[header['standard_value']]
            assay_experiment.standard_relation = record[header['standard_relation']]
            assay_experiment.standard_type = record[header['standard_type']]
            assay_experiment.standard_units = record[header['standard_units']]
            
            assay_experiment.save()
  
    
    def create_assays(self, data, header_dict):
        
        ####find the Chembl_assay_ids in a set - nonredundant lis
        assay_ids = set()
        for row in data:
            assay_ids.add(row[header_dict['assay_chembl_id']])
            
        ####create webResource for chembl assays
        defaults = {
        'name': 'ChEMBL Assays ids',
        'url': 'https://www.ebi.ac.uk/chembl/assay/inspect/$index'
        }
        wr, created = WebResource.objects.get_or_create(slug='chembl_assays', defaults = defaults)

        ####Create an ChemblAssay object (ca) with wiblink
        for chembl_assay_id in assay_ids:
            #print (chembl_assay_id)
            ca, created = ChemblAssay.objects.get_or_create(assay_id=chembl_assay_id)
            wl, created = WebLink.objects.get_or_create(index=chembl_assay_id, web_resource=wr)
            ca.web_links.add(wl)
            ca.save()
        return assay_ids


    def main_func(self, positions, iteration,count,lock):
        #####Create chembl compound link and connect it to the corresponding ligand/cid#####

        print(positions,iteration,count,lock)
        chembl_ids = self.chembl_mol_ids
        list_of_chembl_ids = list(chembl_ids)
        while count.value<len(chembl_ids):
            with lock:
                chembl_ligand = list_of_chembl_ids[count.value]
                count.value +=1 

            if chembl_ligand not in self.chembl_cid_dict.keys():
                cids, not_found = self.find_cid_for_chembl(chembl_ligand)
                if not_found:
                    print('NOT FOUND',chembl_ligand,cids)
                    continue
            else:
                cids = self.chembl_cid_dict[chembl_ligand]
           # cid = self.chembl_cid_dict[chembl_ligand]

            # cids = str(chembl_ligand[1])
            temp = str(cids).split(';') #perhaps we should load all of the CIDs
            #print (temp)
            
            cid = str(temp[0])
            print(count.value)
            #if cid!='3559':
            #    continue
            l = get_or_make_ligand(cid,'PubChem CID') #call the first cid if there are more than one
            #print (cid)
            wl, created = WebLink.objects.get_or_create(index=chembl_ligand, web_resource=self.wr)
            try:
                l.properities.web_links.add(wl)
                
            except AttributeError:
                print (cid) # to do get to the logger 
                continue
            
            ## Check if properities has web resource with pubchem and this cid // otherwise insert
                
            l.properities.save()
            print(cid)#148842 10775772 54477
            try:
                lp = LigandProperities.objects.get(web_links__index = cid, web_links__web_resource__slug = 'pubchem')
            except:
                # NO CID FOR LIGAND! Rare cases where SMILES was used for initial look up
                wl, created = WebLink.objects.get_or_create(index=cid, web_resource=self.wr_pubchem)
                l.properities.web_links.add(wl)
                lp = l.properities
            #print (lp)
            
        ###### Vendor stuff  ######
            cache_dir = ['pubchem', 'cid', 'vendors']
            url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/categories/compound/$index/JSON/'
            vendors = fetch_from_web_api(url, cid, cache_dir)
            
            if vendors:
                for vendor_data in vendors['SourceCategories']['Categories'][0]['Sources'] :
                    lv, created = LigandVendors.objects.get_or_create(slug = slugify(vendor_data['SourceName']))
                    lv.name = vendor_data['SourceName']
                    if 'SourceURL' in vendor_data:
                        lv.url = vendor_data['SourceURL']
                    lv.save()

                    if 'SID' in vendor_data:
                        #print (vendor_data['SID'])
                        lvls = LigandVendorLink.objects.filter(sid = vendor_data['SID'] )
                        if not lvls.exists():
                            lvl = LigandVendorLink()
                            lvl.vendor = lv
                            lvl.lp = l.properities
                            lvl.sid =  vendor_data['SID'] 
                            if 'RegistryID' in vendor_data:
                                lvl.vendor_external_id = vendor_data['RegistryID']
                            if 'SourceRecordURL' in vendor_data:
                                lvl.url = vendor_data['SourceRecordURL']
                            else:
                                continue
                            lvl.save()
 

    def write_dictionary_file(self, chembl_cid_dict):
        with open(self.dictionary_file, 'w') as foo:
#        with open('dict2', 'w') as foo:
            for k in chembl_cid_dict.keys():
                foo.write('{}\t{}\n'.format(k, chembl_cid_dict[k]))

    def find_cid_for_chembl(self, chembl_mol_id):
        cache_dir = ['ebi', 'chembl', 'src_compound_id_all']
        url = 'https://www.ebi.ac.uk/unichem/rest/src_compound_id_all/$index/1/22'
        lig_data = fetch_from_web_api(url, chembl_mol_id, cache_dir)
        # print("Searching for ",chembl_mol_id)
        not_found = False
        cid = False
        temp = []
        if not lig_data:
            #if not successful
            not_found = True
        else:
            try:
                for i, x  in enumerate(lig_data):
                    temp.append(lig_data[i]['src_compound_id'])
                if len(temp) > 1:
                    cid = ';'.join(temp)
                    self.add_cid_to_dict(chembl_mol_id,cid)
                elif len(temp) == 1 and temp[0] != '\n':
                    cid = temp[0] #lig_data[0]['src_compound_id'] 
                    self.add_cid_to_dict(chembl_mol_id,cid)
                else:
                    not_found = True
                    #print (chembl_mol_id)
            except KeyError:
                not_found = True
        if not cid:
            # not found
            cache_dir = ['pubchem', 'chembl', 'compound_name']
            url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/$index/json'
            lig_data = fetch_from_web_api(url, chembl_mol_id, cache_dir)
            if lig_data:
                try:
                    cid = lig_data['PC_Compounds'][0]['id']['id']['cid']
                    self.add_cid_to_dict(chembl_mol_id,cid)
                    not_found = False
                except KeyError:
                    not_found = True

        return cid,not_found

    def call_pubchem_service(self, chembl_cid_dict, notfound_ids_set):
        notfound = set()
        for chembl_id in notfound_ids_set:
            cache_dir = ['pubchem', 'chembl', 'compound_name']
            url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/$index/json'
            lig_data = fetch_from_web_api(url, chembl_id, cache_dir)
            try:
                cid = lig_data['PC_Compounds'][0]['id']['id']['cid']
                #updating the existing dictionary
                chembl_cid_dict[chembl_id] = cid
                self.add_cid_to_dict(chembl_mol_id,cid)
                #.discard(x)
            except KeyError:
                #print (chembl_id)
                notfound.add(chembl_id)


        return chembl_cid_dict, notfound 

    def add_cid_to_dict(self,chembl_id, cid):
        with open(self.dictionary_file, 'a') as myfile:
            myfile.write('{}\t{}\n'.format(chembl_id, cid))
    
    ##call pubchem service to find cids for chembl
    def find_cid(self, chembl_mol_ids, chembl_cid_dict):
        notfound = set()
        for chembl_mol_id in chembl_mol_ids:
            
            if chembl_mol_id not in chembl_cid_dict.keys():
                temp = []
                
                # url = 'https://www.ebi.ac.uk/unichem/rest/src_compound_id_all/'+chembl_mol_id+'/1/22'
                # response = requests.get(url)
                # lig_data = response.json()
                cache_dir = ['ebi', 'chembl', 'src_compound_id_all']
                url = 'https://www.ebi.ac.uk/unichem/rest/src_compound_id_all/$index/1/22'
                lig_data = fetch_from_web_api(url, chembl_mol_id, cache_dir)
                print("Searching for ",chembl_mol_id,len(chembl_mol_ids))
                if not lig_data:
                    #if not successful
                    notfound.add(chembl_mol_id)
                    continue
                try:
                    for i, x  in enumerate(lig_data):
                        temp.append(lig_data[i]['src_compound_id'])
                    if len(temp) > 1:
                        cid = ';'.join(temp)
                        chembl_cid_dict[chembl_mol_id] = cid
                        self.add_cid_to_dict(chembl_mol_id,cid)
                    elif len(temp) == 1 and temp[0] != '\n':
                        cid = temp[0] #lig_data[0]['src_compound_id'] 
                        #updating the existing dictionary
                        chembl_cid_dict[chembl_mol_id] = cid
                        self.add_cid_to_dict(chembl_mol_id,cid)
                    else:
                        notfound.add(chembl_mol_id)
                        #print (chembl_mol_id)
                except KeyError:
                    notfound.add(chembl_mol_id)
                    #print (chembl_mol_id)
            elif chembl_mol_id in chembl_cid_dict:
                continue
            else:
                #raise KeyError:
                print (chembl_mol_id) #to do put it to logfile where it should be
    
        return chembl_cid_dict, notfound # to do perhaps a redundant to have found and chembl_cid_dict.


#    #Read pregenerated dictionary file that contain chemblId\tCID
    def read_dict_file(self, dictionary_file):
        temp = []
        chembl_cid_dict = {}
        with open(dictionary_file, 'r') as fii:
            for line in fii:
                if line[0] != '\n':
                    try:
                        temp = line.rstrip().split()
                        chembl_cid_dict[temp[0]] = temp[1]
                    except IndexError:
                        print(line) #to do record to the log if tehre is some errors
        return chembl_cid_dict
#        return chembl_cid_dict


    ##read pre-generated file and extract the chembl_ids
    def extract_chembl_mol_ids(self, filenames=False):
        #self.logger.info('CREATING ASSAYS IMPORT LIGANDS START')
#        
#        # read source files
#        if not filenames:
#            filenames = os.listdir(self.links_data_dir)
        for filename in filenames:
            if filename=='dictionary.txt':
                continue
            full_file = os.sep.join([self.links_data_dir,filename])
            print (filename,full_file)
            header_dict = {}
            data = []
            chembl_assay_ids = set ()
            chembl_mol_ids = set()

            if filename[-2:]=='gz':
                import gzip, io
                with gzip.open(full_file, "r") as fii:
                    linereader = csv.reader(io.TextIOWrapper(fii, newline=""))
                    #ensure that i got the correct columns - the ones with the correct headers
                    for x,row in enumerate(linereader):
                        if x == 0:
                            for y,column in enumerate(row):
                                header_dict[column] = y
                            #print (header_dict)
                        
                        elif row[0] != '\n':
                            data.append(row)
                            #chembl_mol_ids.add(row[header_dict['molecule_chembl_id']])
                            #chembl_target_ids.add(row[header_dict['target_chembl_id']])
                            #chembl_assay_ids.add(row[header_dict['assay_chembl_id']])
                #print (len(chembl_assay_ids))
            
        #return chembl_mol_ids, chembl_assay_ids
        return data, header_dict