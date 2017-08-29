from django.core.management.base import BaseCommand, CommandError
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

class Command(BaseCommand):
    help = 'Reads source data and creates links to other databases'

    def add_arguments(self, parser):
        parser.add_argument('--filename', action='append', dest='filename',
            help='Filename to import. Can be used multiple times')

    logger = logging.getLogger(__name__)

    # source file directory
    links_data_dir = os.sep.join([settings.DATA_DIR, 'ligand', 'assay'])
    dictionary_file = os.sep.join([settings.DATA_DIR, 'ligand', 'assay', 'dictionary.txt'])
    
    
    
    ##chembl_to_cid = load file
    
    def handle(self, *args, **options):
        if options['filename']:
            filenames = options['filename']
        else:
            filenames = False
        
        # try:
        # def purge_drugs(self):
        #try:
        #    Drugs.objects.all().delete()
        #except Drugs.DoesNotExist:
        #    self.logger.warning('Drugs mod not found: nothing to delete.')
        self.control(filenames)
        #self.create_assays(filenames)
   
    def control (self, filenames):

        #chembl_mol_ids,assay_ids = self.extract_chembl_mol_ids(filenames)
        data,header_dict = self.extract_chembl_mol_ids(filenames) #reading the entry file from the data filder.
        chembl_cid_dict = self.read_dict_file(self.dictionary_file)
        #found_ids_dict = chembl_cid_dict #testing time. else commented.

##### find the chembl_molecular_ids ############
        chembl_mol_ids = set()
        for i in data:
            chembl_mol_ids.add(i[header_dict['molecule_chembl_id']])
            #print (i[header_dict['pchembl_value']])
######## end ##############
            

        ##notfound_ids_set = self.read_notfound(filenames) # testing
        chembl_cid_dict, notfound_ids_set = self.find_cid(chembl_mol_ids, chembl_cid_dict)
        chembl_cid_dict, notfound_ids_set = self.call_pubchem_service(chembl_cid_dict, notfound_ids_set)
        self.write_dictionary_file(chembl_cid_dict)
        ###print (len(chembl_mol_ids))
        self.create_ligand(chembl_cid_dict)
        self.create_assays(data,header_dict)
        self.create_assaysExperiment(data,header_dict)
        
        

    
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


    def create_ligand(self, found_ids_dict):
        #####Create chembl compound link and connect it to the corresponding ligand/cid#####
        defaults = {
        'name': 'ChEMBL_compound_ids',
        'url': 'https://www.ebi.ac.uk/chembl/compound/inspect/$index'
        }
        wr, created = WebResource.objects.get_or_create(slug='chembl_ligand', defaults = defaults)
        wr_pubchem = WebResource.objects.get(slug='pubchem')
        
        
        for chembl_ligand  in found_ids_dict:
            
            cids = str(found_ids_dict[chembl_ligand])
            temp = cids.split(';') #perhaps we should load all of the CIDs
            #print (temp)
            
            cid = str(temp[0])
            #if cid!='3559':
            #    continue
            l = get_or_make_ligand(cid,'PubChem CID') #call the first cid if there are more than one
            #print (cid)
            wl, created = WebLink.objects.get_or_create(index=chembl_ligand, web_resource=wr)
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
                wl, created = WebLink.objects.get_or_create(index=cid, web_resource=wr_pubchem)
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


    def call_pubchem_service(self, chembl_cid_dict, notfound_ids_set):
        notfound = set()
        for chembl_id in notfound_ids_set:
            url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/'+chembl_id+'/json'
            response = requests.get(url)
            lig_data = response.json()
            try:
                cid = lig_data['PC_Compounds'][0]['id']['id']['cid']
                #updating the existing dictionary
                chembl_cid_dict[chembl_id] = cid
                #.discard(x)
            except KeyError:
                #print (chembl_id)
                notfound.add(chembl_id)


        return chembl_cid_dict, notfound 

    
    ##call pubchem service to find cids for chembl
    def find_cid(self, chembl_mol_ids, chembl_cid_dict):
        notfound = set()
        for chembl_mol_id in chembl_mol_ids:
            
            if chembl_mol_id not in chembl_cid_dict.keys():
                temp = []
                
                url = 'https://www.ebi.ac.uk/unichem/rest/src_compound_id_all/'+chembl_mol_id+'/1/22'
                response = requests.get(url)
                lig_data = response.json()
                #print (lig_data)
                try:
                    for i, x  in enumerate(lig_data):
                        temp.append(lig_data[i]['src_compound_id'])
                    if len(temp) > 1:
                        cid = ';'.join(temp)
                        chembl_cid_dict[chembl_mol_id] = cid
                    elif len(temp) == 1 and temp[0] != '\n':
                        cid = temp[0] #lig_data[0]['src_compound_id'] 
                        #updating the existing dictionary
                        chembl_cid_dict[chembl_mol_id] = cid
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
            #print (filename)
            header_dict = {}
            data = []
            chembl_assay_ids = set ()
            chembl_mol_ids = set()
            with open (filename, 'r') as fii:
                
                linereader = csv.reader(fii)
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



     ##for i in notfound_ids_set:
        #    #inchi_key = self.find_inchi_key_from_chembl_id(i)
        #    
        ##
        #
        #
    #def read_notfound(self,filenames):
    #    for filename in filenames:
    #    
    #        ids = set()
    #        with open(filename, 'r') as fii:
    #            for line in fii:
    #                line = line.rstrip()
    #                ids.add(line)
    #                
    #        return ids

#not in use
    #def find_inchi_key_from_chembl_id(self, chembl_mol_id):
    #    url = 'https://www.ebi.ac.uk/chembl/api/data/molecule/'+chembl_mol_id+'.json'
    #    response = requests.get(url)
    #    mol_data = response.json()
    #    try:
    #        inchy_key = (mol_data['molecule_structures']['standard_inchi_key'])
    #    except TypeError:
    #        print (chembl_mol_id)

        #return inchi_key