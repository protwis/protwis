from django.core.management.base import BaseCommand, CommandError
from build.management.commands.base_build import Command as BaseBuild
from django.utils.text import slugify
from django.conf import settings
from django.db import IntegrityError
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
import datetime

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
        parser.add_argument('--test_run', action='store_true', help='Skip this during a test run', default=False)

    logger = logging.getLogger(__name__)

    # source file directory
    links_data_dir = os.sep.join([settings.DATA_DIR, 'ligand_data', 'assay_data'])
    dictionary_file = os.sep.join([settings.DATA_DIR, 'ligand_data', 'assay_data', 'dictionary.txt'])

    wr = WebResource.objects.get(slug='chembl_ligand')
    wr_pubchem = WebResource.objects.get(slug='pubchem')
        
    
    
    def handle(self, *args, **options):
        if options['test_run']:
            print('Skipping in test run')
            return
        if options['filename']:
            filenames = options['filename']
        else:
            filenames = os.listdir(self.links_data_dir)

        self.chembl_cid_dict = self.read_dict_file(self.dictionary_file)
        self.chembl_cid_dict = OrderedDict(self.chembl_cid_dict)

        print("size of current dict",len(self.chembl_cid_dict))
        self.data,self.header_dict = self.load_data(filenames) #reading the entry file from the data filder.

        ##### find the chembl_molecular_ids ############
        self.chembl_mol_ids = set()
        for i in self.data:
            self.chembl_mol_ids.add(i[self.header_dict['molecule_chembl_id']])
        self.chembl_mol_ids = sorted(list(self.chembl_mol_ids))
        ######## end ##############
        print("Distinct ligands in Dataset",len(self.chembl_mol_ids))
        print("Total rows in dataset",len(self.data))

        # Code that counts duplicates
        check = []
        lists = {}
        for i,record in enumerate(self.data):
            target = record[self.header_dict['target_chembl_id']]
            assay = record[self.header_dict['assay_chembl_id']]
            ligand =record[self.header_dict['molecule_chembl_id']]
            common = '%s_%s_%s' % (target,assay,ligand)
            check.append(common)
            if not common in lists:
                lists[common] = []
            lists[common].append(record)
        
        from collections import Counter
        a = dict(Counter(check))
        count = 0
        for k,v in a.items():
            if v>1:
                count += 1
                # print(k)
                # for l in lists[k]:
                #     print(l)
        print('Data rows that seem duplicated',count)
            
        # Load all ligands ( possible to skip if believed to be included or already imported )
        #self.prepare_input(options['proc'], self.chembl_mol_ids, 0)
        # Insert the actual data points
        self.prepare_input(options['proc'], self.data ,1)

        
    def main_func(self, positions, iteration,count,lock):
        #####Create chembl compound link and connect it to the corresponding ligand/cid#####

        if iteration==0:
            # First load makes sure ligands are there
            list_of_chembl_ids = self.chembl_mol_ids
            while count.value<len(list_of_chembl_ids):
                with lock:
                    chembl_ligand = list_of_chembl_ids[count.value]
                    count.value +=1 
                    if count.value % 1000 == 0:
                        print('{} Status {} out of {}'.format(
                        datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d %H:%M:%S'), count.value, len(list_of_chembl_ids)))

                l = Ligand.objects.filter(properities__web_links__web_resource__slug = 'chembl_ligand', properities__web_links__index=chembl_ligand).first()
                if l:
                    cid = l.properities.web_links.filter(web_resource__slug = 'pubchem').first()
                    if cid:
                        cid = cid.index
                    else:
                        l = None
                        # make sure code blow is run

                if not l:
                    # if l already has chembl link, assume all is good.
                    if chembl_ligand not in self.chembl_cid_dict.keys():
                        cids, not_found = self.find_cid_for_chembl(chembl_ligand)
                        if not_found:
                            print('SKIPPED: Could not determine CID',chembl_ligand,cids)
                            continue
                    else:
                        cids = self.chembl_cid_dict[chembl_ligand]

                    temp = str(cids).split(';') #perhaps we should load all of the CIDs
                    cid = str(temp[0])

                    l = get_or_make_ligand(cid,'PubChem CID') #call the first cid if there are more than one
                    if not l:
                        print('SKIPPED: Ligand not found in PubChem', cid)
                        continue

                    if not l.properities.web_links.filter(web_resource__slug = 'pubchem',index = cid).exists():
                        # NO CID FOR LIGAND! Rare cases where SMILES was used for initial look up
                        wl, created = WebLink.objects.get_or_create(index=cid, web_resource=self.wr_pubchem)
                        l.properities.web_links.add(wl)

                    if not l.properities.web_links.filter(web_resource__slug = 'chembl_ligand',index = chembl_ligand).exists():
                        wl, created = WebLink.objects.get_or_create(index=chembl_ligand, web_resource=self.wr)
                        l.properities.web_links.add(wl)
                
                ###### Vendor stuff  ######
                if not len(l.properities.vendors.all()):
                    # If it has some, assume they are all loaded
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

        elif iteration==1:
            # Third load loads the exp (based on ligand/assay)
            header = self.header_dict
            skipped = 0
            non_p = []
            wr_chembl_assays = WebResource.objects.get(slug='chembl_assays')
            while count.value<len(self.data):
                with lock:
                    record = self.data[count.value]
                    count.value +=1 
                    if count.value % 10000 == 0:
                        print('{} Status {} out of {}'.format(
                        datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d %H:%M:%S'), count.value, len(self.data)))


                target = record[header['target_chembl_id']]
                assay_id = record[header['assay_chembl_id']]

                assay, created = ChemblAssay.objects.get_or_create(assay_id=assay_id)
                if created:
                    wl, created = WebLink.objects.get_or_create(index=assay_id, web_resource=wr_chembl_assays)
                    assay.web_links.add(wl)


                ligand =record[header['molecule_chembl_id']]
                p = Protein.objects.filter(web_links__index = target, web_links__web_resource__slug = 'chembl').first()
                if not p:
                    if not target in non_p:
                        non_p.append(target)
                        print('Not found protein!',target)
                    continue

                ls = Ligand.objects.filter(properities__web_links__index=ligand, properities__web_links__web_resource__slug = 'chembl_ligand', canonical=True)
                if not ls.exists():
                    # if no ligand matches this, then ignore -- be sure this works later.
                    skipped += 1
                    continue
                for l in ls:
                    if len(ls)>1:
                        print('issue with canonical! give to munk',l,l.pk,ligand)
                        break
                assay_experiments = AssayExperiment.objects.filter( protein=p, ligand=l, assay=assay)
                
                if assay_experiments.exists():
                    assay_experiment = assay_experiments.get()
                else:
                    assay_experiment = AssayExperiment()
                    assay_experiment.assay = assay
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
                
                try:
                    assay_experiment.save()
                except IntegrityError:
                    assay_experiment = AssayExperiment.objects.get( protein=p, ligand=l, assay=assay)

            print('done, skipped:',skipped)

    def find_cid_for_chembl(self, chembl_mol_id):
        # function to find cid based on chembl
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
                # print("Searching for ",chembl_mol_id,len(chembl_mol_ids))
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

    ##read pre-generated file and extract the chembl_ids
    def load_data(self, filenames=False):
        for filename in filenames:
            if filename=='dictionary.txt':
                continue
            full_file = os.sep.join([self.links_data_dir,filename])
            print (filename,full_file)
            header_dict = OrderedDict()
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

        return data, header_dict