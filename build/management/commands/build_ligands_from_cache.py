from django.core.management.base import BaseCommand, CommandError
from build.management.commands.base_build import Command as BaseBuild
from django.utils.text import slugify
from django.conf import settings
from common.models import WebResource, WebLink
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
import gzip, json

class Command(BaseBuild):
    help = 'Reads ligand_cache and creates the ligands'

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
    export_dir_path = os.sep.join([settings.DATA_DIR, 'ligand_data','raw_ligands'])


    def handle(self, *args, **options):
        if options['filename']:
            filenames = options['filename']
        else:
            filenames = os.listdir(self.export_dir_path)

        self.create_vendors(filenames)


        with gzip.open(os.sep.join([settings.DATA_DIR, 'ligand_data','raw_ligands','ligands.json.gz']), "rb") as f:
            ligand_dump = json.loads(f.read().decode("ascii"))
        print(len(ligand_dump),"ligands to load")
        # self.prepare_input(options['proc'], self.chembl_mol_ids)


    def create_vendors(self,filenames):
        for filename in filenames:
            if filename.startswith("vendors"):
                print(filename)
                with gzip.open(os.sep.join([settings.DATA_DIR, 'ligand_data','raw_ligands',filename]), "rb") as f:
                    d = json.loads(f.read().decode("ascii"))
                create_count = 0
                for v in d:
                    lv, created = LigandVendors.objects.get_or_create(slug = v['slug'])
                    if created:
                        create_count += 1
                        lv.name = v['name']
                        lv.url = v['url']
                        lv.save()
                print(len(d),"vendors",create_count,"vendors created")  

    def main_func(self, positions, iteration,count,lock):
        #####Create chembl compound link and connect it to the corresponding ligand/cid#####

        print(positions,iteration,count,lock)
        list_of_chembl_ids = self.chembl_mol_ids
        while count.value<len(list_of_chembl_ids):
            with lock:
                chembl_ligand = list_of_chembl_ids[count.value]
                count.value +=1 
                if count.value % 1000 == 0:
                    print('{} Status {} out of {}'.format(
                    datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d %H:%M:%S'), count.value, len(list_of_chembl_ids)))

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
            #if cid!='3559':
            #    continue
            l = get_or_make_ligand(cid,'PubChem CID') #call the first cid if there are more than one
            if not l:
                print('Ligand not found in PubChem', cid)
                continue

            if not l.properities.web_links.filter(web_resource__slug = 'pubchem',index = cid).exists():
                # NO CID FOR LIGAND! Rare cases where SMILES was used for initial look up
                wl, created = WebLink.objects.get_or_create(index=cid, web_resource=self.wr_pubchem)
                l.properities.web_links.add(wl)

            if not l.properities.web_links.filter(web_resource__slug = 'chembl_ligand',index = chembl_ligand).exists():
                wl, created = WebLink.objects.get_or_create(index=chembl_ligand, web_resource=self.wr)
                l.properities.web_links.add(wl)
            
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
 