from django.core.management.base import BaseCommand, CommandError
from build.management.commands.base_build import Command as BaseBuild
from django.utils.text import slugify
from django.conf import settings
from common.models import WebResource, WebLink
from ligand.models import Ligand, LigandProperities, LigandRole, LigandType, ChemblAssay, AssayExperiment
from ligand.models import LigandVendorLink, LigandVendors
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
        parser.add_argument('--test_run', action='store_true', help='Skip this during a test run', default=False)
    logger = logging.getLogger(__name__)

    # source file directory
    export_dir_path = os.sep.join([settings.DATA_DIR, 'ligand_data','raw_ligands'])


    def handle(self, *args, **options):

        if options['test_run']:
            print('Skipping in test run')
            return
        if options['filename']:
            filenames = options['filename']
        else:
            filenames = os.listdir(self.export_dir_path)

        self.create_vendors(filenames)

        with gzip.open(os.sep.join([settings.DATA_DIR, 'ligand_data','raw_ligands','ligands.json.gz']), "rb") as f:
            self.ligand_dump = json.loads(f.read().decode("ascii"))
        print(len(self.ligand_dump),"ligands to load")
        self.prepare_input(options['proc'], self.ligand_dump)


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

        # print(positions,iteration,count,lock)
        ligands = self.ligand_dump

        while count.value<len(ligands):
            with lock:
                l = ligands[count.value]
                count.value +=1
                if count.value % 10000 == 0:
                    print('{} Status {} out of {}'.format(
                    datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d %H:%M:%S'), count.value, len(ligands)))

            if 'logp' not in l:
                # temp skip to only use "full" annotated ligands
                continue

            lp = LigandProperities.objects.filter(inchikey=l['inchikey']).first()
            ligand = None

            if lp:
                # Check if inchikey is there
                ligand = Ligand.objects.filter(name=l['name'], properities=lp).prefetch_related('properities__ligand_type','properities__web_links','properities__vendors').first()

            # The name with corresponding inchikey is there, assume all is good and skip.
            # Will add links to make sure they're there.
            if not ligand:
                if lp:
                    print(l['name'],'is there! (but not by name, only inchi')
                    ligand = Ligand()
                    ligand.properities = lp
                    ligand.name = l['name']
                    ligand.canonical = l['canonical']
                    ligand.ambigious_alias = l['ambigious_alias']
                    ligand.save()
                else:
                    # No ligand seems to match by inchikey -- start creating it.
                    # Make LigandProperities first
                    lt, created = LigandType.objects.get_or_create(slug=l['ligand_type__slug'],defaults = {'name':l['ligand_type__name']})
                    lp = LigandProperities()
                    lp.inchikey = l['inchikey']
                    lp.smiles = l['smiles']
                    lp.mw = l['mw']
                    lp.logp = l['logp']
                    lp.rotatable_bonds = l['rotatable_bonds']
                    lp.hacc = l['hacc']
                    lp.hdon = l['hdon']
                    lp.ligand_type = lt

                    lp.save()

                    ligand = Ligand()
                    ligand.properities = lp
                    ligand.name = l['name']
                    ligand.canonical = l['canonical']
                    ligand.ambigious_alias = l['ambigious_alias']
                    ligand.save()


            # create links - impossible to make duplicates so no need to check if there already
            if ligand.properities.web_links.count()<len(l['web_links']):
                for link in l['web_links']:
                    wr = WebResource.objects.get(slug=link['web_resource'])
                    wl, created = WebLink.objects.get_or_create(index=link['index'], web_resource=wr)
                    ligand.properities.web_links.add(wl)

            # create vendors - impossible to make duplicates so no need to check if there already
            if ligand.properities.vendors.count()<len(l['vendors']):
                for link in l['vendors']:
                    lv = LigandVendors.objects.get(slug = link['vendor_slug'])
                    check = LigandVendorLink.objects.filter(sid=link['sid']).exists()
                    if not check:
                        lvl = LigandVendorLink()
                        lvl.sid = link['sid']
                        lvl.vendor = lv
                        lvl.lp = ligand.properities
                        lvl.vendor_external_id = link['vendor_external_id']
                        lvl.url = link['url']
                        lvl.save()
