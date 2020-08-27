from django.core.management.base import BaseCommand, CommandError
from django.conf import settings

from ligand.models import *
import logging
import os
from decimal import Decimal
import gzip, json
import datetime

class Command(BaseCommand):
    help = 'Exports ligands, to lessen load on pubchem when loading'

    logger = logging.getLogger(__name__)

    export_dir_path = os.sep.join([settings.DATA_DIR, 'ligand_data','raw_ligands'])
    if not os.path.exists(export_dir_path):
        os.makedirs(export_dir_path)
    def handle(self, *args, **options):
        functions = [
            'export_vendors',
            'export_ligands',
        ]

        # execute functions
        for f in functions:
            getattr(self, f)()

    def export_vendors(self):
        self.logger.info('EXPORTING LIGAND VENDORS')
        # fetch all vendors
        lvs = LigandVendors.objects.all()
        export_data = []
        for lv in lvs:
            export = {}
            export['slug'] = lv.slug
            export['name'] = lv.name
            export['url'] = lv.url
            export_data.append(export)
        export_file_path_gz = os.sep.join([self.export_dir_path,'vendors.json.gz'])
        with gzip.open(export_file_path_gz, mode="wt") as f:
          json.dump(export_data, f)

    def export_ligands(self):
        self.logger.info('EXPORTING LIGANDS')

        # fetch all ligands
        ls = Ligand.objects.all().prefetch_related('properities__ligand_type','properities__web_links','properities__vendors')
        export_data = []
        number_of_ligands = len(ls)
        for i, l in enumerate(ls):
            if i % 1000 == 0:
                print('{} Status {} out of {}'.format(
                datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d %H:%M:%S'), i, number_of_ligands))
            export = {}
            # Main ligand fields
            export['name'] = l.name
            export['canonical'] = l.canonical
            export['ambigious_alias'] = l.ambigious_alias

            # Properities
            export['smiles'] = l.properities.smiles
            export['inchikey'] = l.properities.inchikey
            if not l.properities.inchikey:
                # If no inchikey, do not dump
                continue
            if l.properities.mw:
                export['mw'] = str(Decimal(l.properities.mw))
            else:
                # If no mw, do not dump, as it needs to be created and attempted to find it
                continue
            if l.properities.logp:
                export['logp'] = str(Decimal(l.properities.logp))
            else:
                # If no mw, do not dump, as it needs to be created and attempted to find it
                continue
            export['rotatable_bonds'] = l.properities.rotatable_bonds
            export['hacc'] = l.properities.hacc
            export['hdon'] = l.properities.hdon
            if l.properities.ligand_type:
                export['ligand_type__slug'] = l.properities.ligand_type.slug
                export['ligand_type__name'] = l.properities.ligand_type.name

            # Links
            export['web_links'] = []
            for wl in l.properities.web_links.all():
                wl_export = {}
                wl_export['web_resource'] = wl.web_resource.slug
                wl_export['index'] = wl.index
                export['web_links'].append(wl_export)

            if not len(export['web_links']):
                # If no weblinks, do not dump, as it needs to be created and attempted to find them
                continue

            # Vendors
            export['vendors'] = []
            for wl in l.properities.vendors.all():
                v_export = {}
                v_export['url'] = wl.url
                v_export['vendor_external_id'] = wl.vendor_external_id
                v_export['sid'] = wl.sid
                v_export['vendor_slug'] = wl.vendor.slug
                export['vendors'].append(v_export)
            export_data.append(export)

        print('total ligands to dump',len(export_data))
        export_file_path_gz = os.sep.join([self.export_dir_path,'ligands.json.gz'])
        with gzip.open(export_file_path_gz, mode="wt") as f:
            json.dump(export_data, f)

        self.logger.info('COMPLETED EXPORTING LIGANDS - Entries:'+str(len(export_data)))
