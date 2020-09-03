from django.core.management.base import BaseCommand, CommandError
from django.conf import settings

from protein.models import Protein, ProteinGProteinPair
from ligand.models import *
from common.models import WebLink, Publication
import pandas as pd

import logging
import os
from decimal import Decimal
import gzip, json
import datetime

class Command(BaseCommand):
    help = 'Exports ligands, to lessen load on pubchem when loading'

    logger = logging.getLogger(__name__)

    export_dir_path = os.sep.join([settings.DATA_DIR, 'ligand_assay'])
    if not os.path.exists(export_dir_path):
        os.makedirs(export_dir_path)
    def handle(self, *args, **options):
        functions = [
            'export_ligand_assay',
        ]

        # execute functions
        for f in functions:
            getattr(self, f)()

    def export_ligand(self, l):
        export = {}
        export['name'] = l.name
        export['canonical'] = l.canonical
        export['ambigious_alias'] = l.ambigious_alias

        # Properities
        export['smiles'] = l.properities.smiles
        export['inchikey'] = l.properities.inchikey
        if not l.properities.inchikey:
            # If no inchikey, do not dump
            export = None
        if l.properities.mw:
            export['mw'] = str(Decimal(l.properities.mw))
        else:
            # If no mw, do not dump, as it needs to be created and attempted to find it
            export = None
        if l.properities.logp:
            export['logp'] = str(Decimal(l.properities.logp))
        else:
            # If no mw, do not dump, as it needs to be created and attempted to find it
            export = None
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
            export = None
        if export:
            return export



    def export_ligand_assay(self):
        self.logger.info('EXPORTING LIGANDS')
        # fetch all ligands
        ls = AssayExperiment.objects.all().prefetch_related('ligand', 'protein', 'publication')
        export_data = []
        print(len(ls))
        number_of_assays = len(ls)
        increment  = 0
        for l in enumerate(ls):
            export = {}
            increment=increment+1
            if increment%10000 == 0:
                print('---status---', increment)
            # Main ligand field
            export['ligand'] = str(l[1].ligand)
            # export['ligand'] = Ligand.objects.filter(name=l[1].ligand).values()[0]
            #ligands data
            export['protein'] = str(l[1].protein)
            # export['protein'] = Protein.objects.filter(entry_name=l[1].protein).values()[0]

            export['assay'] = str(l[1].assay)
            # chembl_assay_len = ChemblAssay.objects.filter(assay_id=l[1].assay).values()
            # if len(chembl_assay_len)>0:
            #     export['assay'] = chembl_assay_len[0]
            export['publication'] = None

            if l[1].publication != None:
                export['publication']= str(l[1].publication.reference)

            # Properities
            export['assay_type'] = l[1].assay_type
            export['assay_description'] = l[1].assay_description
            export['pchembl_value'] = l[1].pchembl_value
            export['published_value'] = None
            if export['published_value']:
                export['published_value'] = float(l[1].published_value)

            export['published_relation'] = l[1].published_relation
            export['published_type'] = l[1].published_type
            export['published_units'] = l[1].published_units

            export['standard_value'] = l[1].standard_value
            export['standard_relation'] = l[1].standard_relation
            export['standard_type'] = l[1].standard_type
            export['standard_units'] = l[1].standard_units

            export['chembl'] = l[1].chembl
            export['smiles'] = l[1].smiles
            export['activity'] = l[1].activity
            export['document_chembl_id'] = l[1].document_chembl_id
            export['cell_line'] = l[1].cell_line
            # print('\n',export)
            export_data.append(export)

        print('total ligands to dump',len(export_data))

        export_file_path_gz = os.sep.join([self.export_dir_path,'ligands.json.gz'])
        with gzip.open(export_file_path_gz, mode="wt") as f:
            json.dump(export_data, f)
        print('path', self.export_dir_path)
        self.logger.info('COMPLETED EXPORTING LIGANDS - Entries:'+str(len(export_data)))
