from django.core.management.base import BaseCommand, CommandError
from django.conf import settings

from ligand.models import Ligand, LigandProperities, LigandRole

import logging
import os
import yaml

class Command(BaseCommand):
    help = 'Exports reference positions for all proteins'

    logger = logging.getLogger(__name__)

    export_dir_path = os.sep.join([settings.DATA_DIR, 'ligand_data'])
    if not os.path.exists(export_dir_path):
        os.makedirs(export_dir_path)
    def handle(self, *args, **options):
        functions = [
            'export_ligands',
        ]
        
        # execute functions
        for f in functions:
            try:
                getattr(self, f)()
            except Exception as msg:
                print(msg)
                self.logger.error(msg)

    def export_ligands(self):
        self.logger.info('EXPORTING LIGANDS')

        # fetch all ligands
        ls = Ligand.objects.all().prefetch_related('properities','properities__web_links')
        export_data = []

        for l in ls:
            export = {}
            export['name'] = l.name
            export['canonical'] = l.canonical
            export['ambigious_alias'] = l.ambigious_alias
            export['smiles'] = l.properities.smiles
            export['inchikey'] = l.properities.inchikey
            if l.properities.ligand_type:
                export['ligand_type__slug'] = l.properities.ligand_type.slug
                export['ligand_type__name'] = l.properities.ligand_type.name
            export['ligand__weblinks'] = []
            for link in l.properities.web_links.all():
                export['ligand__weblinks'].append([link.index,link.web_resource.slug])
            #export['ligand__weblinks'] = l.properities.web_links
            export_data.append(export)

            # open export file for writing
        export_file_path = os.sep.join([self.export_dir_path,'ligands.yaml'])
        with open(export_file_path, 'w') as export_file:
            yaml.dump(export_data, export_file, default_flow_style=False)

        self.logger.info('COMPLETED EXPORTING LIGANDS - Entries:'+str(len(export_data)))