from django.core.management.base import BaseCommand, CommandError
from django.conf import settings

from protein.models import ProteinConformation
from residue.models import Residue

import logging
import os
import yaml

class Command(BaseCommand):
    help = 'Exports reference positions for all proteins'

    logger = logging.getLogger(__name__)

    export_dir_path = os.sep.join([settings.DATA_DIR, 'residue_data', 'reference_positions'])

    def handle(self, *args, **options):
        functions = [
            'export_reference_positions',
        ]
        
        # execute functions
        for f in functions:
            try:
                getattr(self, f)()
            except Exception as msg:
                print(msg)
                self.logger.error(msg)

    def export_reference_positions(self):
        self.logger.info('EXPORTING REFERENCE POSITIONS')

        # fetch all wild-type proteins
        pcs = ProteinConformation.objects.filter(protein__sequence_type__slug='wt', protein__species_id=1,
            state__slug=settings.DEFAULT_PROTEIN_STATE).select_related('protein')

        for pc in pcs:
            # fetch reference residues
            export_data = {}
            rs = Residue.objects.filter(protein_conformation=pc,
                generic_number__label__in=settings.REFERENCE_POSITIONS.values()).select_related('generic_number')
            for r in rs:
                export_data[r.generic_number.label] = r.sequence_number

            # open export file for writing
            export_file_path = os.sep.join([self.export_dir_path, pc.protein.entry_name + '.yaml'])
            with open(export_file_path, 'w') as export_file:
                yaml.dump(export_data, export_file, default_flow_style=False)

            self.logger.info('Exported reference posistions for protein {}'.format(pc.protein.entry_name))

        self.logger.info('COMPLETED EXPORTING REFERENCE POSITIONS')