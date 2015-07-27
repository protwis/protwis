from django.core.management.base import BaseCommand, CommandError
from django.conf import settings

from protein.models import ProteinSet
from structure.models import Structure

import logging

class Command(BaseCommand):
    help = 'Reads source data and creates common database tables'

    logger = logging.getLogger(__name__)

    def handle(self, *args, **options):
        functions = [
            'create_protein_sets',
        ]

        # execute functions
        for f in functions:
            try:
                getattr(self, f)()
            except Exception as msg:
                print(msg)
                self.logger.error(msg)

    def create_protein_sets(self):
        self.logger.info('CREATING PROTEIN SETS')

        # proteins with a structure
        structures = Structure.objects.order_by('protein_conformation__protein__parent__entry_name').distinct(
            'protein_conformation__protein__parent__entry_name')
        if structures:
            ps = ProteinSet.objects.create(name='Structure determined')
            for structure in structures:
                ps.proteins.add(structure.protein_conformation.protein.parent)

        self.logger.info('COMPLETED CREATING PROTEIN SETS')
