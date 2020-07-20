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

        ProteinSet.objects.all().delete()

        class_dict = {'001': 'A', '002': 'B1', '003': 'B2', '004': 'C', '005': 'D1', '006': 'F', '007': 'T', '008': 'Other'}

        # proteins with a structure
        structures = Structure.objects.order_by('protein_conformation__protein__parent__entry_name').distinct(
            'protein_conformation__protein__parent__entry_name')
        if structures:
            ps = ProteinSet.objects.create(name='All') # David's request
            ps_class = {}
            for structure in structures:
                # Grab the class slug
                pc = structure.protein_conformation.protein.parent.family.slug.split("_")[0]
                if pc not in ps_class:
                    ps_class[pc] = ProteinSet.objects.create(name='{}'.format(class_dict[pc])) # David's request

                ps.proteins.add(structure.protein_conformation.protein.parent)
                ps_class[pc].proteins.add(structure.protein_conformation.protein.parent)

        self.logger.info('COMPLETED CREATING PROTEIN SETS')
