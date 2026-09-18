from django.core.management.base import BaseCommand

from table_provider.models import GpcrStructureBrowserTable
from structure.tables.structure_browser_table_query import StructureBrowserTableRows
from common.tools import test_model_updates
import logging
import django.apps


class Command(BaseCommand):
    help = ('Builds table_provider.GpcrStructureBrowserTable, a flat per-structure table '
            'consumed by the structure browser and (round 2) the ligand-interaction '
            'PDB selection page.')

    logger = logging.getLogger(__name__)
    tracker = {}
    all_models = django.apps.apps.get_models()[6:]
    test_model_updates(all_models, tracker, initialize=True)

    def handle(self, *args, **options):
        functions = [
            'build_structure_browser_table',
        ]

        for f in functions:
            try:
                getattr(self, f)()
            except Exception as msg:
                print(msg)
                self.logger.error(msg)
        test_model_updates(self.all_models, self.tracker, check=True)

    def build_structure_browser_table(self):
        self.logger.info('BUILDING STRUCTURE BROWSER TABLE')

        GpcrStructureBrowserTable.objects.all().delete()

        rows = [GpcrStructureBrowserTable(**row) for row in StructureBrowserTableRows()]
        GpcrStructureBrowserTable.objects.bulk_create(rows, batch_size=1000)

        self.logger.info('COMPLETED BUILDING STRUCTURE BROWSER TABLE ({} rows)'.format(len(rows)))
