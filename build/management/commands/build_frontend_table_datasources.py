# Base imports
import logging

# Django imports
import django.apps
from django.core.management.base import BaseCommand

# Table model imports
from table_provider.models import GpcrStructureBrowserTable
from table_provider.models import GpcrStructureStatisticsTable

# Table query immports
from structure.tables.structure_browser_table_query import StructureBrowserTableRows
from structure.tables.structure_coverage_statistics_query import GpcrStructureCoverageStatisticsQuery

# Other imports
from common.tools import test_model_updates


class Command(BaseCommand):
    help = ()

    logger = logging.getLogger(__name__)
    tracker = {}
    all_models = django.apps.apps.get_models()[6:]
    test_model_updates(all_models, tracker, initialize=True)

    def handle(self, *args, **options):
        functions = [
            'build_structure_browser_table',
            'build_structure_statistics_table',
        ]

        for f in functions:
            try:
                getattr(self, f)()
            except Exception as msg:
                print(msg)
                self.logger.error(msg)
        test_model_updates(self.all_models, self.tracker, check=True)

    def build_structure_browser_table(self):
        """ Builds table_provider.GpcrStructureBrowserTable. 
        
        A flat per-structure table consumed by the structure browser and (round 2) the ligand-interaction PDB selection page.
        """
        self.logger.info('BUILDING STRUCTURE BROWSER TABLE')

        GpcrStructureBrowserTable.objects.all().delete()

        rows = [GpcrStructureBrowserTable(**row) for row in StructureBrowserTableRows()]
        GpcrStructureBrowserTable.objects.bulk_create(rows, batch_size=1000)

        self.logger.info('COMPLETED BUILDING STRUCTURE BROWSER TABLE ({} rows)'.format(len(rows)))

    def build_structure_statistics_table(self):
        """ Builds table_provider.GpcrStructureStatisticsTable, a flat per-structure table
            'consumed by the structure statistics page."""
        self.logger.info('BUILDING STRUCTURE STATISTICS TABLE')

        GpcrStructureStatisticsTable.objects.all().delete()

        rows = [ GpcrStructureStatisticsTable(**row) for row in GpcrStructureCoverageStatisticsQuery() ]
        GpcrStructureStatisticsTable.objects.bulk_create(rows, batch_size=10000)

        self.logger.info('COMPLETED BUILDING STRUCTURE STATISTICS TABLE ({} rows)'.format(len(rows)))
