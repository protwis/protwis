from django.core.management.base import BaseCommand

from construct.models import CrystalInfo
from interaction.models import StructureLigandInteraction
from signprot.models import SignprotStructure
from structure.models import Fragment, PdbData, Rotamer, Structure, StructureComplexModel, StructureModel

# (model, field) pairs for every FK pointing at PdbData.
PDB_DATA_REFERRERS = [
    (Structure, 'pdb_data'),
    (StructureModel, 'pdb_data'),
    (StructureComplexModel, 'pdb_data'),
    (Rotamer, 'pdbdata'),
    (Fragment, 'pdbdata'),
    (CrystalInfo, 'pdb_data'),
    (SignprotStructure, 'pdb_data'),
    (StructureLigandInteraction, 'pdb_file'),
]


class Command(BaseCommand):
    help = 'Delete PdbData rows that are no longer referenced by any of the models that point to PdbData. Dry run by default; pass --apply to actually delete.'

    def add_arguments(self, parser):
        parser.add_argument('--apply', action='store_true', help='Actually delete the orphaned rows. Without this flag, only reports counts.')
        parser.add_argument('--batch-size', type=int, default=5000, help='Number of rows to delete per batch when --apply is passed (default: 5000).')

    def handle(self, *args, **options):
        referenced_ids = set()
        for model, field in PDB_DATA_REFERRERS:
            ids = model.objects.exclude(**{field: None}).values_list(field, flat=True).distinct()
            referenced_ids.update(ids)

        total_count = PdbData.objects.count()
        orphaned_qs = PdbData.objects.exclude(id__in=referenced_ids)
        orphaned_count = orphaned_qs.count()

        self.stdout.write(f'Total PdbData rows: {total_count}')
        self.stdout.write(f'Referenced PdbData rows: {len(referenced_ids)}')
        self.stdout.write(f'Orphaned PdbData rows: {orphaned_count}')

        if not options['apply']:
            self.stdout.write(self.style.WARNING('Dry run (no changes made). Re-run with --apply to delete the orphaned rows.'))
            return

        if orphaned_count == 0:
            self.stdout.write(self.style.SUCCESS('Nothing to delete.'))
            return

        batch_size = options['batch_size']
        deleted_total = 0
        while True:
            batch_ids = list(orphaned_qs.order_by('id').values_list('id', flat=True)[:batch_size])
            if not batch_ids:
                break
            PdbData.objects.filter(id__in=batch_ids).delete()
            deleted_total += len(batch_ids)
            self.stdout.write(f'Deleted {deleted_total}/{orphaned_count}')

        self.stdout.write(self.style.SUCCESS(f'Done. Deleted {deleted_total} orphaned PdbData rows.'))
