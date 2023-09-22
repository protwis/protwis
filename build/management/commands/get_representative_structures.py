from django.core.management.base import BaseCommand, CommandError
from django.conf import settings

from structure.models import Structure
from io import StringIO
from Bio.PDB import PDBIO, PDBParser

class Command(BaseCommand):

    def handle(self, *args, **options):

        repr_structs = list(Structure.objects.exclude(structure_type__slug__startswith='af-').order_by('protein_conformation__protein__parent', 'state',
            'resolution').distinct('protein_conformation__protein__parent'))
        print(len(repr_structs))
        print("PDB code\tProtein name\tPreferred_chain")
        for struct in repr_structs:
            #print(struct.pdb_data.pdb)
            #pdb_structure = PDBParser(PERMISSIVE=True).get_structure('name', StringIO(struct.pdb_data.pdb))[0]
            #io = PDBIO()
            #io.set_structure(pdb_structure)
            #io.save("{!s}.pdb".format(struct.pdb_code.index))
            print("{!s}\t{!s}\t{!s}".format(struct.pdb_code.index, struct.protein_conformation.protein.parent.entry_name, struct.preferred_chain))
