from build.management.commands.base_build import Command as BaseBuild

from structure.models import Structure
from structure.functions import StructureBuildCheck


class Command(BaseBuild):

    def add_arguments(self, parser):
        super(Command, self).add_arguments(parser=parser)

    def handle(self, *args, **options):
        sbc = StructureBuildCheck()
        sbc.check_structures()
        structs = Structure.objects.all()
        for s in structs:
            sbc.check_structure_residues(s)
