from django.core.management.base import BaseCommand

from structure.models import Structure
from structure.functions import StructureBuildCheck
from signprot.models import SignprotComplex


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--signprot', help='Test G alpha structures', action='store_true', default=False)

    def handle(self, *args, **options):
        sbc = StructureBuildCheck()
        if not options['signprot']:
            sbc.check_structures()
            structs = Structure.objects.all().exclude(structure_type__slug__startswith='af-')
            sbc.check_duplicate_residues()
            for s in structs:
                sbc.check_segment_ends(s)
            print("Missing segments: ", len(sbc.missing_seg))
            for i in sbc.missing_seg:
                print("Error: Missing segment {} {} has no residue objects. Should have {} to {}".format(i[0],i[1],i[2],i[3]))
            print("Very short helix segments: ", len(sbc.helix_length_error))
            for i in sbc.helix_length_error:
                print("Error: Very short helix segment for {}: {} {}".format(i[0],i[1],i[2]))
            print("Start errors: ", len(sbc.start_error))
            for i in sbc.start_error:
                print("Error: {} {} starts at {} instead of annotated {}".format(i[0],i[1],i[2],i[3]))
            print("End errors: ", len(sbc.end_error))
            for i in sbc.end_error:
                print("Error: {} {} ends at {} instead of annotated {}".format(i[0],i[1],i[2],i[3]))
            print("Residue duplicate errors: ", len(sbc.duplicate_residue_error))
            for i, j in sbc.duplicate_residue_error.items():
                print("Error: {} has duplicate residue for {}".format(i,j))
        else:
            for sc in SignprotComplex.objects.all():
                sbc.check_signprot_struct_residues(sc)
