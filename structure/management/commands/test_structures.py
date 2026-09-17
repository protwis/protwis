from django.core.management.base import BaseCommand

from structure.models import Structure
from structure.functions import StructureBuildCheck
from signprot.models import SignprotComplex


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--signprot', help='Test G alpha structures', action='store_true', default=False)
        parser.add_argument('-s', '--structures', dest='structures', nargs='+', default=None,
            help='Limit all checks to these PDB codes (space-separated), e.g. -s 7XW5 6TPJ')

    def handle(self, *args, **options):
        sbc = StructureBuildCheck()
        pdb_filter = [p.upper() for p in options['structures']] if options['structures'] else None
        if not options['signprot']:
            if pdb_filter:
                sbc.pdbs = [p for p in sbc.pdbs if p in pdb_filter]
            sbc.check_structures()
            structs = Structure.objects.all().exclude(structure_type__slug__startswith='af-')
            if pdb_filter:
                structs = structs.filter(pdb_code__index__in=pdb_filter)
            sbc.check_duplicate_residues(structs)
            for s in structs:
                sbc.check_segment_ends(s)
                sbc.check_residue_mismatches(s)
            print("Missing segments: ", len(sbc.missing_seg))
            p = []
            for i in sbc.missing_seg:
                print("Error: Missing segment {} {} has no residue objects. Should have {} to {}".format(i[0],i[1],i[2],i[3]))
                if str(i[0]) not in p:
                    p.append(str(i[0]))
            print(' '.join(p))
            print("Missing parent segment residues: ", len(sbc.missing_parent_seg))
            for i in sbc.missing_parent_seg:
                print("Error: {} {} — parent protein has no residues for this segment (residue build likely needs to be (re-)run for the parent protein)".format(i[0], i[1]))
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
            print("Structures with >10 consecutive residue mismatches vs parent: ", len(sbc.residue_mismatch_structures))
            if pdb_filter:
                for s in sbc.residue_mismatch_structures:
                    structure_residues, parent_by_gn, max_streak = sbc.residue_mismatch_data[s]
                    sbc.print_residue_mismatch_alignment(s, structure_residues, parent_by_gn, max_streak)
            else:
                print(' '.join(str(s) for s in sbc.residue_mismatch_structures))

            for s in structs:
                sbc.check_ligand_interactions(s)
            print("=== Ligand interaction checks ===")
            print("Ligand count mismatches: ", len(sbc.ligand_count_error))
            for i in sbc.ligand_count_error:
                print("Error: {} has {} ligands in ligands.tsv but {} StructureLigandInteraction objects".format(i[0],i[1],i[2]))
            print("Missing residue-fragment interactions: ", len(sbc.missing_ligand_interaction))
            for i in sbc.missing_ligand_interaction:
                print("Error: {} has no ResidueFragmentInteraction objects for ligand {}".format(i[0],i[1]))
            print("Peptide/protein ligand count mismatches: ", len(sbc.peptide_count_error))
            for i in sbc.peptide_count_error:
                print("Error: {} has {} peptide/protein ligands in ligands.tsv but {} LigandPeptideStructure objects".format(i[0],i[1],i[2]))
            print("Missing peptide residue pairs: ", len(sbc.missing_peptide_residue_pair))
            for i in sbc.missing_peptide_residue_pair:
                print("Error: {} has no InteractingPeptideResiduePair objects for peptide ligand {}".format(i[0],i[1]))
        else:
            sc_qs = SignprotComplex.objects.all()
            if pdb_filter:
                sc_qs = sc_qs.filter(structure__pdb_code__index__in=pdb_filter)
            for sc in sc_qs:
                sbc.check_signprot_struct_residues(sc)
            scs = sc_qs.exclude(structure__structure_type__slug__startswith='af-')
            for i in sbc.missing_seg:
                print("Error: Missing segment {} {} has {} residue objects.".format(i[0],i[1],i[2]))
            sbc.check_duplicate_residues(scs)
            for i, j in sbc.duplicate_residue_error.items():
                print("Error: {} has duplicate residue for {}".format(i,j))
            sbc.duplicate_residue_error = {}
            sbc.check_duplicate_residues(scs, 'display_generic_number__label')
            for i, j in sbc.duplicate_residue_error.items():
                print("Error: {} has duplicate GN for {} {}".format(i,j, [k.display_generic_number for k in j]))
