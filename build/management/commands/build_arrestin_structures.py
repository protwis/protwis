from build.management.commands.base_build import Command as BaseBuild

from protein.models import Protein, ProteinConformation, ProteinSequenceType, ProteinSource, ProteinState
from residue.models import Residue, ResidueNumberingScheme
from signprot.models import SignprotComplex
from structure.functions import create_structure_rotamer, get_pdb_ids, fetch_signprot_data, build_signprot_struct
from structure.models import Rotamer
from common.tools import test_model_updates
from Bio.PDB import PDBParser, PPBuilder, PDBIO, Polypeptide
from io import StringIO
from Bio import pairwise2
import logging
import urllib
import json
import datetime
import django.apps


class Command(BaseBuild):

    #Setting the variables for the test tracking of the model upadates
    tracker = {}
    all_models = django.apps.apps.get_models()[6:]
    test_model_updates(all_models, tracker, initialize=True)

    def add_arguments(self, parser):
        parser.add_argument('-p', '--proc', type=int, action='store', dest='proc', default=1, help='Number of processes to run')
        parser.add_argument("-s", default=False, type=str, action="store", nargs="+", help="PDB codes to build")
        parser.add_argument("--debug", default=False, action="store_true", help="Debug mode")
        parser.add_argument('-u', '--purge', default=False, action="store_true", help="Purge arrestin complex structures from database")


    def handle(self, *args, **options):
        startTime = datetime.datetime.now()
        if options['purge']:
            self.purge()
            self.tracker = {}
            test_model_updates(self.all_models, self.tracker, initialize=True)

        # Complex structures
        self.scs = SignprotComplex.objects.filter(protein__family__slug__startswith="200")
        for sc in self.scs:
            pdbp = PDBParser(PERMISSIVE=True, QUIET=True)
            s = pdbp.get_structure("struct", StringIO(sc.structure.pdb_data.pdb))
            chain = s[0][sc.alpha]
            structure_seq = ''
            nums = []
            for res in chain:
                if "CA" in res and res.id[0]==" ":
                    if sc.structure.pdb_code.index=='5DGY':
                        if not 1012<=res.get_id()[1]<=1361:
                            continue
                    elif sc.structure.pdb_code.index=='5W0P':
                        if not 2012<=res.get_id()[1]<=2362:
                            continue
                    elif sc.structure.pdb_code.index=='4ZWJ':
                        if not 2012<=res.get_id()[1]<=2361:
                            continue
                    nums.append(res.get_id()[1])
                    structure_seq+=Polypeptide.three_to_one(res.get_resname())

            if options['debug']:
                print(sc)
                print('Structure seq:')
                print(structure_seq)

            protein, created = Protein.objects.get_or_create(entry_name=sc.structure.pdb_code.index.lower()+'_arrestin', accession=None,
                                                             name=sc.structure.pdb_code.index.lower()+'_arrestin',
                                                             sequence=structure_seq, family=sc.protein.family, parent=sc.protein,
                                                             residue_numbering_scheme=ResidueNumberingScheme.objects.get(slug='can'),
                                                             sequence_type=ProteinSequenceType.objects.get(slug="mod"), source=ProteinSource.objects.get(name="OTHER"),
                                                             species=sc.protein.species)
            protconf, created = ProteinConformation.objects.get_or_create(protein=protein, state=ProteinState.objects.get(slug='active'))

            pw2 = pairwise2.align.localms(protein.parent.sequence, structure_seq, 3, -4, -3, -1)
            ref_seq, temp_seq = str(pw2[0][0]), str(pw2[0][1])

            parent_residues = Residue.objects.filter(protein_conformation__protein=protein.parent)

            bulked_rotamers = []
            r_c, s_c = 0, 0
            for r, t in zip(ref_seq, temp_seq):
                if options['debug']:
                    print(r,t)
                if t!='-':
                    res = Residue()
                    res.sequence_number = nums[s_c]
                    res.amino_acid = t
                    res.display_generic_number = parent_residues[r_c].display_generic_number
                    res.generic_number = parent_residues[r_c].generic_number
                    res.protein_conformation = protconf
                    res.protein_segment = parent_residues[r_c].protein_segment
                    res.save()
                    rot = create_structure_rotamer(chain[nums[s_c]], res, sc.structure)
                    bulked_rotamers.append(rot)
                    s_c+=1
                r_c+=1
            Rotamer.objects.bulk_create(bulked_rotamers)

        # Non-complex arrestin structures
        arrestin_proteins = Protein.objects.filter(family__slug__startswith="200", accession__isnull=False)
        complex_structures = self.scs.values_list("structure__pdb_code__index", flat=True)
        for a in arrestin_proteins:
            pdb_list = get_pdb_ids(a.accession)
            for pdb in pdb_list:
                if pdb not in complex_structures:
                    try:
                        data = fetch_signprot_data(pdb, a)
                        if data:
                            build_signprot_struct(a, pdb, data)
                    except Exception as msg:
                        self.logger.error("SignprotStructure of {} {} failed\n{}: {}".format(a.entry_name, pdb, type(msg), msg))
        if options["debug"]:
            print(datetime.datetime.now() - startTime)

        test_model_updates(self.all_models, self.tracker, check=True)


    def purge(self):
        Residue.objects.filter(protein_conformation__protein__entry_name__endswith='_arrestin')
        ProteinConformation.objects.filter(protein__entry_name__endswith='_arrestin').delete()
        Protein.objects.filter(entry_name__endswith='_arrestin').delete()
