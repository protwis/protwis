from io import StringIO
from collections import OrderedDict
from django.core.management.base import BaseCommand
from copy import deepcopy
import pprint

import Bio.PDB.Polypeptide as polypeptide
from Bio.PDB import PDBParser
import structure.structural_superposition as sp
from structure.models import Structure, Rotamer
from structure.functions import run_residue_flip
from residue.functions import dgn


class Command(BaseCommand):
    def add_arguments(self, parser):
        pass

    def handle(self, *args, **options):
        ### To be filled
        s = PDBParser(PERMISSIVE=True).get_structure('ref', './Alphafold/7ul3.pdb')[0]['A']
        af = PDBParser(PERMISSIVE=True).get_structure('af', '../../data/protwis/gpcr/structure_data/Alphafold/hrh2_human_inactive.pdb')[0]['A']
        resnums =  []
        gn_dict = OrderedDict()
        ###

        templates = OrderedDict()
        coverage = {}
        for r, g in zip(resnums, gn_dict):
            struct_atoms = [a for a in s[r] if not a.get_id().startswith('H')]
            af_atoms = deepcopy([a for a in af[r] if not a.get_id().startswith('H')])
            superpose = sp.RotamerSuperpose(sorted(struct_atoms), sorted(af_atoms))
            superpose.run()
            rmsd = deepcopy(superpose.rmsd)
            if af_atoms[0].get_parent().get_resname() in ['TYR','PHE','ARG','ASP','GLU']:
                new_af_atoms = run_residue_flip(af_atoms)
                superpose2 = sp.RotamerSuperpose(sorted(struct_atoms), sorted(new_af_atoms))
                superpose2.run()
                if rmsd>superpose2.rmsd:
                    rmsd = round(superpose2.rmsd,3)

            AA = polypeptide.three_to_one(af_atoms[0].get_parent().get_resname())
            this_temp = ['AF', round(rmsd,3)]
            for pdb in Structure.objects.all().exclude(structure_type__slug__startswith='af-'):#gn_dict[g]:
                struct = pdb#Structure.objects.get(pdb_code__index=pdb)
                try:
                    dgn_lab = dgn(g, struct.protein_conformation)
                    rot_obj = Rotamer.objects.filter(structure=struct, residue__amino_acid=AA, residue__display_generic_number__label=dgn_lab)
                except ResidueGenericNumberEquivalent.DoesNotExist:
                    continue
                if len(rot_obj)==0:
                    continue
                rot_obj = self.right_rotamer_select(rot_obj)
                rot = PDBParser().get_structure('rot', StringIO(rot_obj.pdbdata.pdb))[0]
                for c in rot:
                    for res in c:
                        alt_atoms = deepcopy([a for a in res])
                if len(af_atoms)!=len(alt_atoms):
                    continue
                superpose = sp.RotamerSuperpose(sorted(struct_atoms), sorted(alt_atoms))
                superpose.run()
                alt_rmsd = deepcopy(superpose.rmsd)
                if alt_atoms[0].get_parent().get_resname() in ['TYR','PHE','ARG','ASP','GLU']:
                    new_alt_atoms = run_residue_flip(alt_atoms)
                    superpose2 = sp.RotamerSuperpose(sorted(struct_atoms), sorted(new_alt_atoms))
                    superpose2.run()
                    if alt_rmsd>superpose2.rmsd:
                        alt_rmsd = superpose2.rmsd
                if alt_rmsd<rmsd:
                    print(r,g,struct, rmsd, alt_rmsd)
                    if struct not in coverage:
                        coverage[struct] = 1
                    else:
                        coverage[struct]+=1
                if alt_rmsd<this_temp[1]:
                    this_temp = [struct, round(alt_rmsd,3)]
            templates[r] = [res.get_resname(), 'AF', round(rmsd,3), this_temp]
        pprint.pprint(templates)
        coverage = sorted(coverage.items(), key=lambda x: (-x[1]))
        pprint.pprint(coverage)


    def right_rotamer_select(self, rotamer):
        ''' Filter out compound rotamers.
        '''
        if len(rotamer)>1:
            for i in rotamer:
                if not i.pdbdata.pdb.startswith('COMPND'):
                    rotamer = i
                    break
        else:
            rotamer=rotamer[0]
        return rotamer
