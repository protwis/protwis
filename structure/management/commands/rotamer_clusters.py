from io import StringIO
import numpy as np
from django.core.management.base import BaseCommand
import pprint

import Bio.PDB.Polypeptide as polypeptide
from Bio.PDB import PDBParser
import structure.structural_superposition as sp
from structure.models import Rotamer
from residue.models import ResidueGenericNumberEquivalent

# from sklearn.decomposition import TruncatedSVD
from sklearn.cluster import MeanShift

AAs = {'E':9, 'S':6, 'Y':12, 'G':4, 'A':5, 'V':7, 'M':8, 'L':8, 'I':8, 'T':7, 'F':11, 'H':10, 'K':9,
                         'D':8, 'C':6, 'R':11, 'P':7, 'Q':9, 'N':8, 'W':14, '-':0}

class Command(BaseCommand):
    def add_arguments(self, parser):
        pass

    def handle(self, *args, **options):
        out = {}
        for aa in AAs:
            if aa in ['A','-','G']:
                continue
            tms = ['TM1','TM2','TM3','TM4','TM5','TM6','TM7']
            rotamers = Rotamer.objects.filter(residue__amino_acid=aa, residue__display_generic_number__isnull=False, residue__protein_segment__slug__in=tms)
            data = np.array([0,0])

            first = True
            ref_atoms = []
            alt_atoms = []
            print(aa)
            for rot_obj in rotamers:
                if rot_obj.missing_atoms:
                    continue
                chi_angles = np.array([])
                rot = PDBParser(QUIET=True).get_structure('rot', StringIO(rot_obj.pdbdata.pdb))
                rot.atom_to_internal_coordinates()

                for res in rot.get_residues():
                    try:
                        if polypeptide.three_to_one(res.get_resname())!=aa:
                            continue
                    except KeyError:
                        continue
                    for i in range(1,6):
                        chi = res.internal_coord.get_angle('chi{}'.format(i))
                        # print(gn,s,res,chi,i)
                        if chi:
                            chi_angles = np.append(chi_angles, chi)
                    # if len(chi_angles)>2:
                    #     tsvd = TruncatedSVD(n_components=2, random_state=1)
                    #     if len(chi_angles)==3:
                    #         dummy_array = np.array([0,0,0])
                    #     elif len(chi_angles)==4:
                    #         dummy_array = np.array([0,0,0,0])
                    #     else:
                    #         dummy_array = np.array([0,0,0,0,0])
                    #     chi_angles = np.vstack((chi_angles, dummy_array))
                    #     tsvd_out = tsvd.fit_transform(chi_angles)
                    #     chi_angles = tsvd_out[0]
                    if len(chi_angles)==1:
                        chi_angles = np.append(chi_angles, 0)
                    if first:
                        ref_atoms = [a for a in res]
                        if len(ref_atoms)!=AAs[aa]:
                            continue
                        first = False
                        data = chi_angles
                    else:
                        alt_atoms = [a for a in res]
                        sup = sp.RotamerSuperpose(sorted(ref_atoms), sorted(alt_atoms))
                        sup.run()
                if len(ref_atoms)!=len(alt_atoms):
                    continue
                # print(rot_obj.residue, aa, len(ref_atoms), len(alt_atoms), chi_angles)
                
                if len(chi_angles)>0:
                    data = np.vstack((data, chi_angles))
            ms = MeanShift(bin_seeding=True).fit(data)
            labels = ms.labels_
            _ = ms.cluster_centers_
            out[aa] = len(set(labels))
        pprint.pprint(out)
            
                
