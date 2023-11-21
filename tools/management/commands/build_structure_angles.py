from django.core.management.base import BaseCommand
from django.db import connection

import contactnetwork.pdb as pdb
from structure.models import Structure, StructureVectors
from residue.models import Residue
from angles.models import ResidueAngle as Angle
from contactnetwork.models import Distance, distance_scaling_factor
from common.tools import test_model_updates
import Bio.PDB
import copy
import freesasa
import io
import logging
import math
import subprocess
import os
import re
import traceback
import django.apps


import numpy as np
import scipy.stats as stats

from scipy.spatial.transform import Rotation as R
from collections import OrderedDict
from sklearn.decomposition import PCA
from numpy.core.umath_tests import inner1d


from multiprocessing import Queue, Process, Value, Lock

SASA = True
HSE  = True
extra_pca = True
print_pdb = False
GN_only = False
incremental_update = False

# atom name dictionary
# Based on https://github.com/fomightez/structurework/blob/master/spartan_fixer/SPARTAN08_Fixer_standalone.py
residue_atom_names = {
    'ALA': ['N', 'CA', 'C', 'O', 'CB'],
    'ARG': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2'],
    'ASN': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'OD1', 'ND2'],
    'ASP': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'OD1', 'OD2'],
    'CYS': ['N', 'CA', 'C', 'O', 'CB', 'SG'],
    'GLU': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'OE1', 'OE2'],
    'GLN': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'OE1', 'NE2'],
    'GLY': ['N', 'CA', 'C', 'O'],
    'HIS': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2'],
    'ILE': ['N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2', 'CD1'],
    'LEU': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2'],
    'LYS': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'CE', 'NZ'],
    'MET': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'SD', 'CE'],
    'PHE': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
    'PRO': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD'],
    'SER': ['N', 'CA', 'C', 'O', 'CB', 'OG'],
    'THR': ['N', 'CA', 'C', 'O', 'CB', 'OG1', 'CG2'],
    'TRP': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
    'TYR': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH'],
    'VAL': ['N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2']
}
residue_atom_names['YCM'] = residue_atom_names['CYS']
residue_atom_names['CSD'] = residue_atom_names['ALA']
residue_atom_names['TYS'] = residue_atom_names['TYR']
residue_atom_names['SEP'] = residue_atom_names['SER']

# sidechain torsion angles - chi dihedrals
# From https://gist.github.com/lennax/0f5f65ddbfa278713f58
chi_atoms = dict(
    chi1=dict(
        ARG=['N', 'CA', 'CB', 'CG'],
        ASN=['N', 'CA', 'CB', 'CG'],
        ASP=['N', 'CA', 'CB', 'CG'],
        CYS=['N', 'CA', 'CB', 'SG'],
        GLN=['N', 'CA', 'CB', 'CG'],
        GLU=['N', 'CA', 'CB', 'CG'],
        HIS=['N', 'CA', 'CB', 'CG'],
        ILE=['N', 'CA', 'CB', 'CG1'],
        LEU=['N', 'CA', 'CB', 'CG'],
        LYS=['N', 'CA', 'CB', 'CG'],
        MET=['N', 'CA', 'CB', 'CG'],
        PHE=['N', 'CA', 'CB', 'CG'],
        PRO=['N', 'CA', 'CB', 'CG'],
        SER=['N', 'CA', 'CB', 'OG'],
        THR=['N', 'CA', 'CB', 'OG1'],
        TRP=['N', 'CA', 'CB', 'CG'],
        TYR=['N', 'CA', 'CB', 'CG'],
        VAL=['N', 'CA', 'CB', 'CG1'],
        YCM=['N', 'CA', 'CB', 'SG'],
        TYS=['N', 'CA', 'CB', 'CG'],
        SEP=['N', 'CA', 'CB', 'OG'],
    ),
    altchi1=dict(
        VAL=['N', 'CA', 'CB', 'CG2'],
    ),
    chi2=dict(
        ARG=['CA', 'CB', 'CG', 'CD'],
        ASN=['CA', 'CB', 'CG', 'OD1'],
        ASP=['CA', 'CB', 'CG', 'OD1'],
        GLN=['CA', 'CB', 'CG', 'CD'],
        GLU=['CA', 'CB', 'CG', 'CD'],
        HIS=['CA', 'CB', 'CG', 'ND1'],
        ILE=['CA', 'CB', 'CG1', 'CD1'],
        LEU=['CA', 'CB', 'CG', 'CD1'],
        LYS=['CA', 'CB', 'CG', 'CD'],
        MET=['CA', 'CB', 'CG', 'SD'],
        PHE=['CA', 'CB', 'CG', 'CD1'],
        PRO=['CA', 'CB', 'CG', 'CD'],
        TRP=['CA', 'CB', 'CG', 'CD1'],
        TYR=['CA', 'CB', 'CG', 'CD1'],
        TYS=['CA', 'CB', 'CG', 'CD1'],
    ),
    altchi2=dict(
        ASP=['CA', 'CB', 'CG', 'OD2'],
        LEU=['CA', 'CB', 'CG', 'CD2'],
        PHE=['CA', 'CB', 'CG', 'CD2'],
        TYR=['CA', 'CB', 'CG', 'CD2'],
        TYS=['CA', 'CB', 'CG', 'CD2'],
    ),
    chi3=dict(
        ARG=['CB', 'CG', 'CD', 'NE'],
        GLN=['CB', 'CG', 'CD', 'OE1'],
        GLU=['CB', 'CG', 'CD', 'OE1'],
        LYS=['CB', 'CG', 'CD', 'CE'],
        MET=['CB', 'CG', 'SD', 'CE'],
    ),
    chi4=dict(
        ARG=['CG', 'CD', 'NE', 'CZ'],
        LYS=['CG', 'CD', 'CE', 'NZ'],
    ),
    chi5=dict(
        ARG=['CD', 'NE', 'CZ', 'NH1'],
    ),
)

# Empirical values as defined by Tien et al. Plos ONE 2013
maxSASA = {
    "A": 121,
    "C": 148,
    "D": 187,
    "E": 214,
    "F": 228,
    "G":  97,
    "H": 216,
    "I": 195,
    "K": 230,
    "L": 191,
    "M": 203,
    "N": 187,
    "P": 154,
    "Q": 214,
    "R": 265,
    "S": 143,
    "T": 163,
    "V": 165,
    "W": 264,
    "Y": 255,
    "ALA": 121,
    "CYS": 148,
    "ASP": 187,
    "GLU": 214,
    "PHE": 228,
    "GLY":  97,
    "HIS": 216,
    "ILE": 195,
    "LYS": 230,
    "LEU": 191,
    "MET": 203,
    "ASN": 187,
    "PRO": 154,
    "GLN": 214,
    "ARG": 265,
    "SER": 143,
    "THR": 163,
    "VAL": 165,
    "TRP": 264,
    "TYR": 255
}
# TODO adjust to capture actual SASA of modified residue
maxSASA['YCM'] = maxSASA['CYS']
maxSASA['CSD'] = maxSASA['ALA']
maxSASA['TYS'] = maxSASA['TYR']
maxSASA['SEP'] = maxSASA['SER']

# Most outer residue atom
outerAtom = {
    "ALA": 'CB', # endpoint
    "CYS": 'SG', # endpoint
    "ASP": 'CG', # middle point - rotation small effect
    "GLU": 'CD', # middle point - rotation small effect
    "PHE": 'CZ', # endpoint
    "GLY": 'CA', # no sidechain
    "HIS": 'CG', # no sidechain
    "ILE": 'CD1', # outer endpoint
    "LYS": 'NZ', # endpoint
    "LEU": 'CG', # middle point - rotation small effect
    "MET": 'CE', # endpoint
    "ASN": 'CG', # middle point - flippable residue
    "PRO": 'CG', # rigid
    "GLN": 'CD', # middle point - flippable residue
    "ARG": 'CZ', # middle point - rotation small effect
    "SER": 'OG', # endpoint
    "THR": 'OG1', # endpoint donor - capture H-bond change
    "VAL": 'CB', # middle point - rotation small effect
    "TRP": 'CZ3', # second ring - capture rotation
    "TYR": 'OH' # endpoint
}
outerAtom['YCM'] = outerAtom['CYS']
outerAtom['CSD'] = outerAtom['ALA']
outerAtom['TYS'] = outerAtom['TYR']
outerAtom['SEP'] = outerAtom['SER']


class NonHetSelect(Bio.PDB.Select):
    def accept_residue(self, residue):
        return 1 if residue.id[0] == " " else 0

class Command(BaseCommand):

    help = "Command to calculate all angles for residues in each TM helix."
    tracker = {}
    all_models = django.apps.apps.get_models()[6:]
    test_model_updates(all_models, tracker, initialize=True)
    np.set_printoptions(suppress=True)
    logger = logging.getLogger(__name__)

    ###########################################################################
    ############################ Helper  Functions ############################
    ###########################################################################

    processes = 2

    def prepare_input(self, proc, items, iteration=1):
        q = Queue()
        procs = list()
        num_items = len(items)
        num = Value('i', 0)
        lock = Lock()

        if not num_items:
            return False

        # make sure not to use more jobs than proteins (chunk size will be 0, which is not good)
        if proc > num_items:
            proc = num_items

        chunk_size = int(num_items / proc)
        connection.close()
        for i in range(0, proc):
            first = chunk_size * i
            if i == proc - 1:
                last = False
            else:
                last = chunk_size * (i + 1)

            p = Process(target=self.main_func, args=([(first, last), iteration,num,lock]))
            procs.append(p)
            p.start()

        for p in procs:
            p.join()

    def add_arguments(self, parser):
        parser.add_argument('-p', '--proc',
            type=int,
            action='store',
            dest='proc',
            default=2,
            help='Number of processes to run')

    def load_pdb_var(self, pdb_code, var):
        """
        load string of pdb as pdb with a file handle. Would be nicer to do this
        directly, but no such function implemented in Bio PDB
        """
        parser = pdb.PDBParser(QUIET=True)
        with io.StringIO(var) as f:
            return parser.get_structure(pdb_code,f)


    def handle(self, *args, **options):
        if incremental_update:
            done_structures = Angle.objects.values('structure_id').distinct()
            # TODO add filter here for non-processed structures
            self.references = Structure.objects.all().exclude(structure_type__slug__startswith='af-').exclude(id__in=done_structures).prefetch_related('pdb_code','pdb_data','protein_conformation__protein','protein_conformation__state').order_by('protein_conformation__protein')
        else:
            Angle.objects.all().delete()
            Distance.objects.all().delete()
            StructureVectors.objects.all().delete()
            print("All Angle, Distance, and StructureVector data cleaned")
            self.references = Structure.objects.all().exclude(structure_type__slug__startswith='af-').prefetch_related('pdb_code','pdb_data','protein_conformation__protein','protein_conformation__state').order_by('protein_conformation__protein')

        # DEBUG for a specific PDB
        # self.references = Structure.objects.filter(pdb_code__index="2RH1").prefetch_related('pdb_code','pdb_data','protein_conformation__protein','protein_conformation__state').order_by('protein_conformation__protein')

        if 'proc' in options and options['proc']>0:
            self.processes = options['proc']

        print(len(self.references),'structures')
        self.references = list(self.references)
        self.prepare_input(self.processes, self.references)
        test_model_updates(self.all_models, self.tracker, check=True)

    def main_func(self, positions, iteration,count,lock):
        def recurse(entity,slist):
            """
            filter a pdb structure in a recursive way

            entity: the pdb entity, a structure should be given on the top level

            slist: the list of filter criterias, for each level.
            """
            for subenty in entity.get_list():
                if not subenty.id in slist[0]: entity.detach_child(subenty.id)
                elif slist[1:]: recurse(subenty, slist[1:])

        def cal_pseudo_CB(r):
            """
            Calculate pseudo CB for Glycin
            from Bio pdb faq
            """
            a =r['CA'].get_vector()
            n = r['N'].get_vector() - a
            c = r['C'].get_vector() - a
            rot = pdb.rotaxis(-np.pi*120.0/180.0, c)
            b = n.left_multiply(rot) + a
            return b.get_array()

        def pca_line(pca,h, r=0):
            """
            Calculate the pca for h and return the first pc transformed back to
            the original coordinate system
            """
            if ((not r) if pca.fit_transform(h)[0][0] < 0 else r):
                return pca.inverse_transform(np.asarray([[0,0,0],[1,0,0]]))
            else:return pca.inverse_transform(np.asarray([[0,0,0],[-1,0,0]]))

        def calc_angle(b,c):
            """
            Calculate the angle between c, b and the orthogonal projection of b
            to the x axis.
            """
            ba = -b
            bc = c + ba
            ba[:,0] = 0
            #return np.degrees(np.arccos(inner1d(ba, bc) / (np.linalg.norm(ba,axis=1) * np.linalg.norm(bc,axis=1))))

            # Alternative and clockwise angle implementation - angles left/right different value
            ba = ba[:,1:3]
            bc = bc[:,1:3]
            return np.degrees(np.arctan2(ba[:,0]*bc[:,1]-ba[:,1]*bc[:,0], inner1d(ba, bc)))

        def ca_cb_calc(ca,cb,pca):
            """
            Calculate the angles between ca, cb and center axis
            """
            return calc_angle(pca.transform(ca),pca.transform(cb))

        def ca_distance_calc(ca,pca):
            """
            Calculate the smallest distance between the ca and the center axis
            """
            return np.sqrt(np.sum(np.power(pca.transform(ca)[:,1:],2), axis = 1))

        def center_coordinates(h, p, pca):
            """
            Calculate the orthogonal projection of the CA to the helix axis
            which is moved to the mean of seven consecutive amino acids
            """
            h_middle = np.transpose(np.stack((moving_average(h[:,0], 7), moving_average(h[:,1], 7), moving_average(h[:,2], 7))))
            a = np.concatenate((h_middle[(0,0,0),:], h_middle, h_middle[(-1,-1,-1),:]))

            # b = p.transform(h)
            # b[:,1:] = p.transform(a)[:,1:]
            # b = p.inverse_transform(b)

            # TODO cleanup
            helper_lines = a - np.roll(a,-1,0)
            helper_lines[-1] = a[-1] - a[-5]
            helper_lines[-2] = helper_lines[-1]
            helper_lines[-3] = helper_lines[-1]
            helper_lines[-4] = helper_lines[-1]
            helper_lines[0] = a[0] - a[4]
            helper_lines[1] = helper_lines[0]
            helper_lines[2] = helper_lines[0]
            helper_lines[3] = helper_lines[0]
            helper_lines = np.array([ line/np.linalg.norm(line) for line in helper_lines])

            # loop over points
            #c = np.round([ a[idx] + np.dot(h_ca - a[idx], helper_lines[idx]) * helper_lines[idx] for idx, h_ca in enumerate(h)],3)
            c = [ a[idx] + np.dot(h_ca - a[idx], helper_lines[idx]) * helper_lines[idx] for idx, h_ca in enumerate(h)]


            # return b
            return c

        def axes_calc(h,p,pca):
            """
            Calculate the orthogonal projection of the CA to the helix axis
            which is moved to the mean of three consecutive amino acids
            """

            # Running average - over 3 residues at a same time
            #a = (np.roll(np.vstack((h,h[0])),1,axis=0)[:-1] + h + np.roll(np.vstack((h,h[-1])),-1,axis=0)[:-1])/3

            # Running average - over 7 residues at a same time
            h_middle = np.transpose(np.stack((moving_average(h[:,0], 7), moving_average(h[:,1], 7), moving_average(h[:,2], 7))))
            a = np.concatenate((h_middle[(0,0,0),:], h_middle, h_middle[(-1,-1,-1),:]))

            # PCA transform moved to running average
            # b = p.transform(h)
            # b[:,1:] = p.transform(a)[:,1:]
            # b = p.inverse_transform(b)

            # Use reference of 7 average as placement point
            # TODO cleanup
            helper_lines = a - np.roll(a,-1,0)
            helper_lines[-1] = a[-1] - a[-5]
            helper_lines[-2] = helper_lines[-1]
            helper_lines[-3] = helper_lines[-1]
            helper_lines[-4] = helper_lines[-1]
            helper_lines[0] = a[0] - a[4]
            helper_lines[1] = helper_lines[0]
            helper_lines[2] = helper_lines[0]
            helper_lines[3] = helper_lines[0]
            helper_lines = np.array([ line/np.linalg.norm(line) for line in helper_lines])

            # loop over points
            # b = np.round([ a[idx] + np.dot(h_ca - a[idx], helper_lines[idx]) * helper_lines[idx] for idx, h_ca in enumerate(h)],3)
            b = [ a[idx] + np.dot(h_ca - a[idx], helper_lines[idx]) * helper_lines[idx] for idx, h_ca in enumerate(h)]

            # count=0
            # for row in h:
            #     count += 1
            #     print("pseudo orig" + str(count) ,", pos=[", row[0], ",", row[1], ",", row[2] ,"];")
            #
            # count=0
            # for row in b:
            #     count += 1
            #     print("pseudo line" + str(count) ,", pos=[", row[0], ",", row[1], ",", row[2] ,"];")
            #
            # print("hide everything")
            # print("set sphere_scale, .5")
            # print("show spheres, fake*")
            # print("show spheres, orig*")
            # print("show spheres, line*")
            # print("color red, line*")
            # print("color cyan, fake*")
            # exit(0)

            return calc_angle(pca.transform(b),pca.transform(h))

        def moving_average(a, n=3) :
            ret = np.cumsum(a, dtype=float)
            ret[n:] = ret[n:] - ret[:-n]
            return ret[n - 1:] / n

        def set_bfactor(chain,angles):
            """
            simple helper to set the bfactor of all residues by some value of a
            list
            """
            for r,an in zip(chain.get_list(),angles):
                for a in r: a.set_bfactor(an)

        def qgen(x, qset):
            """
            Helper function to slice a list of all residues of a protein of the
            list of the residues of all proteins
            """
            start = False
            for i in range(len(qset)-1,0,-1):
                if not start and qset[i].protein_conformation.protein == x:
                    start = i
                if start and qset[i].protein_conformation.protein != x:
                    if start != len(qset)-1:
                        del qset[start+1:]
                        return qset[i+1:]
                    return qset[i+1:]
            del qset[start+1:]
            return qset


        def calculate_missing_atoms(poly):
            """
            Helper function to calculate missing atoms for all residues in poly
            """
            # loop over each residue
            for r in poly:
                # store missing atoms in residue annotations
                if r.resname in residue_atom_names:
                    r.xtra["MISSING"] = sum([0 if atom in r else 1 for atom in residue_atom_names[r.resname]])

        def calculate_chi_angles(poly):
            """
            Helper function to calculate all chi angles for all residues in poly
            """
            # loop over each residue
            for r in poly:
                # check for all chi_atoms if residue is present -
                chi_angles = [None] * 5
                for chi_index in range(1,6):
                    if r.resname in chi_atoms["chi" + str(chi_index)]:
                        try:
                            atom_list = chi_atoms["chi" + str(chi_index)][r.resname]
                            vec_atoms = [r[a] for a in atom_list]
                        except KeyError:
                            continue

                        vectors = [a.get_vector() for a in vec_atoms]
                        chi_angles[chi_index-1] = round(np.rad2deg(Bio.PDB.vectors.calc_dihedral(*vectors)),3)

                # store in residue annotations
                r.xtra["CHI"] = chi_angles.copy()


        #######################################################################
        ######################### Start of main loop ##########################
        #######################################################################

        failed = []
        dblist = []

        # Get all structures
        #references = Structure.objects.filter(protein_conformation__protein__family__slug__startswith="001").prefetch_related('pdb_code','pdb_data','protein_conformation__protein','protein_conformation__state').order_by('protein_conformation__protein')
        #references = Structure.objects.all().prefetch_related('pdb_code','pdb_data','protein_conformation__protein','protein_conformation__state').order_by('protein_conformation__protein')
        # DEBUG for a specific PDB
        #references = Structure.objects.filter(pdb_code__index="6AK3").prefetch_related('pdb_code','pdb_data','protein_conformation__protein','protein_conformation__state').order_by('protein_conformation__protein')

        # references = list(references)
        references = self.references

        pids = [ref.protein_conformation.protein.id for ref in references]

        qset = Residue.objects.filter(protein_conformation__protein__id__in=pids)
        if GN_only:
            qset = qset.filter(generic_number__label__regex=r'^[1-7]x[0-9]+').order_by('-protein_conformation__protein','-generic_number__label')
        else:
            qset = qset.order_by('-protein_conformation__protein','-generic_number__label')
        qset = list(qset.prefetch_related('generic_number', 'protein_conformation__protein','protein_conformation__state'))

        res_dict = {ref.pdb_code.index:qgen(ref.protein_conformation.protein,qset) for ref in references}

        #######################################################################
        ######################### Start of main loop ##########################
        #######################################################################
        angle_dict = [{},{},{},{}]
        median_dict = [{},{},{},{}]

        #for reference in references:
        while count.value<len(references):
            with lock:
                if count.value<len(references):
                    reference = references[count.value]
                    count.value +=1
                else:
                    break
            preferred_chain = reference.preferred_chain.split(',')[0]
            pdb_code = reference.pdb_code.index
#            print(pdb_code)

            try:
                structure = self.load_pdb_var(pdb_code,reference.pdb_data.pdb)
                pchain = structure[0][preferred_chain]
                state_id = reference.protein_conformation.state.id

                # DSSP
                filename = "{}_temp.pdb".format(pdb_code)
                pdbio = Bio.PDB.PDBIO()
                pdbio.set_structure(pchain)
                pdbio.save(filename, NonHetSelect())
                if os.path.exists("/env/bin/dssp"):
                    dssp = Bio.PDB.DSSP(structure[0], filename, dssp='/env/bin/dssp')
                elif os.path.exists("/env/bin/mkdssp"):
                    dssp = Bio.PDB.DSSP(structure[0], filename, dssp='/env/bin/mkdssp')
                elif os.path.exists("/usr/local/bin/mkdssp"):
                    dssp = Bio.PDB.DSSP(structure[0], filename, dssp='/usr/local/bin/mkdssp')

                # DISABLED STRIDE - selected DSSP 3 over STRIDE
#                try:
#                    if os.path.exists("/env/bin/stride"):
#                       stride = subprocess.Popen(['/env/bin/stride', filename], stdout=subprocess.PIPE)
#                       # Grab SS assignment (ASG) and parse residue (cols 12-15) and SS (cols 25-25)
#                       for line in io.TextIOWrapper(stride.stdout, encoding="utf-8"):
#                           if line.startswith("ASG"):
#                               res_id = int(line[11:15].strip())
#                               res_ss = line[24:25].strip()
#                               # assign to residue
#                               pchain[res_id].xtra["SS_STRIDE"] = res_ss.upper()
#                except OSError:
#                   print(pdb_code, " - STRIDE ERROR - ", e)

                # CLEANUP
                os.remove(filename)

                #######################################################################
                ###################### prepare and evaluate query #####################

                db_reslist = res_dict[pdb_code]
                gn_reslist = []
                tm_reslist = []
                for i in db_reslist:
                    if i.generic_number:
                        gn_reslist.append(i)
                        if re.match(r'^[1-7]x[0-9]+', i.generic_number.label):
                            tm_reslist.append(i)


                full_resdict = {str(r.sequence_number):r for r in db_reslist}

                #######################################################################
                ######################### filter data from db #########################

                def reslist_gen(x):
                    try:
                        while tm_reslist[-1].generic_number.label[0] == x:
                            yield tm_reslist.pop()
                    except IndexError:
                        pass

                # Fix IDs matching for handling PTM-ed residues
                ids_in_pchain = []
                for residue in pchain:
                    if residue.id[1] not in pchain:
                        residue.id = (' ', residue.id[1], ' ')

                # when gdict is not needed the helper can be removed
                db_helper = [[(r,r.sequence_number) for r in reslist_gen(x) if r.sequence_number in pchain] for x in ["1","2","3","4","5","6","7"]]
                gdict = {r[1]:r[0] for hlist in db_helper for r in hlist}
                tm_keys = [r[1] for hlist in db_helper for r in hlist]
                tm_keys_str = [str(i) for i in tm_keys]
                tm_keys_int = [int(i) for i in tm_keys]
                db_tmlist = [[(' ',r[1],' ') for r in sl] for sl in db_helper]
                db_set = set(db_tmlist[0]+db_tmlist[1]+db_tmlist[2]+db_tmlist[3]+db_tmlist[4]+db_tmlist[5]+db_tmlist[6])

                #######################################################################
                ##################### Angles/dihedrals residues #######################

                #for line in pdblines_temp: #Get rid of all odd records
                #polychain = [ residue for residue in pchain if Bio.PDB.Polypeptide.is_aa(residue) and "CA" in residue]
                polychain = [ residue for residue in pchain if (Bio.PDB.Polypeptide.is_aa(residue) or residue.resname in ['YCM','CSD','TYS','SEP']) and "CA" in residue]
                poly = Bio.PDB.Polypeptide.Polypeptide(polychain)

                # Calculate backbone and sidechain dihedrals + missing atoms
                poly.get_phi_psi_list() # backbone dihedrals
                calculate_chi_angles(poly) # calculate chi1-chi5
                calculate_missing_atoms(poly) # calculate missing

                # Possibly only relevant for helices?
                poly.get_theta_list() # angle three consecutive Ca atoms
                poly.get_tau_list() # dihedral four consecutive Ca atoms

                ### DEPRECATED: clean the structure to solely the 7TM bundle
                #recurse(structure, [[0], preferred_chain, db_set])

                ### NEW: clean the structure to all protein residues in DB
                db_fullset = set([(' ',r.sequence_number,' ') for r in db_reslist])
                recurse(structure, [[0], preferred_chain, db_fullset])

                dihedrals = {}
                for r in poly:
                  angle_list = ["PHI", "PSI", "THETA", "TAU", "SS_DSSP", "SS_STRIDE", "CHI", "MISSING"]
                  for angle in angle_list:
                      if angle not in r.xtra:
                          r.xtra[angle] = None

                  # Add outer angle
                  outer = None
                  try:
                      angle_atoms = [r[a].get_vector() for a in ['N','CA', outerAtom[r.resname]]]

                      # use pseudo CB placement when glycine (or in case of missing CB)
                      #if r.resname == 'GLY':
                      if 'CB' not in r:
                          angle_atoms[2] = Bio.PDB.vectors.Vector(*cal_pseudo_CB(r))

                      outer = Bio.PDB.calc_angle(*angle_atoms)
                  except Exception as e:
#                      print(pdb_code, " - OUTER ANGLE ERROR - ", e)
                      outer = None

                  # Add tau (N-Ca-C) backbone angle (in addition to the tau dihedral)
                  tau_angles = None
                  try:
                      angle_atoms = [r[a].get_vector() for a in ['N','CA', 'C']]

                      tau_angles = Bio.PDB.calc_angle(*angle_atoms)
                  except Exception as e:
#                      print(pdb_code, " - TAU ANGLE ERROR - ", e)
                      tau_angles = None

                  dihedrals[str(r.id[1])] = [r.xtra["PHI"], r.xtra["PSI"], r.xtra["THETA"], r.xtra["TAU"], r.xtra["SS_DSSP"], r.xtra["SS_STRIDE"], outer, tau_angles, r.xtra["CHI"], r.xtra["MISSING"]]

                # Extra: remove hydrogens from structure (e.g. 5VRA)
                for residue in structure[0][preferred_chain]:
                    for id in [atom.id for atom in residue if atom.element == "H"]:
                        residue.detach_child(id)

                # List of CA coordinates for all residues with GN for distances
                gn_res_gns = [res.generic_number.label for res in gn_reslist]
                gn_res_ids = [res.sequence_number for res in gn_reslist]

                # Order by GNs
                gns_order = []
                for idx, gn in enumerate(gn_res_gns):
                    part1, part2 = gn.split("x")
                    if part1.isnumeric():
                        multiply1 = 1000 if len(part1)>=2 else 10000
                        multiply2 = 10 if (len(part2))<=2 else 1
                        gns_order.append(int(part1)*multiply1 + int(part2)*multiply2)
                    else:
                        # When not a numerical GN (e.g. D1e1) - just use seq num
                        gns_order.append(gn_res_ids[idx])


                gns_ids_list = [gn_res_ids[key] for key in np.argsort(gns_order)]

                #gn_res_gns = [gn_res_gns[key] for key in np.argsort(gns_order)]

                #gns_ca_list = {resid:pchain[resid]["CA"].get_coord() for resid in gns_ids_list if resid in pchain}
                gns_ca_list = {residue.id[1]:residue["CA"].get_coord() for residue in poly if residue.id[1] in gns_ids_list}
                #gns_cb_list = {resid:np.asarray(pchain[resid]["CB"].get_coord() if "CB" in pchain[resid] else cal_pseudo_CB(pchain[resid]), dtype=float) for resid in gns_ids_list if resid in pchain}
                gns_cb_list = {residue.id[1]:np.asarray(residue["CB"].get_coord() if "CB" in residue else cal_pseudo_CB(residue), dtype=float) for residue in poly if residue.id[1] in gns_ids_list}

                ### AXES through each of the TMs and the TM bundle (center axis)
                hres_list = [np.asarray([pchain[r]["CA"].get_coord() for r in sl], dtype=float) for sl in db_tmlist]
                #h_cb_list = [np.asarray([pchain[r]["CB"].get_coord() if "CB" in pchain[r] else cal_pseudo_CB(pchain[r]) for r in sl], dtype=float) for sl in db_tmlist]
                h_cb_list = [[gns_cb_list[r[1]] for r in sl] for sl in db_tmlist]

                # fast and fancy way to take the average of N consecutive elements
                N = 3

                # NOTE - can result in ragged nested sequences
                hres_three = np.asarray([sum([h[i:-(len(h) % N) or None:N] for i in range(N)])/N for h in hres_list], dtype=object)

                ### PCA - determine axis through center + each transmembrane helix
                helix_pcas = [PCA() for i in range(7)]
                helix_pca_vectors = [pca_line(helix_pcas[i], h,i%2) for i,h in enumerate(hres_three)]

                # Calculate PCA based on the upper (extracellular) half of the GPCR (more stable, except class B)
                pca = PCA()
                pos_list = []
                if extra_pca:
                    minlength = 100
                    for i,h in enumerate(hres_three):
                        if len(h)<minlength:
                            minlength = len(h)

                    if minlength > 6:
                        minlength = 6

                    # create PCA per helix using extracellular half
                    # Exclude the first turn if possible (often still part of loop)
                    for i,h in enumerate(hres_three):
                        if i%2: # reverse directionality of even helices (TM2, TM4, TM6)
                            h = np.flip(h, 0)

                        if len(h)>minlength+2:
                            pos_list.append(pca_line(PCA(), h[2:minlength+2]))
                        else:
                            pos_list.append(pca_line(PCA(), h[0:minlength]))


                    # create fake coordinates along each helix PCA to create center PCA
                    # UGLY hack - should be cleaned up
                    coord_list = []
                    for pos in pos_list:
                        start = pos[0]
                        vector = pos[1]-pos[0]
                        line_points = []
                        for i in range(-45,55):
                            line_points.append(start+i*vector)

                        coord_list.append(line_points)
                    center_vector = pca_line(pca, np.vstack(coord_list))
                else:
                    # Create PCA line through whole helix
                    # NOTE: much less robust with differing TM lengths, bends, kinks, etc.
                    center_vector = pca_line( pca, np.vstack(hres_three))

                # DEBUG print arrow for PyMol
                # a = [str(i) for i in center_vector[0]]
                # b = [str(i) for i in center_vector[1]]
                # print("cgo_arrow [" + a[0] + ", " + a[1] + ", " + a[2] + "], [" + b[0] + ", " + b[1] + ", " + b[2] + "]")

                # Measure level of activation by TM6 tilt
                # Residue most often found at kink start
                # TODO: check for numbering at other classes
                # kink_start = 44 #  general class A number 6x44 seems quite conserved to be the kink start

                # Select all residues before indicated residue
                # lower_tm6 = []
                # kink_start_res = None
                # kink_measure = None
                # for res in db_tmlist[5]:
                #     gnlabel = gdict[res[1]].generic_number.label
                #     if int(gnlabel.replace("6x","")) <= kink_start:
                #         lower_tm6.append(pchain[res]["CA"].get_coord())
                #         if int(gnlabel.replace("6x","")) == kink_start:
                #             kink_start_res = pchain[res]["CA"].get_coord()
                #         if int(gnlabel.replace("6x","")) == 38:
                #             kink_measure = pchain[res]["CA"].get_coord()
                #
                # lower_tm6 = np.asarray(lower_tm6)

                # TM2 intracellular for comparison
                # lower_tm2 = []
                # ref_tm2 = None
                # for res in db_tmlist[1]:
                #     gnlabel = gdict[res[1]].generic_number.label
                #     gn_id = int(gnlabel.replace("2x",""))
                #     if gn_id >= 40 and gn_id <= 50: # Lower well-defined half of TM2
                #         lower_tm2.append(pchain[res]["CA"].get_coord())
                #         if gn_id == 41:
                #             ref_tm2 = pchain[res]["CA"].get_coord()
                # lower_tm2 = np.asarray(lower_tm2)

                # posb_list = []
                # # create PCA per helix using full helix
                # for i,h in enumerate(hres_three):
                #     if i%2: # reverse directionality of even helices (TM2, TM4, TM6)
                #         h = np.flip(h, 0)
                #     posb_list.append(pca_line(PCA(), h))


                # NOTE: Slight variations between the mid membrane residues can have a strong affect on the plane
                # For now just use a single residue to deduce the mid membrane height (next code)
                # if len(kink_measure) == 3:
                #     # Use membrane middle references 1x44 - 2x52 - 4x54
                #     membrane_mid = []
                #     for res in gdict:
                #         if gdict[res].generic_number.label in ["1x44", "2x52", "4x54"]:
                #             membrane_mid.append(pchain[res]["CA"].get_coord())
                #
                #     if len(membrane_mid) == 3:
                #         v1 = membrane_mid[1] - membrane_mid[0]
                #         v2 = membrane_mid[2] - membrane_mid[0]
                #         plane_normal = np.cross(v1 / np.linalg.norm(v1), v2 / np.linalg.norm(v2))
                #         plane_normal = plane_normal / np.linalg.norm(plane_normal)
                #
                #         rayDirection = center_vector[1] - center_vector[0]
                #         ndotu = plane_normal.dot(rayDirection)
                #         w = center_vector[0] - membrane_mid[0]
                #         si = -plane_normal.dot(w) / ndotu
                #         membrane_point = w + si * rayDirection + membrane_mid[0]
                #
                #         # calculate distances to this point
                #         midpoint_distances = np.round([ np.linalg.norm(membrane_point-ca) for helix in hres_list for ca in helix ],3)

                # Find 5x46 (residue at membrane middle)
                membrane_mid = None
                for res in db_tmlist[4]:
                    gnlabel = gdict[res[1]].generic_number.label
                    if gnlabel == "5x46":
                        membrane_mid = pchain[res]["CA"].get_coord()
                        break

                # Calculate distances to the mid of membrane plane
                if len(membrane_mid) == 3:
                    # 1. Find intersect of membrane mid with 7TM axis (project point to plane)
                    membrane_mid_pca = pca.transform([membrane_mid])
                    membrane_mid_pca[0,1:3] = 0 # project onto the same axis
                    membrane_point = pca.inverse_transform(membrane_mid_pca)
                    plane_normal = membrane_point - center_vector[1]
                    plane_normal = plane_normal / np.linalg.norm(plane_normal)

                    # calculate distances to the mid of membrane plane
                    mid_membrane_distances = np.round([ np.dot(ca - membrane_point[0], plane_normal[0]) for helix in hres_list for ca in helix ],3)

                    # calculate distances to the midpoint
                    midpoint_distances = np.round([ np.linalg.norm(membrane_point[0]-ca) for helix in hres_list for ca in helix ],3)


                # TM6 tilt with respect to 7TM bundle axis and plane through 6x44
                # if len(lower_tm6) >= 3 and len(membrane_mid) == 3:
                #     # Take the average of N consecutive elements
                #     tm6_lower_three = sum([lower_tm6[i:-(len(lower_tm6) % N) or None:N] for i in range(N)])/N
                #     if len(tm6_lower_three) > 2:
                #         tm6_pca_vector = pca_line(PCA(), tm6_lower_three, 1)
                #     else:
                #         tm6_pca_vector = pca_line(PCA(), lower_tm6, 1)
                #
                #     # 1. Find intersect of membrane mid with 7TM axis (project point to plane)
                #     membrane_mid_pca = pca.transform([membrane_mid])
                #     membrane_mid_pca[0,1:3] = 0 # project onto the same axis
                #     midpoint = pca.inverse_transform(membrane_mid_pca)

                # Distance TM2-3-4-5-6-7 at height of 6x38 using Mid membrane as plane
                # if len(kink_measure) == 3:
                #     # 1x44 - 2x52 - 4x54
                #     membrane_mid = []
                #     for res in gdict:
                #         if gdict[res].generic_number.label in ["1x44", "2x52", "4x54"]:
                #             membrane_mid.append(pchain[res]["CA"].get_coord())
                #
                #     if len(membrane_mid) == 3:
                #         v1 = membrane_mid[1] - membrane_mid[0]
                #         v2 = membrane_mid[2] - membrane_mid[0]
                #         plane_normal = np.cross(v1 / np.linalg.norm(v1), v2 / np.linalg.norm(v2))
                #         planeNormal = plane_normal / np.linalg.norm(plane_normal)
                #
                #         points = []
                #         for i in [1,2,4,5,6]: # TM number - 1
                #             rayDirection = posb_list[i][1] - posb_list[i][0]
                #             ndotu = planeNormal.dot(rayDirection)
                #             w = posb_list[i][0] - kink_measure
                #             si = -planeNormal.dot(w) / ndotu
                #             intersect = w + si * rayDirection + kink_measure
                #             points.append(intersect)
                #
                #         distance = 0
                #         last = points[len(points)-1]
                #         for point in points:
                #             distance += np.linalg.norm(last - point)
                #             last = point
                #
                #         print("MIDMEM DISTANCE {} {}".format(pdb_code, distance))
                #         #reference.tm6_angle = distance
                #         #reference.save()

                # Distance based on distance pairs
                # if len(kink_measure) == 3:
                #     # 1x50 2x41 3x26 4x42 5x42 6x37 7x49
                #     points = []
                #     for gn in ["1x50", "2x41", "3x26", "4x42", "5x42", "6x37", "7x49"]:
                #         res = [key for (key, value) in gdict.items() if value.generic_number.label == gn]
                #         if len(res) > 0:
                #             points.append(pchain[res[0]]["CA"].get_coord())
                #
                #     if len(points) != 7:
                #         continue
                #
                #     distance = 0
                #     last = points[len(points)-1]
                #     for point in points:
                #         distance += np.linalg.norm(last - point)
                #         last = point
                #
                #     print("PAIRS DISTANCE {} {}".format(pdb_code, distance))


                # Distances between TM points
                # if len(kink_measure) == 3:
                #     minlength = 100
                #     posb_list = []
                #
                #     # create PCA per helix using full helix
                #     for i,h in enumerate(hres_three):
                #         if i%2: # reverse directionality of even helices (TM2, TM4, TM6)
                #             h = np.flip(h, 0)
                #
                #         posb_list.append(pca_line(PCA(), h))
                #
                #     points = []
                #     for i in range(7): # TM number - 1
                #         rayDirection = posb_list[i][1] - posb_list[i][0]
                #         ndotu = planeNormal.dot(rayDirection)
                #         w = pos_list[i][0] - kink_measure
                #         si = -planeNormal.dot(w) / ndotu
                #         intersect = w + si * rayDirection + kink_measure
                #         points.append(intersect)
                #
                #
                #     hstr = ""
                #     dstr = ""
                #     for i in range(len(points)):
                #         for j in range(i+1,len(points)):
                #             hstr += "," + str(i+1) + "x" + str(j+1)
                #             dstr += "," + str(np.linalg.norm(points[i] - points[j]))
                #
                #     print("HEADER {}".format(hstr))
                #     print(pdb_code + dstr)

                        #print("REFERENCE {} {}".format(pdb_code, distance))

                # Distance TM2-3-4-5-6-7 at height of 6x38 using TM bundle axis as plane normal
                # if len(kink_measure) == 3:
                #     planeNormal = center_vector[0]-center_vector[1]
                #
                #     points = []
                #     for i in [1,2,4,5,6]: # TM number - 1
                #         rayDirection = posb_list[i][1] - posb_list[i][0]
                #         ndotu = planeNormal.dot(rayDirection)
                #         w = posb_list[i][0] - kink_measure
                #         si = -planeNormal.dot(w) / ndotu
                #         intersect = w + si * rayDirection + kink_measure
                #         points.append(intersect)
                #
                #     distance = 0
                #     last = points[len(points)-1]
                #     for point in points:
                #         distance += np.linalg.norm(last - point)
                #         last = point
                #
                #     print("BUNDLEAXIS DISTANCE {} {}".format(pdb_code, distance))
                #     #reference.tm6_angle = distance
                #     #reference.save()
                #     #print("REFERENCE {} {}".format(pdb_code, distance))

                # TM6 tilt compared to lower TM2 using pca vectors
#                 if len(lower_tm6) >= 3 and len(lower_tm2) >= 3:
#                      # Take the average of N consecutive elements of TM2
#                      tm2_lower_three = sum([lower_tm2[i:-(len(lower_tm2) % N) or None:N] for i in range(N)])/N
#                      if len(tm2_lower_three) > 2:
#                          tm2_pca_vector = pca_line(PCA(), tm2_lower_three, 1)
#                      else:
#                          tm2_pca_vector = pca_line(PCA(), lower_tm2, 1)
#
#                      # Take the average of N consecutive elements of TM6
#                      tm6_lower_three = sum([lower_tm6[i:-(len(lower_tm6) % N) or None:N] for i in range(N)])/N
#                      if len(tm6_lower_three) > 2:
#                          tm6_pca_vector = pca_line(PCA(), tm6_lower_three, 1)
#                      else:
#                          tm6_pca_vector = pca_line(PCA(), lower_tm6, 1)
#
#                      membrane_mid_pca = pca.transform([membrane_mid])
#                      membrane_mid_pca[0,1:3] = 0 # project onto the same axis
#                      midpoint = pca.inverse_transform(membrane_mid_pca)[0]
#
#
#                      # Project pca vector onto plane to TM2
#                      # Find normal of the plane through Axis and 2x41
#                      v1 = center_vector[1] - midpoint
#                      v2 = ref_tm2 - midpoint
#                      plane_normal = np.cross(v1 / np.linalg.norm(v1), v2 / np.linalg.norm(v2))
#                      plane_normal = plane_normal / np.linalg.norm(plane_normal)
#
#                      # Find projected tm6 angle to plane
#                      displaced_point = tm6_pca_vector[1] - tm6_pca_vector[0]
#                      dist = np.dot(displaced_point, plane_normal)
#                      proj_point = (midpoint+displaced_point) - dist*plane_normal
#                      tm6_tilt_proj = proj_point - midpoint
# #                     tm6_tilt_proj = tm6_tilt_proj/np.linalg.norm(tm6_tilt_proj)
#
#                      # Find projected tm2 angle to plane
#                      displaced_point = tm2_pca_vector[1] - tm2_pca_vector[0]
#                      dist = np.dot(displaced_point, plane_normal)
#                      proj_point = (midpoint[0]+displaced_point) - dist*plane_normal
#                      tm2_tilt_proj = proj_point - midpoint[0]
# #                     tm2_tilt_proj = tm2_tilt_proj/np.linalg.norm(tm2_tilt_proj)
#
#                      # angle between these vectors
#                      tm6_angle = np.arccos(np.dot(tm6_tilt_proj, tm2_tilt_proj)/(np.linalg.norm(tm6_tilt_proj)*np.linalg.norm(tm2_tilt_proj)))
#                      print("REFERENCE {} {}".format(pdb_code, np.degrees(tm6_angle)))
#
#                      # Store as structure property
#                      reference.tm6_angle = np.degrees(tm6_angle)
#                      reference.save()

                # # Angle of 6x38 to 2x41 via midpoint in membrane
                # if len(ref_tm2) == 3 and len(kink_start_res) == 3:
                #     membrane_mid_pca = pca.transform([membrane_mid])
                #     membrane_mid_pca[0,1:3] = 0 # project onto the same axis
                #     midpoint = pca.inverse_transform(membrane_mid_pca)
                #
                #     v1 = ref_tm2 - midpoint[0]
                #     v2 = kink_start_res - midpoint[0]
                #
                #     # angle between these vectors
                #     tm6_angle = np.arccos(np.dot(v1, v2)/(np.linalg.norm(v1)*np.linalg.norm(v2)))
                #     print(np.degrees(tm6_angle))
                #
                #     # Store as structure property
                #     reference.tm6_angle = np.degrees(tm6_angle)
                #     reference.save()

                # TM6 tilt compared to lower TM2 using pca vectors
                # if len(lower_tm6) >= 3 and len(lower_tm2) >= 3:
                #      # Take the average of N consecutive elements of TM2
                #      tm2_lower_three = sum([lower_tm2[i:-(len(lower_tm2) % N) or None:N] for i in range(N)])/N
                #      if len(tm2_lower_three) > 2:
                #          tm2_pca_vector = pca_line(PCA(), tm2_lower_three, 1)
                #      else:
                #          tm2_pca_vector = pca_line(PCA(), lower_tm2, 1)
                #
                #      # Take the average of N consecutive elements of TM6
                #      tm6_lower_three = sum([lower_tm6[i:-(len(lower_tm6) % N) or None:N] for i in range(N)])/N
                #      if len(tm6_lower_three) > 2:
                #          tm6_pca_vector = pca_line(PCA(), tm6_lower_three, 1)
                #      else:
                #          tm6_pca_vector = pca_line(PCA(), lower_tm6, 1)
                #
                #      # angle between these vectors
                #      tm6_angle = np.arccos(np.clip(np.dot(tm6_pca_vector[1]-tm6_pca_vector[0], tm2_pca_vector[1]-tm2_pca_vector[0]), -1.0, 1.0))
                #      print(tm6_angle)
                #
                #      # Store as structure property
                #      reference.tm6_angle = np.degrees(tm6_angle)
                #      reference.save()


                # TM6 tilt with respect to 7TM bundle axis and plane through 6x44
                # if len(lower_tm6) >= 3 and len(membrane_mid) == 3:
                #     # Take the average of N consecutive elements
                #     tm6_lower_three = sum([lower_tm6[i:-(len(lower_tm6) % N) or None:N] for i in range(N)])/N
                #     if len(tm6_lower_three) > 2:
                #         tm6_pca_vector = pca_line(PCA(), tm6_lower_three, 1)
                #     else:
                #         tm6_pca_vector = pca_line(PCA(), lower_tm6, 1)
                #
                #     # 1. Find intersect of membrane mid with 7TM axis (project point to plane)
                #     membrane_mid_pca = pca.transform([membrane_mid])
                #     membrane_mid_pca[0,1:3] = 0 # project onto the same axis
                #     midpoint = pca.inverse_transform(membrane_mid_pca)
                #
                #     # 2. Find normal of plane through origin, kink start and project kink start
                #     #    A) Find projected point of kink_start
                #     kink_start_res_pca = pca.transform([kink_start_res])
                #     kink_start_res_pca[0,1:3] = 0 # project onto the same axis
                #     kink_start_proj = pca.inverse_transform(kink_start_res_pca)
                #
                #     #    B) Find normal of the new plane through kink start
                #     v1 = kink_start_res - center_vector[0]
                #     v2 = kink_start_proj - center_vector[0]
                #     plane_normal = np.cross(v1 / np.linalg.norm(v1), v2 / np.linalg.norm(v2))[0]
                #     plane_normal = plane_normal / np.linalg.norm(plane_normal)
                #
                #     #    C) Find projected tm6 angle to plane
                #     displaced_point = tm6_pca_vector[1] - tm6_pca_vector[0]
                #     dist = np.dot(displaced_point, plane_normal)
                #     proj_point = (center_vector[0]+displaced_point) - dist*plane_normal
                #     tm6_tilt_proj = proj_point - center_vector[0]
                #
                #     # Calculate angle between vectors
                #     tm6_angle = np.arccos(np.dot(tm6_tilt_proj,center_vector[1]-center_vector[0])/(np.linalg.norm(tm6_tilt_proj)*np.linalg.norm(center_vector[1]-center_vector[0])))
                #
                #     # Check change in distance for projected angle point on tm axis
                #     # distance increased? -> negative angle - distance decreased -> positive angle
                #     distance_kink_start = np.linalg.norm(kink_start_res - center_vector[0])
                #
                #     # VERIFY: not sure if correct
                #     dist = np.dot(tm6_tilt_proj, center_vector[1] - center_vector[0])
                #     proj_proj_tilt_point = (center_vector[0]+tm6_tilt_proj) - dist*(center_vector[1] - center_vector[0])
                #     distance_tilt = np.linalg.norm(kink_start_res - proj_proj_tilt_point)
                #     if (distance_tilt-distance_kink_start) > 0:
                #         tm6_angle = -1*tm6_angle
                #
                #     # Store as structure property
                #     reference.tm6_angle = np.degrees(tm6_angle)
                #     reference.save()

                # TM6 tilt with respect to 7TM bundle axis
                # if len(lower_tm6) >= 3:
                #     # Take the average of N consecutive elements
                #     tm6_lower_three = sum([lower_tm6[i:-(len(lower_tm6) % N) or None:N] for i in range(N)])/N
                #     if len(tm6_lower_three) > 2:
                #         tm6_pca_vector = pca_line(PCA(), tm6_lower_three, 1)
                #     else:
                #         tm6_pca_vector = pca_line(PCA(), lower_tm6, 1)
                #
                #     # Calculate angle between vectors
                #     tm6_angle = np.arccos(np.clip(np.dot(tm6_pca_vector[1]-tm6_pca_vector[0], center_vector[1]-center_vector[0]), -1.0, 1.0))
                #
                #     cen_tm6 = center_vector[0]+(tm6_pca_vector[1]-tm6_pca_vector[0])
                #     print("pseudoatom center1, pos=[{},{},{}]".format(center_vector[0][0],center_vector[0][1],center_vector[0][2]))
                #     print("pseudoatom center2, pos=[{},{},{}]".format(center_vector[1][0],center_vector[1][1],center_vector[1][2]))
                #     print("pseudoatom cen_tm6, pos=[{},{},{}]".format(cen_tm6[0],cen_tm6[1],cen_tm6[2]))
                #
                #     # Check distance - closer to axis - negative angle - further away - positive angle
                #     tm6_inward = ca_distance_calc(np.asarray([tm6_pca_vector[1]]),pca) - ca_distance_calc(np.asarray([tm6_pca_vector[0]]),pca)
                #     print(np.degrees(tm6_angle))
                #
                #     if tm6_inward < 0:
                #         tm6_angle = -1*tm6_angle
                #
                #     # Store as structure property
                #     reference.tm6_angle = np.degrees(tm6_angle)
                #     reference.save()

                c_vector = np.array2string(center_vector[0] - center_vector[1], separator=',')
                translation = np.array2string(-1*center_vector[0], separator=',')

                StructureVectors.objects.filter(structure = reference).all().delete()
                sv = StructureVectors(structure = reference, translation = str(translation), center_axis = str(c_vector))
                sv.save()

                # TODO:
                # FIX RESIDUE ORDER
                # FIX requirement checking

                ### DISTANCES - moved here from cube
                # REMOVE OLD distances
                # Distance.objects.filter(structure=reference).all().delete()

                # Perpendicular projection of Ca onto helical PCA
                h_center_list = np.concatenate([center_coordinates(h,p,pca) for h,p in zip(hres_list,helix_pcas)])
                gns_center_list = dict(zip(tm_keys_int, h_center_list))

                # New rotation angle
                # Angle between normal from center axis to 1x46 and normal from helix axis to CA
                key_tm1 = gn_res_ids[gn_res_gns.index("1x46")]
                ref_tm1 = gns_center_list[key_tm1]
                # print(gns_order)
                # print(np.argsort(gns_order))
                # print(gn_res_gns)
                # print(gn_res_ids)
                # print(key_tm1)

                # Project 1x46 onto center axis
                axis_vector = (center_vector[0] - center_vector[1])/np.linalg.norm(center_vector[0] - center_vector[1])
                center_tm1 = center_vector[0] + np.dot(ref_tm1 - center_vector[0], axis_vector) * axis_vector
                tm1_vector = (center_tm1 - ref_tm1)/np.linalg.norm(center_tm1 - ref_tm1)

                # Calculate CA to helix center vectors
                ca_center_vectors = {resid:(gns_ca_list[resid] - gns_center_list[resid])/np.linalg.norm(gns_ca_list[resid] - gns_center_list[resid]) for resid in gns_center_list }

                # Calculate rotation angles
                rotation_angles = { resid:np.rad2deg(np.arccos(np.dot(tm1_vector, ca_center))) for resid, ca_center in ca_center_vectors.items() }

                # Rotate tm1_vector by 90 degrees and then check the angles again  - if angle gets larger -> 360 - angle otherwise ok
                rotation_vector = np.radians(90) * axis_vector
                rotation = R.from_rotvec(rotation_vector)
                rotated_tm1_vector = rotation.apply(tm1_vector)
                rotation_angles_ref = { resid:np.rad2deg(np.arccos(np.dot(rotated_tm1_vector, ca_center))) for resid, ca_center in ca_center_vectors.items() }
                # Make key a string to match with other dictionaries
                rotation_angles = {str(resid):(round(rotation_angles[resid],3) if rotation_angles_ref[resid] - rotation_angles[resid] < 0 else round(360 - rotation_angles[resid],3)) for resid in rotation_angles }
                # Making the rotation angle compliant with the other angles (-180 to 180 degrees)
                rotation_angles = {resid:rotation_angles[resid]-180 for resid in rotation_angles }
                # print(pdb_code, "1x45", rotation_angles[str(gn_res_ids[gn_res_gns.index("1x45")])], "and 1x47", rotation_angles[str(gn_res_ids[gn_res_gns.index("1x47")])])


                # print("pseudo center, pos=[", ref_tm1[0], ",", ref_tm1[1], ",", ref_tm1[2] ,"];")
                # print("pseudo ca, pos=[", gns_ca_list[key_tm1][0], ",", gns_ca_list[key_tm1][1], ",", gns_ca_list[key_tm1][2] ,"];")
                # print("pseudo mid, pos=[", center_tm1[0], ",", center_tm1[1], ",", center_tm1[2] ,"];")
                # print(rotation_angles[key_tm1])

                # triangular matrix for distances
                up_ind = np.triu_indices(len(gns_ca_list), 1)
                bulk_distances = []

                for i1, i2 in zip(up_ind[0], up_ind[1]):
                    key1 = gns_ids_list[i1]
                    key2 = gns_ids_list[i2]
                    res1 = full_resdict[str(key1)]
                    res2 = full_resdict[str(key2)]

                    ca_dist = int(np.linalg.norm(gns_ca_list[key1] - gns_ca_list[key2])*distance_scaling_factor)
                    cb_dist = int(np.linalg.norm(gns_cb_list[key1] - gns_cb_list[key2])*distance_scaling_factor)
                    center_dist = None
                    if key1 in gns_center_list and key2 in gns_center_list:
                        center_dist = int(np.linalg.norm(gns_center_list[key1] - gns_center_list[key2])*distance_scaling_factor)

                    # residues in gn_reslist, structure in structure
                    distance = Distance(distance = ca_dist, distance_cb = cb_dist, distance_helix_center = center_dist, res1=res1, res2=res2, gn1=res1.generic_number.label, gn2=res2.generic_number.label, gns_pair='_'.join([res1.generic_number.label, res2.generic_number.label]), structure=reference)
                    bulk_distances.append(distance)

                # Bulk insert
                Distance.objects.bulk_create(bulk_distances, batch_size=5000)

                ### ANGLES
                # Center axis to helix axis to CA
                a_angle = np.concatenate([axes_calc(h,p,pca) for h,p in zip(hres_list,helix_pcas)]).round(3)

                # Center axis to CA to CB
                b_angle = np.concatenate([ca_cb_calc(ca,cb,pca) for ca,cb in zip(hres_list,h_cb_list)]).round(3)

                # Distance from center axis to CA
                core_distances = np.concatenate([ca_distance_calc(ca,pca) for ca in hres_list]).round(3)

                ### freeSASA (only for TM bundle)
                # SASA calculations - results per atom
                clean_structure = self.load_pdb_var(pdb_code,reference.pdb_data.pdb)
                clean_pchain = clean_structure[0][preferred_chain]

                # PTM residues give an FreeSASA error - remove
                db_fullset = set([(' ',r.sequence_number,' ') for r in db_reslist])
                recurse(clean_structure, [[0], preferred_chain, db_fullset])
                # Remove hydrogens from structure (e.g. 5VRA)
                for residue in clean_structure[0][preferred_chain]:
                    for id in [atom.id for atom in residue if atom.element == "H"]:
                        residue.detach_child(id)

                res, trash = freesasa.calcBioPDB(clean_structure)

                # create results dictionary per residue
                asa_list = {}
                rsa_list = {}
                atomlist = list(clean_pchain.get_atoms())
                for i in range(res.nAtoms()):
                    resnum = str(atomlist[i].get_parent().id[1])
                    if resnum not in asa_list:
                        asa_list[resnum] = 0
                        rsa_list[resnum] = 0

                    resname = atomlist[i].get_parent().get_resname()
                    if resname in maxSASA:
                        rsa_list[resnum] += res.atomArea(i)/maxSASA[resname]*100
                    else:
                        rsa_list[resnum] = None

                    asa_list[resnum] += res.atomArea(i)

                # correct for N/C-term exposure
                for i in rsa_list:
                    if (rsa_list[i]>100):
                        rsa_list[i] = 100

                ### Half-sphere exposure (HSE)
                hse = pdb.HSExposure.HSExposureCB(structure[0][preferred_chain])

                # x[1] contains HSE - 0 outer half, 1 - inner half, 2 - ?
                hselist = dict([ (str(x[0].id[1]), x[1][0]) if x[1][0] > 0 else (str(x[0].id[1]), 0) for x in hse ])

                # Few checks
                if GN_only:
                    if len(pchain) != len(a_angle):
                        raise Exception("\033[91mLength mismatch a-angles " + pdb_code + "\033[0m")

                        if len(pchain) != len(b_angle):
                            raise Exception("\033[91mLength mismatch b-angles " + pdb_code + "\033[0m")

                ### Collect all data in database list
                #print(a_angle) # only TM
                #print(b_angle) # only TM
                #print(asa_list) # only TM
                #print(hselist) # only TM
                #print(dihedrals) # HUSK: contains full protein!

                ### PCA space can be upside down - in that case invert the results
                # Check rotation of 1x49 - 1x50
                inversion_ref = -1
                for res in tm_keys:
                    inversion_ref += 1
                    if gdict[res].generic_number.label == "1x49":
                        break

                signed_diff = (a_angle[inversion_ref + 1] - a_angle[inversion_ref] + 540 ) % 360 - 180
                if signed_diff > 0:
#                     print("{} Rotating the wrong way {}".format(pdb_code, signed_diff))
                     a_angle = -1*a_angle
                     b_angle = -1*b_angle
#                else:
#                    print("{} Rotating the right way  {}".format(pdb_code, signed_diff))

                # tm_keys_str = [str(i) for i in tm_keys]
                a_angle = dict(zip(tm_keys_str, a_angle))
                b_angle = dict(zip(tm_keys_str, b_angle))
                core_distances = dict(zip(tm_keys_str, core_distances))
                midpoint_distances = dict(zip(tm_keys_str, midpoint_distances))
                mid_membrane_distances = dict(zip(tm_keys_str, mid_membrane_distances))

                # Correct for missing values
                for res in polychain:
                    residue_id = str(res.id[1])

                    if not residue_id in a_angle:
                        a_angle[residue_id] = None
                    if not residue_id in b_angle:
                        b_angle[residue_id] = None
                    if not residue_id in core_distances:
                        core_distances[residue_id] = None
                    if not residue_id in midpoint_distances:
                        midpoint_distances[residue_id] = None
                    if not residue_id in mid_membrane_distances:
                        mid_membrane_distances[residue_id] = None
                    if not residue_id in rsa_list:
                        rsa_list[residue_id] = None
                    if not residue_id in hselist:
                        hselist[residue_id] = None
                    if not residue_id in dihedrals:
                        dihedrals[residue_id] = None
                    if not residue_id in asa_list:
                        asa_list[residue_id] = None
                    if not residue_id in rotation_angles:
                        rotation_angles[residue_id] = None

                #for res, angle1, angle2, distance, midpoint_distance, mid_membrane_distance in zip(pchain, a_angle, b_angle, core_distances, midpoint_distances, mid_membrane_distances):
                for res in polychain:
                    residue_id = str(res.id[1])

                    # structure, residue, A-angle, B-angle, RSA, HSE, "PHI", "PSI", "THETA", "TAU", "SS_DSSP", "SS_STRIDE", "OUTER", "TAU_ANGLE", "CHI", "MISSING", "ASA", "DISTANCE", "ROTATION_ANGLE"
                    if residue_id in full_resdict:
                        dblist.append([reference, full_resdict[residue_id], a_angle[residue_id], b_angle[residue_id], \
                            rsa_list[residue_id], \
                            hselist[residue_id]] + \
                            dihedrals[residue_id] + \
                            [asa_list[residue_id], core_distances[residue_id], midpoint_distances[residue_id], mid_membrane_distances[residue_id], rotation_angles[residue_id]])
            except Exception as e:
                print(pdb_code, " - ERROR - ", e)
                failed.append(pdb_code)
                # DEBUGGING
                if True:
                    traceback.print_exc()
#                    exit(0)
                continue

#        for i in range(4):
#            for key in angle_dict[i]:
#                sortlist = np.array(angle_dict[i][key])
#                median_dict[i][key] = np.median(sortlist)

#        for i, res in enumerate(dblist):
#            g = res[0]
#            a = res[1]
#
#            templist = copy.copy(angle_dict[res[4]][g.generic_number.label])
#            del templist[templist.index(a)]

#            std_test = abs(np.average(templist) - int(a))/np.std(templist)
#            std_len  = len(templist) - 1
#            std = stats.t.cdf(std_test, df=std_len)
#            dblist[i].append(0.501 if np.isnan(std) else std)

        # structure, residue, A-angle, B-angle, RSA, HSE, "PHI", "PSI", "THETA", "TAU", "SS_DSSP", "SS_STRIDE", "OUTER", "TAU_ANGLE", "CHI", "MISSING", "ASA", "DISTANCE", "ROTATION_ANGLE"
        object_list = []
        for ref,res,a1,a2,rsa,hse,phi,psi,theta,tau,ss_dssp,ss_stride,outer,tau_angle,chi_angles,missing,asa,distance,midpoint_distance,mid_membrane_distance,rotation_angle in dblist:
            try:
                if asa != None:
                    asa = round(asa,1)
                if outer != None:
                    outer = round(np.rad2deg(outer),3)
                if phi != None:
                    phi = round(np.rad2deg(phi),3)
                if psi != None:
                    psi = round(np.rad2deg(psi),3)
                if rsa != None:
                    rsa = round(rsa,1)
                if theta != None:
                    theta = round(np.rad2deg(theta),3)
                if tau != None:
                    tau = round(np.rad2deg(tau),3)
                if tau_angle != None:
                    tau_angle = round(np.rad2deg(tau_angle),3)
                object_list.append(Angle(residue=res, a_angle=a1, b_angle=a2, structure=ref, sasa=asa, rsa=rsa, hse=hse, phi=phi, psi=psi, theta=theta, tau=tau, tau_angle=tau_angle, chi1=chi_angles[0], chi2=chi_angles[1], chi3=chi_angles[2], chi4=chi_angles[3], chi5=chi_angles[4], missing_atoms=missing, ss_dssp=ss_dssp, ss_stride=ss_stride, outer_angle=outer, core_distance=distance, mid_distance=midpoint_distance, midplane_distance=mid_membrane_distance, rotation_angle=rotation_angle))
            except Exception as e:
                print(e)
                print([ref,res,a1,a2,rsa,hse,phi,psi,theta,tau,ss_dssp,ss_stride,outer,tau_angle,asa,distance,midpoint_distance,mid_membrane_distance])

        print("created list")
        print(len(object_list))

        # Store the results
        # faster than updating: deleting and recreating
        Angle.objects.bulk_create(object_list,batch_size=5000)
