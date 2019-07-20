from django.core.management.base import BaseCommand
from django.db import connection

import contactnetwork.pdb as pdb
from structure.models import Structure, StructureVectors
from residue.models import Residue
from angles.models import ResidueAngle as Angle

import Bio.PDB
import copy
import freesasa
import io
import logging
import math
import subprocess
import os

import numpy as np
import scipy.stats as stats

from collections import OrderedDict
from sklearn.decomposition import PCA
from numpy.core.umath_tests import inner1d


from multiprocessing import Queue, Process, Value, Lock

TMNUM = 7

SASA = True
HSE  = True
extra_pca = True
print_pdb = False

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


class NonHetSelect(Bio.PDB.Select):
    def accept_residue(self, residue):
        return 1 if residue.id[0] == " " else 0

class Command(BaseCommand):

    help = "Command to calculate all angles for residues in each TM helix."

    np.set_printoptions(suppress=True)
    logger = logging.getLogger(__name__)

    ###########################################################################
    ############################ Helper  Functions ############################
    ###########################################################################

    processes = 4

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

    def load_pdb_var(self, pdb_code, var):
        """
        load string of pdb as pdb with a file handle. Would be nicer to do this
        directly, but no such function implemented in Bio PDB
        """
        parser = pdb.PDBParser(QUIET=True)
        with io.StringIO(var) as f:
            return parser.get_structure(pdb_code,f)


    def handle(self, *args, **options):

        Angle.objects.all().delete()
        self.references = Structure.objects.all().exclude(refined=True).prefetch_related('pdb_code','pdb_data','protein_conformation__protein','protein_conformation__state').order_by('protein_conformation__protein')
        print(len(self.references),'structures')
        self.references = list(self.references)
        self.prepare_input(self.processes, self.references)

    def main_func(self, positions, iteration,count,lock):
        print('main_func')
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

        def axes_calc(h,p,pca):
            """
            Calculate the orthogonal projection of the CA to the helix axis
            which is moved to the mean of three consecutive amino acids
            """
            a = (np.roll(np.vstack((h,h[0])),1,axis=0)[:-1] + h + np.roll(np.vstack((h,h[-1])),-1,axis=0)[:-1])/3
            b = p.transform(h)
            b[:,1:] = p.transform(a)[:,1:]
            b = p.inverse_transform(b)
            return calc_angle(pca.transform(b),pca.transform(h))

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

        #######################################################################
        ######################### Start of main loop ##########################
        #######################################################################

        failed = []
        dblist = []

        failed = []

        # Get all structures
        #references = Structure.objects.filter(protein_conformation__protein__family__slug__startswith="001").exclude(refined=True).prefetch_related('pdb_code','pdb_data','protein_conformation__protein','protein_conformation__state').order_by('protein_conformation__protein')
        #references = Structure.objects.all().exclude(refined=True).prefetch_related('pdb_code','pdb_data','protein_conformation__protein','protein_conformation__state').order_by('protein_conformation__protein')
        # DEBUG for a specific PDB
        #references = Structure.objects.filter(pdb_code__index="6AK3").exclude(refined=True).prefetch_related('pdb_code','pdb_data','protein_conformation__protein','protein_conformation__state').order_by('protein_conformation__protein')

        # references = list(references)
        references = self.references

        pids = [ref.protein_conformation.protein.id for ref in references]

        qset = Residue.objects.filter(protein_conformation__protein__id__in=pids)
        qset = qset.filter(generic_number__label__regex=r'^[1-7]x[0-9]+').order_by('-protein_conformation__protein','-generic_number__label')
        qset = list(qset.prefetch_related('generic_number', 'protein_conformation__protein','protein_conformation__state'))

        res_dict = {ref.pdb_code.index:qgen(ref.protein_conformation.protein,qset) for ref in references}

        # clean structure vectors table
        StructureVectors.objects.all().delete()

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
            print(pdb_code)

            try:
#            if True:
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
                if os.path.exists("/env/bin/mkdssp"):
                    dssp = Bio.PDB.DSSP(structure[0], filename, dssp='/env/bin/mkdssp')

                # STRIDE
                try:
                    if os.path.exists("/env/bin/stride"):
                       stride = subprocess.Popen(['/env/bin/stride', filename], stdout=subprocess.PIPE)
                       # Grab SS assignment (ASG) and parse residue (cols 12-15) and SS (cols 25-25)
                       for line in io.TextIOWrapper(stride.stdout, encoding="utf-8"):
                           if line.startswith("ASG"):
                               res_id = int(line[11:15].strip())
                               res_ss = line[24:25].strip()
                               # assign to residue
                               pchain[res_id].xtra["SS_STRIDE"] = res_ss.upper()
                except OSError:
                   print(pdb_code, " - STRIDE ERROR - ", e)

                # CLEANUP
                os.remove(filename)

                #######################################################################
                ###################### prepare and evaluate query #####################

                db_reslist = res_dict[pdb_code]
                #print(db_reslist)

                #######################################################################
                ######################### filter data from db #########################

                def reslist_gen(x):
                    try:
                        while db_reslist[-1].generic_number.label[0] == x:
                            yield db_reslist.pop()
                    except IndexError:
                        pass

                # when gdict is not needed the helper can be removed
                db_helper = [[(r,r.sequence_number) for r in reslist_gen(x) if r.sequence_number in pchain] for x in ["1","2","3","4","5","6","7"]]
                gdict = {r[1]:r[0] for hlist in db_helper for r in hlist}
                db_tmlist = [[(' ',r[1],' ') for r in sl] for sl in db_helper]
                db_set = set(db_tmlist[0]+db_tmlist[1]+db_tmlist[2]+db_tmlist[3]+db_tmlist[4]+db_tmlist[5]+db_tmlist[6])

                #######################################################################
                ##################### Angles/dihedrals residues #######################

                polychain = [ residue for residue in pchain if Bio.PDB.Polypeptide.is_aa(residue) and "CA" in residue]
                poly = Bio.PDB.Polypeptide.Polypeptide(polychain)
                poly.get_phi_psi_list() # backbone dihedrals

                #poly.get_theta_list() # angle three consecutive Ca atoms
                #poly.get_tau_list() # dihedral four consecutive Ca atoms

                # TODO: extend with Chi1-5?
                # https://gist.github.com/lennax/0f5f65ddbfa278713f58
                # Definition http://www.ccp14.ac.uk/ccp/web-mirrors/garlic/garlic/commands/dihedrals.html
                # http://biopython.org/DIST/docs/api/Bio.PDB.Polypeptide-pysrc.html#Polypeptide.get_phi_psi_list

                ### clean the structure to solely the 7TM bundle
                recurse(structure, [[0], preferred_chain, db_set])
                poly.get_theta_list() # angle three consecutive Ca atoms
                poly.get_tau_list() # dihedral four consecutive Ca atoms
                dihedrals = {}
                for r in poly:
                  angle_list = ["PHI", "PSI", "THETA", "TAU", "SS_DSSP", "SS_STRIDE"]
                  for angle in angle_list:
                      if angle not in r.xtra:
                          r.xtra[angle] = None

                  # Add outer angle
                  outer = None
                  try:
                      angle_atoms = [r[a].get_vector() for a in ['N','CA', outerAtom[r.resname]]]

                      # use pseudo CB placement when glycine
                      if r.resname == 'GLY':
                          angle_atoms[2] = Bio.PDB.vectors.Vector(*cal_pseudo_CB(r))

                      outer = Bio.PDB.calc_angle(*angle_atoms)
                  except Exception as e:
#                      print(pdb_code, " - ANGLE ERROR - ", e)
                      outer = None

                  dihedrals[r.id[1]] = [r.xtra["PHI"], r.xtra["PSI"], r.xtra["THETA"], r.xtra["TAU"], r.xtra["SS_DSSP"], r.xtra["SS_STRIDE"], outer]

                # Extra: remove hydrogens from structure (e.g. 5VRA)
                for residue in structure[0][preferred_chain]:
                    for id in [atom.id for atom in residue if atom.element == "H"]:
                        residue.detach_child(id)

                ### AXES through each of the TMs and the TM bundle (center axis)
                hres_list = [np.asarray([pchain[r]["CA"].get_coord() for r in sl], dtype=float) for sl in db_tmlist]
                h_cb_list = [np.asarray([pchain[r]["CB"].get_coord() if "CB" in pchain[r] else cal_pseudo_CB(pchain[r]) for r in sl], dtype=float) for sl in db_tmlist]

                # fast and fancy way to take the average of N consecutive elements
                N = 3
                hres_three = np.asarray([sum([h[i:-(len(h) % N) or None:N] for i in range(N)])/N for h in hres_list])

                ### PCA - determine axis through center + each transmembrane helix
                helix_pcas = [PCA() for i in range(7)]
                helix_pca_vectors = [pca_line(helix_pcas[i], h,i%2) for i,h in enumerate(hres_three)]

                # Calculate PCA based on the upper (extracellular) half of the GPCR (more stable, except class B)
                pca = PCA()
                if extra_pca:
                    minlength = 100
                    for i,h in enumerate(hres_three):
                        if len(h)<minlength:
                            minlength = len(h)

                    if minlength > 6:
                        minlength = 6

                    # create PCA per helix using extracellular half
                    # Exlude the first turn if possible (often still part of loop)
                    pos_list = []
                    for i,h in enumerate(hres_three):
                        if i%2: # reverse directionality of even helices (TM2, TM4, TM6)
                            h = np.flip(h, 0)

                        if len(h)>minlength+2:
                            pos_list.append(pca_line(PCA(), h[2:minlength+2]))
                        else:
                            pos_list.append(pca_line(PCA(), h[0:minlength]))


                    # create fake coordinates along each helix PCA to create center PCA
                    # UGLY hack - should be cleand up
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
                    # Less robust with differing TM lengths
                    center_vector = pca_line( pca, np.vstack(hres_three))

                # DEBUG print arrow for PyMol
                # a = [str(i) for i in center_vector[0]]
                # b = [str(i) for i in center_vector[1]]
                # print("cgo_arrow [" + a[0] + ", " + a[1] + ", " + a[2] + "], [" + b[0] + ", " + b[1] + ", " + b[2] + "]")

                ### ANGLES
                # Center axis to helix axis to CA
                a_angle = np.concatenate([axes_calc(h,p,pca) for h,p in zip(hres_list,helix_pcas)]).round(3)

                # Center axis to CA to CB
                b_angle = np.concatenate([ca_cb_calc(ca,cb,pca) for ca,cb in zip(hres_list,h_cb_list)]).round(3)

                # Distance from center axis to CA
                core_distance = np.concatenate([ca_distance_calc(ca,pca) for ca in hres_list]).round(3)

                # STORE STRUCTURE REFERENCES
                # center axis
                c_vector = np.array2string(center_vector[0] - center_vector[1], separator=',')
                translation = np.array2string(-1*center_vector[0], separator=',')

                sv = StructureVectors(structure = reference, translation = str(translation), center_axis = str(c_vector))
                sv.save()

                ### freeSASA (only for TM bundle)
                # SASA calculations - results per atom
                res, trash = freesasa.calcBioPDB(structure)

                # create results dictionary per residue
                asa_list = {}
                rsa_list = {}
                atomlist = list(pchain.get_atoms())
                for i in range(res.nAtoms()):
                    resnum = atomlist[i].get_parent().id[1]
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
                hselist = dict([ (x[0].id[1], x[1][0]) if x[1][0] > 0 else (x[0].id[1], 0) for x in hse ])

                # Few checks
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

                # Correct for missing values
                for res in pchain:
                    residue_id = res.id[1]
                    if not residue_id in rsa_list:
                        rsa_list[residue_id] = None
                    if not residue_id in hselist:
                        hselist[residue_id] = None
                    if not residue_id in dihedrals:
                        dihedrals[residue_id] = None
                    if not residue_id in asa_list:
                        asa_list[residue_id] = None


                for res, angle1, angle2, distance in zip(pchain, a_angle, b_angle, core_distance):
                    residue_id = res.id[1]
                    # structure, residue, A-angle, B-angle, RSA, HSE, "PHI", "PSI", "THETA", "TAU", "SS_DSSP", "SS_STRIDE", "OUTER", "ASA", "DISTANCE"
                    dblist.append([reference, gdict[residue_id], angle1, angle2, \
                        rsa_list[residue_id], \
                        hselist[residue_id]] + \
                        dihedrals[residue_id] + \
                        [asa_list[residue_id], distance])

            except Exception as e:
#            else:
                print(pdb_code, " - ERROR - ", e)
                failed.append(pdb_code)
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

        # structure, residue, A-angle, B-angle, RSA, HSE, "PHI", "PSI", "THETA", "TAU", "SS_DSSP", "SS_STRIDE", "OUTER", "ASA", "DISTANCE"
        object_list = []
        for ref,res,a1,a2,rsa,hse,phi,psi,theta,tau,ss_dssp,ss_stride,outer,asa,distance in dblist:
            try:
                if phi != None:
                    phi = round(np.rad2deg(phi),3)
                if psi != None:
                    psi = round(np.rad2deg(psi),3)
                if theta != None:
                    theta = round(np.rad2deg(theta),3)
                if tau != None:
                    tau = round(np.rad2deg(tau),3)
                if outer != None:
                    outer = round(np.rad2deg(outer),3)
                object_list.append(Angle(residue=res, a_angle=a1, b_angle=a2, structure=ref, sasa=round(asa,1), rsa=round(rsa,1), hse=hse, phi=phi, psi=psi, theta=theta, tau=tau, ss_dssp=ss_dssp, ss_stride=ss_stride, outer_angle=outer, core_distance=distance))
            except Exception as e:
                print([ref,res,a1,a2,rsa,hse,phi,psi,theta,tau,ss_dssp,ss_stride,outer,asa,distance])

        print("created list")
        print(len(object_list))

        # Store the results
        # faster than updating: deleting and recreating
        Angle.objects.bulk_create(object_list,batch_size=5000)
