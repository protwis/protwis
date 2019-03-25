from django.core.management.base import BaseCommand

import contactnetwork.pdb as pdb
from structure.models import Structure
from residue.models import Residue
from angles.models import ResidueAngle as Angle

import Bio.PDB
import copy
import freesasa
import io
import logging

import numpy as np
import scipy.stats as stats

from collections import OrderedDict
from sklearn.decomposition import PCA
from numpy.core.umath_tests import inner1d

TMNUM = 7

SASA = True
HSE  = True
extra_pca = True
print_pdb = False

class Command(BaseCommand):

    help = "Command to calculate all angles for residues in each TM helix."

    np.set_printoptions(suppress=True)
    logger = logging.getLogger(__name__)

    ###########################################################################
    ############################ Helper  Functions ############################
    ###########################################################################

    def load_pdb_var(self, pdb_code, var):
        """
        load string of pdb as pdb with a file handle. Would be nicer to do this
        directly, but no such function implemented in Bio PDB
        """
        parser = pdb.PDBParser(QUIET=True)
        with io.StringIO(var) as f:
            return parser.get_structure(pdb_code,f)

    def handle(self, *args, **options):
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
            return np.degrees(np.arccos(inner1d(ba, bc) / (np.linalg.norm(ba,axis=1) * np.linalg.norm(bc,axis=1))))

        def ca_cb_calc(ca,cb,pca):
            """
            Calcuate the angles between ca, cb and center axis
            """
            return calc_angle(pca.transform(ca),pca.transform(cb))

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
        references = Structure.objects.all().exclude(refined=True).prefetch_related('pdb_code','pdb_data','protein_conformation__protein','protein_conformation__state').order_by('protein_conformation__protein')
        # DEBUG for a specific PDB
        #references = Structure.objects.filter(pdb_code__index="6H7O").filter(protein_conformation__protein__family__slug__startswith="001").exclude(refined=True).prefetch_related('pdb_code','pdb_data','protein_conformation__protein','protein_conformation__state').order_by('protein_conformation__protein')

        references = list(references)

        pids = [ref.protein_conformation.protein.id for ref in references]

        qset = Residue.objects.filter(protein_conformation__protein__id__in=pids)
        qset = qset.filter(generic_number__label__regex=r'^[1-7]x[0-9]+').order_by('-protein_conformation__protein','-generic_number__label')
        qset = list(qset.prefetch_related('generic_number', 'protein_conformation__protein','protein_conformation__state'))

        res_dict = {ref.pdb_code.index:qgen(ref.protein_conformation.protein,qset) for ref in references}

        #######################################################################
        ######################### Start of main loop ##########################
        #######################################################################

        angle_dict = [{},{},{},{}]
        median_dict = [{},{},{},{}]

        for reference in references:
            preferred_chain = reference.preferred_chain.split(',')[0]
            pdb_code = reference.pdb_code.index
            print(pdb_code)

            try:
            #if True:
                structure = self.load_pdb_var(pdb_code,reference.pdb_data.pdb)
                pchain = structure[0][preferred_chain]
                state_id = reference.protein_conformation.state.id

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
                #db_tmlist = [[(' ',r.sequence_number,' ') for r in reslist_gen(x) if r.sequence_number in pchain and r.sequence_number < 1000] for x in ["1","2","3","4","5","6","7"]]
                # db_helper = [[(r,r.sequence_number) for r in reslist_gen(x) if r.sequence_number in pchain and r.sequence_number < 1000] for x in ["1","2","3","4","5","6","7"]]
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
                  angle_list = ["PHI", "PSI", "THETA", "TAU"]
                  for angle in angle_list:
                      if angle not in r.xtra:
                          r.xtra[angle] = None
                  dihedrals[r.id[1]] = [r.xtra["PHI"], r.xtra["PSI"], r.xtra["THETA"], r.xtra["TAU"]]

                # Extra: remove hydrogens from structure (e.g. 5VRA)
                for residue in structure[0][preferred_chain]:
                    for id in [atom.id for atom in residue if atom.element == "H"]:
                        residue.detach_child(id)

                ### AXES through each of the TMs and the TM bundle (center axis)
                hres_list = [np.asarray([pchain[r]["CA"].get_coord() for r in sl], dtype=float) for sl in db_tmlist]
                h_cb_list = [np.asarray([pchain[r]["CB"].get_coord() if "CB" in pchain[r] else cal_pseudo_CB(pchain[r]) for r in sl], dtype=float) for sl in db_tmlist]
                #print(hres_list)
                # fast and fancy way to take the average of N consecutive elements
                N = 3
                hres_three = np.asarray([sum([h[i:-(len(h) % N) or None:N] for i in range(N)])/N for h in hres_list])
                #print(hres_three)
                ### PCA - determine axis through center + each transmembrane helix
                helix_pcas = [PCA() for i in range(7)]
                helix_pca_vectors = [pca_line(helix_pcas[i], h,i%2) for i,h in enumerate(hres_three)]

                # Calculate PCA based on the upper (extracellular) half of the GPCR (more stable)
                pca = PCA()
                # DISABLED FOR NOW - not stable at all
                if False and extra_pca:
                    helices_mn = np.asarray([np.mean(h, axis=0) for h in hres_three])

                    # TODO: investigate and make more robust
                    # FOR EXAMPLE: skip the first residue + focus on middle membrane and up residues (knowledge-based not structure-dependent)
                    #pos_list = np.asarray([pca_line(PCA(), h[:len(h)//2:(-(i%2) or 1)]) for i,h in enumerate(hres_three)])
                    pos_list = []
                    for i,h in enumerate(hres_three):
                        if len(h)>6:
                            pos_list.append(pca_line(PCA(), h[:len(h)//2:(-(i%2) or 1)]))
                        else:
                            pos_list.append(pca_line(PCA(), h))

                    pos_list = np.asarray(pos_list)
                    pos_list = pos_list - (np.mean(pos_list,axis=1)-helices_mn).reshape(-1,1,3)

                    # print cgo_arrows for each TM helix
#                    for pos in pos_list:
#                        a = [str(i) for i in pos[0]]
#                        b = [str(i) for i in pos[1]]
#                        print("cgo_arrow [" + a[0] + ", " + a[1] + ", " + a[2] + "], [" + b[0] + ", " + b[1] + ", " + b[2] + "]")

                    # TODO store center_vector + vector to reference residue
                    center_vector = pca_line(pca, np.vstack(pos2_list))

                    # average vector -> center vector ok, but PCA not initalized
#                    center_point = [0,0,0]
#                    vector_average = [0,0,0]
#                    for pos in pos_list:
#                        center_point += pos[0]/len(pos_list)
#                        vector_average += (pos[1]-pos[0])/len(pos_list)

#                    center_vector = [center_point, center_point + vector_average ]
                else:
                    center_vector = pca_line(pca, np.vstack(hres_three))

#                a = [str(i) for i in center_vector[0]]
#                b = [str(i) for i in center_vector[1]]
#                print("cgo_arrow [" + a[0] + ", " + a[1] + ", " + a[2] + "], [" + b[0] + ", " + b[1] + ", " + b[2] + "]")

                ### ANGLES
                # Center axis to helix axis to CA
                a_angle = np.concatenate([axes_calc(h,p,pca) for h,p in zip(hres_list,helix_pcas)]).round(3)

                # Center axis to CA to CB
                b_angle = np.concatenate([ca_cb_calc(ca,cb,pca) for ca,cb in zip(hres_list,h_cb_list)]).round(3)


                ### freeSASA (only for TM bundle)
                # SASA calculations - results per atom
                res, trash = freesasa.calcBioPDB(structure)

                # create results dictionary per residue
                asa_list = {}
                atomlist = list(pchain.get_atoms())
                for i in range(res.nAtoms()):
                    resnum = atomlist[i].get_parent().id[1]
                    if resnum not in asa_list:
                        asa_list[resnum] = 0
                    asa_list[resnum] += res.atomArea(i)

                ### Half-sphere exposure (HSE)
                hse = pdb.HSExposure.HSExposureCB(structure[0])
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
                for res, angle1, angle2 in zip(pchain, a_angle, b_angle):
                    residue_id = res.id[1]
                    # structure, residue, A-angle, B-angle, ASA, HSE, "PHI", "PSI", "THETA", "TAU"
                    dblist.append([reference, gdict[residue_id], angle1, angle2, asa_list[residue_id], hselist[residue_id]] + dihedrals[residue_id])

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

        # structure, residue, A-angle, B-angle, ASA, HSE, "PHI", "PSI", "THETA", "TAU"
        object_list = []
        for ref,res,a1,a2,asa,hse,phi,psi,theta,tau in dblist:
            try:
                if phi != None:
                    phi = round(np.rad2deg(phi),3)
                if psi != None:
                    psi = round(np.rad2deg(psi),3)
                if theta != None:
                    theta = round(np.rad2deg(theta),3)
                if tau != None:
                    tau = round(np.rad2deg(tau),3)
                object_list.append(Angle(residue=res, a_angle=a1, b_angle=a2, structure=ref, sasa=round(asa,1), hse=hse, phi=phi, psi=psi, theta=theta, tau=tau))
            except Exception as e:
                print([ref,res,a1,a2,asa,hse,phi,psi,theta,tau])
        #object_list = [Angle(residue=res, a_angle=a1, b_angle=a2, structure=ref, sasa=round(asa,1), hse=hse, phi=round(np.rad2deg(phi),3), psi=round(np.rad2deg(psi),3), theta=round(np.rad2deg(theta),3), tau=round(np.rad2deg(tau),3)) for ref,res,a1,a2,asa,hse,phi,psi,theta,tau in dblist]
        print("created list")
        print(len(object_list))

        # Store the results
        # faster than updating: deleting and recreating
        Angle.objects.all().delete()
        Angle.objects.bulk_create(object_list,batch_size=5000)
