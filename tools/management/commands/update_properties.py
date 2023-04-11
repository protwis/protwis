from django.core.management.base import BaseCommand

import copy

import contactnetwork.pdb as pdb

from structure.models import Structure
from residue.models import Residue
from angles.models import Angle
import logging

import numpy as np
from sklearn.decomposition import PCA
from numpy.core.umath_tests import inner1d
import io
import freesasa
import scipy.stats as stats

def load_pdb_var(pdb_code, var):
    """
    load string of pdb as pdb with a file handle. Would be nicer to do this
    directly, but no such function implemented in Bio PDB
    """
    parser = pdb.PDBParser(QUIET=True)
    with io.StringIO(var) as f:
        return parser.get_structure(pdb_code,f)

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



class Command(BaseCommand):

    help = "Command to calculate an axis for a TM helix."

    np.set_printoptions(suppress=True)
    logger = logging.getLogger(__name__)

    def handle(self, *args, **options):
        ###########################################################################
        ############################ Helper  Functions ############################
        ###########################################################################

        dblist = []
        extra_pca = True

        ###########################################################################
        ############################ Helper  Functions ############################
        ###########################################################################

        failed = []

        # get preferred chain for PDB-code
        references = Structure.objects.filter(protein_conformation__protein__family__slug__startswith="001").exclude(structure_type__slug__startswith='af-').prefetch_related('pdb_code','pdb_data','protein_conformation__protein','protein_conformation__state').order_by('protein_conformation__protein')
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


            try:

                print(pdb_code)

                structure = load_pdb_var(pdb_code,reference.pdb_data.pdb)
                pchain = structure[0][preferred_chain]
                state_id = reference.protein_conformation.state.id

                #######################################################################
                ###################### prepare and evaluate query #####################

                db_reslist = res_dict[pdb_code]

                #######################################################################
                ######################### filter data from db #########################

                def reslist_gen(x, rlist):
                    try:
                        while rlist[-1].generic_number.label[0] == x:
                            yield rlist.pop()
                    except IndexError:
                        pass

                # when gdict is not needed the helper can be removed
                #db_tmlist = [[(' ',r.sequence_number,' ') for r in reslist_gen(x) if r.sequence_number in pchain and r.sequence_number < 1000] for x in ["1","2","3","4","5","6","7"]]
                db_helper = [[(r,r.sequence_number) for r in reslist_gen(x,db_reslist) if r.sequence_number in pchain and r.sequence_number < 1000] for x in ["1","2","3","4","5","6","7"]]
                gdict = {r[1]:r[0] for hlist in db_helper for r in hlist}
                db_tmlist = [[(' ',r[1],' ') for r in sl] for sl in db_helper]
                db_set = set(db_tmlist[0]+db_tmlist[1]+db_tmlist[2]+db_tmlist[3]+db_tmlist[4]+db_tmlist[5]+db_tmlist[6])

                #######################################################################
                ############################# filter  pdb #############################

                recurse(structure, [[0], preferred_chain, db_set])

                #######################################################################
                ############### Calculate the axes through the helices ################
                #######################################################################
                N = 3

                hres_list = [np.asarray([pchain[r]["CA"].get_coord() for r in sl], dtype=float) for sl in db_tmlist]
                h_cb_list = [np.asarray([pchain[r]["CB"].get_coord() if "CB" in pchain[r] else cal_pseudo_CB(pchain[r]) for r in sl], dtype=float) for sl in db_tmlist]

                # fast and fancy way to take the average of N consecutive elements
                hres_three = np.asarray([sum([h[i:-(len(h) % N) or None:N] for i in range(N)])/N for h in hres_list])

                #######################################################################
                ################################# PCA #################################
                #######################################################################

                helix_pcas = [PCA() for i in range(7)]
                helix_pca_vectors = [pca_line(helix_pcas[i], h,i%2) for i,h in enumerate(hres_three)]

                # extracellular part
                if extra_pca:
                    helices_mn = np.asarray([np.mean(h, axis=0) for h in hres_three])
                    pos_list = np.asarray([pca_line(PCA(), h[:len(h)//2:(-(i%2) or 1)]) for i,h in enumerate(hres_three)])
                    pos_list = pos_list - (np.mean(pos_list,axis=1)-helices_mn).reshape(-1,1,3)

                    pca = PCA()
                    center_vector = pca_line(pca, np.vstack(pos_list))
                else:
                    pca = PCA()
                    center_vector = pca_line(pca, np.vstack(hres_three))

                #######################################################################
                ################################ Angles ###############################
                #######################################################################

                ########################### Axis to CA to CB ##########################

                b_angle = np.concatenate([ca_cb_calc(ca,cb,pca) for ca,cb in zip(hres_list,h_cb_list)]).round(3)

                #set_bfactor(pchain,angle)

                ######################### Axis to Axis to CA ##########################

                a_angle = np.concatenate([axes_calc(h,p,pca) for h,p in zip(hres_list,helix_pcas)]).round(3)

                #set_bfactor(pchain,angle2)

                ############## SASA

                res, trash = freesasa.calcBioPDB(structure)

                asa_list = []
                oldnum   = -1
                atomlist = list(pchain.get_atoms())
                for i in range(res.nAtoms()):
                    resnum = atomlist[i].get_parent().id[1]
                    if resnum == oldnum:
                        asa_list[-1] += res.atomArea(i)
                    else:
                        asa_list.append(res.atomArea(i))
                        oldnum = resnum

                ################ HSE


                hse = pdb.HSExposure.HSExposureCB(structure[0])
                hselist = [x[1][1] if x[1][1] > 0 else 0 for x in hse ]


                ################ dblist gen
                if len(pchain) != len(hselist):
                    raise Exception("\033[91mLength mismatch hse " + pdb_code + "\033[0m")

                if len(pchain) != len(asa_list):
                    raise Exception("\033[91mLength mismatch sasa " + pdb_code + "\033[0m")

                if np.isnan(np.sum(asa_list)):
                    raise Exception("\033[91mNAN sasa " + pdb_code + "\033[0m")

                if np.isnan(np.sum(hselist)):
                    raise Exception("\033[91mNAN hse " + pdb_code + "\033[0m")


                for res,a1,a2,asa,hse in zip(pchain,a_angle,b_angle,asa_list,hselist):
                    dblist.append([gdict[res.id[1]], a1, a2, reference,state_id-1,asa,hse])
                    if gdict[res.id[1]].generic_number.label not in angle_dict[state_id-1]:
                        angle_dict[state_id-1][gdict[res.id[1]].generic_number.label] = [round(a1,3)]
                    else:
                        angle_dict[state_id-1][gdict[res.id[1]].generic_number.label].append(round(a1,3))


            except Exception as e:
                print("ERROR!!", pdb_code, e)
                failed.append(pdb_code)
                continue

        for i in range(4):
            for key in angle_dict[i]:
                sortlist = np.array(angle_dict[i][key])
                median_dict[i][key] = np.median(sortlist)

        for i, res in enumerate(dblist):

            g = res[0]
            a = res[1]

            templist = copy.copy(angle_dict[res[4]][g.generic_number.label])
            del templist[templist.index(a)]

            std_test = abs(np.average(templist) - int(a))/np.std(templist)
            std_len  = len(templist) - 1
            std = stats.t.cdf(std_test, df=std_len)
            dblist[i].append(0.501 if np.isnan(std) else std)



        dblist = [Angle(residue=g, diff_med=round(abs(median_dict[i][g.generic_number.label]-a1),3), angle=a1, b_angle=a2, structure=ref, sasa=round(asa,3), hse=hse, sign_med=round(sig,3)) for g,a1,a2,ref,i,asa,hse,sig in dblist]

        print("created list")
        print(len(dblist))

        # faster than updating: deleting and recreating
        Angle.objects.all().delete()
        Angle.objects.bulk_create(dblist,batch_size=5000)
