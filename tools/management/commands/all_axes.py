from django.core.management.base import BaseCommand#, CommandError
#from django.core.management import call_command
#from django.conf import settings
#from django.db import connection

import contactnetwork.pdb as pdb

from structure.models import Structure
from residue.models import Residue
#from angles.model import Angle

import logging
import numpy as np
from sklearn.decomposition import PCA
from numpy.core.umath_tests import inner1d
import io
import freesasa
import pickle
TMNUM = 7

SASA = False
HSE  = True
extra_pca = True
print_pdb = False

class Command(BaseCommand):

    help = "Command to calculate an axis for a TM helix."

    np.set_printoptions(suppress=True)
    logger = logging.getLogger(__name__)

    ###########################################################################
    ############################ Helper  Functions ############################
    ###########################################################################

    def load_pdb_var(self,pdb_code, var):
        """
        load string of pdb as pdb with a file handle. Would be nicer to do this
        directly, but no such function implemented in Bio PDB
        """
        parser = pdb.PDBParser(QUIET=True)
        with io.StringIO(var) as f:
            return parser.get_structure(pdb_code,f)

    def write_cgo_arrow_pml(self, pdb_code, name, pos_list):
        """
        function to write a pymol script to automatically create cgo arrows for
        a list of positions
        """
        with open("pymol_output/"+pdb_code + name + ".pml", "w") as ps:
            ps.write("run cgo_arrow.py\n")
            for i, p in enumerate(pos_list):
                ps.write("cgo_arrow  " + str(list(p[0])) +", "+ str(list(p[1])) + ", name="+pdb_code + name + str(i) +"\n")

    def save_pdb(self, strct, name):
        """
        save a pdb structure as file
        """
        io1 = pdb.PDBIO()
        io1.set_structure(strct)
        io1.save("pymol_output/"+name)

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
                return pca.inverse_transform(np.asarray([[-20,0,0],[20,0,0]]))
            else:return pca.inverse_transform(np.asarray([[20,0,0],[-20,0,0]]))

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

        def qgen(x):
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

        failed = []

        # get preferred chain for PDB-code
        references = Structure.objects.filter(protein_conformation__protein__family__slug__startswith="001").exclude(refined=True).prefetch_related('pdb_code','pdb_data','protein_conformation').order_by('protein_conformation__protein')
        references = list(references)

        pids = [ref.protein_conformation.protein.id for ref in references]

        qset = Residue.objects.filter(protein_conformation__protein__id__in=pids)
        qset = qset.filter(generic_number__label__regex=r'^[1-7]x[0-9]+').order_by('-protein_conformation__protein','-generic_number__label')
        qset = list(qset.prefetch_related('generic_number', 'protein_conformation'))

        res_dict = {ref.pdb_code.index:qgen(ref.protein_conformation.protein) for ref in references}

        #######################################################################
        ######################### Start of main loop ##########################
        #######################################################################

        for reference in references:

            preferred_chain = reference.preferred_chain.split(',')[0]
            pdb_code = reference.pdb_code.index
            state_id = reference.protein_conformation.state.id

            try:

                print(pdb_code)

                structure = self.load_pdb_var(pdb_code,reference.pdb_data.pdb)
                pchain = structure[0][preferred_chain]

                #######################################################################
                ###################### prepare and evaluate query #####################

                db_reslist = res_dict[pdb_code]

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
                db_helper = [[(r.generic_number.label,r.sequence_number) for r in reslist_gen(x) if r.sequence_number in pchain and r.sequence_number < 1000] for x in ["1","2","3","4","5","6","7"]]
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
                [pca_line(helix_pcas[i], h,i%2) for i,h in enumerate(hres_three)]

                # extracellular part
                if extra_pca:
                    helices_mn = np.asarray([np.mean(h, axis=0) for h in hres_three])
                    pos_list = np.asarray([pca_line(PCA(), h[:len(h)//2:(-(i%2) or 1)]) for i,h in enumerate(hres_three)])
                    pos_list = pos_list - (np.mean(pos_list,axis=1)-helices_mn).reshape(-1,1,3)

                    pca = PCA()
                    pca_line(pca, np.vstack(pos_list))
                else:
                    pca = PCA()
                    pca_line(pca, np.vstack(hres_three))

                #######################################################################
                ################################ Angles ###############################
                #######################################################################

                ########################### Axis to CA to CB ##########################

                angle = np.concatenate([ca_cb_calc(ca,cb,pca) for ca,cb in zip(hres_list,h_cb_list)])

                set_bfactor(pchain,angle)

                if print_pdb:
                    self.save_pdb(structure, pdb_code+'angle_colored_ca_cb.pdb')

                ######################### Axis to Axis to CA ##########################

                angle2 = np.concatenate([axes_calc(h,p,pca) for h,p in zip(hres_list,helix_pcas)])

                set_bfactor(pchain,angle2)
                if print_pdb:
                    self.save_pdb(structure, pdb_code+'angle_colored_axes.pdb')

                ################################ SASA #################################
                if SASA:
                    pdbstruct = freesasa.Structure("pymol_output/" + pdb_code+'angle_colored_axes.pdb')
                    res = freesasa.calc(pdbstruct)

                    asa_list = []
                    oldnum   = -1
                    for i in range(res.nAtoms()):
                        resnum = pdbstruct.residueNumber(i)
                        if resnum == oldnum:
                            asa_list[-1] += res.atomArea(i)
                        else:
                            asa_list.append(res.atomArea(i))
                            oldnum = resnum

                    set_bfactor(pchain,asa_list)
                    if print_pdb:
                        self.save_pdb(structure, pdb_code+'asa_colored.pdb')

                ################################# HSE #################################
                if HSE:
                    hse = pdb.HSExposure.HSExposureCB(structure[0])
                    [[a.set_bfactor(x[1][1]) for a in x[0]] for x in hse]

                    if print_pdb:
                        self.save_pdb(structure, pdb_code+'hsea_colored.pdb')

                ############################### pickle ################################
                if HSE and SASA:
                    reslist = []
                    grslist = []
                    hse     = []
                    for r in pchain:
                        reslist.append(r.id[1])
                        grslist.append(gdict[r.id[1]])
                        hse.append(r["CA"].get_bfactor())
                    with open('pymol_output/'+ pdb_code + '_measures.pickle', 'wb') as handle:
                        pickle.dump((np.array(reslist),grslist, np.array(asa_list),np.array(hse),angle,angle2, state_id), handle)

                #Angle.objects.bulk_create([Angle(residue=gdict[res.id[1]], angle=res["CA"].get_bfactor(), structure=reference) for res in pchain])

            except Exception as e:
                print("ERROR!!", pdb_code, e)
                failed.append(pdb_code)
                continue

        print(len(failed), "of", len(references), "failed:", failed)
