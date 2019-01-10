from django.core.management.base import BaseCommand#, CommandError
#from django.core.management import call_command
#from django.conf import settings
#from django.db import connection

import freesasa
import contactnetwork.pdb as pdb

from structure.models import Structure
from residue.models import Residue

import logging, os
import numpy as np
from sklearn.decomposition import PCA
from numpy.core.umath_tests import inner1d
import io

TMNUM = 7
SASA  = True
HSE   = True

class Command(BaseCommand):

    help = "Command to calculate an axis for a TM helix."
    
    np.set_printoptions(suppress=True)
    logger = logging.getLogger(__name__)

    def add_arguments(self, parser):
        parser.add_argument('pdb_code')
    
    ###########################################################################
    ############################ Helper  Functions ############################
    ###########################################################################

    def load_pdb_var(self,pdb_code, var):
        """
        load string of pdb as pdb with a file handle. Would be nicer to do this
        directly, but no such function implemented in Bio PDB
        """
        parser = pdb.PDBParser()
        with io.StringIO(var) as f:
            return parser.get_structure(pdb_code,f)


    def save_pseudo(self, chainlist, pname):

        pseudopdb = pdb.Structure.Structure(pname)
        pseudopdb.add(pdb.Model.Model(0))
        
        for hi, h in enumerate(chainlist):
            pseudopdb[0].add((pdb.Chain.Chain(str(hi+1))))
            for j, r in enumerate(h):
                res = pdb.Residue.Residue((" ",j," "),"X",j)
                res.add(pdb.Atom.Atom("CA",r,0,0,"X","PSO",0,"U"))
                pseudopdb[0][str(hi+1)].add(res)
        
        io1 = pdb.PDBIO()
        io1.set_structure(pseudopdb)
        io1.save(pname + '.pdb')
        
    def write_cgo_arrow_pml(self, pdb_code, name, pos_list):
        """
        function to write a pymol script to automatically create cgo arrows for
        a list of positions
        """
        with open(pdb_code + name + ".pml", "w") as ps:
            ps.write("run cgo_arrow.py\n")
            for i, p in enumerate(pos_list):
                ps.write("cgo_arrow  " + str(list(p[0])) +", "+ str(list(p[1])) + ", name="+pdb_code + name + str(i) +"\n")

    def save_pdb(self,strct, name):
        """
        save a pdb structure as file
        """
        io1 = pdb.PDBIO()
        io1.set_structure(strct)
        io1.save(name)

    def handle(self, *args, **options):
        # grab PDB
        pdb_code = options.get('pdb_code', None).upper()
        
        reference = Structure.objects.get(pdb_code__index=pdb_code)
        # get preferred chain for PDB structure
        preferred_chain = reference.preferred_chain.split(',')[0]
        # read pdb structure (from RCSB) using Biopython
        structure = self.load_pdb_var(pdb_code,reference.pdb_data.pdb)

        # grab residues with the generic numbering for this structure
        db_reslist = list(Residue.objects.filter(generic_number__label__regex=r'^[1-7]x[0-9]+').filter(protein_conformation__protein=reference.protein_conformation.protein).prefetch_related('generic_number').order_by('-generic_number__label'))
        
        #######################################################################
        ############################# filter  pdb #############################
        
        os.chdir("pymol_output")
        pchain = structure[0][preferred_chain]
        
        def reslist_gen(x):
            try:
                while db_reslist[-1].generic_number.label[0] == x:
                    yield db_reslist.pop()
            except IndexError:
                pass
        
        db_tmlist = [[(' ',r.sequence_number,' ') for r in reslist_gen(x) if r.sequence_number in pchain and r.sequence_number < 1000] for x in ["1","2","3","4","5","6","7"]]
        db_set = set(db_tmlist[0]+db_tmlist[1]+db_tmlist[2]+db_tmlist[3]+db_tmlist[4]+db_tmlist[5]+db_tmlist[6])
        
        def recurse(entity,slist):
            for subenty in entity.get_list():
                if not subenty.id in slist[0]: entity.detach_child(subenty.id)
                elif slist[1:]: recurse(subenty, slist[1:])
        
        recurse(structure, [[0], preferred_chain, db_set])
        
        #######################################################################
        ############### Calculate the axes through the helices ################
        #######################################################################
        N = 3
        
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
        
        hres_list = [np.asarray([pchain[r]["CA"].get_coord() for r in sl], dtype=float) for sl in db_tmlist]
        h_cb_list = [np.asarray([pchain[r]["CB"].get_coord() if "CB" in pchain[r] else cal_pseudo_CB(pchain[r]) for r in sl], dtype=float) for sl in db_tmlist]
        
        # fast and fancy way to take the average of N consecutive elements
        hres_three = np.asarray([sum([h[i:-(len(h) % N) or None:N] for i in range(N)])/N for h in hres_list])
        helices_mn = np.asarray([np.mean(h, axis=0) for h in hres_three ])
        self.save_pseudo(hres_three, pdb_code+"helper")
        
        #######################################################################
        ################################# PCA #################################
        #######################################################################
        
        def pca_line(pca,h, r=0):
            if ((not r) if pca.fit_transform(h)[0][0] < 0 else r):
                return pca.inverse_transform(np.asarray([[-20,0,0],[20,0,0]]))
            else:return pca.inverse_transform(np.asarray([[20,0,0],[-20,0,0]]))  
        
        helix_pcas = [PCA() for i in range(7)]
        pos_list = np.asarray([pca_line(helix_pcas[i], h,i%2) for i,h in enumerate(hres_three)])
        self.write_cgo_arrow_pml(pdb_code, "pca",pos_list)
        
        pos_list = np.mean(pos_list,axis=0)
        self.write_cgo_arrow_pml(pdb_code, "pca_mean",[pos_list])
        
        pca = PCA()
        pos_list = pca_line(pca, np.vstack(hres_three))
        self.write_cgo_arrow_pml(pdb_code, "pca_all",[pos_list])
        
        pos_list = np.asarray([pca_line(PCA(), h[:len(h)//2:(-(i%2) or 1)]) for i,h in enumerate(hres_three)])
        pos_list = pos_list - (np.mean(pos_list,axis=1)-helices_mn).reshape(-1,1,3)
        self.write_cgo_arrow_pml(pdb_code, "pca_extra",pos_list)
        self.write_cgo_arrow_pml(pdb_code, "pca_extra_mean",[np.mean(pos_list,axis=0)])
        
        pca_extra = PCA()
        pos_list = pca_line(pca_extra, np.vstack(pos_list))
        self.write_cgo_arrow_pml(pdb_code, "pca_extra_pca",[pos_list])
        
        #######################################################################
        ################################ Angles ###############################
        #######################################################################
        
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
        
        centerpca = pca
        
        ########################### Axis to CA to CB ##########################

        angle = np.concatenate([ca_cb_calc(ca,cb,centerpca) for ca,cb in zip(hres_list,h_cb_list)])
        set_bfactor(pchain,angle)
        self.save_pdb(structure, pdb_code+'angle_colored_ca_cb.pdb')
        
        ######################### Axis to Axis to CA ##########################
        
        angle2 = np.concatenate([axes_calc(h,p,centerpca) for h,p in zip(hres_list,helix_pcas)])
        set_bfactor(pchain,angle2)

        self.save_pdb(structure, pdb_code+'angle_colored_axes.pdb')
        
        #######################################################################
        ############################ HSE and SASA #############################
        #######################################################################
        
        ################################ SASA #################################
        
        if SASA:
            pdbstruct = freesasa.Structure(pdb_code+'angle_colored_axes.pdb')
            res = freesasa.calc(pdbstruct)
            asa_list = []
            oldnum = -1
            for i in range(res.nAtoms()):
                resnum = pdbstruct.residueNumber(i)
                if resnum == oldnum:
                    asa_list[-1] += res.atomArea(i)
                else:
                    asa_list.append(res.atomArea(i))
                    oldnum = resnum
            
            set_bfactor(pchain,asa_list)
            self.save_pdb(structure, pdb_code+'asa_colored.pdb')
        
        ################################# HSE #################################
        if HSE:
            exp_cb = pdb.HSExposure.HSExposureCB(structure[0])
            
            [[a.set_bfactor(x[1][1]) for a in x[0]] for x in exp_cb]
            
            exp_cb = [a["CA"].get_bfactor() for a in pchain.get_list()]
            
            self.save_pdb(structure, pdb_code+'hsea_colored.pdb')

        
        
        
        
        
        
        
        
        
        
        
        
        