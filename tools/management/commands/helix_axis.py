from django.core.management.base import BaseCommand#, CommandError
#from django.core.management import call_command
#from django.conf import settings
#from django.db import connection

import contactnetwork.pdb as pdb

from structure.models import Structure
from protein.models import Protein #, ProteinSegment
from residue.models import Residue

#from Bio.PDB import PDBParser, PPBuilder, parse_pdb_header 

import logging, os #, json, os
from time import time
import numpy as np
from sklearn import linear_model
from sklearn.decomposition import PCA

class Command(BaseCommand):

    help = "Command to calculate an axis for a TM helix."
    
    np.set_printoptions(suppress=True)
    logger = logging.getLogger(__name__)

    def add_arguments(self, parser):
        parser.add_argument('pdb_code')
        
        
    def load_pdb(self,pdb_code, no_save=False):
        if no_save:
            return pdb.pdb_get_structure(pdb_code)
        else:
            if not os.path.exists("pdbfiles/pdb" + pdb_code.lower() + ".ent"):
                pdbl = pdb.PDBList()
                pdbl.retrieve_pdb_file(pdb_code, file_format="pdb", pdir="pdbfiles")
                
            parser = pdb.PDBParser()
            return parser.get_structure(pdb_code, "pdbfiles/pdb" + pdb_code.lower() + ".ent")
      
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
        
    def write_cgo_arrow_pml(self, name, pos_list):
        with open(name + ".pml", "w") as ps:
            ps.write("run cgo_arrow.py\n")
            for i, p in enumerate(pos_list):
                ps.write("cgo_arrow  " + str(list(p[0])) +", "+ str(list(p[1])) + ", name="+ name + str(i) +"\n")

    def handle(self, *args, **options):
        
        begin = time()
        
        # grab PDB
        pdb_code = options.get('pdb_code', None).upper()
        
        t1 = time()
        # read pdb structure (from RCSB) using Biopython
        structure = self.load_pdb(pdb_code)
        t2 = time()
        print("PDB file loaded:", t2-t1)
        
        # get preferred chain for PDB-code
        reference = Structure.objects.get(pdb_code__index=pdb_code)
        preferred_chain = reference.preferred_chain.split(',')[0]

        # grab residues with the generic numbering for this structure
        # TODO: one query?
        p = Protein.objects.get(protein=reference.protein_conformation.protein)
        db_reslist = Residue.objects.filter(protein_conformation__protein=p).prefetch_related( 'generic_number')
        
        start = time()
        print("Database structures fetched:",start-t2)
        
        #######################################################################
        ############################# filter  pdb #############################
        
        db_tmlist = [[] for i in range(7)]
        db_set    = set()
        for r in db_reslist:
            if r.generic_number and r.generic_number.label[:2] in ["1x","2x","3x","4x","5x","6x","7x"]:
                db_tmlist[int(r.generic_number.label[0])-1].append(r)
                db_set.add(r.sequence_number)
        
        #better if no tmlist
        #db_set = {r.sequence_number for r in db_reslist if r.generic_number and  r.generic_number.label[:2] in ["1x","2x","3x","4x","5x","6x","7x"]} # r.protein_segment.label in [1,2,3,4,5,6,7] } #

        def recurse(entity,slist):
            if slist:
                for subenty in entity.get_list():
                    if not slist[0](subenty): entity.detach_child(subenty.get_id())
                    else: recurse(subenty, slist[1:])

        recurse(structure,[lambda entity: True if entity.get_id() == 0 else False, #model
                           lambda entity: True if entity.get_id() == preferred_chain else False, #chain
                           lambda entity: True if entity.get_id()[1] in db_set else False, # residue
                           lambda entity: True if entity.get_id() == "CA" else False])

        io1 = pdb.PDBIO()
        io1.set_structure(structure)
        io1.save('outS.pdb')
        
        pchain = structure[0][preferred_chain]
        
        end = time()
        print("pdb filtered and out:", end-start)
        
        #######################################################################
        ############### Calculate the axes through the helices ################
        #######################################################################
        N = 3
        
        # TODO: "in" inefficient?
        hres_list = [np.asarray([pchain[r.sequence_number]["CA"].get_coord() for r in sl if r.sequence_number in pchain], dtype=float) for sl in db_tmlist]

        # fast and fancy way to take the average of N consecutive elements
        hres_three = np.asarray([sum([h[i:-(len(h) % N) or None:N] for i in range(N)])/N for h in hres_list])
        helix_mean = np.asarray([np.mean(h, axis=0) for h in hres_three ])
        self.save_pseudo(hres_three, "helper")
        
        #######################################################################
        ################################ Mean #################################
        
        pos_list = [(np.average(tm[-6:-3],axis=0), np.average(tm[3:6],axis=0)) for tm in hres_list]
        self.write_cgo_arrow_pml("no_first_three",pos_list)

        vmean = np.asarray([sum( np.sum(h[i+2:] - l,axis=0) for i,l in enumerate(h[:-2])) for h in hres_three])
        vmean /= vmean[:,0].reshape(-1,1)

        pos_list = [(h[0],h[0]+5*v) for h,v in zip(hres_three,vmean)]
        self.write_cgo_arrow_pml("vvv",pos_list)
        
        #######################################################################
        ################################## LR #################################
        
        def lr_line(regr,s):
            regr.fit(s[:,:2],s[:,2].reshape(-1,1))
            return regr.predict(s[[0,-1],:2])

        regr = linear_model.LinearRegression()
        
        pos_list = [np.hstack((s[[0,-1],:2],lr_line(regr, s))) for s in hres_three]
        self.write_cgo_arrow_pml("LR",pos_list)
        
        pos_list = [np.hstack((s[[0,-1],:2],lr_line(regr, s))) for s in [np.vstack(hres_three)]]
        self.write_cgo_arrow_pml("LR_all",pos_list)
        
        #######################################################################
        ################################# PCA #################################
        
        pca = PCA()
        #def pca_line(pca,h, r=0): return (pca.inverse_transform(np.asarray([[-20,0,0],[20,0,0]])) if ((not r) if pca.fit_transform(h)[0][0] < 0 else r) else pca.inverse_transform(np.asarray([[20,0,0],[-20,0,0]])))
        def pca_line(pca,h, r=0):
            if ((not r) if pca.fit_transform(h)[0][0] < 0 else r):
                return pca.inverse_transform(np.asarray([[-20,0,0],[20,0,0]]))
            else:return pca.inverse_transform(np.asarray([[20,0,0],[-20,0,0]]))        
        
        pos_list = np.asarray([pca_line(pca, h,i%2) for i,h in enumerate(hres_three)])
        self.write_cgo_arrow_pml("pca",pos_list)
        comb = pos_list[:,1] - pos_list[:,0]
        print("TM:",repr(comb / np.linalg.norm(comb,axis=1).reshape(-1,1)),sep="\n")
        
        pos_list = np.mean(pos_list,axis=0)
        self.write_cgo_arrow_pml("pca_mean",[pos_list])
        comb = pos_list[1] - pos_list[0]
        print("Mean:",repr(comb / np.linalg.norm(comb)),sep="\n")
        
        pos_list = pca_line(pca, np.vstack(hres_three))
        self.write_cgo_arrow_pml("pca_all",[pos_list])
        comb = pos_list[1] - pos_list[0]
        print("All:",repr(comb / np.linalg.norm(comb)),sep="\n")
        
        pos_list = np.asarray([pca_line(pca, h[:len(h)//2:-(-(i%2) or 1)]) for i,h in enumerate(hres_three)])
        self.write_cgo_arrow_pml("pca_intra",pos_list)
        
        pos_list = np.asarray([pca_line(pca, h[:len(h)//2:(-(i%2) or 1)]) for i,h in enumerate(hres_three)])
        self.write_cgo_arrow_pml("pca_extra",pos_list)
        
        pos_list = pos_list - (np.mean(pos_list,axis=1)-helix_mean).reshape(-1,1,3)
        self.write_cgo_arrow_pml("pca_extra2",pos_list)
        
        self.write_cgo_arrow_pml("pca_extra2_mean",[np.mean(pos_list,axis=0)])
        
        print("Total runtime:", time()-begin)
        # Use db numbering for figuring out the rigid half
