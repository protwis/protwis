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
        with open(name, "w") as ps:
            ps.write("run cgo_arrow.py\n")
            for i, p in enumerate(pos_list):
                ps.write("pseudoatom " +str(i) + "_bot, pos=" + str(list(p[0])) + "\n" \
                         "pseudoatom " +str(i) + "_top, pos=" + str(list(p[1])) + "\n" \
                         "cgo_arrow  " +str(i) + "_bot, "   + str(i)+"_top"     + "\n" )

    def handle(self, *args, **options):
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
        
        #######################################
        
        db_tmlist = [[] for i in range(7)]
        db_set    = set()
        for r in db_reslist:
            if r.generic_number and r.generic_number.label[:2] in ["1x","2x","3x","4x","5x","6x","7x"]:
                db_tmlist[int(r.generic_number.label[0])-1].append(r)
                db_set.add(r.sequence_number)
        
        #better if no tmlist
        #db_set = {r.sequence_number for r in db_reslist if r.generic_number and  r.generic_number.label[:2] in ["1x","2x","3x","4x","5x","6x","7x"]} # r.protein_segment.label in [1,2,3,4,5,6,7] } #
        
        class MySelect2(pdb.Select):
            def accept_atom(self, atom):
                return (True if (atom.get_name() == "CA") else False)
            
            def accept_residue(self, residue):
                return (True if residue.get_id()[1] in db_set else False)
            
            def accept_chain(self, chain):
                return (True if chain.get_id() == preferred_chain else False)
            
            def accept_model(self, model):
                return (True if model.get_id() == 0 else False)
        
        select = MySelect2()
        for model in structure:
            if not select.accept_model(model):
                structure.detach_child(model.get_id())
                continue
            for chain in model.get_list():
                if not select.accept_chain(chain):
                    model.detach_parent(chain.get_id())
                    continue
                for residue in chain.get_unpacked_list():
                    if not select.accept_residue(residue):
                        chain.detach_child(residue.get_id())
                        continue
                    for atom in residue.get_unpacked_list():
                        if not select.accept_atom(atom):
                            residue.detach_child(atom.get_id())

        io1 = pdb.PDBIO()
        io1.set_structure(structure)
        io1.save('outS.pdb')
        
        pchain = structure[0][preferred_chain]
        
        end = time()
        print("pdb filtered and out:", end-start)
        
        #######################################################################
        # Calculate the axss through the helices
        
        # TODO: "in" inefficient?
        hres_list = [np.asarray([pchain[r.sequence_number]["CA"].get_coord() for r in sl if r.sequence_number in pchain], dtype=float) for sl in db_tmlist]

        N = 3
        # fast and fancy way to take the average of three consecutive elements
        hres_three = np.asarray([sum([h[i:-(len(h) % N) or None:N] for i in range(N)])/N for h in hres_list])
        self.save_pseudo(hres_three, "helper")    
        
        pos_list = [(np.average(tm[-6:-3],axis=0), np.average(tm[3:6],axis=0)) for tm in hres_list]
        self.write_cgo_arrow_pml("no_first_three.pml",pos_list)

        vmean = np.asarray([sum( np.sum(h[i+2:] - l,axis=0) for i,l in enumerate(h[:-2])) for h in hres_three])
        vmean /= vmean[:,0].reshape(-1,1)

        yyy = [(h[0],h[0]+5*v) for h,v in zip(hres_three,vmean)]
        self.write_cgo_arrow_pml("vvv.pml",yyy)
        
        
        def lr_line(regr,s):
            regr.fit(s[:,:2],s[:,2].reshape(-1,1))

            return regr.predict(s[[0,-1],:2])

        regr = linear_model.LinearRegression()
        
        pos_list = [np.hstack((s[[0,-1],:2],lr_line(regr, s))) for s in hres_three]
        self.write_cgo_arrow_pml("LR.pml",pos_list)
        
        pos_list = [np.hstack((s[[0,-1],:2],lr_line(regr, s))) for s in [np.vstack(hres_three)]]
        self.write_cgo_arrow_pml("LR_all.pml",pos_list)
        
        pca = PCA()
        def pca_line(pca,h):
            pca.fit(h)
            return pca.inverse_transform(np.asarray([[20,0,0],[-20,0,0]]))
        
        pos_list = [pca_line(pca, h) for h in hres_three]
        self.write_cgo_arrow_pml("pca.pml",pos_list)
        
        pos_list = [pca_line(pca, h) for h in [np.vstack(hres_three)]]
        self.write_cgo_arrow_pml("pca_all.pml",pos_list)
