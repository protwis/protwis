from django.core.management.base import BaseCommand#, CommandError
#from django.core.management import call_command
#from django.conf import settings
#from django.db import connection

import contactnetwork.pdb as pdb

from structure.models import Structure
from residue.models import Residue

import logging, os
import numpy as np
from sklearn.decomposition import PCA
from numpy.core.umath_tests import inner1d
import io
import freesasa
import pickle
TMNUM = 7

SASA = False
HSE  = True

class Command(BaseCommand):

    help = "Command to calculate an axis for a TM helix."
    
    np.set_printoptions(suppress=True)
    logger = logging.getLogger(__name__)

#    def add_arguments(self, parser):
#        parser.add_argument('pdb_code')
    
    ###########################################################################
    ############################ Helper  Functions ############################
    ###########################################################################
        
    def load_pdb(self,pdb_code, no_save=False):
        if no_save:
            return pdb.pdb_get_structure(pdb_code)
        else:
            if not os.path.exists("pdbfiles/pdb" + pdb_code.lower() + ".ent"):
                pdbl = pdb.PDBList()
                pdbl.retrieve_pdb_file(pdb_code, file_format="pdb", pdir="pdbfiles")
                
            parser = pdb.PDBParser(QUIET=True)
            return parser.get_structure(pdb_code, "pdbfiles/pdb" + pdb_code.lower() + ".ent")

    def load_pdb_var(self,pdb_code, var):
        parser = pdb.PDBParser(QUIET=True)
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
        io1.save("pymol_output/"+ pname + '.pdb')
        
    def write_cgo_arrow_pml(self, pdb_code, name, pos_list):
        with open("pymol_output/"+pdb_code + name + ".pml", "w") as ps:
            ps.write("run cgo_arrow.py\n")
            for i, p in enumerate(pos_list):
                ps.write("cgo_arrow  " + str(list(p[0])) +", "+ str(list(p[1])) + ", name="+pdb_code + name + str(i) +"\n")

    def save_pdb(self, strct, name):
        io1 = pdb.PDBIO()
        io1.set_structure(strct)
        io1.save("pymol_output/"+name)

    def handle(self, *args, **options):
        
        def recurse(entity,slist):
            for subenty in entity.get_list():
                if not subenty.id in slist[0]: entity.detach_child(subenty.id)
                elif slist[1:]: recurse(subenty, slist[1:])
        
        def pca_line(pca,h, r=0):
            if ((not r) if pca.fit_transform(h)[0][0] < 0 else r):
                return pca.inverse_transform(np.asarray([[-20,0,0],[20,0,0]]))
            else:return pca.inverse_transform(np.asarray([[20,0,0],[-20,0,0]])) 
        
        def calc_angle(b,c):
            ba = -b
            bc = c + ba
            ba[:,0] = 0
            return np.degrees(np.arccos(inner1d(ba, bc) / (np.linalg.norm(ba,axis=1) * np.linalg.norm(bc,axis=1))))
        
        def ca_cb_calc(i,pca):
            fin = np.isfinite(h_cb_list[i][:,0])
            return calc_angle(pca.transform(hres_list[i][fin]),pca.transform(h_cb_list[i][fin]))
        
        def axes_calc(i,pca_list,pca):
            p = pca_list[i]
            h = hres_list[i]
            a = (np.roll(np.vstack((h,h[0])),1,axis=0)[:-1] + h + np.roll(np.vstack((h,h[-1])),-1,axis=0)[:-1])/3
            b = p.transform(h)
            b[:,1:] = p.transform(a)[:,1:]
            b = p.inverse_transform(b)
            return calc_angle(pca.transform(b),pca.transform(h))
        
        def set_bfactor(chain,angles):
            for r,an in zip(chain.get_list(),angles):
                for a in r: a.set_bfactor(an)

        
        failed = []
        
        # get preferred chain for PDB-code
        references = Structure.objects.filter(protein_conformation__protein__family__slug__startswith="001").prefetch_related('pdb_code','pdb_data')
        
        for reference in references:
            
            preferred_chain = reference.preferred_chain.split(',')[0]
            pdb_code = reference.pdb_code.index
            state_id = reference.protein_conformation.state.id
            try:
                print(pdb_code)
                if "refined" in pdb_code:
                    failed.append(pdb_code)
                    continue
                #structure = self.load_pdb(pdb_code)
                structure = self.load_pdb_var(pdb_code,reference.pdb_data.pdb)
                # grab residues with the generic numbering for this structure
                db_reslist = list(Residue.objects.exclude(generic_number__isnull=True).filter(protein_conformation__protein=reference.protein_conformation.protein).prefetch_related('generic_number'))
                
                #######################################################################
                ############################# filter  pdb #############################
                
                db_tmlist = [[] for i in range(TMNUM)]
                db_genlst = [[] for i in range(TMNUM)]
                gdict  = {}
                db_set = set()
                for r in db_reslist:
                    if r.generic_number.label[:2] in ["1x","2x","3x","4x","5x","6x","7x"] and r.sequence_number in structure[0][preferred_chain]:
                        if not r.sequence_number > 1000:
                            db_tmlist[int(r.generic_number.label[0])-1].append(r.sequence_number)
                            db_genlst[int(r.generic_number.label[0])-1].append(int(r.generic_number.label[2:]))
                            gdict[r.sequence_number] = r.generic_number.label
                            db_set.add((' ',r.sequence_number,' '))
        
                recurse(structure, [[0], preferred_chain, db_set])
                pchain = structure[0][preferred_chain]
                
                #######################################################################
                ############### Calculate the axes through the helices ################
                #######################################################################
                N = 3
                
                hres_list = [np.asarray([pchain[r]["CA"].get_coord() for r in sl], dtype=float) for sl in db_tmlist]
                h_cb_list = [np.asarray([pchain[r]["CB"].get_coord() if "CB" in pchain[r] else np.array([None,None,None]) for r in sl], dtype=float) for sl in db_tmlist]

                # fast and fancy way to take the average of N consecutive elements
                hres_three = np.asarray([sum([h[i:-(len(h) % N) or None:N] for i in range(N)])/N for h in hres_list])
                helices_mn = np.asarray([np.mean(h, axis=0) for h in hres_three ])
                self.save_pseudo(hres_three, pdb_code+"helper")
                
                #######################################################################
                ################################# PCA #################################
                #######################################################################
                
                helix_pcas = [PCA() for i in range(7)]
                pos_list = np.asarray([pca_line(helix_pcas[i], h,i%2) for i,h in enumerate(hres_three)])
                #self.write_cgo_arrow_pml(pdb_code, "pca",pos_list)
                
                pos_list = np.mean(pos_list,axis=0)
                #self.write_cgo_arrow_pml(pdb_code, "pca_mean",[pos_list])
                
                pca = PCA()
                pos_list = pca_line(pca, np.vstack(hres_three))
                #self.write_cgo_arrow_pml(pdb_code, "pca_all",[pos_list])
                
                pos_list = np.asarray([pca_line(PCA(), h[:len(h)//2:(-(i%2) or 1)]) for i,h in enumerate(hres_three)])
                pos_list = pos_list - (np.mean(pos_list,axis=1)-helices_mn).reshape(-1,1,3)
                #self.write_cgo_arrow_pml(pdb_code, "pca_extra",pos_list)
                #self.write_cgo_arrow_pml(pdb_code, "pca_extra_mean",[np.mean(pos_list,axis=0)])
                
                pca_extra = PCA()
                pos_list = pca_line(pca_extra, np.vstack(pos_list))
                #self.write_cgo_arrow_pml(pdb_code, "pca_extra_pca",[pos_list])
                
                #######################################################################
                ################################ Angles ###############################
                #######################################################################
                
                ########################### Axis to CA to CB ##########################
                centerpca = pca
        
                i = 0
                pl = 0
                plist = pchain.get_unpacked_list()
                for ind in np.where(np.isnan(np.concatenate(h_cb_list)[:,0]))[0]:
                    while ind - pl >= len(h_cb_list[i]):
                        pl += len(h_cb_list[i])
                        i +=1
                    t = plist[ind]
                    n = t['N'].get_vector()
                    c = t['C'].get_vector()
                    ca = t['CA'].get_vector()
                    n = n - ca
                    c = c - ca
                    rot = pdb.rotaxis(-np.pi*120.0/180.0, c)
                    cb_at_origin = n.left_multiply(rot)
                    cb = cb_at_origin + ca
                    h_cb_list[i][ind-pl] = cb.get_array()
                angle = np.concatenate([ca_cb_calc(i,centerpca) for i in range(TMNUM)])
                
                set_bfactor(pchain,angle)
                
                self.save_pdb(structure, pdb_code+'angle_colored_ca_cb.pdb')
                
                ######################### Axis to Axis to CA ##########################
                
                angle2 = np.concatenate([axes_calc(i,helix_pcas,centerpca) for i in range(TMNUM)])
                
                set_bfactor(pchain,angle2)
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
                    self.save_pdb(structure, pdb_code+'asa_colored.pdb')
                
                ################################# HSE #################################
                if HSE:
                    hse = pdb.HSExposure.HSExposureCB(structure[0])
                    [[a.set_bfactor(x[1][1]) for a in x[0]] for x in hse]
                    
                    self.save_pdb(structure, pdb_code+'hsea_colored.pdb')
                
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
        
            except Exception as e:
                print("ERROR!!", pdb_code, e)
                failed.append(pdb_code)
                continue
            
        print(len(failed), "of", len(references), "failed:", failed)