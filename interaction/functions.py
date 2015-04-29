#import urllib2
import urllib
import subprocess
import time
import os.path
import sys
import getopt

from Bio.PDB import *
import openbabel
#sudo apt-get install e-openbabel
#sudo python3 -m pip install openbabel

import pybel
import yaml

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

import re

import os
from collections import Counter
import numpy as np
import collections
from math import pi,degrees
from operator import itemgetter, attrgetter, methodcaller

import getopt
import sys

AA = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D',
     'CYS':'C', 'GLN':'Q', 'GLU':'E', 'GLY':'G',
     'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K',
     'MET':'M', 'PHE':'F', 'PRO':'P', 'SER':'S',
     'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V'}

AROMATIC = {'TYR','TRP','PHE','HIS'}

CHARGEDAA = {'ARG','LYS','ASP','GLU'} #skip ,'HIS'

module_dir = os.path.dirname(__file__)
projectdir = module_dir + '/temp/'
if not os.path.exists(projectdir):
        os.makedirs(projectdir)
        os.makedirs(projectdir+'/temp/')
ignore_het = ['NA','W'] #ignore sodium and water


radius = 4.5
hydrophob_radius = 4.5
ignore_het = ['NA','W'] #ignore sodium and water




def fetch_pdb(id):
  url = 'http://www.rcsb.org/pdb/files/%s.pdb' % id
  return urllib.urlopen(url).read()


def check_unique_ligand_mol(filename):
    f_in = open(filename, 'r')
    tempstr = ''
    check = []
    #print filename
    ligandid = 0
    chainid = 0
    for line in f_in:
        if line.startswith('HETATM'): 
            temp = line.split()
            m = re.match("(\d+\.\d{2})([\d\.]+)",temp[8]) ### need to fix bad PDB formatting where col4 and col5 are put together for some reason -- usually seen when the id is +1000
            if (m):
                temp.extend(temp[9])
                temp[9] = m.group(2)
                temp[8] = m.group(1)
            m = re.match("(\w)(\d+)",temp[4]) ### need to fix bad PDB formatting where col4 and col5 are put together for some reason -- usually seen when the id is +1000
            if (m):
                temp.extend(temp[10])

                temp[10] = temp[9]
                temp[9] = temp[8]
                temp[8] = temp[7]
                temp[7] = temp[6]
                temp[6] = temp[5]
                
                temp[4] = m.group(1)
                temp[5] = m.group(2)

            if (temp[5]!=ligandid and ligandid!=0) or (temp[4]!=chainid and chainid!=0): continue

            ligandid = temp[5]
            chainid = temp[4]

        tempstr += line
    #print tempstr
    f_in.close();
    f=open(filename,'w')
    f.write(tempstr)
    f.close();


def check_pdb():
    if not os.path.exists(projectdir+'pdbs/'):
        os.makedirs(projectdir+'pdbs/')


    if not os.path.isfile(projectdir+'pdbs/'+pdbname+'.pdb'):
        pdbfile = fetch_pdb(pdbname)
        temp_path = projectdir+'pdbs/'+pdbname+'.pdb'
        f=open(temp_path,'w')
        f.write(pdbfile)
        f.close();


def checkdirs():
    if not os.path.isfile(projectdir+'pdbs/'+pdbname+'.pdb'):
        pdbfile = fetch_pdb(pdbname)
        temp_path = projectdir+'pdbs/'+pdbname+'.pdb'
        f=open(temp_path,'w')
        f.write(pdbfile)
        f.close();

    directory = projectdir + 'results/'+pdbname+'/interaction'
    if not os.path.exists(directory):
        os.makedirs(directory)
    directory = projectdir + 'results/'+pdbname+'/ligand'
    if not os.path.exists(directory):
        os.makedirs(directory)
    directory = projectdir + 'results/'+pdbname+'/output'
    if not os.path.exists(directory):
        os.makedirs(directory)
    directory = projectdir + 'results/'+pdbname+'/png'
    if not os.path.exists(directory):
        os.makedirs(directory)
    directory = projectdir + 'results/'+pdbname+'/fragments'
    if not os.path.exists(directory):
        os.makedirs(directory)

def fragment_library(ligand,atomvector,atomname,residuenr,typeinteraction):
    #print "Make fragment pdb file for ligand:",ligand,"atom vector",atomvector,"atomname",atomname,"residuenr from protein", residuenr

    ligand_pdb = projectdir+'results/'+pdbname+'/ligand/'+ligand+'_'+pdbname+'.pdb'
    #print "Look in",ligand_pdb
    mol = pybel.readfile("pdb", ligand_pdb).next()
    mol.removeh()
    listofvectors = []
    for atom in mol:
        distance = (Vector(getattr(atom,'coords'))-atomvector).norm()
        if distance>0.1:
            continue
        #print "Parent:",getattr(atom,'type'),getattr(atom,'idx') ,Vector(getattr(atom,'coords'))
        listofvectors.append(Vector(getattr(atom ,'coords')))
        for neighbour_atom in openbabel.OBAtomAtomIter(atom.OBAtom):
           #print neighbour_atom.GetAtomicNum()
           neighbor = pybel.Atom(neighbour_atom)
           #print "Neighbour:",neighbour_atom.GetType(),Vector(getattr(neighbor,'coords'))
           listofvectors.append(Vector(getattr(neighbor ,'coords')))
           for neighbour_atom2 in openbabel.OBAtomAtomIter(neighbour_atom):
           #print neighbour_atom.GetAtomicNum()
               neighbor2 = pybel.Atom(neighbour_atom2)
               #print "Neighbour2:",neighbour_atom2.GetType(),Vector(getattr(neighbor2,'coords'))
               listofvectors.append(Vector(getattr(neighbor2,'coords')))
    #print "vectors:",listofvectors

    pdbfile = projectdir+'pdbs/'+pdbname+'.pdb'


    f_in = open(pdbfile, 'r')
    tempstr = ''
    for line in f_in:
        if line.startswith('HETATM'): 
            temp = line.split()

            m = re.match("(\d+\.\d{2})([\d\.]+)",temp[8]) ### need to fix bad PDB formatting where col4 and col5 are put together for some reason -- usually seen when the id is +1000
            if (m):
                temp.extend(temp[9])
                temp[9] = m.group(2)
                temp[8] = m.group(1)
            m = re.match("(\w)(\d+)",temp[4]) ### need to fix bad PDB formatting where col4 and col5 are put together for some reason -- usually seen when the id is +1000
            if (m):
                temp.extend(temp[10])

                temp[10] = temp[9]
                temp[9] = temp[8]
                temp[8] = temp[7]
                temp[7] = temp[6]
                temp[6] = temp[5]
                
                temp[4] = m.group(1)
                temp[5] = m.group(2)
            
            atomvector = Vector(temp[6],temp[7],temp[8])
            skip = 1
            for targetvector in listofvectors:
                distance = (targetvector-atomvector).norm()
                if distance<0.1:
                    #print "FOUND!"
                    skip = 0
            if skip==1: continue
        elif line.startswith('ATOM'): 
            temp = line.split()
            m = re.match("(\d+\.\d{2})([\d\.]+)",temp[8]) ### need to fix bad PDB formatting where col4 and col5 are put together for some reason -- usually seen when the id is +1000
            if (m):
                temp.extend(temp[9])
                temp[9] = m.group(2)
                temp[8] = m.group(1)
            m = re.match("(\w)(\d+)",temp[4]) ### need to fix bad PDB formatting where col4 and col5 are put together for some reason -- usually seen when the id is +1000
            if (m):
                temp.extend(temp[10])

                temp[10] = temp[9]
                temp[9] = temp[8]
                temp[8] = temp[7]
                temp[7] = temp[6]
                temp[6] = temp[5]
                
                temp[4] = m.group(1)
                temp[5] = m.group(2)
            if temp[5]!=residuenr:
                continue
            residuename = temp[3]
            chain = temp[4]
        else:
            continue #ignore all other lines

        tempstr += line


    filename = projectdir + 'results/'+pdbname+'/fragments/'+pdbname+"_"+ligand+"_"+residuename+residuenr+chain+"_"+atomname+"_"+typeinteraction+".pdb"
    #print tempstr
    f_in.close();
    f=open(filename,'w')
    f.write(tempstr)
    f.close();

def fragment_library_aromatic(ligand,atomvectors,residuenr,ringnr):
    #print "Make aromatic fragment pdb file for ligand:",ligand,"atom vectors",atomvectors,"residuenr from protein", residuenr

    pdbfile = projectdir+'pdbs/'+pdbname+'.pdb'
    residuename = ''

    f_in = open(pdbfile, 'r')
    tempstr = ''
    for line in f_in:
        if line.startswith('HETATM'): 
            temp = line.split()
            m = re.match("(\d+\.\d{2})([\d\.]+)",temp[8]) ### need to fix bad PDB formatting where col4 and col5 are put together for some reason -- usually seen when the id is +1000
            if (m):
                temp.extend(temp[9])
                temp[9] = m.group(2)
                temp[8] = m.group(1)
            m = re.match("(\w)(\d+)",temp[4]) ### need to fix bad PDB formatting where col4 and col5 are put together for some reason -- usually seen when the id is +1000
            if (m):
                temp.extend(temp[10])

                temp[10] = temp[9]
                temp[9] = temp[8]
                temp[8] = temp[7]
                temp[7] = temp[6]
                temp[6] = temp[5]
                
                temp[4] = m.group(1)
                temp[5] = m.group(2)
            
            atomvector = Vector(temp[6],temp[7],temp[8])
            skip = 1
            for targetvector in atomvectors:
                distance = (targetvector-atomvector).norm()
                if distance<0.1:
                    #print "FOUND!"
                    skip = 0
            if skip==1: continue
        elif line.startswith('ATOM'): 
            temp = line.split()
            m = re.match("(\d+\.\d{2})([\d\.]+)",temp[8]) ### need to fix bad PDB formatting where col4 and col5 are put together for some reason -- usually seen when the id is +1000
            if (m):
                temp.extend(temp[9])
                temp[9] = m.group(2)
                temp[8] = m.group(1)
            m = re.match("(\w)(\d+)",temp[4]) ### need to fix bad PDB formatting where col4 and col5 are put together for some reason -- usually seen when the id is +1000
            if (m):
                temp.extend(temp[10])

                temp[10] = temp[9]
                temp[9] = temp[8]
                temp[8] = temp[7]
                temp[7] = temp[6]
                temp[6] = temp[5]
                
                temp[4] = m.group(1)
                temp[5] = m.group(2)
            if temp[5]!=residuenr:
                continue
            residuename = temp[3]
        else:
            continue #ignore all other lines
        

        tempstr += line


    filename = projectdir + 'results/'+pdbname+'/fragments/'+pdbname+"_"+ligand+"_"+residuename+str(residuenr)+"_aromatic_"+str(ringnr)+".pdb"
    #print tempstr
    f_in.close();
    f=open(filename,'w')
    f.write(tempstr)
    f.close();


def create_ligands_and_poseview():
    
    class HetSelect(Select):
        def accept_residue(self, residue):
            if residue.get_resname().strip()==HETNAM:
                return 1
            else:
                return 0

    p=PDBParser(QUIET=True)
    s = p.get_structure(pdbname, projectdir+'pdbs/'+pdbname+'.pdb') #Disable warnings
    hetflag_done = {}
    for model in s:
        for chain in model:
            for residue in chain:
                hetresname = residue.get_resname()
                hetflag = residue.get_full_id()[3][0].strip() #catch residues with hetflag

                hetflag = hetflag.replace("H_","").strip()
                #hetflag = hetflag.replace("W","")

                if hetflag and hetflag not in ignore_het: 
                    if not hetflag in hetflag_done: 

                        hetflag_done[hetflag] = 1
                        HETNAM = hetflag

                        temp_path = projectdir+'results/'+pdbname+'/ligand/'+HETNAM+'_'+pdbname+'.sdf'

                        ligand_pdb = projectdir+'results/'+pdbname+'/ligand/'+HETNAM+'_'+pdbname+'.pdb'
                        ligand_sdf = projectdir+'results/'+pdbname+'/ligand/'+HETNAM+'_'+pdbname+'.sdf'
                        ligand_inchi = projectdir+'results/'+pdbname+'/ligand/'+HETNAM+'_'+pdbname+'.inchi'
                        ligand_poseview = projectdir+'results/'+pdbname+'/png/'+pdbname+'_'+HETNAM+'.png'
                        ligand_png = projectdir+'results/'+pdbname+'/png/'+HETNAM+'.png'


                        if not os.path.isfile(ligand_pdb) or 1==1: #if sdf not made, make it #Always make them for now
                            io = PDBIO()
                            io.set_structure(s)
                            io.save(ligand_pdb,HetSelect())


                            check_unique_ligand_mol(ligand_pdb)


                            if len(list(pybel.readfile("pdb", ligand_pdb)))==0:
                                continue

                            obConversion = openbabel.OBConversion()
                            obConversion.SetInAndOutFormats("pdb", "inchi")
                            obConversion.SetOptions("K", obConversion.OUTOPTIONS)
                            mol = openbabel.OBMol()
                            obConversion.ReadFile(mol, ligand_pdb)   # Open Babel will uncompress automatically
                            obConversion.WriteFile(mol, ligand_inchi)
                            inchikey = obConversion.WriteString(mol)

                            inchikeys[HETNAM] = inchikey.strip()

                            #smiles[HETNAM] = smile

                            smiles[HETNAM] = pybel.readfile("pdb", ligand_pdb).next().write("smi").split("\t")[0]
 

                            mol = pybel.readfile("pdb", ligand_pdb).next()
                            mol.OBMol.AddHydrogens(False, True, 7.4)
                            mol.write("pdb", ligand_pdb,overwrite=True)

                            obConversion = openbabel.OBConversion()
                            obConversion.SetInAndOutFormats("pdb", "sdf")
                            mol = openbabel.OBMol()
                            obConversion.ReadFile(mol, ligand_pdb)   # Open Babel will uncompress automatically

                            obConversion.WriteFile(mol, ligand_sdf)
                           
                        
                        if not os.path.isfile(ligand_png):  #if png of ligand not made, make it
                            m = Chem.MolFromMolFile(ligand_sdf)
                            Draw.MolToFile(m,ligand_png)


                        if not os.path.isfile(ligand_poseview) and 1==2:  #if interaction png not made, make it #SKIP poseview stuff
                            cmd = "poseview -l "+ligand_sdf+" -p "+projectdir+"pdbs/"+pdbname+".pdb -o "+ligand_poseview

                            #print('Running cmd ' + cmd)
                            proc = subprocess.Popen([cmd], stdout=subprocess.PIPE, shell=True)
                            while proc.poll() is None:
                                time.sleep(1)
                            
                            #(out, err) = proc.communicate()
                        else:
                            #print "Already made Poseview:",pdbname+"_"+HETNAM+".png"
                            continue
    #print "Done "+str(len(hetflag_done))

def addresiduestoligand(ligand,pdb,residuelist):
    temp_path = projectdir+'pdbs/'+pdb+'.pdb'
    f_in = open(temp_path, 'r')
    inserstr = ''
    check = []
    #print filename
    ligandid = 0
    chainid = 0
    for line in f_in:
        if line.startswith('ATOM'): 
            temp = line.split()
            m = re.match("(\w)(\d+)",temp[4]) ### need to fix bad PDB formatting where col4 and col5 are put together for some reason -- usually seen when the id is +1000
            if (m):
                temp[4] = m.group(1)
                temp[5] = m.group(2)

            aaname =  temp[3]+temp[5]+temp[4] 

            if aaname in residuelist:
                #print aaname
                inserstr += line
    #print inserstr
    f_in.close();

    #ligands/'+hetflag+'_'+pdbname+".pdb")

    temp_path = projectdir+'results/'+pdbname+'/ligand/'+ligand+'_'+pdb+'.pdb'
    f_in = open(temp_path, 'r')
    tempstr = ''
    inserted = 0
    for line in f_in:
        if line.startswith('ATOM'): 
            temp = line.split()
            if temp[2]=='H': continue #skip hydrogen in model

        if (line.startswith('CONECT') or line.startswith('MASTER')  or line.startswith('END')) and inserted==0:
            tempstr += inserstr 
            inserted = 1
        tempstr += line
    #print tempstr


    #print tempstr
    f_in.close();
    f=open(projectdir+'results/'+pdbname+'/interaction/'+pdb+'_'+ligand+'.pdb','w')
    f.write(tempstr)
    f.close();


def get_ring_from_aa(residueid):

    class AAselect(Select):
        def accept_residue(self, residue):
            #print residue.get_full_id()[3][1],residueid
            if str(residue.get_full_id()[3][1])==residueid:
                return 1
            else:
                return 0
    ptemp=PDBParser(QUIET=True) #disable warnings
    stemp = ptemp.get_structure(pdbname, projectdir+'pdbs/'+pdbname+'.pdb')
    temp_aa_id = residueid

    io = PDBIO()
    io.set_structure(stemp)
    io.save(projectdir+'temp/'+residueid+'.pdb',AAselect())

    mol = pybel.readfile("pdb", projectdir+'temp/'+residueid+'.pdb').next()
    #print hetflag
    rings = getattr(mol,"OBMol").GetSSSR()
    ringlist = []
    for ring in rings:
        center = Vector(0.0, 0.0, 0.0)
        members = ring.Size()
        if ring.IsAromatic():
            atomlist = []
            atomnames = []
            for atom in mol:
                if ring.IsMember( atom.OBAtom): 
                    a_vector = Vector(getattr(atom,'coords'))
                    center += a_vector
                    atomlist.append(atom.idx)
                    atomnames.append(getattr(atom,'type'))
            center = center/members
            normal = center-a_vector #vector in plane
            ringlist.append([atomlist,center,normal,atomnames])
    return ringlist

def get_hydrogen_from_aa(residueid):

    class AAselect(Select):
        def accept_residue(self, residue):
            #print residue.get_full_id()[3][1],residueid
            if str(residue.get_full_id()[3][1])==residueid:
                return 1
            else:
                return 0
    ptemp=PDBParser(QUIET=True)
    stemp = ptemp.get_structure(pdbname, projectdir+'pdbs/'+pdbname+'.pdb')
    temp_aa_id = residueid

    io = PDBIO()
    io.set_structure(stemp)
    io.save(projectdir+'temp/'+residueid+'.pdb',AAselect())

    mol = pybel.readfile("pdb", projectdir+'temp/'+residueid+'.pdb').next()


    mol.OBMol.AddHydrogens(False, True, 7.4)
    #print hetflag
    donors = []
    for atom in mol:
        if getattr(atom,'OBAtom').IsHbondDonor():
            chargevector = Vector(getattr(atom,'coords'))
            #print getattr(atom,'type')," is Donor",chargevector
            temphatoms = []
            for neighbor in pybel.ob.OBAtomAtomIter(atom.OBAtom):
                neighbor = pybel.Atom(neighbor)
                if getattr(neighbor,'type')=="H":
                    #print "neighbor Atom",getattr(neighbor,'type'),"Coords:",getattr(neighbor,'coords')
                    temphatoms.append(Vector(getattr(neighbor,'coords')))

            donors.append([getattr(atom,'type'),chargevector,temphatoms])
        
        if getattr(atom,'OBAtom').IsHbondAcceptor():
            chargevector = Vector(getattr(atom,'coords'))

    return donors


def build_ligand_info():
    count_atom_ligand = {}
    p=PDBParser(QUIET=True)
    s = p.get_structure(pdbname, projectdir+'pdbs/'+pdbname+'.pdb')
    for model in s:
        for chain in model:
            for residue in chain:
                hetresname = residue.get_resname()
                hetflag = residue.get_full_id()[3][0].strip() #catch residues with hetflag

                hetflag = hetflag.replace("H_","").strip()
                #hetflag = hetflag.replace("W","")

                if hetflag and hetflag not in ignore_het: 
                    #if goodhet!='' and hetflag!=goodhet and "H_"+goodhet!=hetflag: continue ### Only look at the ligand that has an image from poseview made for it.
                    if not hetflag in hetlist: 

                        if len(list(pybel.readfile("pdb", projectdir+'results/'+pdbname+'/ligand/'+hetflag+'_'+pdbname+'.pdb')))==0:
                            #This ligand has no molecules
                            continue
                                

                        hetlist[hetflag] = []
                        ligand_charged[hetflag] = []
                        ligand_donors[hetflag] = []
                        count_atom_ligand[hetflag] = 0



                        mol = pybel.readfile("pdb", projectdir+'results/'+pdbname+'/ligand/'+hetflag+'_'+pdbname+".pdb").next()
                        #print "LIGAND",hetflag
                        

                        rings = getattr(mol,"OBMol").GetSSSR()


                        #http://python.zirael.org/e-openbabel4.html
                        ringlist = []
                        for ring in rings:
                            center = Vector(0.0, 0.0, 0.0)
                            members = ring.Size()
                            if ring.IsAromatic():
                                #print "Found an aromatic ring"
                                atomlist = []
                                atomnames = []
                                vectorlist = []
                                for atom in mol:
                                    if ring.IsMember( atom.OBAtom): 
                                        #print atom.idx,getattr(atom,'type'), ring.IsMember( atom.OBAtom)
                                        a_vector = Vector(getattr(atom,'coords'))
                                        center += a_vector
                                        atomlist.append(atom.idx)
                                        vectorlist.append(a_vector)
                                        atomnames.append(getattr(atom,'type'))
                                center = center/members
                                normal = center-a_vector #vector in plane
                                ringlist.append([atomlist,center,normal,atomnames,vectorlist])

                        ligand_rings[hetflag] = ringlist


                        for atom in mol:
                            #print "Atom",getattr(atom,'type'),"Coords:",getattr(atom,'coords'),"FormalCharge:",getattr(atom,'formalcharge'),"PartialCharge",getattr(atom,'partialcharge')
                            if getattr(atom,'formalcharge')!=0:
                                chargevector = Vector(getattr(atom,'coords'))
                                ligand_charged[hetflag].append([getattr(atom,'type'),chargevector,getattr(atom,'formalcharge')])
                            if getattr(atom,'OBAtom').IsCarboxylOxygen():
                                chargevector = Vector(getattr(atom,'coords'))
                                #print getattr(atom,'type')," is CarboxylOxygen",chargevector
                                ligand_charged[hetflag].append([getattr(atom,'type'),chargevector,-1])
                            if getattr(atom,'OBAtom').IsHbondDonor():
                                chargevector = Vector(getattr(atom,'coords'))
                                #print getattr(atom,'type')," is Donor",chargevector
                                temphatoms = []
                                for neighbor in pybel.ob.OBAtomAtomIter(atom.OBAtom):
                                    neighbor = pybel.Atom(neighbor)
                                    if getattr(neighbor,'type')=="H":
                                        #print "neighbor Atom",getattr(neighbor,'type'),"Coords:",getattr(neighbor,'coords')
                                        temphatoms.append(Vector(getattr(neighbor,'coords')))

                                ligand_donors[hetflag].append([getattr(atom,'type'),chargevector,temphatoms])
                            
                            if getattr(atom,'OBAtom').IsHbondAcceptor():
                                chargevector = Vector(getattr(atom,'coords'))
                                #print getattr(atom,'type')," is Acceptor",chargevector
                                #ligand_charged[hetflag].append([getattr(atom,'type'),chargevector,-1])


                        #Function to get ligand centers to maybe skip some residues
                        check = 0
                        center = Vector(0.0, 0.0, 0.0)

                        for atom in residue:
                            if check==0 and hetflag in ligand_atoms: continue #skip when there are many of same ligand
                            het_atom = atom.name

                            check = 1

                            atom_vector = atom.get_vector()
                            center += atom_vector
                            hetlist[hetflag].append([hetresname, het_atom,atom_vector])

                            if not hetflag in ligand_atoms: ligand_atoms[hetflag] = [] #make the ligand_atoms ready
                            ligand_atoms[hetflag].append([count_atom_ligand[hetflag],atom_vector,het_atom])
                            count_atom_ligand[hetflag] += 1
                        center = center / count_atom_ligand[hetflag]
                        ligandcenter[hetflag] = [center,count_atom_ligand[hetflag]]




#LOOP OVER RECEPTOR AND FIND INTERACTIONS
def find_interactions():
    global count_calcs, count_skips
    count_atom = 0
    count_skips = 0
    count_calcs = 0
    p=PDBParser(QUIET=True)
    s = p.get_structure(pdbname, projectdir+'pdbs/'+pdbname+'.pdb')
    for model in s:
        for chain in model:
            chainid = chain.get_id()
            for residue in chain:
                aa_resname=residue.get_resname()
                aa_seqid=str(residue.get_full_id()[3][1])
                hetflagtest=str(residue.get_full_id()[3][0]).strip()
                aaname = aa_resname+aa_seqid+chainid


                hetflagtest = hetflagtest.replace("H_","")
                #hetflagtest = hetflagtest.replace("W","")

                if hetflagtest: continue #residue is a hetnam
                if hetflagtest in hetlist: continue #residue is a hetnam
                #print "Looking at ",aa_resname,aa_seqid
                countresidue = count_atom
                #print aaname
                for hetflag, atomlist in hetlist.iteritems(): #could probably make a check here to see if this residue was anywhere near the ligand, otherwise skip the check per atom
                    #print aa_resname
                    ca = residue['CA'].get_vector()
                    if (ca-ligandcenter[hetflag][0]).norm()>ligandcenter[hetflag][1]:
                        #print "skipping"
                        count_skips += 1
                        continue

                    count_atom = countresidue
                    sum = 0
                    hydrophobic_count = 0

                    # if goodhet!='' and hetflag!=goodhet and "H_"+goodhet!=hetflag: continue ### Only look at the ligand that has an image from poseview made for it.
                    tempdistance = radius   

                    for atom in atomlist:
                        hetresname = atom[0]
                        het_atom = atom[1]
                        het_vector = atom[2]
                        hydrophobic_check = 1

                        aaatomlist = []
                        for atom in residue:
                            count_atom += 1
                            aa_vector = atom.get_vector()
                            aa_atom = atom.name
                            aaatomlist.append([count_atom,aa_vector,aa_atom])

                            d=(het_vector-aa_vector)
                            count_calcs += 1
                            if d.norm()<radius:
                                if not hetflag in results: 
                                    results[hetflag] = {}
                                    summary_results[hetflag] = {'score':[],'hbond':[],'hbondplus':[], 'hbond_confirmed' : [],'waals':[],'aromatic':[],'aromaticplus':[],'hydrophobic':[]}
                                if not aaname in results[hetflag]: 
                                    results[hetflag][aaname] = []
                                results[hetflag][aaname].append([het_atom,aa_atom,round(d.norm(),2),het_vector,aa_vector,aa_seqid])
                                tempdistance = round(d.norm(),2)
                                sum += 1
                            if het_atom[0]=='C' and aa_atom[0]=='C' and d.norm()<hydrophob_radius and hydrophobic_check: #if both are carbon then we are making a hydrophic interaction
                                hydrophobic_count +=1
                                hydrophobic_check = 0

                    if hydrophobic_count>2: #min 3 c-c interactions
                        summary_results[hetflag]['hydrophobic'].append([aaname,hydrophobic_count])

                    if sum>5 and aa_resname in AROMATIC:
                        #print aaatomlist
                        #print "Need to analyse aromatic ring in ",aaname#, get_ring_atoms(aaatomlist)
                        #aaring = get_ring_atoms(aaatomlist)
                        aaring = get_ring_from_aa(aa_seqid)
                        if not aaring:
                            #print "Could not find aromatic ring in",aaname
                            continue
                        aaring = aaring[0]
                        center = aaring[1]
                        count = 0
                        for ring in ligand_rings[hetflag]:
                            #print ring
                            count += 1
                            angle = Vector.angle(center-ring[1],ring[2]) #take vector from two centers, and compare against vector from center to outer point -- this will give the perpendicular angel.
                            angle2 = Vector.angle(center-ring[1],aaring[2]) #take vector from two centers, and compare against vector from center to outer point -- this will give the perpendicular angel.
                            #print "angleaa",aaring[2],"anglelig",ring[2]
                            angle_degrees = [round(degrees(angle),1),round(degrees(angle2),1)]
                            distance = (center-ring[1]).norm()
                            #print "Ring #",count,"Distance:",round(distance,2), "Angle:",angle_degrees
                            if distance<5: #poseview uses <5
                                #print "Ring #",count,"Distance:",round(distance,2), "Angle:",round(angle_degrees,2)
                                summary_results[hetflag]['aromatic'].append([aaname,count,round(distance,2),angle_degrees])

                                fragment_library_aromatic(hetflag,ring[4],aa_seqid,count)

                        for charged in ligand_charged[hetflag]:
                            distance = (center-charged[1]).norm()
                            if distance<4.2 and charged[2]>0: ### needs max 4.2 distance to make aromatic+
                                print "Ring #",count,"Distance:",round(distance,2), "Angle:",round(angle_degrees,2)
                                summary_results[hetflag]['aromaticplus'].append([aaname,count,round(distance,2),charged])
                    #print aaname, ligand_rings,hetflag
                    if sum>2 and aa_resname in CHARGEDAA and ligand_rings[hetflag]:
                        #print "check for charged AA to aromatic rings!",aa_resname,hetflag

                        for atom in residue:
                            aa_vector = atom.get_vector()
                            aa_atom = atom.name
                            for ring in ligand_rings[hetflag]:
                                d=(ring[2]-aa_vector).norm()
                                #if d<10: print "aa_atom",aa_atom,aaname,"distance to a ring",d,hetflag,aa_resname


def analyze_interactions():
    for ligand,result in results.iteritems():

        #print "AA close to ligands ("+ligand+"): ",list(result.keys())
        #print "Results for"+ligand
        sortedresults = []
        ligscore = 0
        for residue,interaction in result.iteritems():
            #print residue[0:2]
            sum = 0
            score = 0
            hbond = []
            hbondplus = []
            type = 'waals'
            for entry in interaction:
                hbondconfirmed = []
                if entry[2]<3.3:
                    
                    #print "Likely H-Bond",entry


                    if entry[0][0] == 'C' or entry[1][0] == 'C': 
                        #print "No Hbond possible",entry
                        continue #If either atom is C then no hydrogen bonding
                    #summary_results[ligand]['hbond'].append([residue,entry])

                    aa_donors = get_hydrogen_from_aa(entry[5])
                    hydrogenmatch = 0
                    for donor in aa_donors:
                        d = (donor[1]-entry[4]).norm()
                        if d<0.5:
                            #print 'found donor',residue,d,entry,donor
                            hydrogens = donor[2]
                            for hydrogen in hydrogens:
                                hydrogenvector = hydrogen-donor[1]
                                bindingvector = entry[3]-hydrogen
                                angle = round(degrees(Vector.angle(hydrogenvector,bindingvector)),2)
                                distance = round(bindingvector.norm(),2)
                                #print "RESDONOR",residue,"From ligand",entry[0],"To AA",entry[1],"HydrogenCheck angle",angle,"Distance from hydrogen to acceptor",distance
                                if distance>2.5:
                                    #print "Too far away"
                                    continue
                                if angle>60:
                                    #print "Bad angle"
                                    continue
                                hydrogenmatch = 1
                                hbondconfirmed.append(["D",entry[0],entry[1],angle,distance])


                    #print "aadonors:",aa_donors

                    for donor in ligand_donors[ligand]:
                        #dummy = Bio.PDB.Vector(0.0, 0.0, 0.0)
                        d = (donor[1]-entry[3]).norm()
                        #print charged,d,residue,entry
                        if d<0.5:
                            #print 'found donor',residue,d,entry,donor
                            hydrogens = donor[2]
                            for hydrogen in hydrogens:
                                hydrogenvector = hydrogen-donor[1]
                                bindingvector = entry[4]-hydrogen
                                angle = round(degrees(Vector.angle(hydrogenvector,bindingvector)),2)
                                distance = round(bindingvector.norm(),2)
                                #print "LIGDONOR",residue,"From ligand",entry[0],"To AA",entry[1],"HydrogenCheck angle",angle,"Distance from hydrogen to acceptor",distance
                                if distance>2.5:
                                    #print "Too far away"
                                    continue
                                if angle>60:
                                    #print "Bad angle"
                                    continue
                                hydrogenmatch = 1
                                hbondconfirmed.append(["A",entry[0],entry[1],angle,distance])
                                
                            #chargedcheck = 1
                        # elif entry[0][0]=='O' and charged[0][0]=='O' and d<2.5:
                        #     print "found O close to O that is charged -- defining as COO-"
                        #     chargedcheck = 1

                    

                    chargedcheck = 0
                    for charged in ligand_charged[ligand]:
                        #dummy = Bio.PDB.Vector(0.0, 0.0, 0.0)
                        d = (charged[1]-entry[3]).norm()
                        #print charged,d,residue,entry
                        if d<0.5:
                            #print 'found charge',residue,d,entry
                            chargedcheck = 1
                            hydrogenmatch = 0 #Replace previous match!

                        # elif entry[0][0]=='O' and charged[0][0]=='O' and d<2.5:
                        #     print "found O close to O that is charged -- defining as COO-"
                        #     chargedcheck = 1

                    if residue[0:3] in CHARGEDAA:
                        #print "check for hbondplus!",residue,entry
                        #Need to check which atoms, but for now assume charged
                        chargedcheck = 1
                        hydrogenmatch = 0 #Replace previous match!



                    if hydrogenmatch:
                        #hbondconfirmed.append(entry)
                        found = 0
                        for x in summary_results[ligand]['hbond_confirmed']:
                            if residue==x[0]:
                                #print "Already key there",residue
                                key = summary_results[ligand]['hbond_confirmed'].index(x)
                                summary_results[ligand]['hbond_confirmed'][key][1].extend(hbondconfirmed)
                                found = 1

                        #if any(residue == x[0] for x in summary_results[ligand]['hbond_confirmed']):
                            
                            #print any(residue == x[0] for x in summary_results[ligand]['hbond_confirmed'])

                        if found==0: summary_results[ligand]['hbond_confirmed'].append([residue,hbondconfirmed])
                        if chargedcheck: 
                            type = 'hbondplus'
                            hbondplus.append(entry)

                        fragment_library(ligand,entry[3],entry[0],entry[5],'HB')

                    elif chargedcheck:
                        type = 'hbondplus'
                        hbondplus.append(entry)
                        fragment_library(ligand,entry[3],entry[0],entry[5],'HBC')
                    else:
                        type = 'hbond'
                        hbond.append(entry)
                        fragment_library(ligand,entry[3],entry[0],entry[5],'HB')


                    entry[3] = ''

                if (entry[2]<4.5): 
                    sum += 1
                    score += 4.5-entry[2]
            score = round(score,2)

            if type == 'waals' and score>2: #mainly no hbond detected
                summary_results[ligand]['waals'].append([residue,score,sum])
            elif type == 'hbond':
                summary_results[ligand]['hbond'].append([residue,score,sum,hbond])
            elif type == 'hbondplus':
                summary_results[ligand]['hbondplus'].append([residue,score,sum,hbondplus])
            # elif type == 'hbond_confirmed':
            #     summary_results[ligand]['hbond_confirmed'].append([residue,score,sum,hbondconfirmed])

            ligscore += score



            #print "Total <4 (score is combined diff from 4)",sum,"score",score
            sortedresults.append([residue,score,sum,hbond,type])
        
        summary_results[ligand]['score'].append([ligscore])
        summary_results[ligand]['inchikey'] = inchikeys[ligand]
        summary_results[ligand]['smiles'] = smiles[ligand]
        #print ligand,"Ligand score:"+str(ligscore) 

        sortedresults = sorted(sortedresults, key=itemgetter(1), reverse=True)  

def pretty_results():
    for ligand,result in summary_results.iteritems():
        output = ''
        bindingresidues = []
        #output += "Results for "+str(ligand)+"\n"
        for type, typelist in result.iteritems():
            if type == 'waals': continue
            output += type+"\n"
            if type == 'waals': 
                typelist = sorted(typelist, key=itemgetter(2), reverse=True)  
            if type == 'hydrophobic': 
                typelist = sorted(typelist, key=itemgetter(1), reverse=True)  
            for entry in typelist:
                if type!='score': bindingresidues.append(entry[0])
                if type == 'hbond':
                    output += '\t'.join(map(str, entry[0:1]))+'\n'
                    for bond in entry[3]:
                        output += '\t'.join(map(str, bond[0:3]))+'\n'
                elif type == 'hbondplus':
                    output += '\t'.join(map(str, entry[0:1]))+'\n'
                    for bond in entry[3]:
                        output += '\t'.join(map(str, bond[0:3]))+'\n'
                elif type == 'hbond_confirmed':
                    output += '\t'.join(map(str, entry[0:1]))+'\n'
                    for bond in entry[1]:
                        output += '\t'.join(map(str, bond))+'\n'
                else:
                    #print entry
                    output += '\t'.join(map(str, entry))+'\n'
        #print output
        # temp_path = projectdir+'results/'+pdbname+'/output/'+pdbname+'_'+ligand.replace("H_","")+'.txt'
        # #print "writing to ",temp_path
        # f=open(temp_path,'w')
        # f.write(output)
        # f.close();
        temp_path = projectdir+'results/'+pdbname+'/output/'+pdbname+'_'+ligand.replace("H_","")+'.yaml'
        #print "writing to ",temp_path
        yaml.dump(result, open(temp_path,'w'))

        addresiduestoligand(ligand,pdbname,bindingresidues)
    #print bindingresidues

    #print "Total atom-pair distance calculations",count_calcs,"Skipped AA",count_skips


def calculate_interactions(pdb):
    global pdbname,hetlist,ligand_atoms,ligand_charged,ligandcenter,ligand_rings,ligand_donors,results,sortedresults,summary_results,inchikeys,smiles

    hetlist = {}
    ligand_atoms = {}
    ligand_charged = {}
    ligandcenter = {}
    ligand_rings = {}
    ligand_donors = {}
    results = {}
    sortedresults = {}
    summary_results = {}
    inchikeys = {}
    smiles = {}

    pdbname = pdb
    #print "checking ",pdbname
    check_pdb()
    checkdirs()
    create_ligands_and_poseview()
    build_ligand_info()
    find_interactions()
    analyze_interactions()

    pretty_results()

def main(argv): 
    pdbname = ''                             
    try:                                
        opts, args = getopt.getopt(argv, "p:", ["pdb"]) 
    except getopt.GetoptError:           
        print "Remember PDB name -p "                      
        sys.exit(2)

    for opt, arg in opts:                
        if opt in ("-p"):      
            pdbname = arg    

    if not pdbname:
        print "Remember PDB name -p "                      
        sys.exit(2)
    #print "looking for",pdbname

    calculate_interactions(pdbname)


if __name__ == "__main__":
    main(sys.argv[1:])
    #pdbname = '1F88'
    #calculate_interactions(pdbname)
