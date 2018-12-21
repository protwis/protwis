#import urllib2
import urllib
import subprocess
import time
import os.path
import sys
import getopt

from Bio.PDB import *
import openbabel

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
from math import pi, degrees
from operator import itemgetter, attrgetter, methodcaller

import getopt
import sys
import shutil

AA = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
      'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
      'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
      'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
      'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}

HBD = {'H', 'K', 'N', 'Q', 'R', 'S', 'T', 'W', 'Y'}
HBA = {'D', 'E', 'H', 'N', 'Q', 'S', 'T', 'Y'}
NEGATIVE = {'D', 'E'}
POSITIVE = {'H', 'K', 'R'}

AROMATIC = {'TYR', 'TRP', 'PHE', 'HIS'}

CHARGEDAA = {'ARG', 'LYS', 'ASP', 'GLU'}  # skip ,'HIS'

HYDROPHOBIC_AA = {'A', 'C', 'F', 'I', 'L', 'M', 'P', 'V', 'W', 'Y'}

projectdir = '/tmp/interactions/'
if not os.path.exists(projectdir):
    os.makedirs(projectdir)
    os.chmod(projectdir, 0o777)
tempdir = projectdir + 'temp/'
if not os.path.exists(tempdir):
    os.makedirs(tempdir)
    os.chmod(tempdir, 0o777)
ignore_het = ['NA', 'W']  # ignore sodium and water


radius = 5
hydrophob_radius = 4.5
ignore_het = ['NA', 'W']  # ignore sodium and water


debug = False


def fetch_pdb(id):
    url = 'http://www.rcsb.org/pdb/files/%s.pdb' % id
    return urllib.urlopen(url).read()


def check_unique_ligand_mol(filename):
    # check that only HETATM are exported to file
    f_in = open(filename, 'r')
    tempstr = ''
    check = []
    ligandid = 0
    chainid = 0
    for line in f_in:
        if line.startswith('HETATM'):
            residue_number = line[22:26]
            chain = line[21]

            if (residue_number != ligandid and ligandid != 0) or (chain != chainid and chainid != 0):
                continue

            ligandid = residue_number
            chainid = chain

        tempstr += line

    f_in.close()
    f = open(filename, 'w')
    f.write(tempstr)
    f.close()


def check_pdb():
    # check if PDB is there, otherwise fetch
    if not os.path.exists(projectdir + 'pdbs/'):
        os.makedirs(projectdir + 'pdbs/')

    if not os.path.isfile(projectdir + 'pdbs/' + pdbname + '.pdb'):
        pdbfile = fetch_pdb(pdbname)
        temp_path = projectdir + 'pdbs/' + pdbname + '.pdb'
        f = open(temp_path, 'w')
        f.write(pdbfile)
        f.close()


def checkdirs():
    # check that dirs are there and have right permissions
    directory = projectdir + 'results/' + pdbname
    if os.path.exists(directory):
        shutil.rmtree(directory)

    directory = projectdir + 'results/' + pdbname + '/interaction'
    if not os.path.exists(directory):
        os.makedirs(directory)
        os.chmod(directory, 0o777)
    directory = projectdir + 'results/' + pdbname + '/ligand'
    if not os.path.exists(directory):
        os.makedirs(directory)
        os.chmod(directory, 0o777)
    directory = projectdir + 'results/' + pdbname + '/output'
    if not os.path.exists(directory):
        os.makedirs(directory)
        os.chmod(directory, 0o777)
    directory = projectdir + 'results/' + pdbname + '/png'
    if not os.path.exists(directory):
        os.makedirs(directory)
        os.chmod(directory, 0o777)
    directory = projectdir + 'results/' + pdbname + '/fragments'
    if not os.path.exists(directory):
        os.makedirs(directory)
        os.chmod(directory, 0o777)


def find_ligand_full_names():
    pdbfile = projectdir + 'pdbs/' + pdbname + '.pdb'
    residuename = ''

    f_in = open(pdbfile, 'r')
    d = {}

    for line in f_in:
        if line.startswith('HETSYN'):
            # need to fix bad PDB formatting where col4 and col5 are put
            # together for some reason -- usually seen when the id is +1000
            m = re.match("HETSYN[\s]+([\w]{3})[\s]+(.+)", line)
            if (m):
                d[m.group(1)] = m.group(2).strip()
    return d


def fragment_library(ligand, atomvector, atomname, residuenr, chain, typeinteraction):
    #if debug:
        #print "Make fragment pdb file for ligand:", ligand, "atom vector", atomvector, "atomname", atomname, "residuenr from protein", residuenr, typeinteraction, 'chain', chain
    residuename = 'unknown'
    ligand_pdb = projectdir + 'results/' + pdbname + \
        '/ligand/' + ligand + '_' + pdbname + '.pdb'
    mol = pybel.readfile("pdb", ligand_pdb).next()
    mol.removeh()
    listofvectors = []
    chain = chain.strip()
    if atomvector is not None:
        for atom in mol:
            distance = (Vector(getattr(atom, 'coords')) - atomvector).norm()
            if distance > 0.1:
                continue
            # print "Parent:",getattr(atom,'type'),getattr(atom,'idx')
            # ,Vector(getattr(atom,'coords'))
            listofvectors.append(Vector(getattr(atom, 'coords')))
            for neighbour_atom in openbabel.OBAtomAtomIter(atom.OBAtom):
                # print neighbour_atom.GetAtomicNum()
                neighbor = pybel.Atom(neighbour_atom)
                # print
                # "Neighbour:",neighbour_atom.GetType(),Vector(getattr(neighbor,'coords'))
                listofvectors.append(Vector(getattr(neighbor, 'coords')))
                for neighbour_atom2 in openbabel.OBAtomAtomIter(neighbour_atom):
                    # print neighbour_atom.GetAtomicNum()
                    neighbor2 = pybel.Atom(neighbour_atom2)
                    # print
                    # "Neighbour2:",neighbour_atom2.GetType(),Vector(getattr(neighbor2,'coords'))
                    listofvectors.append(Vector(getattr(neighbor2, 'coords')))
        #if debug:
            #print "vectors:", listofvectors

    pdbfile = projectdir + 'pdbs/' + pdbname + '.pdb'

    f_in = open(pdbfile, 'r')
    tempstr = ''
    for line in f_in:
        if line.startswith('HETATM'):
            atomvector = Vector(line[30:38], line[38:46], line[46:54])
            residue_number = line[22:26]
            tempchain = line[21]
            skip = 1
            for targetvector in listofvectors:
                distance = (targetvector - atomvector).norm()
                if distance < 0.1:
                    # print "FOUND!"
                    skip = 0
            if skip == 1:
                continue
        elif line.startswith('ATOM'):

            residue_number = line[22:26].strip()
            tempchain = line[21].strip()
            if residue_number != residuenr:
                continue
            if tempchain != chain:
                continue
            residuenr = residue_number
            chain = tempchain
            residuename = line[17:20].strip()
        else:
            continue  # ignore all other lines

        tempstr += line

    filename = projectdir + 'results/' + pdbname + '/fragments/' + pdbname + "_" + ligand + \
        "_" + residuename + residuenr + chain + "_" + \
        atomname + "_" + typeinteraction + ".pdb"
    # if debug:
        # print filename
    f_in.close()
    f = open(filename, 'w')
    f.write(tempstr)
    f.close()
    mol = pybel.readfile("pdb", filename).next()
    mol.write("pdb", filename, overwrite=True)

    return filename


def fragment_library_aromatic(ligand, atomvectors, residuenr, chain, ringnr):
    # print "Make aromatic fragment pdb file for ligand:",ligand,"atom
    # vectors",atomvectors,"residuenr from protein", residuenr
    chain = chain.strip()
    pdbfile = projectdir + 'pdbs/' + pdbname + '.pdb'
    residuename = ''

    f_in = open(pdbfile, 'r')
    tempstr = ''
    for line in f_in:
        if line.startswith('HETATM'):
            atomvector = Vector(line[30:38], line[38:46], line[46:54])

            skip = 1
            for targetvector in atomvectors:
                distance = (targetvector - atomvector).norm()
                if distance < 0.1:
                    # print "FOUND!"
                    skip = 0
            if skip == 1:
                continue
        elif line.startswith('ATOM'):

            residue_number = line[22:26].strip()
            tempchain = line[21].strip()

            if residue_number != residuenr:
                continue
            if tempchain != chain:
                continue
            residuename = line[17:20].strip()
            chain = tempchain
        else:
            continue  # ignore all other lines

        tempstr += line

    filename = projectdir + 'results/' + pdbname + '/fragments/' + pdbname + "_" + ligand + \
        "_" + residuename + str(residuenr) + chain + \
        "_aromatic_" + str(ringnr) + ".pdb"
    # print tempstr
    f_in.close()
    f = open(filename, 'w')
    f.write(tempstr)
    f.close()
    return filename


def create_ligands_and_poseview():

    class HetSelect(Select):

        def accept_residue(self, residue):
            if residue.get_resname().strip() == HETNAM:
                return 1
            else:
                return 0

    class ClassSelect(Select):

        def accept_residue(self, residue):
            if residue.get_parent().id == peptideligand:
                return 1
            else:
                return 0

    p = PDBParser(QUIET=True)
    s = p.get_structure(pdbname, projectdir + 'pdbs/' +
                        pdbname + '.pdb')  # Disable warnings
    hetflag_done = {}
    for model in s:
        for chain in model:
            for residue in chain:
                hetresname = residue.get_resname()
                # catch residues with hetflag
                hetflag = residue.get_full_id()[3][0].strip()

                hetflag = hetflag.replace("H_", "").strip()
                #hetflag = hetflag.replace("W","")
                #print(hetflag)
                if peptideligand and chain.id==peptideligand:
                    hetflag= 'pep'

                if peptideligand and chain.id!=peptideligand:
                    continue

                if hetflag and hetflag not in ignore_het:
                    if not hetflag in hetflag_done:

                        hetflag_done[hetflag] = 1
                        HETNAM = hetflag

                        temp_path = projectdir + 'results/' + pdbname + \
                            '/ligand/' + HETNAM + '_' + pdbname + '.sdf'

                        ligand_pdb = projectdir + 'results/' + pdbname + \
                            '/ligand/' + HETNAM + '_' + pdbname + '.pdb'
                        ligand_sdf = projectdir + 'results/' + pdbname + \
                            '/ligand/' + HETNAM + '_' + pdbname + '.sdf'
                        ligand_inchi = projectdir + 'results/' + pdbname + \
                            '/ligand/' + HETNAM + '_' + pdbname + '.inchi'
                        ligand_poseview = projectdir + 'results/' + \
                            pdbname + '/png/' + pdbname + '_' + HETNAM + '.png'
                        ligand_png = projectdir + 'results/' + pdbname + '/png/' + HETNAM + '.png'

                        # if sdf not made, make it #Always make them for now
                        if not os.path.isfile(ligand_pdb) or 1 == 1:
                            io = PDBIO()
                            io.set_structure(s)
                            if peptideligand and chain.id==peptideligand:
                                io.save(ligand_pdb, ClassSelect())
                            else:
                                io.save(ligand_pdb, HetSelect())

                            check_unique_ligand_mol(ligand_pdb)

                            if len(list(pybel.readfile("pdb", ligand_pdb))) == 0:
                                continue

                            obConversion = openbabel.OBConversion()
                            obConversion.SetInAndOutFormats("pdb", "inchi")
                            obConversion.SetOptions(
                                "K", obConversion.OUTOPTIONS)
                            mol = openbabel.OBMol()
                            # Open Babel will uncompress automatically
                            obConversion.ReadFile(mol, ligand_pdb)
                            obConversion.WriteFile(mol, ligand_inchi)
                            inchikey = obConversion.WriteString(mol)

                            inchikeys[HETNAM] = inchikey.strip()

                            #smiles[HETNAM] = smile

                            smiles[HETNAM] = pybel.readfile(
                                "pdb", ligand_pdb).next().write("smi").split("\t")[0]

                            mol = pybel.readfile("pdb", ligand_pdb).next()
                            mol.OBMol.AddHydrogens(False, True, 7.4)
                            mol.write("pdb", ligand_pdb, overwrite=True)

                            obConversion = openbabel.OBConversion()
                            obConversion.SetInAndOutFormats("pdb", "sdf")
                            mol = openbabel.OBMol()
                            # Open Babel will uncompress automatically
                            obConversion.ReadFile(mol, ligand_pdb)

                            obConversion.WriteFile(mol, ligand_sdf)

                        # if png of ligand not made, make it
                        if not os.path.isfile(ligand_png):
                            m = Chem.MolFromMolFile(ligand_sdf)
                            # Draw.MolToFile(m,ligand_png)

                        # if interaction png not made, make it #SKIP poseview
                        # stuff
                        if not os.path.isfile(ligand_poseview) and 1 == 2:
                            cmd = "poseview -l " + ligand_sdf + " -p " + projectdir + \
                                "pdbs/" + pdbname + ".pdb -o " + ligand_poseview

                            #print('Running cmd ' + cmd)
                            proc = subprocess.Popen(
                                [cmd], stdout=subprocess.PIPE, shell=True)
                            while proc.poll() is None:
                                time.sleep(1)

                            #(out, err) = proc.communicate()
                        else:
                            # print "Already made
                            # Poseview:",pdbname+"_"+HETNAM+".png"
                            continue
    # print "Done "+str(len(hetflag_done))


def addresiduestoligand(ligand, pdb, residuelist):
    temp_path = projectdir + 'pdbs/' + pdb + '.pdb'
    f_in = open(temp_path, 'r')
    inserstr = ''
    check = []
    # print filename
    ligandid = 0
    chainid = 0
    for line in f_in:
        if line.startswith('ATOM'):
            temp = line.split()
            # need to fix bad PDB formatting where col4 and col5 are put
            # together for some reason -- usually seen when the id is +1000
            m = re.match("(\w)(\d+)", temp[4])
            if (m):
                temp[4] = m.group(1)
                temp[5] = m.group(2)

            aaname = temp[3] + temp[5] + temp[4]

            if aaname in residuelist:
                # print aaname
                inserstr += line
    # print inserstr
    f_in.close()

    # ligands/'+hetflag+'_'+pdbname+".pdb")

    temp_path = projectdir + 'results/' + pdbname + \
        '/ligand/' + ligand + '_' + pdb + '.pdb'
    f_in = open(temp_path, 'r')
    tempstr = ''
    inserted = 0
    for line in f_in:
        if line.startswith('ATOM'):
            temp = line.split()
            if temp[2] == 'H':
                continue  # skip hydrogen in model

        if (line.startswith('CONECT') or line.startswith('MASTER') or line.startswith('END')) and inserted == 0:
            tempstr += inserstr
            inserted = 1
        tempstr += line
    # print tempstr

    # print tempstr
    f_in.close()
    f = open(projectdir + 'results/' + pdbname +
             '/interaction/' + pdb + '_' + ligand + '.pdb', 'w')
    f.write(tempstr)
    f.close()


def get_ring_from_aa(residueid):

    class AAselect(Select):

        def accept_residue(self, residue):
            # print residue.get_full_id()[3][1],residueid
            if str(residue.get_full_id()[3][1]) == residueid:
                return 1
            else:
                return 0
    ptemp = PDBParser(QUIET=True)  # disable warnings
    stemp = ptemp.get_structure(
        pdbname, projectdir + 'pdbs/' + pdbname + '.pdb')
    temp_aa_id = residueid

    io = PDBIO()
    io.set_structure(stemp)
    io.save(projectdir + 'temp/' + residueid + '.pdb', AAselect())

    mol = pybel.readfile("pdb", projectdir + 'temp/' +
                         residueid + '.pdb').next()
    # print hetflag
    rings = getattr(mol, "OBMol").GetSSSR()
    ringlist = []
    for ring in rings:
        center = Vector(0.0, 0.0, 0.0)
        members = ring.Size()
        if ring.IsAromatic():
            atomlist = []
            atomnames = []
            atomvectors = []
            for atom in mol:
                if ring.IsMember(atom.OBAtom):
                    a_vector = Vector(getattr(atom, 'coords'))
                    center += a_vector
                    atomlist.append(atom.idx)
                    atomvectors.append(a_vector)
                    atomnames.append(getattr(atom, 'type'))
            center = center / members
            normal = center - a_vector  # vector in plane
            normal1 = center - atomvectors[0]
            normal2 = center - atomvectors[2]
            normal = Vector(np.cross([normal1[0],normal1[1],normal1[2]],[normal2[0],normal2[1],normal2[2]]))
            ringlist.append([atomlist, center, normal, atomnames, atomvectors])
    return ringlist


def get_hydrogen_from_aa(residueid):

    class AAselect(Select):

        def accept_residue(self, residue):
            # print residue.get_full_id()[3][1],residueid
            if str(residue.get_full_id()[3][1]) == residueid:
                return 1
            else:
                return 0
    ptemp = PDBParser(QUIET=True)
    stemp = ptemp.get_structure(
        pdbname, projectdir + 'pdbs/' + pdbname + '.pdb')
    temp_aa_id = residueid

    io = PDBIO()
    io.set_structure(stemp)
    io.save(projectdir + 'temp/' + residueid + '.pdb', AAselect())

    mol = pybel.readfile("pdb", projectdir + 'temp/' +
                         residueid + '.pdb').next()

    mol.OBMol.AddHydrogens(False, True, 7.4)
    # print hetflag
    donors = []
    for atom in mol:
        if getattr(atom, 'OBAtom').IsHbondDonor():
            chargevector = Vector(getattr(atom, 'coords'))
            # print getattr(atom,'type')," is Donor",chargevector
            temphatoms = []
            for neighbor in pybel.ob.OBAtomAtomIter(atom.OBAtom):
                neighbor = pybel.Atom(neighbor)
                if getattr(neighbor, 'type') == "H":
                    # print "neighbor
                    # Atom",getattr(neighbor,'type'),"Coords:",getattr(neighbor,'coords')
                    temphatoms.append(Vector(getattr(neighbor, 'coords')))

            donors.append([getattr(atom, 'type'), chargevector, temphatoms,getattr(atom, 'OBAtom').IsHbondAcceptor()])

        if getattr(atom, 'OBAtom').IsHbondAcceptor():
            chargevector = Vector(getattr(atom, 'coords'))
            #print getattr(atom, 'type'),chargevector,'acceptor!'

    return donors


def build_ligand_info():
    count_atom_ligand = {}
    p = PDBParser(QUIET=True)
    s = p.get_structure(pdbname, projectdir + 'pdbs/' + pdbname + '.pdb')
    for model in s:
        for chain in model:
            for residue in chain:
                hetresname = residue.get_resname()
                # catch residues with hetflag
                hetflag = residue.get_full_id()[3][0].strip()

                hetflag = hetflag.replace("H_", "").strip()
                #hetflag = hetflag.replace("W","")

                if peptideligand and chain.id==peptideligand:
                    hetflag= 'pep'
                if peptideligand and chain.id!=peptideligand:
                    continue

                if hetflag and hetflag not in ignore_het:
                    # if goodhet!='' and hetflag!=goodhet and
                    # "H_"+goodhet!=hetflag: continue ### Only look at the
                    # ligand that has an image from poseview made for it.
                    if  hetflag not in hetlist or (peptideligand and chain.id==peptideligand):

                        if len(list(pybel.readfile("pdb", projectdir + 'results/' + pdbname + '/ligand/' + hetflag + '_' + pdbname + '.pdb'))) == 0:
                            # This ligand has no molecules
                            # print('no info for',hetflag)
                            continue

                        if hetflag not in hetlist: #do not recreate for peptides
                            hetlist[hetflag] = []
                            ligand_charged[hetflag] = []
                            ligand_donors[hetflag] = []
                            ligand_acceptors[hetflag] = []
                            count_atom_ligand[hetflag] = 0

                            mol = pybel.readfile(
                                "pdb", projectdir + 'results/' + pdbname + '/ligand/' + hetflag + '_' + pdbname + ".pdb").next()
                            # print "LIGAND",hetflag

                            rings = getattr(mol, "OBMol").GetSSSR()

                            # http://python.zirael.org/e-openbabel4.html
                            ringlist = []
                            for ring in rings:
                                center = Vector(0.0, 0.0, 0.0)
                                members = ring.Size()
                                if ring.IsAromatic():
                                    # print "Found an aromatic ring"
                                    atomlist = []
                                    atomnames = []
                                    vectorlist = []
                                    for atom in mol:
                                        if ring.IsMember(atom.OBAtom):
                                            # print atom.idx,getattr(atom,'type'),
                                            # ring.IsMember( atom.OBAtom)
                                            a_vector = Vector(
                                                getattr(atom, 'coords'))
                                            center += a_vector
                                            atomlist.append(atom.idx)
                                            vectorlist.append(a_vector)
                                            atomnames.append(getattr(atom, 'type'))
                                    center = center / members
                                    normal = center - a_vector  # vector in plane
                                    #print center - vectorlist[0],center - vectorlist[2]
                                    normal1 = center - vectorlist[0]
                                    normal2 = center - vectorlist[2]
                                    normal = Vector(np.cross([normal1[0],normal1[1],normal1[2]],[normal2[0],normal2[1],normal2[2]]))
                                    ringlist.append(
                                        [atomlist, center, normal, atomnames, vectorlist])

                            ligand_rings[hetflag] = ringlist

                            for atom in mol:
                                #print "Atom",getattr(atom,'type'),"Coords:",getattr(atom,'coords'),"FormalCharge:",getattr(atom,'formalcharge'),"PartialCharge",getattr(atom,'partialcharge')
                                if getattr(atom, 'formalcharge') != 0:
                                    chargevector = Vector(getattr(atom, 'coords'))
                                    ligand_charged[hetflag].append(
                                        [getattr(atom, 'type'), chargevector, getattr(atom, 'formalcharge')])
                                if getattr(atom, 'OBAtom').IsCarboxylOxygen():
                                    chargevector = Vector(getattr(atom, 'coords'))
                                    # print getattr(atom,'type')," is
                                    # CarboxylOxygen",chargevector
                                    ligand_charged[hetflag].append(
                                        [getattr(atom, 'type'), chargevector, -1])
                                if getattr(atom, 'OBAtom').IsHbondDonor():
                                    chargevector = Vector(getattr(atom, 'coords'))
                                    # print getattr(atom,'type')," is
                                    # Donor",chargevector
                                    temphatoms = []
                                    for neighbor in pybel.ob.OBAtomAtomIter(atom.OBAtom):
                                        neighbor = pybel.Atom(neighbor)
                                        if getattr(neighbor, 'type') == "H":
                                            # print "neighbor
                                            # Atom",getattr(neighbor,'type'),"Coords:",getattr(neighbor,'coords')
                                            temphatoms.append(
                                                Vector(getattr(neighbor, 'coords')))

                                    ligand_donors[hetflag].append(
                                        [getattr(atom, 'type'), chargevector, temphatoms])

                                if getattr(atom, 'OBAtom').IsHbondAcceptor():
                                    chargevector = Vector(getattr(atom, 'coords'))
                                   # print getattr(atom,'type')," is Acceptor",chargevector
                                    ligand_acceptors[hetflag].append([getattr(atom, 'type'), chargevector])
                                    # ligand_charged[hetflag].append([getattr(atom,'type'),chargevector,-1])

                        # Function to get ligand centers to maybe skip some
                        # residues
                        check = 0
                        center = Vector(0.0, 0.0, 0.0)

                        if peptideligand and chain.id==peptideligand:

                            if hetflag in ligandcenter:
                                center = ligandcenter[hetflag][2]

                            for atom in residue:
                                het_atom = atom.name
                                atom_vector = atom.get_vector()
                                center += atom_vector
                                hetlist[hetflag].append(
                                    [hetresname, het_atom, atom_vector])

                                if not hetflag in ligand_atoms:
                                    # make the ligand_atoms ready
                                    ligand_atoms[hetflag] = []
                                ligand_atoms[hetflag].append(
                                    [count_atom_ligand[hetflag], atom_vector, het_atom])
                                count_atom_ligand[hetflag] += 1
                                ligandcenter[hetflag] = [center, count_atom_ligand[hetflag]]

                        else:

                            for atom in residue:
                                if check == 0 and hetflag in ligand_atoms:
                                    continue  # skip when there are many of same ligand
                                het_atom = atom.name

                                check = 1

                                atom_vector = atom.get_vector()
                                center += atom_vector
                                hetlist[hetflag].append(
                                    [hetresname, het_atom, atom_vector])

                                if not hetflag in ligand_atoms:
                                    # make the ligand_atoms ready
                                    ligand_atoms[hetflag] = []
                                ligand_atoms[hetflag].append(
                                    [count_atom_ligand[hetflag], atom_vector, het_atom])
                                count_atom_ligand[hetflag] += 1


                        center2 = center / count_atom_ligand[hetflag]
                        ligandcenter[hetflag] = [
                            center2, count_atom_ligand[hetflag],center]

def remove_hyd(aa,ligand):
    templist = []
    for res in new_results[ligand]['interactions']:
        #print res[0],res[2],aa
        if res[0]==aa and (res[2]=='HYD' or res[2]=='hyd'):
            continue
        else:
            templist.append(res)
    new_results[ligand]['interactions'] = templist

def check_other_aromatic(aa,ligand,info):
    templist = []
    check = True
    for res in new_results[ligand]['interactions']:
        #print res[0],res[2],aa
        if res[0]==aa and res[4]=='aromatic':
            #if the new aromatic interaction has a center-center distance greater than the old one, keep old.
            if info['Distance']>res[6]['Distance']:
                templist.append(res)
                check = False #Do not add the new one.
            else: #if not, delete the old one, as the new is better.
                check = True #add the new one
                continue
        else:
            templist.append(res)
    new_results[ligand]['interactions'] = templist
    return check

# LOOP OVER RECEPTOR AND FIND INTERACTIONS
def find_interactions():
    global count_calcs, count_skips
    count_atom = 0
    count_skips = 0
    count_calcs = 0
    p = PDBParser(QUIET=True)
    s = p.get_structure(pdbname, projectdir + 'pdbs/' + pdbname + '.pdb')
    for model in s:
        for chain in model:
            chainid = chain.get_id()

            if peptideligand and chainid==peptideligand:
                continue
            for residue in chain:
                aa_resname = residue.get_resname()
                aa_seqid = str(residue.get_full_id()[3][1])
                hetflagtest = str(residue.get_full_id()[3][0]).strip()
                aaname = aa_resname + aa_seqid + chainid

                hetflagtest = hetflagtest.replace("H_", "")
                #hetflagtest = hetflagtest.replace("W","")

                if hetflagtest:
                    continue  # residue is a hetnam
                if hetflagtest in hetlist:
                    continue  # residue is a hetnam
                # print "Looking at ",aa_resname,aa_seqid,chainid
                countresidue = count_atom
                # print aaname
                # could probably make a check here to see if this residue was
                # anywhere near the ligand, otherwise skip the check per atom
                for hetflag, atomlist in hetlist.iteritems():
                    if not 'CA' in residue:  # prevent errors
                        continue

                    ca = residue['CA'].get_vector()
                    if (ca - ligandcenter[hetflag][0]).norm() > ligandcenter[hetflag][1]:
                        # print "skipping"
                        count_skips += 1
                        continue

                    count_atom = countresidue
                    sum = 0
                    hydrophobic_count = 0
                    accesible_check = 0

                    # if goodhet!='' and hetflag!=goodhet and
                    # "H_"+goodhet!=hetflag: continue ### Only look at the
                    # ligand that has an image from poseview made for it.
                    tempdistance = radius

                    for atom in atomlist:
                        #print(hetflag,atom)
                        hetresname = atom[0]
                        het_atom = atom[1]
                        het_vector = atom[2]
                        hydrophobic_check = 1

                        aaatomlist = []
                        for atom in residue:
                            count_atom += 1
                            aa_vector = atom.get_vector()
                            aa_atom = atom.name
                            aa_atom_type = atom.element
                            aaatomlist.append([count_atom, aa_vector, aa_atom])

                            d = (het_vector - aa_vector)
                            count_calcs += 1
                            if d.norm() < radius:
                                if not hetflag in results:
                                    results[hetflag] = {}
                                    summary_results[hetflag] = {'score': [], 'hbond': [], 'hbondplus': [],
                                                                'hbond_confirmed': [], 'aromatic': [],'aromaticff': [],
                                                                'ionaromatic': [], 'aromaticion': [], 'aromaticef': [],
                                                                'aromaticfe': [], 'hydrophobic': [], 'waals': [], 'accessible':[]}
                                    new_results[hetflag] = {'interactions':[]}
                                if not aaname in results[hetflag]:
                                    results[hetflag][aaname] = []
                                if not (het_atom[0] == 'H' or aa_atom[0] == 'H' or aa_atom_type=='H'):
                                    #print(aa_atom_type)
                                    results[hetflag][aaname].append([het_atom, aa_atom, round(
                                        d.norm(), 2), het_vector, aa_vector, aa_seqid, chainid])
                                    tempdistance = round(d.norm(), 2)
                                    sum += 1
                            # if both are carbon then we are making a hydrophic
                            # interaction
                            if het_atom[0] == 'C' and aa_atom[0] == 'C' and d.norm() < hydrophob_radius and hydrophobic_check:
                                hydrophobic_count += 1
                                hydrophobic_check = 0

                            if d.norm() < 5 and (aa_atom!='C' and aa_atom!='O' and aa_atom!='N'):
                                #print(aa_atom)
                                accesible_check = 1

                    if accesible_check: #if accessible!
                        summary_results[hetflag]['accessible'].append(
                            [aaname])

                        fragment_file = fragment_library(hetflag, None, '',
                                         aa_seqid, chainid, 'access')

                        new_results[hetflag]['interactions'].append([aaname,fragment_file,'acc','accessible','hidden',''])

                    if hydrophobic_count > 2 and AA[aaname[0:3]] in HYDROPHOBIC_AA:  # min 3 c-c interactions
                        summary_results[hetflag]['hydrophobic'].append(
                            [aaname, hydrophobic_count])

                        fragment_file = fragment_library(hetflag, None, '',
                                         aa_seqid, chainid, 'hydrop')

                        new_results[hetflag]['interactions'].append([aaname,fragment_file,'hyd','hydrophobic','hydrophobic',''])


                    if sum > 1 and aa_resname in AROMATIC:
                        # if debug:
                            # , get_ring_atoms(aaatomlist)
                            # print "Need to analyse aromatic ring in ", aaname
                        aarings = get_ring_from_aa(aa_seqid)
                        if not aarings:
                            # print "Could not find aromatic ring in",aaname
                            continue
                        #print "amount of rings in AA",len(aarings)
                        for aaring in aarings:
                            #aaring = aaring[0]  # res_ring
                            center = aaring[1]
                            count = 0
                            #print "AARING",aaring
                            for ring in ligand_rings[hetflag]:
                                # print ring

                                shortest_center_het_ring_to_res_atom = 10
                                shortest_center_aa_ring_to_het_atom = 10
                                # print aaring[4]
                                # print ring[4]
                                for a in aaring[4]:
                                    if (ring[1] - a).norm() < shortest_center_het_ring_to_res_atom:
                                        shortest_center_het_ring_to_res_atom = (ring[1] - a).norm()

                                for a in ring[4]:
                                    if (center - a).norm() < shortest_center_aa_ring_to_het_atom:
                                        shortest_center_aa_ring_to_het_atom = (center - a).norm()

                                count += 1
                                # take vector from two centers, and compare against
                                # vector from center to outer point -- this will
                                # give the perpendicular angel.
                                angle = Vector.angle(center - ring[1], ring[2]) #aacenter to ring center vs ring normal
                                # take vector from two centers, and compare against
                                # vector from center to outer point -- this will
                                # give the perpendicular angel.
                                angle2 = Vector.angle(center - ring[1], aaring[2]) #aacenter to ring center vs AA normal

                                angle3 = Vector.angle(ring[2], aaring[2]) #two normal vectors against eachother
                                #print "angleaa",aaring[2],"anglelig",ring[2]
                                angle_degrees = [
                                    round(degrees(angle), 1), round(degrees(angle2), 1), round(degrees(angle3), 1)]
                                distance = (center - ring[1]).norm()
                                #if debug:
                                    #print aaname,"Ring #", count, "Distance:", round(distance, 2), "Angle:", angle_degrees, 'Shortest res->ligcenter', shortest_center_het_ring_to_res_atom, 'Shortest lig->rescenter', shortest_center_aa_ring_to_het_atom
                                if distance < 5 and (angle_degrees[2]<20 or abs(angle_degrees[2]-180)<20):  # poseview uses <5
                                    # print "Ring
                                    # #",count,"Distance:",round(distance,2),
                                    # "Angle:",round(angle_degrees,2)
                                    summary_results[hetflag]['aromatic'].append(
                                        [aaname, count, round(distance, 2), angle_degrees])

                                    fragment_file = fragment_library_aromatic(
                                        hetflag, ring[4], aa_seqid, chainid, count)

                                    if debug:
                                        print aaname,"F2F Ring #", count, "Distance:", round(distance, 2), "Angle:", angle_degrees, 'Shortest res->ligcenter', round(shortest_center_het_ring_to_res_atom,2), 'Shortest lig->rescenter', round(shortest_center_aa_ring_to_het_atom,2)
                                    if check_other_aromatic(aaname,hetflag,{'Distance':round(distance, 2),'Angles':angle_degrees}):
                                        new_results[hetflag]['interactions'].append([aaname,fragment_file,'aro_ff','aromatic (face-to-face)','aromatic','none',{'Distance':round(distance, 2),'ResAtom to center':round(shortest_center_het_ring_to_res_atom,2),'LigAtom to center': round(shortest_center_aa_ring_to_het_atom,2),'Angles':angle_degrees}])
                                        remove_hyd(aaname,hetflag)

                                # need to be careful for edge-edge
                                elif (shortest_center_aa_ring_to_het_atom < 4.5) and abs(angle_degrees[0]-90)<30 and abs(angle_degrees[2]-90)<30:
                                    summary_results[hetflag]['aromaticfe'].append(
                                        [aaname, count, round(distance, 2), angle_degrees])

                                    fragment_file = fragment_library_aromatic(
                                        hetflag, ring[4], aa_seqid, chainid, count)

                                    if debug:
                                        print aaname,"FE Ring #", count, "Distance:", round(distance, 2), "Angle:", angle_degrees, 'Shortest res->ligcenter', round(shortest_center_het_ring_to_res_atom,2), 'Shortest lig->rescenter', round(shortest_center_aa_ring_to_het_atom,2)
                                    if check_other_aromatic(aaname,hetflag,{'Distance':round(distance, 2),'Angles':angle_degrees}):
                                        new_results[hetflag]['interactions'].append([aaname,fragment_file,'aro_fe_protein','aromatic (face-to-edge)','aromatic','protein',{'Distance':round(distance, 2),'ResAtom to center':round(shortest_center_het_ring_to_res_atom,2),'LigAtom to center': round(shortest_center_aa_ring_to_het_atom,2),'Angles':angle_degrees}])
                                        remove_hyd(aaname,hetflag)
                                # need to be careful for edge-edge
                                elif (shortest_center_het_ring_to_res_atom < 4.5) and abs(angle_degrees[1]-90)<30 and abs(angle_degrees[2]-90)<30:
                                    summary_results[hetflag]['aromaticef'].append(
                                        [aaname, count, round(distance, 2), angle_degrees])

                                    fragment_file = fragment_library_aromatic(
                                        hetflag, ring[4], aa_seqid, chainid, count)


                                    if debug:
                                        print aaname,"EF Ring #", count, "Distance:", round(distance, 2), "Angle:", angle_degrees, 'Shortest res->ligcenter', round(shortest_center_het_ring_to_res_atom,2), 'Shortest lig->rescenter', round(shortest_center_aa_ring_to_het_atom,2)
                                    if check_other_aromatic(aaname,hetflag,{'Distance':round(distance, 2),'Angles':angle_degrees}):
                                        new_results[hetflag]['interactions'].append([aaname,fragment_file,'aro_ef_protein','aromatic (edge-to-face)','aromatic','protein',{'Distance':round(distance, 2),'ResAtom to center':round(shortest_center_het_ring_to_res_atom,2),'LigAtom to center': round(shortest_center_aa_ring_to_het_atom,2),'Angles':angle_degrees}])
                                        remove_hyd(aaname,hetflag)
                            for charged in ligand_charged[hetflag]:
                                distance = (center - charged[1]).norm()
                                # needs max 4.2 distance to make aromatic+
                                if distance < 4.2 and charged[2] > 0:
                                    if debug:
                                        print "Ring #", count, "Distance:", round(distance, 2), "Angle:", round(angle_degrees, 2)
                                    summary_results[hetflag]['aromaticion'].append(
                                        [aaname, count, round(distance, 2), charged])

                                    #FIXME fragment file
                                    new_results[hetflag]['interactions'].append([aaname,'','aro_ion_protein','aromatic (pi-cation)','aromatic','protein',{'Distance':round(distance, 2)}])
                                    remove_hyd(aaname,hetflag)

                    if sum > 2 and aa_resname in CHARGEDAA and ligand_rings[hetflag]:
                        # print "check for charged AA to aromatic
                        # rings!",aa_resname,hetflag

                        for atom in residue:
                            aa_vector = atom.get_vector()
                            aa_atom = atom.name
                            for ring in ligand_rings[hetflag]:
                                d = (ring[2] - aa_vector).norm()
                                # if d<10: print
                                # "aa_atom",aa_atom,aaname,"distance to a
                                # ring",d,hetflag,aa_resname


def analyze_interactions():
    for ligand, result in results.iteritems():

        # print "AA close to ligands ("+ligand+"): ",list(result.keys())
        # print "Results for"+ligand
        sortedresults = []
        ligscore = 0
        for residue, interaction in result.iteritems():
            sum = 0
            score = 0
            hbond = []
            hbondplus = []
            type = 'waals'
            for entry in interaction:
                hbondconfirmed = []
                if entry[2] <= 3.5:

                    # print(entry)
                    # if debug:
                    #     print "Likely H-Bond", entry

                    if entry[0][0] == 'C' or entry[1][0] == 'C':
                        continue  # If either atom is C then no hydrogen bonding

                    # if entry[1] == 'N': #if residue atom is N, then it is backbone!
                    #     print('backbone interaction!')

                    aa_donors = get_hydrogen_from_aa(entry[5])
                    hydrogenmatch = 0
                    res_is_acceptor = False
                    res_is_donor = False
                    for donor in aa_donors:
                        d = (donor[1] - entry[4]).norm()
                        if d < 0.5:
                            #print 'found donor in residue',residue,entry,donor
                            hydrogens = donor[2]
                            res_is_acceptor = donor[3]
                            res_is_donor = True
                            for hydrogen in hydrogens:
                                hydrogenvector = hydrogen - donor[1]
                                bindingvector = entry[3] - hydrogen
                                angle = round(degrees(Vector.angle(
                                    hydrogenvector, bindingvector)), 2)
                                distance = round(bindingvector.norm(), 2)
                                # print "RESDONOR",residue,"From
                                # ligand",entry[0],"To
                                # AA",entry[1],"HydrogenCheck
                                # angle",angle,"Distance from hydrogen to
                                # acceptor",distance
                                if distance > 2.5:
                                    # print "Too far away"
                                    continue
                                if angle > 60:
                                    # print "Bad angle"
                                    continue
                                hydrogenmatch = 1
                                hbondconfirmed.append(
                                    ["D", entry[0], entry[1], angle, distance])

                    # print "aadonors:",aa_donors

                    found_donor = 0
                    for donor in ligand_donors[ligand]:
                        d = (donor[1] - entry[3]).norm()
                        # print charged,d,residue,entry
                        if d < 0.5:
                            found_donor = 1
                            hydrogens = donor[2]
                            for hydrogen in hydrogens:
                                hydrogenvector = hydrogen - donor[1]
                                bindingvector = entry[4] - hydrogen
                                angle = round(degrees(Vector.angle(
                                    hydrogenvector, bindingvector)), 2)
                                distance = round(bindingvector.norm(), 2)
                                # print "LIGDONOR",residue,"From
                                # ligand",entry[0],"To
                                # AA",entry[1],"HydrogenCheck
                                # angle",angle,"Distance from hydrogen to
                                # acceptor",distance
                                if distance > 2.5:
                                    # print "Too far away"
                                    continue
                                if angle > 60:
                                    # print "Bad angle"
                                    continue
                                hydrogenmatch = 1
                                hbondconfirmed.append(
                                    ["A", entry[0], entry[1], angle, distance])

                    found_acceptor = 0
                    for acceptor in ligand_acceptors[ligand]:
                        d = (acceptor[1] - entry[3]).norm()
                            # print charged,d,residue,entry
                        if d < 0.5:
                            found_acceptor = 1
                            if found_donor==0 and res_is_donor:
                                hydrogenmatch = 1
                                hbondconfirmed.append(['D']) #set residue as donor
                                #print 'found acceptor which is not donor',residue,entry[0],acceptor

                    if not found_acceptor and found_donor and res_is_acceptor:
                        hydrogenmatch = 1
                        hbondconfirmed.append(['A']) #set residue as acceptor
                        #print 'donor which is not acceptor',residue,entry[0]

                    if found_acceptor and found_donor:
                        if res_is_donor and not res_is_acceptor:
                            hydrogenmatch = 1
                            hbondconfirmed.append(['D'])
                        elif not res_is_donor and res_is_acceptor:
                            hydrogenmatch = 1
                            hbondconfirmed.append(['A'])
                        else:
                            pass
                        #print 'can be both donor and acceptor'

                    chargedcheck = 0
                    charge_value = 0
                    res_charge_value = 0
                    doublechargecheck = 0
                    for charged in ligand_charged[ligand]:
                        d = (charged[1] - entry[3]).norm()
                        if d < 0.5:
                            # print 'found charge',residue,d,entry
                            chargedcheck = 1
                            hydrogenmatch = 0  # Replace previous match!
                            charge_value = charged[2]


                    if residue[0:3] in CHARGEDAA:
                        # print "check for hbondplus!",residue,entry
                        # Need to check which atoms, but for now assume charged
                        if chargedcheck:
                            doublechargecheck = 1
                        chargedcheck = 1
                        hydrogenmatch = 0  # Replace previous match!

                        if AA[residue[0:3]] in POSITIVE:
                            res_charge_value = 1
                        elif AA[residue[0:3]] in NEGATIVE:
                            res_charge_value = -1


                    if entry[1] == 'N': #backbone connection!
                        fragment_file = fragment_library(ligand, entry[3], entry[
                                         0], entry[5], entry[6], 'HB_backbone')
                        new_results[ligand]['interactions'].append([residue,fragment_file,'polar_backbone','polar (hydrogen bond with backbone)','polar','protein',entry[0],entry[1],entry[2]])
                        remove_hyd(residue,ligand)
                    elif entry[1] == 'O': #backbone connection!
                        fragment_file = fragment_library(ligand, entry[3], entry[
                                         0], entry[5], entry[6], 'HB_backbone')
                        new_results[ligand]['interactions'].append([residue,fragment_file,'polar_backbone','polar (hydrogen bond with backbone)','polar','protein',entry[0],entry[1],entry[2]])
                        remove_hyd(residue,ligand)
                    elif hydrogenmatch:
                        found = 0

                        fragment_file = fragment_library(ligand, entry[3], entry[
                                         0], entry[5], entry[6], 'HB')

                        for x in summary_results[ligand]['hbond_confirmed']:
                            if residue == x[0]:
                                # print "Already key there",residue
                                key = summary_results[ligand][
                                    'hbond_confirmed'].index(x)
                                summary_results[ligand]['hbond_confirmed'][
                                    key][1].extend(hbondconfirmed)
                                found = 1

                        if hbondconfirmed[0][0]=="D":
                            new_results[ligand]['interactions'].append([residue,fragment_file,'polar_donor_protein','polar (hydrogen bond)','polar','protein',entry[0],entry[1],entry[2]])
                            remove_hyd(residue,ligand)
                        if hbondconfirmed[0][0]=="A":
                            new_results[ligand]['interactions'].append([residue,fragment_file,'polar_acceptor_protein','polar (hydrogen bond)','polar','protein',entry[0],entry[1],entry[2]])
                            remove_hyd(residue,ligand)

                        if found == 0:
                            summary_results[ligand]['hbond_confirmed'].append(
                                [residue, hbondconfirmed])
                        if chargedcheck:
                            type = 'hbondplus'
                            hbondplus.append(entry)



                    elif chargedcheck:
                        type = 'hbondplus'
                        hbondplus.append(entry)
                        fragment_file = fragment_library(ligand, entry[3], entry[
                                         0], entry[5], entry[6], 'HBC')

                        remove_hyd(residue,ligand)
                        if doublechargecheck:
                            if (res_charge_value>0):
                                new_results[ligand]['interactions'].append([residue,fragment_file,'polar_double_pos_protein','polar (charge-charge)','polar','',entry[0],entry[1],entry[2]])
                            elif (res_charge_value<0):
                                new_results[ligand]['interactions'].append([residue,fragment_file,'polar_double_neg_protein','polar (charge-charge)','polar','',entry[0],entry[1],entry[2]])
                        elif (charge_value>0):
                            new_results[ligand]['interactions'].append([residue,fragment_file,'polar_pos_ligand','polar (charge-assisted hydrogen bond)','polar','ligand',entry[0],entry[1],entry[2]])
                        elif (charge_value<0):
                            new_results[ligand]['interactions'].append([residue,fragment_file,'polar_neg_ligand','polar (charge-assisted hydrogen bond)','polar','ligand',entry[0],entry[1],entry[2]])
                        else:
                            if (res_charge_value>0):
                                new_results[ligand]['interactions'].append([residue,fragment_file,'polar_pos_protein','polar (charge-assisted hydrogen bond)','polar','protein',entry[0],entry[1],entry[2]])
                            elif (res_charge_value<0):
                                new_results[ligand]['interactions'].append([residue,fragment_file,'polar_neg_protein','polar (charge-assisted hydrogen bond)','polar','protein',entry[0],entry[1],entry[2]])
                            else:
                                new_results[ligand]['interactions'].append([residue,fragment_file,'polar_unknown_protein','polar (charge-assisted hydrogen bond)','polar','protein',entry[0],entry[1],entry[2]])

                    else:
                        type = 'hbond'
                        hbond.append(entry)
                        fragment_file = fragment_library(ligand, entry[3], entry[
                                         0], entry[5], entry[6], 'HB')
                        new_results[ligand]['interactions'].append([residue,fragment_file,'polar_unspecified','polar (hydrogen bond)','polar','',entry[0],entry[1],entry[2]])
                        remove_hyd(residue,ligand)
                    #print type,hbondconfirmed
                    entry[3] = ''

                if (entry[2] < 4.5):
                    sum += 1
                    score += 4.5 - entry[2]
            score = round(score, 2)

            if type == 'waals' and score > 2:  # mainly no hbond detected
                summary_results[ligand]['waals'].append([residue, score, sum])
            elif type == 'hbond':
                summary_results[ligand]['hbond'].append(
                    [residue, score, sum, hbond])
            elif type == 'hbondplus':
                summary_results[ligand]['hbondplus'].append(
                    [residue, score, sum, hbondplus])
            # elif type == 'hbond_confirmed':
            #     summary_results[ligand]['hbond_confirmed'].append([residue,score,sum,hbondconfirmed])

            ligscore += score

            # print "Total <4 (score is combined diff from
            # 4)",sum,"score",score
            sortedresults.append([residue, score, sum, hbond, type])

        summary_results[ligand]['score'].append([ligscore])
        summary_results[ligand]['inchikey'] = inchikeys[ligand]
        summary_results[ligand]['smiles'] = smiles[ligand]
        new_results[ligand]['score'] = ligscore
        new_results[ligand]['inchikey'] = inchikeys[ligand]
        new_results[ligand]['smiles'] = smiles[ligand]
        if ligand in hetlist_display:
            summary_results[ligand]['prettyname'] = hetlist_display[ligand]
            new_results[ligand]['prettyname'] = hetlist_display[ligand]

        # print ligand,"Ligand score:"+str(ligscore)

        sortedresults = sorted(sortedresults, key=itemgetter(1), reverse=True)


def pretty_results():
    for ligand, result in summary_results.iteritems():
        output = ''
        bindingresidues = []
        #output += "Results for "+str(ligand)+"\n"
        for type, typelist in result.iteritems():
            if type == 'waals':
                continue
            output += type + "\n"
            if type == 'waals':
                typelist = sorted(typelist, key=itemgetter(2), reverse=True)
            if type == 'hydrophobic':
                typelist = sorted(typelist, key=itemgetter(1), reverse=True)
            for entry in typelist:
                if type != 'score':
                    bindingresidues.append(entry[0])
                if type == 'hbond':
                    output += '\t'.join(map(str, entry[0:1])) + '\n'
                    for bond in entry[3]:
                        output += '\t'.join(map(str, bond[0:3])) + '\n'
                elif type == 'hbondplus':
                    output += '\t'.join(map(str, entry[0:1])) + '\n'
                    for bond in entry[3]:
                        output += '\t'.join(map(str, bond[0:3])) + '\n'
                elif type == 'hbond_confirmed':
                    output += '\t'.join(map(str, entry[0:1])) + '\n'
                    for bond in entry[1]:
                        output += '\t'.join(map(str, bond)) + '\n'
                else:
                    # print entry
                    output += '\t'.join(map(str, entry)) + '\n'

        temp_path = projectdir + 'results/' + pdbname + '/output/' + \
            pdbname + '_' + ligand.replace("H_", "") + '.yaml'

        # yaml.dump(result, open(temp_path, 'w'))
        yaml.dump(new_results[ligand], open(temp_path, 'w'))
        if debug:
            print ligand,'\n',open(temp_path,'r').read()

        addresiduestoligand(ligand, pdbname, bindingresidues)


def calculate_interactions(pdb, session=None, peptide=None):
    global pdbname, hetlist, hetlist_display, ligand_atoms, ligand_charged, ligandcenter, ligand_rings, ligand_donors, ligand_acceptors, results, sortedresults, summary_results, inchikeys, smiles, projectdir, new_results, peptideligand

    hetlist = {}
    hetlist_display = {}
    ligand_atoms = {}
    ligand_charged = {}
    ligandcenter = {}
    ligand_rings = {}
    ligand_donors = {}
    ligand_acceptors = {}
    results = {}
    sortedresults = {}
    summary_results = {}
    new_results = {}
    inchikeys = {}
    smiles = {}
    peptideligand = peptide
    if not session:
        pdbname = pdb
        # print "checking normal ",pdbname
        check_pdb()
        checkdirs()
        hetlist_display = find_ligand_full_names()
        create_ligands_and_poseview()
        build_ligand_info()
        find_interactions()
        analyze_interactions()
        pretty_results()
    else:
        pdbname = pdb
        projectdir = '/tmp/interactions/' + session + "/"
        checkdirs()
        hetlist_display = find_ligand_full_names()
        create_ligands_and_poseview()
        build_ligand_info()
        find_interactions()
        analyze_interactions()
        pretty_results()


def main(argv):
    pdbname = ''
    try:
        # print 'ARGV      :', argv
        opts, args = getopt.getopt(argv, "p:s:c:", ["pdb"])
    except getopt.GetoptError as err:
        print "Remember PDB name -p "
        print err
        sys.exit(2)

    session = None
    peptide = None
    for opt, arg in opts:
        if opt in ("-p"):
            pdbname = arg
        elif opt in ("-s"):
            session = arg
        elif opt in ("-c"):
            peptide = arg

    if not pdbname:
        print "Remember PDB name -p "
        sys.exit(2)

    if session:
        calculate_interactions(pdbname, session, peptide=peptide)
    else:
        calculate_interactions(pdbname, peptide=peptide)


if __name__ == "__main__":
    main(sys.argv[1:])
    #pdbname = '1F88'
    # calculate_interactions(pdbname)
