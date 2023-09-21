from django.shortcuts import render, redirect
from django.http import HttpResponse, HttpResponseRedirect
from django.conf import settings
from django.db.models import Count, Sum, Avg, Q
from django.utils.text import slugify
from django.conf import settings

from interaction.models import ResidueFragmentInteraction, StructureLigandInteraction, ResidueFragmentInteractionType
from interaction.forms import PDBform
from ligand.models import Ligand, LigandType, LigandRole
from structure.models import Structure, PdbData, Rotamer, Fragment, StructureModel, StructureComplexModel, StructureExtraProteins, StructureVectors, StructureModelRMSD, StructureModelpLDDT, StructureAFScores
from structure.assign_generic_numbers_gpcr import GenericNumbering
from protein.models import Protein, ProteinSegment
from residue.models import Residue, ResidueGenericNumberEquivalent, ResidueNumberingScheme
from common.models import WebResource, WebLink
from common.tools import fetch_from_web_api
from common.diagrams_gpcr import DrawHelixBox, DrawSnakePlot
from common.selection import Selection, SelectionItem
from common import definitions
from common.views import AbsTargetSelection
from contactnetwork.models import Interaction

import os
import yaml
from operator import itemgetter
from datetime import datetime
import re
import json
import logging
from subprocess import call
import urllib
import collections
from collections import OrderedDict, defaultdict
from io import BytesIO
from Bio.PDB import PDBIO, PDBParser, Select, Vector
import xlsxwriter

#NEW IMPORTS
import shutil
import requests
import numpy as np
from math import degrees, atan2, cos, sin, pi
from rdkit import Chem
from rdkit import RDConfig
from rdkit.Chem import AllChem
from rdkit.Chem import MolFromPDBFile
from rdkit.Chem import ChemicalFeatures
#END NEW IMPORTS


AA = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
      'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
      'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
      'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
      'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}

#IMPLEMENTING THE LEGACY FUNCTIONS HERE!

# Van Der Waals interactions: to define if two atoms are interacting via Van Der Waals
# you can calculate the distance between these two atoms , which should be
# around 3 to 6 Angstrom, then check if their Van Der Waals radius (PeriodicTable)
# is overlapping, is enough to match the distance

HBD = {'H', 'K', 'N', 'Q', 'R', 'S', 'T', 'W', 'Y'}
HBA = {'D', 'E', 'H', 'N', 'Q', 'S', 'T', 'Y'}
NEGATIVE = {'D', 'E'}
POSITIVE = {'H', 'K', 'R'}
cation_atoms =  ['NZ', 'CZ', 'NE', 'NH1', 'NH2']
AROMATIC = {'TYR', 'TRP', 'PHE', 'HIS'}
CHARGEDAA = {'ARG', 'LYS', 'ASP', 'GLU'}
HYDROPHOBIC_AA = {'A', 'C', 'F', 'I', 'L', 'M', 'P', 'V', 'W', 'Y'}
ignore_het = ['NA', 'W']  # ignore sodium and water
radius = 5
hydrophob_radius = 4.5
pdb_dir = os.sep.join([settings.DATA_DIR, 'structure_data', 'pdbs'])

#RETURN THE DICTIONARY RESULTS
def runcalculation_2022(pdbname, peptide="", file_input=False):
    output = calculate_interactions(pdbname, None, peptide, file_input)
    return output

#RETURN THE DICTIONARY RESULTS
def runusercalculation_2022(filename, session):
    output = calculate_interactions(filename, session, None)
    return output

def calculate_interactions(pdb, session=None, peptide=None, file_input=False):
    # REMEMBER TO GET THE RETURNS FROM ALL THE BELOW FUNCTIONS
    hetlist = {}
    ligand_atoms = {}
    ligand_charged = {}
    ligandcenter = {}
    ligand_rings = {}
    ligand_donors = {}
    ligand_acceptors = {}
    results = {}
    sortedresults = []
    summary_results = {}
    new_results = {}
    projectdir = '/tmp/interactions/'
    tempdir = projectdir + 'temp/'
    if not os.path.exists(tempdir):
        os.makedirs(tempdir)
        # os.chmod(tempdir, 0o777)
    if not file_input:
        pdb_location = projectdir + 'pdbs/' + pdb + '.pdb'
    else:
        #/protwis/data/protwis/gpcr/af_arman/fpr2_human-6242-rank0
        complex_name = pdb.split('/')[-1].split('-rank')[0]
        model_name = pdb.split('/')[-1]
        pdb_location = projectdir + 'pdbs/' + complex_name + '/' + model_name + '.pdb'
    if not session:
        # print('Checking PDB')
        check_pdb(projectdir, pdb, file_input)
        # print('Checking Dirs')
        checkdirs(projectdir, pdb, file_input)
        # print('Finding interacting ligand')
        hetlist_display = find_interacting_ligand(pdb_location, pdb, file_input)
        # Defining a shared parser
        parser = PDBParser(QUIET=True)
        if file_input:
            pdb = complex_name
        scroller = parser.get_structure(pdb, pdb_location)
        # print('Creating ligand and poseview')
        create_ligands_and_poseview(hetlist_display, scroller, projectdir, pdb, peptide) #ignore_het (should be global), inchikeys, smiles (should not be used)
        # print('Building ligand info')
        hetlist, ligand_charged, ligand_donors, ligand_atoms, ligand_acceptors, ligandcenter, ligand_rings = build_ligand_info(
                                                                                                                scroller, hetlist_display,
                                                                                                                projectdir, pdb, peptide, hetlist,
                                                                                                                ligand_atoms, ligand_charged, ligand_donors,
                                                                                                                ligand_acceptors, ligandcenter, ligand_rings)
        # print('Finding interactions')
        summary_results, new_results, results = find_interactions(
                                                    scroller, projectdir, pdb, peptide,
                                                    hetlist, ligandcenter, radius, summary_results,
                                                    new_results, results, hydrophob_radius, ligand_rings, ligand_charged, pdb_location)
        # print('Analyzing interactions')
        summary_results, new_results, sortedresults = analyze_interactions(
                                                        projectdir, pdb, results, ligand_donors,
                                                        ligand_acceptors, ligand_charged, new_results,
                                                        summary_results, hetlist_display, sortedresults, pdb_location)
        # print('Making pretty results')
        pretty_results(projectdir, pdb, summary_results, pdb_location)

    else:
        projectdir = projectdir + session + "/"
        checkdirs(projectdir, pdb, file_input)
        hetlist_display = find_interacting_ligand(pdb_location, pdb, file_input)
        # Defining a shared parser
        parser = PDBParser(QUIET=True)
        scroller = parser.get_structure(pdb, pdb_location)
        create_ligands_and_poseview(hetlist_display, scroller, projectdir, pdb, peptide) #ignore_het (should be global), inchikeys, smiles (should not be used)
        hetlist, ligand_charged, ligand_donors, ligand_atoms, ligand_acceptors, ligandcenter, ligand_rings = build_ligand_info(
                                                                                                                scroller, hetlist_display,
                                                                                                                projectdir, pdb, peptide, hetlist,
                                                                                                                ligand_atoms, ligand_charged, ligand_donors,
                                                                                                                ligand_acceptors, ligandcenter, ligand_rings)
        summary_results, new_results, results = find_interactions(
                                                    scroller, projectdir, pdb,
                                                    peptide, hetlist, ligandcenter,
                                                    radius, summary_results, new_results,
                                                    results, hydrophob_radius, ligand_rings, ligand_charged, pdb_location)
        summary_results, new_results, sortedresults = analyze_interactions(
                                                        projectdir, pdb, results, ligand_donors,
                                                        ligand_acceptors, ligand_charged, new_results,
                                                        summary_results, hetlist_display, sortedresults, pdb_location)
        pretty_results(projectdir, pdb, summary_results, pdb_location)
    return new_results


def check_pdb(projectdir, pdb, file_input):  #CAN WE HAVE THE PDB AS A VAR AND NOT A FILE?
    # check if PDB is there, otherwise fetch

    if not os.path.exists(projectdir + 'pdbs/'):
        os.makedirs(projectdir + 'pdbs/')
    if not file_input:
        if not os.path.isfile(projectdir + 'pdbs/' + pdb + '.pdb'):
            url = 'https://www.rcsb.org/pdb/files/%s.pdb' % pdb
            # pdbfile = urllib.request.urlopen(url).read()
            pdbfile = requests.get(url)
            if ("404 Not Found" in pdbfile.text or pdb in ['7F1T', '7XBX']) and pdb+'.pdb' in os.listdir(pdb_dir):
                with open(os.sep.join([pdb_dir, pdb+'.pdb']), 'r') as f:
                    pdbfile = f.read()
            else:
                pdbfile = pdbfile.text

            # output_pdb = pdbfile.decode('utf-8').split('\n')
            temp_path = projectdir + 'pdbs/' + pdb + '.pdb'
            with open(temp_path, "w") as f:
                f.write(pdbfile)
        else:
            with open(projectdir + 'pdbs/' + pdb + '.pdb', 'r') as f:
                pdbfile = f.read()
            if ("404 Not Found" in pdbfile or pdb in ['7F1T', '7XBX']) and pdb+'.pdb' in os.listdir(pdb_dir):
                with open(os.sep.join([pdb_dir, pdb+'.pdb']), 'r') as f:
                    pdbfile = f.read()
                temp_path = projectdir + 'pdbs/' + pdb + '.pdb'
                with open(temp_path, "w") as f:
                    f.write(pdbfile)
    else:
        #/protwis/data/protwis/gpcr/af_arman/npy4r_rat-1521-rank0
        complex_name = pdb.split('/')[-1].split('-rank')[0]
        model_name = pdb.split('/')[-1]
        with open(pdb +'.pdb', 'r') as f:
            pdbfile = f.read()
        # output_pdb = pdbfile.decode('utf-8').split('\n')
        temp_path = projectdir + 'pdbs/' + complex_name + '/' + model_name + '.pdb'
        if not os.path.exists(projectdir + 'pdbs/' + complex_name +'/'):
            os.makedirs(projectdir + 'pdbs/' + complex_name + '/')
        with open(temp_path, "w") as f:
            f.write(pdbfile)
    # return output_pdb

def checkdirs(projectdir, pdb, file_input): # DO WE NEED TO HAVE THIS DATA STORE IN TMP FILES?
    # check that dirs are there and have right permissions
    if not file_input:
        directory = projectdir + 'results/' + pdb
        if os.path.exists(directory):
            shutil.rmtree(directory)
        directory = projectdir + 'results/' + pdb + '/interaction'
        if not os.path.exists(directory):
            os.makedirs(directory)
        directory = projectdir + 'results/' + pdb + '/ligand' #not really necessary
        if not os.path.exists(directory):
            os.makedirs(directory)
        directory = projectdir + 'results/' + pdb + '/output'
        if not os.path.exists(directory):
            os.makedirs(directory)
        directory = projectdir + 'results/' + pdb + '/fragments' #not really necessary
        if not os.path.exists(directory):
            os.makedirs(directory)
        directory = projectdir + 'temp/'
        if not os.path.exists(directory):
            os.makedirs(directory)
    else:
        complex_name = pdb.split('/')[-1].split('-rank')[0]
        directory = projectdir + 'results/' + complex_name
        if os.path.exists(directory):
            shutil.rmtree(directory)
        directory = projectdir + 'results/' + complex_name + '/interaction'
        if not os.path.exists(directory):
            os.makedirs(directory)
        directory = projectdir + 'results/' + complex_name + '/ligand' #not really necessary
        if not os.path.exists(directory):
            os.makedirs(directory)
        directory = projectdir + 'results/' + complex_name + '/output'
        if not os.path.exists(directory):
            os.makedirs(directory)
        directory = projectdir + 'results/' + complex_name + '/fragments' #not really necessary
        if not os.path.exists(directory):
            os.makedirs(directory)
        directory = projectdir + 'temp/'
        if not os.path.exists(directory):
            os.makedirs(directory)


def find_interacting_ligand(pdb_location, pdb, file_input):
    #Compare these names to the ones in the database
    if not file_input:
        db_ligs = list(StructureLigandInteraction.objects.filter(structure_id__pdb_code_id__index=pdb.upper()).values_list('pdb_reference', flat=True))
    else:
        receptor =  pdb.split('/')[-1].split('-r')[0].replace('-','_')
        # lig_id = pdb.split('/')[-1].split('-')[1]
        code = '_'.join(['AFM', receptor]).upper()
        db_ligs = list(StructureLigandInteraction.objects.filter(structure_id__pdb_code_id__index=code).values_list('pdb_reference', flat=True))
    f_in = open(pdb_location, 'r')
    d = {}
    for lig in db_ligs:
        d[lig] = ''
        if lig == 'pep':
            continue
        else:
            for line in f_in:
                if lig in line:
                    if line.startswith('HETSYN'):
                        try:
                            m = re.match("HETSYN[\s]+([\w]{3})[\s]+(.+)", line)
                            d[m.group(1)] = m.group(2).strip()
                        except AttributeError:
                            #apply exception to ligand with multiline names1
                            m = re.match("HETSYN[\s]+([\w]{1})[\s]+(.+)", line)
                            code = m.group(2).split(' ')[0]
                            d[code] += m.group(2).split(' ')[1].strip()
                    else:
                        d[lig] = ''
                # if m.group(1) in db_ligs:
                #     d[m.group(1)] = m.group(2).strip()
    return d

def accept_residue(residue, hetflag, peptide=None):
    if residue.get_parent().id == peptide:
        return 1
    elif residue.get_resname().strip() == hetflag:
        return 1
    else:
        return 0

def create_ligands_and_poseview(ligand_het, scroller, projectdir, pdb, peptide=None):

    class HetSelect(Select):
        @staticmethod
        def accept_residue(residue):
            if residue.get_resname().strip() == hetflag:
                return 1
            else:
                return 0

    class ClassSelect(Select):
        @staticmethod
        def accept_residue(residue):
            if residue.get_parent().id == peptide:
                return 1
            else:
                return 0

    for model in scroller:
        for chain in model:
            for residue in chain:
                # catch residues with hetflag
                hetflag = residue.get_full_id()[3][0].strip()
                hetflag = hetflag.replace("H_", "").strip()
                if peptide and chain.id==peptide:
                    hetflag= 'pep'

                if hetflag in ligand_het.keys():
                    ligand_pdb = projectdir + 'results/' + pdb + '/ligand/' + hetflag + '_' + pdb + '.pdb'

                    # if sdf not made, make it #Always make them for now
                    if not os.path.isfile(ligand_pdb):
                        io = PDBIO()
                        io.set_structure(scroller)
                        if peptide and chain.id==peptide:
                            io.save(ligand_pdb, ClassSelect())
                        else:
                            io.save(ligand_pdb, HetSelect())
                        # check_unique_ligand_mol(ligand_pdb)
                        if MolFromPDBFile(ligand_pdb) == 0:
                            continue
                    else:
                        continue

def check_unique_ligand_mol(filename): # IS THIS NEEDED?
    # check that only HETATM are exported to file
    f_in = open(filename, 'r')
    tempstr = ''
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

def get_sdf_ligand_from_cache(comp_id):
    #CHECK IF THE SDF IS CACHED
    url = 'https://files.rcsb.org/ligands/download/$index'
    cache_dir = ["pdbe", 'sdf_models']
    comp_id += '_model.sdf'
    data = fetch_from_web_api(url, comp_id, cache_dir, raw=True)
    mol = AllChem.MolFromMolBlock(data)
    # AllChem.AssignStereochemistryFrom3D(mol) #FOR PYTHON 3
    return mol

def isRingAromatic(mol, bondRing):
    for ring_id in bondRing:
        if not mol.GetBondWithIdx(ring_id).GetIsAromatic():
            return False
    return True

def build_ligand_info(scroller, lig_het, projectdir, pdb, peptide, hetlist, ligand_atoms, ligand_charged, ligand_donors, ligand_acceptors, ligandcenter, ligand_rings):
    count_atom_ligand = {}

    for model in scroller:
        for chain in model:
            for residue in chain:
                hetresname = residue.get_resname()
                # catch residues with hetflag
                hetflag = residue.get_full_id()[3][0].strip()
                hetflag = hetflag.replace("H_", "").strip()

                if peptide and chain.id==peptide:
                    hetflag= 'pep'
                if peptide and chain.id!=peptide:
                    continue

                # REMEMBER TO PARSE ONLY THE ACTUAL LIGAND
                if hetflag in lig_het.keys():
                    if (hetflag not in hetlist) or (chain.id==peptide):
                        if MolFromPDBFile(projectdir + 'results/' + pdb + '/ligand/' + hetflag + '_' + pdb + '.pdb') == 0:
                            # This ligand has no molecules
                            # print('no info for',hetflag)
                            continue

                        if hetflag not in hetlist: #do not recreate for peptides
                            hetlist[hetflag] = []
                            ligand_charged[hetflag] = []
                            ligand_donors[hetflag] = []
                            ligand_acceptors[hetflag] = []
                            count_atom_ligand[hetflag] = 0
                            mol2 = MolFromPDBFile(projectdir + 'results/' + pdb + '/ligand/' + hetflag + '_' + pdb + ".pdb")
                            if not mol2:
                                mol2 = MolFromPDBFile(projectdir + 'results/' + pdb + '/ligand/' + hetflag + '_' + pdb + ".pdb", sanitize=False)
                            hetflag_sdf = get_sdf_ligand_from_cache(hetflag)
                            try:
                                mol2 = AllChem.AssignBondOrdersFromTemplate(refmol=hetflag_sdf, mol=mol2)
                            except ValueError:
                                try:
                                    smiles = list(Ligand.objects.filter(pdbe=hetflag).values_list('smiles', flat=True))[0]
                                    refmol = AllChem.MolFromSmiles(smiles)
                                    mol2 = AllChem.AssignBondOrdersFromTemplate(refmol=refmol, mol=mol2)
                                except:
                                    pass
                            mol2 = Chem.AddHs(mol2)
                            rings = Chem.rdmolops.GetSSSR(mol2)
                            ringlist = []
                            if len(rings) > 0:
                                counter = -1
                                ri = mol2.GetRingInfo()
                                for ring in ri.BondRings():
                                    counter +=1
                                    center = Vector(0.0, 0.0, 0.0)
                                    members = len(ring)
                                    if isRingAromatic(mol2, ring):
                                        atomlist = []
                                        vectorlist = []
                                        for i, atom in enumerate(mol2.GetAtoms()):
                                            if atom.GetIdx() in ri.AtomRings()[counter]:
                                                positions = mol2.GetConformer().GetAtomPosition(i)
                                                a_vector = Vector(positions)
                                                center += a_vector
                                                atomlist.append(atom.GetIdx())
                                                vectorlist.append(a_vector)
                                        center = center / members
                                        normal1 = center - vectorlist[0]
                                        normal2 = center - vectorlist[2]
                                        # normal = Vector(np.cross(normal1,normal2
                                        normal = Vector(np.cross([normal1[0],normal1[1],normal1[2]],[normal2[0],normal2[1],normal2[2]]))
                                        ringlist.append([atomlist, center, normal, vectorlist])

                            ligand_rings[hetflag] = ringlist
                            fdefName = os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef')
                            factory = ChemicalFeatures.BuildFeatureFactory(fdefName)
                            feats = factory.GetFeaturesForMol(mol2)
                            for i, atom in enumerate(mol2.GetAtoms()):
                                if atom.GetFormalCharge() != 0:
                                    positions = mol2.GetConformer().GetAtomPosition(i)
                                    chargevector = Vector(positions)
                                    ligand_charged[hetflag].append([chargevector, atom.GetFormalCharge()])
                            for feat in feats:
                                if feat.GetFamily() == 'Acceptor':
                                    positions = feat.GetPos()
                                    chargevector = Vector(positions)
                                    ligand_acceptors[hetflag].append(chargevector)
                                if feat.GetFamily() == 'Donor':
                                    positions = feat.GetPos()
                                    chargevector = Vector(positions)
                                    temphatoms = []
                                    atom_id = feat.GetAtomIds()[0]
                                    atom = mol2.GetAtomWithIdx(atom_id)
                                    for j, neighbour_atom in enumerate(atom.GetNeighbors()):
                                        if neighbour_atom.GetSymbol() == 'H':
                                            positions_hs = mol2.GetConformer().GetAtomPosition(j)
                                            temphatoms.append(Vector(positions_hs))
                                    ligand_donors[hetflag].append([chargevector, temphatoms])

                        # Function to get ligand centers to maybe skip some residues
                        check = 0
                        center = Vector(0.0, 0.0, 0.0)
                        if peptide and chain.id==peptide:
                            if hetflag in ligandcenter:
                                center = ligandcenter[hetflag][2]
                            for atom in residue:
                                het_atom = atom.name
                                atom_vector = atom.get_vector()
                                center += atom_vector
                                hetlist[hetflag].append([hetresname, het_atom, atom_vector])
                                if hetflag not in ligand_atoms:
                                    # make the ligand_atoms ready
                                    ligand_atoms[hetflag] = []
                                ligand_atoms[hetflag].append([count_atom_ligand[hetflag], atom_vector, het_atom])
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
                                hetlist[hetflag].append([hetresname, het_atom, atom_vector])
                                if hetflag not in ligand_atoms:
                                    # make the ligand_atoms ready
                                    ligand_atoms[hetflag] = []
                                ligand_atoms[hetflag].append([count_atom_ligand[hetflag], atom_vector, het_atom])
                                count_atom_ligand[hetflag] += 1
                        center2 = center / count_atom_ligand[hetflag]
                        ligandcenter[hetflag] = [center2, count_atom_ligand[hetflag],center]

    return hetlist, ligand_charged, ligand_donors, ligand_atoms, ligand_acceptors, ligandcenter, ligand_rings

# LOOP OVER RECEPTOR AND FIND INTERACTIONS
def find_interactions(scroller, projectdir, pdb, peptide, hetlist, ligandcenter, radius, summary_results, new_results, results, hydrophob_radius, ligand_rings, ligand_charged, pdb_location):
    count_atom = 0
    count_skips = 0
    count_calcs = 0
    # waals = []
    for heteroatom in ligandcenter.keys():
        for model in scroller:
            for chain in model:
                chainid = chain.get_id()
                if peptide and chainid==peptide:
                    continue
                for residue in chain:
                    aa_resname = residue.get_resname()
                    aa_seqid = str(residue.get_full_id()[3][1])
                    aaname = aa_resname + aa_seqid + chainid
                    hetflagtest = str(residue.get_full_id()[3][0]).strip()

                    if hetflagtest or 'CA' not in residue:
                        continue  # residue is a hetnam

                    # could probably make a check here to see if this residue was
                    # anywhere near the ligand, otherwise skip the check per atom
                    ca = residue['CA'].get_vector()
                    if (ca - ligandcenter[heteroatom][0]).norm() > ligandcenter[heteroatom][1]:
                        count_skips += 1
                        continue
                    for hetflag, atomlist in hetlist.items():
                        sum_data = 0
                        hydrophobic_count = 0
                        accesible = False
                        for atom_het in atomlist:
                            het_atom = atom_het[1]
                            het_vector = atom_het[2]
                            aaatomlist = []
                            for atom in residue:
                                count_atom += 1
                                aa_vector = atom.get_vector()
                                aa_atom = atom.name
                                aa_atom_type = atom.element
                                aaatomlist.append([count_atom, aa_vector, aa_atom])
                                distance = (het_vector - aa_vector)
                                count_calcs += 1
                                # if distance.norm() < 6:
                                #     waals.append([aaname])
                                if distance.norm() < radius:
                                    if hetflag not in results:
                                        results[hetflag] = {}
                                        summary_results[hetflag] = {'score': [], 'hbond': [], 'hbondplus': [], 'pistack': [],
                                                                    'hbond_confirmed': [], 'aromatic': [],'aromaticff': [],
                                                                    'ionaromatic': [], 'aromaticion': [], 'aromaticef': [],
                                                                    'aromaticfe': [], 'hydrophobic': [], 'waals': [], 'accessible':[]}
                                        new_results[hetflag] = {'interactions':[]}
                                    if aaname not in results[hetflag]:
                                        results[hetflag][aaname] = []
                                    if (het_atom[0] != 'H') or (aa_atom[0] != 'H') or (aa_atom_type != 'H'):
                                        tempdistance = round(distance.norm(), 2)
                                        results[hetflag][aaname].append([het_atom, aa_atom, tempdistance, het_vector, aa_vector, aa_seqid, chainid])
                                        sum_data += 1
                                # if both are carbon then we are making a hydrophic interaction
                                if (het_atom[0] == 'C') and (aa_atom[0] == 'C') and (distance.norm() < hydrophob_radius):
                                    hydrophobic_count += 1

                                # If within 5 angstrom and not a backbone atom (name C, O, N), then indicate as a residue in vicinity of the ligand
                                if (distance.norm() < radius) and (aa_atom not in ['C', 'O', 'N']):
                                    accesible = True
                        fragment_file = ''
                        # if hetflag in summary_results.keys():
                        #     summary_results[hetflag]['waals'] = waals
                        #     for amino_acid in waals:
                        #         new_results[hetflag]['interactions'].append([amino_acid[0],fragment_file,'waals','accessible','waals',''])
                        if accesible: #if accessible!)
                            summary_results[hetflag]['accessible'].append([aaname])
                            fragment_file = fragment_library(projectdir, pdb, hetflag, None, '', aa_seqid, chainid, 'access', pdb_location)
                            new_results[hetflag]['interactions'].append([aaname,fragment_file,'acc','accessible','hidden',''])
                        if hydrophobic_count > 2 and AA[aaname[0:3]] in HYDROPHOBIC_AA:  # min 3 c-c interactions
                            summary_results[hetflag]['hydrophobic'].append([aaname, hydrophobic_count])
                            fragment_file = fragment_library(projectdir, pdb, hetflag, None, '', aa_seqid, chainid, 'hydrop', pdb_location)
                            new_results[hetflag]['interactions'].append([aaname,fragment_file,'hyd','hydrophobic','hydrophobic',''])
                        if sum_data > 1 and aa_resname in AROMATIC:
                            aarings = get_ring_from_aa(scroller, projectdir, aa_seqid, residue)
                            #aarings.append([atomlist, center, normal, vectorlist])
                            if not aarings:
                                continue
                            for aaring in aarings:
                                center = aaring[1]
                                count = 0
                                for ring in ligand_rings[hetflag]:
                                    shortest_center_het_ring_to_res_atom = 10
                                    shortest_center_aa_ring_to_het_atom = 10
                                    for a in aaring[3]:
                                        if (ring[1] - a).norm() < shortest_center_het_ring_to_res_atom:
                                            shortest_center_het_ring_to_res_atom = (ring[1] - a).norm()
                                    for a in ring[3]:
                                        if (center - a).norm() < shortest_center_aa_ring_to_het_atom:
                                            shortest_center_aa_ring_to_het_atom = (center - a).norm()
                                    count += 1
                                    # take vector from two centers, and compare against
                                    # vector from center to outer point -- this will
                                    # give the perpendicular angle.
                                    angle = Vector.angle(center - ring[1], ring[2]) #aacenter to ring center vs ring normal
                                    # take vector from two centers, and compare against
                                    # vector from center to outer point -- this will
                                    # give the perpendicular angle.
                                    angle2 = Vector.angle(center - ring[1], aaring[2]) #aacenter to ring center vs AA normal
                                    angle3 = Vector.angle(ring[2], aaring[2]) #two normal vectors against eachother
                                    angle_degrees = [round(degrees(angle), 1), round(degrees(angle2), 1), round(degrees(angle3), 1)]
                                    distance = (center - ring[1]).norm()
                                    if distance < 5 and (angle_degrees[2]<20 or abs(angle_degrees[2]-180)<20):  # poseview uses <5
                                        summary_results[hetflag]['aromatic'].append([aaname, count, round(distance, 2), angle_degrees])
                                        fragment_file = fragment_library_aromatic(projectdir, pdb, hetflag, ring[3], aa_seqid, chainid, count, pdb_location)
                                        if check_other_aromatic(aaname, hetflag, {'Distance':round(distance, 2),'Angles':angle_degrees}, new_results):
                                            new_results[hetflag]['interactions'].append([aaname,fragment_file,
                                                                                         'aro_ff','aromatic (face-to-face)','aromatic','none',
                                                                                         {'Distance':round(distance, 2),'ResAtom to center':round(shortest_center_het_ring_to_res_atom,2),
                                                                                         'LigAtom to center': round(shortest_center_aa_ring_to_het_atom,2),'Angles':angle_degrees}])
                                            remove_hyd(aaname, hetflag, new_results)
                                    # need to be careful for edge-edge
                                    elif (shortest_center_aa_ring_to_het_atom < 4.5) and abs(angle_degrees[0]-90)<30 and abs(angle_degrees[2]-90)<30:
                                        summary_results[hetflag]['aromaticfe'].append([aaname, count, round(distance, 2), angle_degrees])
                                        fragment_file = fragment_library_aromatic(projectdir, pdb, hetflag, ring[3], aa_seqid, chainid, count, pdb_location)
                                        if check_other_aromatic(aaname, hetflag, {'Distance':round(distance, 2),'Angles':angle_degrees}, new_results):
                                            new_results[hetflag]['interactions'].append([aaname,fragment_file,
                                                                                         'aro_fe_protein','aromatic (face-to-edge)','aromatic','protein',
                                                                                         {'Distance':round(distance, 2),'ResAtom to center':round(shortest_center_het_ring_to_res_atom,2),
                                                                                         'LigAtom to center': round(shortest_center_aa_ring_to_het_atom,2),'Angles':angle_degrees}])
                                            remove_hyd(aaname, hetflag, new_results)
                                    # need to be careful for edge-edge
                                    elif (shortest_center_het_ring_to_res_atom < 4.5) and abs(angle_degrees[1]-90)<30 and abs(angle_degrees[2]-90)<30:
                                        summary_results[hetflag]['aromaticef'].append([aaname, count, round(distance, 2), angle_degrees])
                                        fragment_file = fragment_library_aromatic(projectdir, pdb, hetflag, ring[3], aa_seqid, chainid, count, pdb_location)
                                        if check_other_aromatic(aaname, hetflag, {'Distance':round(distance, 2),'Angles':angle_degrees}, new_results):
                                            new_results[hetflag]['interactions'].append([aaname,fragment_file,
                                                                                         'aro_ef_protein','aromatic (edge-to-face)','aromatic','protein',
                                                                                         {'Distance':round(distance, 2),'ResAtom to center':round(shortest_center_het_ring_to_res_atom,2),
                                                                                         'LigAtom to center': round(shortest_center_aa_ring_to_het_atom,2),'Angles':angle_degrees}])
                                            remove_hyd(aaname, hetflag, new_results)
                                for charged in ligand_charged[hetflag]:
                                    distance = (center - charged[0]).norm()
                                    # needs max 4.2 distance to make aromatic+
                                    if distance < 4.2 and charged[1] > 0:
                                        summary_results[hetflag]['aromaticion'].append([aaname, count, round(distance, 2), charged])
                                        #FIXME fragment file
                                        new_results[hetflag]['interactions'].append([aaname,'','aro_ion_protein','aromatic (pi-cation)','aromatic','protein',{'Distance':round(distance, 2)}])
                                        remove_hyd(aaname, hetflag, new_results)
                        #Calculate PiStack interactions
                        if (aa_resname in ['ARG', 'LYS']) and ligand_rings[hetflag]:
                            for atom in residue:
                                aa_vector = atom.get_vector()
                                aa_atom = atom.name
                                #First: find the receptor atom that checks
                                if aa_atom in cation_atoms:
                                    #Second: calculate the center of the aromatic ring (ligand)
                                    for ring in ligand_rings[hetflag]:
                                        ring_center = ring[1]
                                        #Third: calculate the distance of the residue atom to the center of the ring
                                        distance = (ring_center - aa_vector).norm()
                                        #and check distance is less than 6.6
                                        if distance < 6.6:
                                            #Fourth: calculate the angle between the atom and the perpendicular center
                                            perp_vector = ring[2]
                                            angle = degrees(Vector.angle(perp_vector, aa_vector))
                                            if angle <= 30:
                                                summary_results[hetflag]['pistack'].append([aaname])
                                                new_results[hetflag]['interactions'].append([aaname,'','aro_ion_protein','aromatic (pi-cation)','aromatic','protein',{'Distance':round(distance, 2)}])
    return summary_results, new_results, results

def get_ring_from_aa(scroller, projectdir, residueid, residue):

    class AAselect(Select):
        def accept_residue(self, residue):
            if str(residue.get_full_id()[3][1]) == residueid:
                return 1
            else:
                return 0

    io = PDBIO()
    io.set_structure(scroller)
    io.save(projectdir + 'temp/' + residueid + '.pdb', AAselect())
    mol = MolFromPDBFile(projectdir + 'temp/' + residueid + '.pdb')
    mol = Chem.AddHs(mol)
    # ANALYZING AROMATIC RINGS
    rings = Chem.rdmolops.GetSSSR(mol)
    ringlist = []
    if len(rings) > 0:
        counter = -1
        ri = mol.GetRingInfo()
        for ring in ri.BondRings():
            counter +=1
            center = Vector(0.0, 0.0, 0.0)
            members = len(ring)
            if isRingAromatic(mol, ring):
                atomlist = []
                vectorlist = []
                for i, atom in enumerate(mol.GetAtoms()):
                    if atom.GetIdx() in ri.AtomRings()[counter]:
                        positions = mol.GetConformer().GetAtomPosition(i)
                        a_vector = Vector(positions)
                        center += a_vector
                        atomlist.append(atom.GetIdx())
                        vectorlist.append(a_vector)
                center = center / members
                normal1 = center - vectorlist[0]
                normal2 = center - vectorlist[2]
                normal = Vector(np.cross([normal1[0],normal1[1],normal1[2]],[normal2[0],normal2[1],normal2[2]]))
                ringlist.append([atomlist, center, normal, vectorlist])
    return ringlist

def remove_hyd(aa, ligand, new_results):
    templist = []
    for res in new_results[ligand]['interactions']:
        if (res[0]==aa) and (res[2] in ['HYD','hyd']):
            continue
        else:
            templist.append(res)
    new_results[ligand]['interactions'] = templist

def check_other_aromatic(aa, ligand, info, new_results):
    templist = []
    check = True
    for res in new_results[ligand]['interactions']:
        if (res[0]==aa) and (res[4]=='aromatic'):
            #if the new aromatic interaction has a center-center distance greater than the old one, keep old.
            if info['Distance'] > res[6]['Distance']:
                templist.append(res)
                check = False #Do not add the new one.
            else: #if not, delete the old one, as the new is better.
                check = True #add the new one
                continue
        else:
            templist.append(res)
    new_results[ligand]['interactions'] = templist
    return check

def get_hydrogen_from_aa(projectdir, pdb, residueid, pdb_location):

    class AAselect(Select):

        def accept_residue(self, residue):
            # print residue.get_full_id()[3][1],residueid
            if str(residue.get_full_id()[3][1]) == residueid:
                return 1
            else:
                return 0
    ptemp = PDBParser(QUIET=True)
    stemp = ptemp.get_structure(pdb, pdb_location)

    io = PDBIO()
    io.set_structure(stemp)
    io.save(projectdir + 'temp/' + residueid + '.pdb', AAselect())

    mol = MolFromPDBFile(projectdir + 'temp/' +residueid + '.pdb')
    mol = Chem.AddHs(mol)

    fdefName = os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef')
    factory = ChemicalFeatures.BuildFeatureFactory(fdefName)
    feats = factory.GetFeaturesForMol(mol)
    donors = []
    for feat in feats:
        acceptor = False
        #FINDING H BOND ACCEPTORS
        if feat.GetFamily() == 'Acceptor':
            positions = feat.GetPos()
            chargevector = Vector(positions)
            acceptor = True
            temphatoms = []
            atom_id = feat.GetAtomIds()[0]
            atom = mol.GetAtomWithIdx(atom_id)
            for j, neighbour_atom in enumerate(atom.GetNeighbors()):
                if neighbour_atom.GetSymbol() == 'H':
                    positions_hs = mol.GetConformer().GetAtomPosition(j)
                    temphatoms.append(Vector(positions_hs))
        if feat.GetFamily() == 'Donor':
            positions = feat.GetPos()
            chargevector = Vector(positions)
            temphatoms = []
            atom_id = feat.GetAtomIds()[0]
            atom = mol.GetAtomWithIdx(atom_id)
            for j, neighbour_atom in enumerate(atom.GetNeighbors()):
                if neighbour_atom.GetSymbol() == 'H':
                    positions_hs = mol.GetConformer().GetAtomPosition(j)
                    temphatoms.append(Vector(positions_hs))
        donors.append([chargevector, temphatoms, acceptor])
    return donors

def fragment_library(projectdir, pdb, ligand, atomvector, atomname, residuenr, chain, typeinteraction, pdbfile):
    #if debug:
        #print "Make fragment pdb file for ligand:", ligand, "atom vector", atomvector, "atomname", atomname, "residuenr from protein", residuenr, typeinteraction, 'chain', chain
    residuename = 'unknown'
    ligand_pdb = projectdir + 'results/' + pdb + '/ligand/' + ligand + '_' + pdb + '.pdb'
    mol2 = MolFromPDBFile(ligand_pdb, removeHs=True)
    listofvectors = []
    chain = chain.strip()
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

    filename = projectdir + 'results/' + pdb + '/fragments/' + pdb + "_" + ligand + \
        "_" + residuename + residuenr + chain + "_" + atomname + "_" + typeinteraction + ".pdb"
    f_in.close()
    f = open(filename, 'w')
    f.write(tempstr)
    f.close()
    try:
        mol2 = MolFromPDBFile(filename)
        Chem.MolToPDBFile(mol2, filename)
    except:
        mol2 = MolFromPDBFile(filename, sanitize=False)
        Chem.MolToPDBFile(mol2, filename)

    return filename

def fragment_library_aromatic(projectdir, pdb, ligand, atomvectors, residuenr, chain, ringnr, pdbfile):
    chain = chain.strip()
    # pdbfile = projectdir + 'pdbs/' + pdb + '.pdb'
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
    filename = projectdir + 'results/' + pdb + '/fragments/' + pdb + "_" + ligand + \
        "_" + residuename + str(residuenr) + chain + "_aromatic_" + str(ringnr) + ".pdb"
    f_in.close()
    f = open(filename, 'w')
    f.write(tempstr)
    f.close()
    return filename

def analyze_interactions(projectdir, pdb, results, ligand_donors, ligand_acceptors, ligand_charged, new_results, summary_results, hetlist_display, sortedresults, pdb_location):
    for ligand, result in results.items():
        ligscore = 0
        fragment_file = ''
        for residue, interaction in result.items():
            sum_data = 0
            score = 0
            hbond = []
            hbondplus = []
            interaction_type = 'waals'
            for entry in interaction:
                hbondconfirmed = []
                if (entry[2] <= 3.5):
                    if entry[0][0] == 'C' or entry[1][0] == 'C':
                        continue  # If either atom is C then no hydrogen bonding
                    aa_donors = get_hydrogen_from_aa(projectdir, pdb, entry[5], pdb_location)
                    hydrogenmatch = False
                    res_is_acceptor = False
                    res_is_donor = False
                    found_donor = False
                    found_acceptor = False
                    chargedcheck = False
                    charge_value = 0
                    res_charge_value = False
                    doublechargecheck = False
                    for donor in aa_donors:
                        d = (donor[0] - entry[4]).norm()
                        if d < 0.5:
                            hydrogens = donor[1]
                            res_is_acceptor = donor[2]
                            res_is_donor = True
                            for hydrogen in hydrogens:
                                hydrogenvector = hydrogen - donor[0]
                                bindingvector = entry[3] - hydrogen
                                angle = round(degrees(Vector.angle(hydrogenvector, bindingvector)), 2)
                                distance = round(bindingvector.norm(), 2)
                                if distance > 2.5:
                                    # "Too far away"
                                    continue
                                if angle > 60:
                                    # "Bad angle"
                                    continue
                                hydrogenmatch = True
                                hbondconfirmed.append(["D", entry[0], entry[1], angle, distance])
                    #Checking donor status
                    for donor in ligand_donors[ligand]:
                        d = (donor[0] - entry[3]).norm()
                        if d < 0.5:
                            found_donor = True
                            hydrogens = donor[1]
                            for hydrogen in hydrogens:
                                hydrogenvector = hydrogen - donor[0]
                                bindingvector = entry[4] - hydrogen
                                angle = round(degrees(Vector.angle(hydrogenvector, bindingvector)), 2)
                                distance = round(bindingvector.norm(), 2)
                                if distance > 2.5:
                                    # "Too far away"
                                    continue
                                if angle > 60:
                                    # "Bad angle"
                                    continue
                                hydrogenmatch = True
                                hbondconfirmed.append(["A", entry[0], entry[1], angle, distance])
                    #Checking acceptor status
                    for acceptor in ligand_acceptors[ligand]:
                        d = (acceptor - entry[3]).norm()
                        if d < 0.5:
                            found_acceptor = 1
                            if not found_donor and res_is_donor:
                                hydrogenmatch = True
                                hbondconfirmed.append(['D']) #set residue as donor

                    if not found_acceptor and found_donor and res_is_acceptor:
                        hydrogenmatch = True
                        hbondconfirmed.append(['A']) #set residue as acceptor
                    if found_acceptor and found_donor:
                        if res_is_donor and not res_is_acceptor:
                            hydrogenmatch = True
                            hbondconfirmed.append(['D'])
                        elif not res_is_donor and res_is_acceptor:
                            hydrogenmatch = True
                            hbondconfirmed.append(['A'])
                        else:
                            pass
                    #Checking charged statuses
                    for charged in ligand_charged[ligand]:
                        d = (charged[0] - entry[3]).norm()
                        if d < 0.5:
                            chargedcheck = True
                            hydrogenmatch = False  # Replace previous match!
                            charge_value = charged[1]
                    if residue[0:3] in CHARGEDAA:
                        if chargedcheck:
                            doublechargecheck = True
                        chargedcheck = True
                        hydrogenmatch = False  # Replace previous match!

                        if AA[residue[0:3]] in POSITIVE:
                            res_charge_value = 1
                        elif AA[residue[0:3]] in NEGATIVE:
                            res_charge_value = -1
                    if entry[1] == 'N': #backbone connection!
                        fragment_file = fragment_library(projectdir, pdb, ligand, entry[3], entry[0], entry[5], entry[6], 'HB_backbone', pdb_location)
                        new_results[ligand]['interactions'].append([residue,fragment_file,
                                                                    'polar_backbone','polar (hydrogen bond with backbone)',
                                                                    'polar','protein',entry[0],entry[1],entry[2]])
                        remove_hyd(residue, ligand, new_results)
                    elif entry[1] == 'O': #backbone connection!
                        fragment_file = fragment_library(projectdir, pdb, ligand, entry[3], entry[0], entry[5], entry[6], 'HB_backbone', pdb_location)
                        new_results[ligand]['interactions'].append([residue,fragment_file,
                                                                    'polar_backbone','polar (hydrogen bond with backbone)',
                                                                    'polar','protein',entry[0],entry[1],entry[2]])
                        remove_hyd(residue, ligand, new_results)
                    elif hydrogenmatch:
                        found = False
                        fragment_file = fragment_library(projectdir, pdb, ligand, entry[3], entry[0], entry[5], entry[6], 'HB', pdb_location)
                        for x in summary_results[ligand]['hbond_confirmed']:
                            if residue == x[0]:
                                # print "Already key there",residue
                                key = summary_results[ligand]['hbond_confirmed'].index(x)
                                summary_results[ligand]['hbond_confirmed'][key][1].extend(hbondconfirmed)
                                found = True
                        if hbondconfirmed[0][0]=="D":
                            new_results[ligand]['interactions'].append([residue,fragment_file,
                                                                        'polar_donor_protein','polar (hydrogen bond)',
                                                                        'polar','protein',entry[0],entry[1],entry[2]])
                            remove_hyd(residue, ligand, new_results)
                        if hbondconfirmed[0][0]=="A":
                            new_results[ligand]['interactions'].append([residue,fragment_file,
                                                                        'polar_acceptor_protein','polar (hydrogen bond)',
                                                                        'polar','protein',entry[0],entry[1],entry[2]])
                            remove_hyd(residue, ligand, new_results)
                        if not found:
                            summary_results[ligand]['hbond_confirmed'].append([residue, hbondconfirmed])
                        if chargedcheck:
                            interaction_type = 'hbondplus'
                            hbondplus.append(entry)
                    elif chargedcheck:
                        interaction_type = 'hbondplus'
                        hbondplus.append(entry)
                        fragment_file = fragment_library(projectdir, pdb, ligand, entry[3], entry[0], entry[5], entry[6], 'HBC', pdb_location)
                        remove_hyd(residue, ligand, new_results)
                        if doublechargecheck:
                            if (res_charge_value>0):
                                new_results[ligand]['interactions'].append([residue, fragment_file,
                                                                            'polar_double_pos_protein','polar (charge-charge)',
                                                                            'polar','',entry[0],entry[1],entry[2]])
                            elif (res_charge_value<0):
                                new_results[ligand]['interactions'].append([residue, fragment_file,
                                                                            'polar_double_neg_protein','polar (charge-charge)',
                                                                            'polar','',entry[0],entry[1],entry[2]])
                        elif (charge_value>0):
                            new_results[ligand]['interactions'].append([residue, fragment_file,
                                                                        'polar_pos_ligand','polar (charge-assisted hydrogen bond)',
                                                                        'polar','ligand',entry[0],entry[1],entry[2]])
                        elif (charge_value<0):
                            new_results[ligand]['interactions'].append([residue, fragment_file,
                                                                        'polar_neg_ligand','polar (charge-assisted hydrogen bond)',
                                                                        'polar','ligand',entry[0],entry[1],entry[2]])
                        else:
                            if (res_charge_value>0):
                                new_results[ligand]['interactions'].append([residue, fragment_file,
                                                                            'polar_pos_protein','polar (charge-assisted hydrogen bond)',
                                                                            'polar','protein',entry[0],entry[1],entry[2]])
                            elif (res_charge_value<0):
                                new_results[ligand]['interactions'].append([residue, fragment_file,
                                                                            'polar_neg_protein','polar (charge-assisted hydrogen bond)',
                                                                            'polar','protein',entry[0],entry[1],entry[2]])
                            else:
                                new_results[ligand]['interactions'].append([residue, fragment_file,
                                                                            'polar_unknown_protein','polar (charge-assisted hydrogen bond)',
                                                                            'polar','protein',entry[0],entry[1],entry[2]])
                    else:
                        interaction_type = 'hbond'
                        hbond.append(entry)
                        fragment_file = fragment_library(projectdir, pdb, ligand, entry[3], entry[0], entry[5], entry[6], 'HB', pdb_location)
                        new_results[ligand]['interactions'].append([residue, fragment_file,
                                                                    'polar_unspecified','polar (hydrogen bond)',
                                                                    'polar','',entry[0],entry[1],entry[2]])
                        remove_hyd(residue, ligand, new_results)
                    #print type,hbondconfirmed
                    entry[3] = ''

                if (entry[2] < 4.5):
                    sum_data += 1
                    score += 4.5 - entry[2]

            score = round(score, 2)
            if interaction_type == 'waals' and score > 1:  # mainly no hbond detected
                summary_results[ligand]['waals'].append([residue, score, sum_data])
                new_results[ligand]['interactions'].append([residue, fragment_file,
                                                            'Van der Waals', 'Van der Waals', 'waals',
                                                            '', entry[0], entry[1], entry[2]])
            elif interaction_type == 'hbond':
                summary_results[ligand]['hbond'].append([residue, score, sum_data, hbond])
            elif interaction_type == 'hbondplus':
                summary_results[ligand]['hbondplus'].append([residue, score, sum_data, hbondplus])
            # elif interaction_type == 'hbond_confirmed':
            #     summary_results[ligand]['hbond_confirmed'].append([residue,score,sum,hbondconfirmed])
            ligscore += score
            # print "Total <4 (score is combined diff from 4)",sum,"score",score
            sortedresults.append([residue, score, sum, hbond, type])

        summary_results[ligand]['score'].append([ligscore])
        new_results[ligand]['score'] = ligscore
        if ligand in hetlist_display:
            summary_results[ligand]['prettyname'] = hetlist_display[ligand]
            new_results[ligand]['prettyname'] = hetlist_display[ligand]

        # print ligand,"Ligand score:"+str(ligscore)

        sortedresults = sorted(sortedresults, key=itemgetter(1), reverse=True)

    return summary_results, new_results, sortedresults

def addresiduestoligand(projectdir, ligand, pdb, residuelist, pdbfile):
    # temp_path = projectdir + 'pdbs/' + pdb + '.pdb'
    f_in = open(pdbfile, 'r')
    inserstr = ''
    for line in f_in:
        if line.startswith('ATOM'):
            temp = line.split()
            m = re.match("(\w)(\d+)", temp[4])
            if (m):
                temp[4] = m.group(1)
                temp[5] = m.group(2)
            aaname = temp[3] + temp[5] + temp[4]
            if aaname in residuelist:
                inserstr += line
    f_in.close()

    temp_path = projectdir + 'results/' + pdb + '/ligand/' + ligand + '_' + pdb + '.pdb'
    f_in = open(temp_path, 'r')
    tempstr = ''
    inserted = 0
    for line in f_in:
        if line.startswith('ATOM'):
            temp = line.split()
            if temp[2] == 'H':
                continue

        if (line.startswith('CONECT') or line.startswith('MASTER') or line.startswith('END')) and inserted == 0:
            tempstr += inserstr
            inserted = 1
        tempstr += line
    f_in.close()

    f = open(projectdir + 'results/' + pdb + '/interaction/' + pdb + '_' + ligand + '.pdb', 'w')
    f.write(tempstr)
    f.close()

def pretty_results(projectdir, pdb, summary_results, pdbfile):
    for ligand, result in summary_results.items():
        bindingresidues = []
        for interaction_type, typelist in result.items():
            if interaction_type == 'waals':
                typelist = sorted(typelist, key=itemgetter(2), reverse=True)
            if interaction_type == 'hydrophobic':
                typelist = sorted(typelist, key=itemgetter(1), reverse=True)
            for entry in typelist:
                if interaction_type not in ['score', 'prettyname']:
                    bindingresidues.append(entry[0])
        bindingresidues = list(set(bindingresidues))
        addresiduestoligand(projectdir, ligand, pdb, bindingresidues, pdbfile)

####### END IMPLEMENTATION OF LEGACY FUNCTIONS ####

def regexaa(aa):
    aaPattern = re.compile(r'^(\w{3})(\d+)([\w\s]+)$')
    # Splits the string into AA number CHAIN : LEU339A => ('LEU', '339', 'A')
    result = aaPattern.search(aa)
    if result:
        result = result.groups()
        aa = AA[result[0]]
        number = result[1]
        chain = result[2]
        return aa, number, chain
    else:
        aaPattern = re.compile(r'^(\w{3})(\d+)$')
        result = aaPattern.search(aa)
        if result:
            result = result.groups()
            aa = AA[result[0]]
            number = result[1]
            chain = ' '
            return aa, number, chain
        else:
            return None, None, None


class InteractionSelection(AbsTargetSelection):

    # Left panel
    step = 1
    number_of_steps = 1
    docs = 'generic_numbering.html'  # FIXME

    # description = 'Select receptors to index by searching or browsing in the middle column. You can select entire' \
    #     + ' receptor families and/or individual receptors.\n\nSelected receptors will appear in the right column,' \
    #     + ' where you can edit the list.\n\nSelect which numbering schemes to use in the middle column.\n\nOnce you' \
    #     + ' have selected all your receptors, click the green button.'

    description = 'Select the structure of interest by using the dropdown in the middle. The selection if viewed to the right and the interactions will be loaded immediately.'

    # Middle section
    numbering_schemes = False
    filters = False
    search = False
    title = "Select a structure based on PDB-code"

    template_name = 'interaction/interactionselection.html'

    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', False),
    ])

    # Buttons
    buttons = {
        'continue': {
            'label': 'Show interactions',
            'onclick': 'submitupload()',
            'color': 'success',
        }
    }

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)

        context['structures'] = ResidueFragmentInteraction.objects.values('structure_ligand_pair__structure__pdb_code__index', 'structure_ligand_pair__structure__protein_conformation__protein__parent__entry_name').annotate(
            num_ligands=Count('structure_ligand_pair', distinct=True), num_interactions=Count('pk', distinct=True)).order_by('structure_ligand_pair__structure__pdb_code__index')
        context['structure_groups'] = sorted(set([ structure['structure_ligand_pair__structure__pdb_code__index'][0] for structure in context['structures'] ]))
        context['form'] = PDBform()
        return context


def StructureDetails(request, pdbname):
    """
    Show structure details
    """
    pdbname = pdbname
    structures = ResidueFragmentInteraction.objects.values('structure_ligand_pair__ligand__name', 'structure_ligand_pair__pdb_reference', 'structure_ligand_pair__annotated').filter(
        structure_ligand_pair__structure__pdb_code__index=pdbname).annotate(numRes=Count('pk', distinct=True)).order_by('-numRes')
    resn_list = ''


    main_ligand = []
    for structure in structures:
        if structure['structure_ligand_pair__annotated']:
            resn_list += ",\"" + structure['structure_ligand_pair__pdb_reference'] + "\""
            main_ligand.append(structure['structure_ligand_pair__pdb_reference'])

    crystal = Structure.objects.get(pdb_code__index=pdbname)
    if pdbname.startswith('AFM'):
        p = Protein.objects.get(id=crystal.protein_conformation.protein.id)
    else:
        p = Protein.objects.get(protein=crystal.protein_conformation.protein)

    residuelist = Residue.objects.filter(protein_conformation__protein=p).prefetch_related('protein_segment','display_generic_number','generic_number')
    lookup = {}

    residues_lookup = {}
    for r in residuelist:
        if r.generic_number:
            lookup[r.generic_number.label] = r.sequence_number
            residues_lookup[r.sequence_number] = r.amino_acid +str(r.sequence_number)+ " "+ r.generic_number.label

    residues = ResidueFragmentInteraction.objects.filter(
        structure_ligand_pair__structure__pdb_code__index=pdbname, structure_ligand_pair__annotated=True).exclude(interaction_type__type ='hidden').order_by('rotamer__residue__sequence_number')
    residues_browser = []
    ligands = []
    display_res = []
    main_ligand_full = []
    residue_table_list = []
    for residue in residues:
        key = residue.interaction_type.name
        aa = residue.rotamer.residue.amino_acid
        pos = residue.rotamer.residue.sequence_number
        wt_pos = -1

        if residue.rotamer.residue.generic_number and residue.rotamer.residue.generic_number.label in lookup:
            residue_table_list.append(residue.rotamer.residue.generic_number.label)
            wt_pos = lookup[residue.rotamer.residue.generic_number.label]

        if residue.rotamer.residue.protein_segment:
            segment = residue.rotamer.residue.protein_segment.slug
        else:
            segment = ''
        if residue.rotamer.residue.display_generic_number:
            display = residue.rotamer.residue.display_generic_number.label
        else:
            display = ''
        ligand = residue.structure_ligand_pair.ligand.name
        display_res.append(str(pos))
        residues_browser.append({'type': key, 'aa': aa, 'ligand': ligand,
                                 'pos': pos, 'wt_pos': wt_pos, 'gpcrdb': display, 'segment': segment})

        if pos not in residues_lookup:
            residues_lookup[pos] = aa + str(pos) + " " +display + " interaction " + key
        else:
            residues_lookup[pos] += " interaction " + key

        if ligand not in ligands:
            ligands.append(ligand)
            main_ligand_full.append(ligand)
    display_res = ' or '.join(display_res)
    # RESIDUE TABLE
    segments = ProteinSegment.objects.all().filter().prefetch_related()

    proteins = [p]

    numbering_schemes_selection = [settings.DEFAULT_NUMBERING_SCHEME]
    numbering_schemes_selection.append(p.residue_numbering_scheme.slug)

    numbering_schemes = ResidueNumberingScheme.objects.filter(
        slug__in=numbering_schemes_selection).all()
    default_scheme = numbering_schemes.get(
        slug=settings.DEFAULT_NUMBERING_SCHEME)
    data = OrderedDict()

    for segment in segments:
        data[segment.slug] = OrderedDict()
        residues = Residue.objects.filter(protein_segment=segment,  protein_conformation__protein=p,
                                          generic_number__label__in=residue_table_list).prefetch_related('protein_conformation__protein',
                                                                                                         'protein_conformation__state', 'protein_segment',
                                                                                                         'generic_number', 'display_generic_number', 'generic_number__scheme',
                                                                                                         'alternative_generic_numbers__scheme')
        for scheme in numbering_schemes:
            if scheme == default_scheme and scheme.slug == settings.DEFAULT_NUMBERING_SCHEME:
                for pos in list(set([x.generic_number.label for x in residues if x.protein_segment == segment])):
                    data[segment.slug][pos] = {
                        scheme.slug: pos, 'seq': ['-'] * len(proteins)}
            elif scheme == default_scheme:
                for pos in list(set([x.generic_number.label for x in residues if x.protein_segment == segment])):
                    data[segment.slug][pos] = {
                        scheme.slug: pos, 'seq': ['-'] * len(proteins)}

        for residue in residues:
            alternatives = residue.alternative_generic_numbers.all()
            pos = residue.generic_number
            for alternative in alternatives:
                if alternative.scheme not in numbering_schemes:
                    continue
                scheme = alternative.scheme
                if default_scheme.slug == settings.DEFAULT_NUMBERING_SCHEME:
                    pos = residue.generic_number
                    if scheme == pos.scheme:
                        data[segment.slug][pos.label]['seq'][proteins.index(
                            residue.protein_conformation.protein)] = str(residue)
                    else:
                        if scheme.slug not in data[segment.slug][pos.label].keys():
                            data[segment.slug][pos.label][
                                scheme.slug] = alternative.label
                        if alternative.label not in data[segment.slug][pos.label][scheme.slug]:
                            data[segment.slug][pos.label][
                                scheme.slug] += " " + alternative.label
                        data[segment.slug][pos.label]['seq'][proteins.index(
                            residue.protein_conformation.protein)] = str(residue)
                else:
                    if scheme.slug not in data[segment.slug][pos.label].keys():
                        data[segment.slug][pos.label][
                            scheme.slug] = alternative.label
                    if alternative.label not in data[segment.slug][pos.label][scheme.slug]:
                        data[segment.slug][pos.label][
                            scheme.slug] += " " + alternative.label
                    data[segment.slug][pos.label]['seq'][proteins.index(
                        residue.protein_conformation.protein)] = str(residue)

    # Preparing the dictionary of list of lists. Dealing with tripple nested
    # dictionary in django templates is a nightmare
    flattened_data = OrderedDict.fromkeys([x.slug for x in segments], [])
    for s in iter(flattened_data):
        flattened_data[s] = [[data[s][x][y.slug]
                              for y in numbering_schemes] + data[s][x]['seq'] for x in sorted(data[s])]

    context = {}
    context['header'] = zip([x.short_name for x in numbering_schemes] + [x.name for x in proteins], [x.name for x in numbering_schemes] + [
                            x.name for x in proteins], [x.name for x in numbering_schemes] + [x.entry_name for x in proteins])
    context['segments'] = [x.slug for x in segments if len(data[x.slug])]
    context['data'] = flattened_data
    context['number_of_schemes'] = len(numbering_schemes)

    residuelist = Residue.objects.filter(protein_conformation__protein=p).prefetch_related('protein_segment','display_generic_number','generic_number')
    HelixBox = DrawHelixBox(
                residuelist, p.get_protein_class(), str(p), nobuttons=1)
    if not pdbname.startswith('AFM'):
        SnakePlot = DrawSnakePlot(
                    residuelist, p.get_protein_class(), str(p), nobuttons=1)
    else:
        SnakePlot = []
    #adjusting main_ligand and main_ligand_full
    if len(main_ligand) == 0:
        multiple_ligands = False
        main_ligand = "None"
        main_ligand_full = "None"
    elif len(main_ligand) == 1:
        multiple_ligands = False
        main_ligand = main_ligand[0]
        main_ligand_full = main_ligand_full[0]
    else:
        multiple_ligands = True

    return render(request, 'interaction/structure.html', {'pdbname': pdbname, 'structures': structures,
                                                          'crystal': crystal, 'protein': p, 'helixbox' : HelixBox, 'snakeplot': SnakePlot, 'residues': residues_browser, 'residues_lookup': residues_lookup, 'display_res': display_res, 'annotated_resn':
                                                          resn_list, 'ligands': ligands,'main_ligand' : main_ligand,'main_ligand_full' : main_ligand_full, 'data': context['data'],
                                                          'header': context['header'], 'segments': context['segments'],
                                                          'number_of_schemes': len(numbering_schemes), "multiple_ligands": multiple_ligands})


def remove_duplicate_dicts(dict_list):
    unique_dicts = []
    seen = set()

    for d in dict_list:
        # Convert the dictionary to a string to make it hashable
        dict_str = json.dumps(d, sort_keys=True)

        # If the string is not in the set, add it back to the list and mark it as seen
        if dict_str not in seen:
            seen.add(dict_str)
            unique_dicts.append(d)

    return unique_dicts

def find_dict_index(dict_list, target_dict):
    for i, d in enumerate(dict_list):
        if d == target_dict:
            return i
    return -1  # return -1 if the dictionary is not found

def sort_and_update(push_gpcr, gpcr, push_gprot, gprot, interactions):
  for key in push_gpcr:
      gpcr[key]['interactions'] = ', '.join(push_gpcr[key])
  for key in push_gprot:
      gprot[key]['interactions'] = ', '.join(push_gprot[key])

  matching_dict = {}
  for record in interactions:
      if record['innerIndex'] not in matching_dict.keys():
          matching_dict[record['innerIndex']] = []
      if record['outerIndex'] not in matching_dict[record['innerIndex']]:
          matching_dict[record['innerIndex']].append(record['outerIndex'])

  outer_to_inner = {}
  for key, value in matching_dict.items():
      for idx in value:
          if idx not in outer_to_inner.keys():
              outer_to_inner[idx] = []
              outer_to_inner[idx].append(key)
  sorted_outer = {k: v for k, v in sorted(outer_to_inner.items(), key=lambda item: len(item[1]), reverse=True)}

  return gpcr, gprot, sorted_outer


def ComplexDetails(request, pdbname):
    """
    Show structure details
    """
    pdbname = pdbname

    structures = Interaction.objects.values('interacting_pair__referenced_structure__protein_conformation__protein__entry_name', 'interacting_pair__referenced_structure__signprot_complex__protein__entry_name').filter(
        interacting_pair__referenced_structure__pdb_code__index=pdbname).annotate(numRes=Count('pk', distinct=True)).order_by('-numRes')
    resn_list = ''

    crystal = Structure.objects.get(pdb_code__index=pdbname)
    if crystal.structure_type.slug.startswith('af-'):
        p = Protein.objects.get(id=crystal.protein_conformation.protein.id)
    else:
        p = Protein.objects.get(protein=crystal.protein_conformation.protein)

    residuelist = Residue.objects.filter(protein_conformation__protein=p).prefetch_related('protein_segment','display_generic_number','generic_number')
    lookup = {}

    residues_lookup = {}
    for r in residuelist:
        if r.generic_number:
            lookup[r.generic_number.label] = r.sequence_number
            residues_lookup[r.sequence_number] = r.amino_acid +str(r.sequence_number)+ " "+ r.generic_number.label

    interactions = Interaction.objects.filter(
        interacting_pair__referenced_structure__pdb_code__index=pdbname).exclude(interaction_type ='hidden')

    residues_browser = []
    gpcr_residues = []
    gprot_residues = []
    display_res = []
    residue_table_list = []

    for residue in interactions:
        type = residue.interaction_type
        gpcr_aa = residue.interacting_pair.res1.amino_acid
        gprot_aa = residue.interacting_pair.res2.amino_acid
        gpcr_pos = residue.interacting_pair.res1.sequence_number
        gprot_pos = residue.interacting_pair.res2.sequence_number
        segment = residue.interacting_pair.res1.protein_segment.slug
        try:
            gpcr_grn = residue.interacting_pair.res1.display_generic_number.label
        except AttributeError:
            gpcr_grn = '-'
        try:
            gprot_grn = residue.interacting_pair.res2.display_generic_number.label
        except AttributeError:
            gprot_grn = '-'
        # gpcr_grn = residue.interacting_pair.res1.generic_number.label
        # gprot_grn = residue.interacting_pair.res2.generic_number.label

        display_res.append(str(gpcr_pos))
        display_res.append(str(gprot_pos))

        residues_browser.append({'type': type, 'gpcr_aa': gpcr_aa, 'gprot_aa': gprot_aa,
                                 'gpcr_pos': gpcr_pos, 'gprot_pos': gprot_pos,
                                 'gpcr_grn': gpcr_grn, 'gprot_grn': gprot_grn, 'segment': segment})

    residues_browser = remove_duplicate_dicts(residues_browser)
    display_res = ' or '.join(display_res)

    # RESIDUE TABLE
    segments = ProteinSegment.objects.all().filter().prefetch_related()
    proteins = [p]
    numbering_schemes_selection = [settings.DEFAULT_NUMBERING_SCHEME]
    numbering_schemes_selection.append(p.residue_numbering_scheme.slug)

    numbering_schemes = ResidueNumberingScheme.objects.filter(slug__in=numbering_schemes_selection).all()
    default_scheme = numbering_schemes.get(slug=settings.DEFAULT_NUMBERING_SCHEME)

    data = OrderedDict()

    for segment in segments:
        data[segment.slug] = OrderedDict()
        residues = Residue.objects.filter(protein_segment=segment,  protein_conformation__protein=p,
                                          generic_number__label__in=residue_table_list).prefetch_related('protein_conformation__protein',
                                                                                                         'protein_conformation__state', 'protein_segment',
                                                                                                         'generic_number', 'display_generic_number', 'generic_number__scheme',
                                                                                                         'alternative_generic_numbers__scheme')
        for scheme in numbering_schemes:
            if scheme == default_scheme and scheme.slug == settings.DEFAULT_NUMBERING_SCHEME:
                for pos in list(set([x.generic_number.label for x in residues if x.protein_segment == segment])):
                    data[segment.slug][pos] = {
                        scheme.slug: pos, 'seq': ['-'] * len(proteins)}
            elif scheme == default_scheme:
                for pos in list(set([x.generic_number.label for x in residues if x.protein_segment == segment])):
                    data[segment.slug][pos] = {
                        scheme.slug: pos, 'seq': ['-'] * len(proteins)}

        for residue in residues:
            alternatives = residue.alternative_generic_numbers.all()
            pos = residue.generic_number
            for alternative in alternatives:
                if alternative.scheme not in numbering_schemes:
                    continue
                scheme = alternative.scheme
                if default_scheme.slug == settings.DEFAULT_NUMBERING_SCHEME:
                    pos = residue.generic_number
                    if scheme == pos.scheme:
                        data[segment.slug][pos.label]['seq'][proteins.index(
                            residue.protein_conformation.protein)] = str(residue)
                    else:
                        if scheme.slug not in data[segment.slug][pos.label].keys():
                            data[segment.slug][pos.label][
                                scheme.slug] = alternative.label
                        if alternative.label not in data[segment.slug][pos.label][scheme.slug]:
                            data[segment.slug][pos.label][
                                scheme.slug] += " " + alternative.label
                        data[segment.slug][pos.label]['seq'][proteins.index(
                            residue.protein_conformation.protein)] = str(residue)
                else:
                    if scheme.slug not in data[segment.slug][pos.label].keys():
                        data[segment.slug][pos.label][
                            scheme.slug] = alternative.label
                    if alternative.label not in data[segment.slug][pos.label][scheme.slug]:
                        data[segment.slug][pos.label][
                            scheme.slug] += " " + alternative.label
                    data[segment.slug][pos.label]['seq'][proteins.index(
                        residue.protein_conformation.protein)] = str(residue)

    # Preparing the dictionary of list of lists. Dealing with tripple nested
    # dictionary in django templates is a nightmare
    flattened_data = OrderedDict.fromkeys([x.slug for x in segments], [])
    for s in iter(flattened_data):
        flattened_data[s] = [[data[s][x][y.slug]
                              for y in numbering_schemes] + data[s][x]['seq'] for x in sorted(data[s])]

    context = {}
    context['header'] = zip([x.short_name for x in numbering_schemes] + [x.name for x in proteins], [x.name for x in numbering_schemes] + [
                            x.name for x in proteins], [x.name for x in numbering_schemes] + [x.entry_name for x in proteins])
    context['segments'] = [x.slug for x in segments if len(data[x.slug])]
    context['data'] = flattened_data
    context['number_of_schemes'] = len(numbering_schemes)

    residuelist = Residue.objects.filter(protein_conformation__protein=p).prefetch_related('protein_segment','display_generic_number','generic_number')
    HelixBox = DrawHelixBox(residuelist, p.get_protein_class(), str(p), nobuttons=1)

    ### Implementing same code for calculating the interactions plot
    model = Structure.objects.get(pdb_code__index=pdbname)
    if model.structure_type.slug == 'af-signprot':
        scores = StructureAFScores.objects.get(structure=model)
    else:
        scores = StructureAFScores()
    #Need to build the plDDT colors
    model_plddt = StructureModelpLDDT.objects.filter(structure=model)
    residues_plddt = {}
    for item in model_plddt:
        residues_plddt[item.residue.id] = [item.residue, item.pLDDT]

    ### Gathering interaction info and structuring JS data
    interactions = Interaction.objects.filter(interacting_pair__referenced_structure=model).prefetch_related('interacting_pair__res1', 'interacting_pair__res2')

    gpcr_aminoacids = []
    gprot_aminoacids = []
    protein_interactions = []
    protein_interactions_strict = []
    conversion = {'aromatic': 'Aromatic',
                  'hydrophobic': 'Hydrophobic',
                  'ionic': 'Ionic',
                  'polar': 'Polar',
                  'van-der-waals': 'Van der waals'}

    for pair in interactions:
        try:
            gn1 = pair.interacting_pair.res1.display_generic_number.label
        except AttributeError:
            gn1 = '-'

        try:
            gn2 = pair.interacting_pair.res2.display_generic_number.label
        except AttributeError:
            gn2 = '-'

        gpcr = {'aminoAcid': pair.interacting_pair.res1.amino_acid,
                'segment': pair.interacting_pair.res1.protein_segment.slug,
                'generic_number': gn1,
                'sequence_number': pair.interacting_pair.res1.sequence_number,
                'interaction_level': pair.interaction_level
                }
        gprot = {'aminoAcid': pair.interacting_pair.res2.amino_acid,
                'segment': pair.interacting_pair.res2.protein_segment.slug,
                'generic_number': gn2,
                'sequence_number': pair.interacting_pair.res2.sequence_number,
                'interaction_level': pair.interaction_level
                }

        gpcr_aminoacids.append(gpcr)
        gprot_aminoacids.append(gprot)

    gpcr_aminoacids_strict = [record for record in gpcr_aminoacids if record['interaction_level'] == 1]
    gprot_aminoacids_strict = [record for record in gprot_aminoacids if record['interaction_level'] == 1]

    gpcr_aminoacids = [{key: value for key, value in d.items() if key != 'interaction_level'} for d in gpcr_aminoacids]
    gprot_aminoacids = [{key: value for key, value in d.items() if key != 'interaction_level'} for d in gprot_aminoacids]
    gpcr_aminoacids_strict = [{key: value for key, value in d.items() if key != 'interaction_level'} for d in gpcr_aminoacids_strict]
    gprot_aminoacids_strict = [{key: value for key, value in d.items() if key != 'interaction_level'} for d in gprot_aminoacids_strict]

    gpcr_aminoacids_strict = remove_duplicate_dicts(gpcr_aminoacids_strict)
    gprot_aminoacids_strict = remove_duplicate_dicts(gprot_aminoacids_strict)
    gpcr_aminoacids = remove_duplicate_dicts(gpcr_aminoacids)
    gprot_aminoacids = remove_duplicate_dicts(gprot_aminoacids)

    segments_order = ['TM1','ICL1', 'TM2', 'ICL2', 'TM3', 'ICL3', 'TM4', 'TM5', 'TM6', 'TM7', 'ICL4', 'H8', 'C-term']
    gprot_segments = ['G.HN','G.hns1','G.S1','G.s1h1','G.H1','G.h1ha','H.HA','H.hahb','H.HB','H.hbhc','H.HC','H.hchd','H.HD','H.hdhe','H.HE','H.hehf','H.HF','G.hfs2','G.S2','G.s2s3','G.S3','G.s3h2','G.H2','G.h2s4','G.S4','G.s4h3','G.H3','G.h3s5','G.S5','G.s5hg','G.HG','G.hgh4','G.H4','G.h4s6','G.S6','G.s6h5','G.H5']
    # Create a dictionary that maps segments to their positions in the custom order
    order_gpcr = {segment: index for index, segment in enumerate(segments_order)}
    order_gprot = {segment: index for index, segment in enumerate(gprot_segments)}

    # Sort the list of dictionaries based on the custom order
    # gpcr_aminoacids = sorted(gpcr_aminoacids, key=lambda x: order_gpcr.get(x['segment'], float('inf')))
    # gprot_aminoacids = sorted(gprot_aminoacids, key=lambda x: x['generic_number'])
    gprot_aminoacids = sorted(gprot_aminoacids, key=lambda x: (order_gprot.get(x['segment'], 9999), int(x['generic_number'].split('.')[-1])))
    gprot_aminoacids_strict = sorted(gprot_aminoacids_strict, key=lambda x: (order_gprot.get(x['segment'], 9999), int(x['generic_number'].split('.')[-1])))

    to_push_gpcr = {}
    to_push_gprot = {}
    to_push_gpcr_strict = {}
    to_push_gprot_strict = {}
    for pair in interactions:
        if pair.atomname_residue1 in ['C', 'CA', 'N', 'O']:
            chain_res1 = 'Main'
        else:
            chain_res1 = 'Side'
        if pair.atomname_residue2 in ['C', 'CA', 'N', 'O']:
            chain_res2 = 'Main'
        else:
            chain_res2 = 'Side'
        try:
            gn1 = pair.interacting_pair.res1.display_generic_number.label
        except AttributeError:
            gn1 = '-'
        try:
            gn2 = pair.interacting_pair.res2.display_generic_number.label
        except AttributeError:
            gn2 = '-'
        gpcr = {'aminoAcid': pair.interacting_pair.res1.amino_acid,
                'segment': pair.interacting_pair.res1.protein_segment.slug,
                'generic_number': gn1,
                'sequence_number': pair.interacting_pair.res1.sequence_number
                }
        gprot = {'aminoAcid': pair.interacting_pair.res2.amino_acid,
                'segment': pair.interacting_pair.res2.protein_segment.slug,
                'generic_number': gn2,
                'sequence_number': pair.interacting_pair.res2.sequence_number
                }

        gpcr_index = find_dict_index(gpcr_aminoacids, gpcr)
        gprot_index = find_dict_index(gprot_aminoacids, gprot)
        protein_interactions.append({'innerIndex': gprot_index, 'outerIndex': gpcr_index, 'type': conversion[pair.interaction_type], 'innerChain': chain_res2, 'outerChain': chain_res1, 'interaction_level': pair.interaction_level})

        if gpcr_index not in to_push_gpcr.keys():
            to_push_gpcr[gpcr_index] = []
        if gprot_index not in to_push_gprot.keys():
            to_push_gprot[gprot_index] = []

        if (conversion[pair.interaction_type]) not in to_push_gpcr[gpcr_index]:
            to_push_gpcr[gpcr_index].append(conversion[pair.interaction_type])
        if (conversion[pair.interaction_type]) not in to_push_gprot[gprot_index]:
            to_push_gprot[gprot_index].append(conversion[pair.interaction_type])

        if pair.interaction_level == 1:
            gpcr_index_strict = find_dict_index(gpcr_aminoacids_strict, gpcr)
            gprot_index_strict = find_dict_index(gprot_aminoacids_strict, gprot)
            protein_interactions_strict.append({'innerIndex': gprot_index_strict, 'outerIndex': gpcr_index_strict, 'type': conversion[pair.interaction_type], 'innerChain': chain_res2, 'outerChain': chain_res1, 'interaction_level': pair.interaction_level})
            ### Copy logic for the strict ones
            if gpcr_index_strict not in to_push_gpcr_strict.keys():
                to_push_gpcr_strict[gpcr_index_strict] = []
            if gprot_index_strict not in to_push_gprot_strict.keys():
                to_push_gprot_strict[gprot_index_strict] = []

            if (conversion[pair.interaction_type]) not in to_push_gpcr_strict[gpcr_index_strict]:
                to_push_gpcr_strict[gpcr_index_strict].append(conversion[pair.interaction_type])
            if (conversion[pair.interaction_type]) not in to_push_gprot_strict[gprot_index_strict]:
                to_push_gprot_strict[gprot_index_strict].append(conversion[pair.interaction_type])


    protein_interactions = remove_duplicate_dicts(protein_interactions)
    protein_interactions_strict = remove_duplicate_dicts(protein_interactions_strict)

    gpcr_aminoacids, gprot_aminoacids, matching_dict = sort_and_update(to_push_gpcr, gpcr_aminoacids, to_push_gprot, gprot_aminoacids, protein_interactions)
    gpcr_aminoacids_strict, gprot_aminoacids_strict, matching_dict_strict = sort_and_update(to_push_gpcr_strict, gpcr_aminoacids_strict, to_push_gprot_strict, gprot_aminoacids_strict, protein_interactions_strict)

    return render(request, 'interaction/structure.html', {'pdbname': pdbname, 'structures': structures,
                                                          'crystal': crystal, 'protein': p, 'helixbox' : HelixBox,
                                                          'residues': residues_browser, 'residues_lookup': residues_lookup,
                                                          'display_res': display_res, 'data': context['data'],
                                                          'header': context['header'], 'segments': context['segments'],
                                                          'number_of_schemes': context['number_of_schemes'], 'complex': True,
                                                          'scores': scores, 'refined': True,
                                                          'outer': json.dumps(gpcr_aminoacids),
                                                          'inner': json.dumps(gprot_aminoacids),
                                                          'interactions': json.dumps(protein_interactions),
                                                          'outer_strict': json.dumps(gpcr_aminoacids_strict),
                                                          'inner_strict': json.dumps(gprot_aminoacids_strict),
                                                          'interactions_strict': json.dumps(protein_interactions_strict),
                                                          'conversion_dict': json.dumps(matching_dict),
                                                          'conversion_dict_strict': json.dumps(matching_dict_strict)})


def list_structures(request):
    form = PDBform()
    #structures = ResidueFragmentInteraction.objects.distinct('structure_ligand_pair__structure').all()
    structures = ResidueFragmentInteraction.objects.exclude(interaction_type__slug='acc').values('structure_ligand_pair__structure__pdb_code__index', 'structure_ligand_pair__structure__protein_conformation__protein__parent__entry_name').annotate(
        num_ligands=Count('structure_ligand_pair', distinct=True), num_interactions=Count('pk', distinct=True)).order_by('structure_ligand_pair__structure__protein_conformation__protein__parent__entry_name')
    #structures = ResidueFragmentInteraction.objects.values('structure_ligand_pair__structure__pdb_code__index','structure_ligand_pair__structure').annotate(Count('structure_ligand_pair__ligand'))
    # print(structures.count())
    genes = {}
    countligands = {}
    totalligands = 0
    totalinteractions = 0
    totaltopinteractions = 0
    for structure in structures:
        # print(structure)
        if structure['structure_ligand_pair__structure__protein_conformation__protein__parent__entry_name'] not in genes:
            genes[structure[
                'structure_ligand_pair__structure__protein_conformation__protein__parent__entry_name']] = 1
        totalligands += structure['num_ligands']
        totalinteractions += structure['num_interactions']
        ligands = ResidueFragmentInteraction.objects.exclude(interaction_type__slug='acc').values('structure_ligand_pair__ligand__name').filter(structure_ligand_pair__structure__pdb_code__index=structure[
            'structure_ligand_pair__structure__pdb_code__index']).annotate(numRes=Count('pk', distinct=True)).order_by('-numRes')
        for ligand in ligands:
            totaltopinteractions += ligand['numRes']
            if ligand['structure_ligand_pair__ligand__name'] not in countligands:
                countligands[ligand['structure_ligand_pair__ligand__name']] = 1
            break

        # print(structure.structure_ligand_pair.structure.pdb_code.index)
        # print(structure.numRes)
    #objects = Model.objects.filter(id__in=object_ids)
    #context = {}
    # print('Structures with ligand information:' + str(structures.count()))
    # print('Distinct genes:' + str(len(genes)))
    #print('ligands:' + str(totalligands))
    # print('interactions:' + str(totalinteractions))
    # print('interactions from top ligands:' + str(totaltopinteractions))
    # print('Distinct ligands as top ligand:' + str(len(countligands)))

    return render(request, 'interaction/list.html', {'form': form, 'structures': structures})


def crystal(request):
    pdbname = request.GET.get('pdb')
    form = PDBform()
    structures = ResidueFragmentInteraction.objects.values('structure_ligand_pair__ligand__name').filter(
        structure_ligand_pair__structure__pdb_code__index=pdbname).annotate(numRes=Count('pk', distinct=True)).order_by('-numRes')
    crystal = Structure.objects.get(pdb_code__index=pdbname)
    p = Protein.objects.get(protein=crystal.protein_conformation.protein)
    residues = ResidueFragmentInteraction.objects.filter(
        structure_ligand_pair__structure__pdb_code__index=pdbname).order_by('rotamer__residue__sequence_number')
    # print("residues", residues)
    return render(request, 'interaction/crystal.html', {'form': form, 'pdbname': pdbname, 'structures': structures, 'crystal': crystal, 'protein': p, 'residues': residues})


def view(request):
    pdbname = request.GET.get('pdb')
    form = PDBform()
    structures = ResidueFragmentInteraction.objects.values('structure_ligand_pair__ligand__name').filter(
        structure_ligand_pair__structure__pdb_code__index=pdbname).annotate(numRes=Count('pk', distinct=True)).order_by('-numRes')
    return render(request, 'interaction/view.html', {'form': form, 'pdbname': pdbname, 'structures': structures})


def ligand(request):
    pdbname = request.GET.get('pdb')
    ligand = request.GET.get('ligand')
    form = PDBform()
    fragments = ResidueFragmentInteraction.objects.filter(structure_ligand_pair__structure__pdb_code__index=pdbname).filter(
        structure_ligand_pair__ligand__name=ligand).order_by('interaction_type')
    return render(request, 'interaction/ligand.html', {'form': form, 'pdbname': pdbname, 'ligand': ligand, 'fragments': fragments})


def fragment(request):
    pdbname = request.GET.get('pdb')
    ligand = request.GET.get('ligand')
    fragment = request.GET.get('fragment')
    form = PDBform()
    fragments = ResidueFragmentInteraction.objects.get(id=fragment)
    return render(request, 'interaction/fragment.html', {'form': form, 'pdbname': pdbname, 'ligand': ligand, 'fragmentid': fragment, 'fragments': fragments})

def check_residue(protein, pos, aa):
    residue = Residue.objects.filter(
        protein_conformation=protein, sequence_number=pos)
    if residue.exists():
        if len(residue) > 1:
            residue = residue[0]
        else:
            residue = Residue.objects.get(
                protein_conformation=protein, sequence_number=pos)
        if residue.amino_acid != aa:
            residue.amino_acid = aa
            residue.save()
    else:
        residue, created = Residue.objects.get_or_create(
            protein_conformation=protein, sequence_number=pos, amino_acid=aa)
        # continue #SKIP THESE -- mostly fusion residues that aren't mapped
        # yet.
    return residue


def extract_fragment_rotamer(f, residue, structure, ligand):
    if os.path.isfile(f):
        f_in = open(f, 'r')
        rotamer_pdb = ''
        fragment_pdb = ''
        for line in f_in:
            if line.startswith('HETATM') or line.startswith('CONECT') or line.startswith('MASTER') or line.startswith('END'):
                fragment_pdb += line
            elif line.startswith('ATOM'):
                rotamer_pdb += line
            else:
                fragment_pdb += line
                rotamer_pdb += line
        f_in.close()

        rotamer_data, created = PdbData.objects.get_or_create(pdb=rotamer_pdb)
        try:
            rotamer = Rotamer.objects.get(residue=residue, structure=structure)
        except Rotamer.DoesNotExist:
            rotamer, _ = Rotamer.objects.get_or_create(residue=residue, structure=structure, pdbdata=rotamer_data)

        fragment_data, _ = PdbData.objects.get_or_create(pdb=fragment_pdb)
        fragment, created = Fragment.objects.get_or_create(ligand=ligand, structure=structure, pdbdata=fragment_data, residue=residue)
    else:
        #quit("Could not find " + residue)
        return None, None

    return fragment, rotamer

# NOTE: this function is solely used by the sitesearch functionality
def calculate(request, redirect=None):
    # convert identified interactions to residue features and add them to the session
    # numbers in lists represent the interaction "hierarchy", i.e. if a residue has more than one
    # interaction,
    interaction_name_dict = {
        'polar_double_neg_protein': [1, 'neg'],
        'polar_double_pos_protein': [1, 'neg'],
        'polar_pos_protein': [2, 'pos'],
        'polar_neg_protein': [3, 'neg'],
        'polar_neg_ligand': [4, 'hbd'],
        'polar_pos_ligand': [5, 'hba'],
        'polar_unknown_protein': [5, 'charge'],
        'polar_donor_protein': [6, 'hbd'],
        'polar_acceptor_protein': [7, 'hba'],
        'polar_unspecified': [8, 'hb'],
        'aro_ff': [9, 'ar'],
        'aro_ef_protein': [10, 'ar'],
        'aro_fe_protein':  [11, 'ar'],
        'aro_ion_protein':  [12, 'pos'],
        'aro_ion_ligand':  [12, 'ar'],
    }

    if request.method == 'POST':
        form = PDBform(request.POST, request.FILES)
        if form.is_valid():

            pdbname = form.cleaned_data['pdbname'].strip()
            results = ''

            if not request.session.exists(request.session.session_key):
                request.session.create()
            session_key = request.session.session_key

            module_dirs = []
            module_dir = '/tmp/interactions'
            module_dirs.append(module_dir)
            module_dir = os.sep.join([module_dir, session_key])
            module_dirs.append(module_dir)
            module_dirs.append(os.sep.join([module_dir, 'pdbs']))
            module_dirs.append(os.sep.join([module_dir, 'temp']))

            # create dirs and set permissions (needed on some systems)
            for mdir in module_dirs:
                os.makedirs(mdir, exist_ok=True)
                # os.chmod(mdir, 0o777)

            if 'file' in request.FILES:
                pdbdata = request.FILES['file']
                pdbname = os.path.splitext(str(pdbdata))[0]
                pdbname = pdbname.replace("_","")
                with open(module_dir + '/pdbs/' + str(pdbdata).replace("_",""), 'wb+') as destination:
                    for chunk in pdbdata.chunks():
                        destination.write(chunk)

                temp_path = module_dir + '/pdbs/' + str(pdbdata).replace("_","")
                pdbdata = open(temp_path, 'r').read()
                results = runusercalculation_2022(pdbname, session_key)

            else:
                pdbname = form.cleaned_data['pdbname'].strip()
                temp_path = module_dir + '/pdbs/' + pdbname + '.pdb'

                if not os.path.isfile(temp_path):
                    url = 'http://www.rcsb.org/pdb/files/%s.pdb' % pdbname
                    pdbdata = urllib.request.urlopen(url).read().decode('utf-8')
                    f = open(temp_path, 'w')
                    f.write(pdbdata)
                    f.close()
                else:
                    pdbdata = open(temp_path, 'r').read()
                results = runusercalculation_2022(pdbname, session_key)

            # MAPPING GPCRdb numbering onto pdb.
            generic_numbering = GenericNumbering(temp_path,top_results=1, blastdb=os.sep.join([settings.STATICFILES_DIRS[0], 'blast', 'protwis_gpcr_blastdb']))
            out_struct = generic_numbering.assign_generic_numbers()
            structure_residues = generic_numbering.residues
            prot_id_list = generic_numbering.prot_id_list
            segments = {}

            generic_ids = []
            generic_number = []
            previous_seg = 'N-term'

            #Get segments built correctly for non-aligned residues
            for c, res in structure_residues.items():
                for i, r in sorted(res.items()):  # sort to be able to assign loops
                    if r.gpcrdb:
                        if r.gpcrdb[0] == '-':
                            # fix stefan - for bulge
                            r.gpcrdb = r.gpcrdb[1:] + "1"
                        r.gpcrdb = str(r.gpcrdb).replace('.', 'x')
                        generic_number.append(r.gpcrdb)
                    if r.gpcrdb_id:
                        generic_ids.append(r.gpcrdb_id)
                    if r.segment:
                        if not r.segment in segments:
                            segments[r.segment] = {}
                        segments[r.segment][r.number] = [r.display, r.name, r.gpcrdb,r.residue_record]
                        previous_seg = r.segment
                    else:  # if no segment assigned by blast
                        if previous_seg in ['N-term', 'ICL1', 'ECL1', 'ICL2', 'ECL2', 'ICL3', 'ICL3', 'C-term']:
                            if not previous_seg in segments:
                                segments[previous_seg] = {}
                            segments[previous_seg][r.number] = ['', r.name, '',r.residue_record]
                        else:
                            if previous_seg == 'TM1':
                                previous_seg = 'ICL1'
                            elif previous_seg == 'TM2':
                                previous_seg = 'ECL1'
                            elif previous_seg == 'TM3':
                                previous_seg = 'ICL2'
                            elif previous_seg == 'TM4':
                                previous_seg = 'ECL2'
                            elif previous_seg == 'TM5':
                                previous_seg = 'ICL3'
                            elif previous_seg == 'TM6':
                                previous_seg = 'ECL3'
                            elif previous_seg == 'TM7':
                                previous_seg = 'C-term'
                            elif previous_seg == 'H8':
                                previous_seg = 'C-term'

                            if not previous_seg in segments:
                                segments[previous_seg] = {}
                            segments[previous_seg][r.number] = ['', r.name, '',r.residue_record]

            residue_list = []
            for seg, reslist in segments.items():
                for seq_number, v in sorted(reslist.items()):
                    if v[3]: #if blast assigned a residue, then use it. Otherwise just make something empty
                        r = v[3]
                    else:
                        r = Residue()
                    r.sequence_number = seq_number
                    r.segment_slug = seg
                    r.amino_acid = v[1]
                    residue_list.append(r)

            xtal = {}
            hetsyn = {}
            hetsyn_reverse = {}
            for line in pdbdata.splitlines():
                if line.startswith('HETSYN'):
                    # need to fix bad PDB formatting where col4 and col5 are
                    # put together for some reason -- usually seen when the id
                    # is +1000
                    m = re.match("HETSYN[\s]+([\w]{3})[\s]+(.+)", line)
                    if (m):
                        hetsyn[m.group(2).strip()] = m.group(1).upper()
                        hetsyn_reverse[m.group(1)] = m.group(2).strip().upper()
                if line.startswith('HETNAM'):
                    # need to fix bad PDB formatting where col4 and col5 are
                    # put together for some reason -- usually seen when the id
                    # is +1000
                    m = re.match("HETNAM[\s]+([\w]{3})[\s]+(.+)", line)
                    if (m):
                        hetsyn[m.group(2).strip()] = m.group(1).upper()
                        hetsyn_reverse[m.group(1)] = m.group(2).strip().upper()
                if line.startswith('REVDAT   1'):
                    xtal['publication_date'] = line[13:22]
                    xtal['pdb_code'] = line[23:27]
                if line.startswith('JRNL        PMID'):
                    xtal['pubmed_id'] = line[19:].strip()
                if line.startswith('JRNL        DOI'):
                    xtal['doi_id'] = line[19:].strip()
                if line.startswith('REMARK   2 RESOLUTION.'):
                    xtal['resolution'] = line[22:].strip()

            simple = collections.OrderedDict()
            simple_generic_number = collections.OrderedDict()
            residues_browser = []
            residue_table_list = []
            mainligand = ''

            #Switching to new results value
            for ligand in results.keys():
                ligand_score = round(results[ligand]['score'])
                # select top hit
                if mainligand == '':
                    mainligand = ligand
                simple[ligand] = {'score': ligand_score}
                simple_generic_number[ligand] = {'score': ligand_score}
                for interaction in results[ligand]['interactions']:
                    aa, pos, chain = regexaa(interaction[0])
                    if int(pos) in structure_residues[chain]:
                        r = structure_residues[chain][int(pos)]
                        display = r.display
                        segment = r.segment
                        generic = r.gpcrdb
                        if generic != "":
                            residue_table_list.append(generic)
                            if generic not in simple_generic_number[ligand]:
                                simple_generic_number[ligand][generic] = []
                            simple_generic_number[ligand][generic].append(interaction[2])
                    else:
                        display = ''
                        segment = ''

                    if interaction[0] in simple[ligand]:
                        simple[ligand][interaction[0]].append(interaction[2])
                    else:
                        simple[ligand][interaction[0]] = [interaction[2]]

                    residues_browser.append({'type': interaction[3], 'aa': aa, 'ligand': ligand, 'pos': pos, 'gpcrdb': display, 'segment': segment, 'slug':interaction[2]})

            # RESIDUE TABLE
            segments = ProteinSegment.objects.all().filter().prefetch_related()
            proteins = []
            protein_list = Protein.objects.filter(pk__in=prot_id_list)
            numbering_schemes_selection = [settings.DEFAULT_NUMBERING_SCHEME]
            for p in protein_list:
                proteins.append(p)
                if p.residue_numbering_scheme.slug not in numbering_schemes_selection:
                    numbering_schemes_selection.append(p.residue_numbering_scheme.slug)

            numbering_schemes = ResidueNumberingScheme.objects.filter(slug__in=numbering_schemes_selection).all()
            default_scheme = numbering_schemes.get(slug=settings.DEFAULT_NUMBERING_SCHEME)
            data = OrderedDict()

            for segment in segments:
                data[segment.slug] = OrderedDict()
                residues = Residue.objects.filter(protein_segment=segment,  protein_conformation__protein__in=proteins,
                                                  generic_number__label__in=residue_table_list).prefetch_related('protein_conformation__protein',
                                                                                                                 'protein_conformation__state', 'protein_segment',
                                                                                                                 'generic_number', 'display_generic_number', 'generic_number__scheme',
                                                                                                                 'alternative_generic_numbers__scheme')
                for scheme in numbering_schemes:
                    if scheme == default_scheme and scheme.slug == settings.DEFAULT_NUMBERING_SCHEME:
                        for pos in list(set([x.generic_number.label for x in residues if x.protein_segment == segment])):
                            data[segment.slug][pos] = {
                                scheme.slug: pos, 'seq': ['-'] * len(proteins)}
                    elif scheme == default_scheme:
                        for pos in list(set([x.generic_number.label for x in residues if x.protein_segment == segment])):
                            data[segment.slug][pos] = {
                                scheme.slug: pos, 'seq': ['-'] * len(proteins)}

                for residue in residues:
                    alternatives = residue.alternative_generic_numbers.all()
                    pos = residue.generic_number
                    for alternative in alternatives:
                        if alternative.scheme not in numbering_schemes:
                            continue
                        scheme = alternative.scheme
                        if default_scheme.slug == settings.DEFAULT_NUMBERING_SCHEME:
                            pos = residue.generic_number
                            if scheme == pos.scheme:
                                data[segment.slug][pos.label]['seq'][proteins.index(
                                    residue.protein_conformation.protein)] = str(residue)
                            else:
                                if scheme.slug not in data[segment.slug][pos.label].keys():
                                    data[segment.slug][pos.label][
                                        scheme.slug] = alternative.label
                                if alternative.label not in data[segment.slug][pos.label][scheme.slug]:
                                    data[segment.slug][pos.label][
                                        scheme.slug] += " " + alternative.label
                                data[segment.slug][pos.label]['seq'][proteins.index(
                                    residue.protein_conformation.protein)] = str(residue)
                        else:
                            if scheme.slug not in data[segment.slug][pos.label].keys():
                                data[segment.slug][pos.label][
                                    scheme.slug] = alternative.label
                            if alternative.label not in data[segment.slug][pos.label][scheme.slug]:
                                data[segment.slug][pos.label][
                                    scheme.slug] += " " + alternative.label
                            data[segment.slug][pos.label]['seq'][proteins.index(
                                residue.protein_conformation.protein)] = str(residue)

            # Preparing the dictionary of list of lists. Dealing with tripple
            # nested dictionary in django templates is a nightmare
            flattened_data = OrderedDict.fromkeys([x.slug for x in segments], [])
            for s in iter(flattened_data):
                flattened_data[s] = [[data[s][x][y.slug] for y in numbering_schemes] + data[s][x]['seq'] for x in sorted(data[s])]

            context = {}
            context['header'] = zip([x.short_name for x in numbering_schemes] + [x.name for x in proteins], [x.name for x in numbering_schemes] +
                                    [x.name for x in proteins], [x.name for x in numbering_schemes] + [x.entry_name for x in proteins])
            context['segments'] = [x.slug for x in segments if len(data[x.slug])]
            context['data'] = flattened_data
            context['number_of_schemes'] = len(numbering_schemes)

            if redirect:
                # get simple selection from session
                simple_selection = request.session.get('selection', False)
                # create full selection and import simple selection (if it
                # exists)
                selection = Selection()
                if simple_selection:
                    selection.importer(simple_selection)

                interaction_counter = 0
                for gn, interactions in simple_generic_number[mainligand].items():
                    if gn != 'score' and gn != 0.0:  # FIXME leave these out when dict is created
                        feature = False
                        for interaction in interactions:
                            if interaction in interaction_name_dict:
                                if (not feature
                                    or interaction_name_dict[interaction][0] < interaction_name_dict[feature][0]):
                                    feature = interaction

                        if not feature:
                            continue

                        # get residue number equivalent object
                        rne = ResidueGenericNumberEquivalent.objects.get(label=gn, scheme__slug='gpcrdba')

                        # create a selection item
                        properties = {
                            'feature': interaction_name_dict[feature][1],
                            'amino_acids': ','.join(definitions.AMINO_ACID_GROUPS_OLD[interaction_name_dict[feature][1]])
                        }
                        selection_item = SelectionItem('site_residue', rne, properties)

                        # add to selection
                        selection.add('segments', 'site_residue', selection_item)

                        # update the minimum match count for the active group
                        interaction_counter += 1
                        selection.site_residue_groups[selection.active_site_residue_group - 1][0] = interaction_counter

                # export simple selection that can be serialized
                simple_selection = selection.exporter()

                # add simple selection to session
                request.session['selection'] = simple_selection

                # re-direct to segment selection (with the extracted interactions already selected)
                return HttpResponseRedirect(redirect)
            else:
                # Only relevant when not redirecting - moved here
                HelixBox = DrawHelixBox(residue_list, 'Class A', str('test'), nobuttons=1)
                SnakePlot = DrawSnakePlot(residue_list, 'Class A', str('test'), nobuttons=1)

                return {'result': "Looking at " + pdbname, 'outputs': results,
                                                                    'simple': simple, 'simple_generic_number': simple_generic_number, 'xtal': xtal, 'pdbname': pdbname, 'mainligand': mainligand, 'residues': residues_browser,
                                                                    'HelixBox': HelixBox, 'SnakePlot': SnakePlot, 'data': context['data'],
                                                                    'header': context['header'], 'segments': context['segments'], 'number_of_schemes': len(numbering_schemes), 'proteins': proteins}

        else:
            print(form.errors)
            return HttpResponse("Error with form ")
    else:
        return HttpResponse("Ooops how did you get here?")


def download(request):
    pdbname = request.GET.get('pdb')
    ligand = request.GET.get('ligand')

    session = request.GET.get('session')

    if session:
        session = request.session.session_key
        pdbdata = open('/tmp/interactions/' + session + '/results/' + pdbname +
                       '/interaction/' + pdbname + '_' + ligand + '.pdb', 'r').read()
        response = HttpResponse(pdbdata, content_type='text/plain')
    else:
        # here's the needed fix
        try:
            pair = StructureLigandInteraction.objects.filter(structure__pdb_code__index=pdbname).filter(
                Q(ligand__inchikey=ligand) | Q(ligand__name=ligand)).exclude(pdb_file__isnull=True).get()
            response = HttpResponse(pair.pdb_file.pdb, content_type='text/plain')
        except:
            return redirect("/interaction/")
    return response

def excel(request, slug, **response_kwargs):
    if ('session' in response_kwargs):
        session = request.session.session_key

        mypath = '/tmp/interactions/' + session + '/pdbs/' + slug + '.pdb'

        generic_numbering = GenericNumbering(mypath)
        out_struct = generic_numbering.assign_generic_numbers()
        structure_residues = generic_numbering.residues
        results = runusercalculation_2022(slug, session)

        data = []
        ligand = list(results.keys())[0]
        for interaction in results[ligand]['interactions']:
            aa, pos, chain = regexaa(interaction[0])
            if int(pos) in structure_residues[chain]:
                r = structure_residues[chain][int(pos)]
                display = r.display
                segment = r.segment
            else:
                display = ''
                segment = ''
            row = {}
            row['Sequence Number'] = pos
            row['Amino Acid'] = aa
            row['Generic Number'] = display
            row['Segment'] = segment
            row['Interaction'] = interaction[3]
            row['Interaction Slug'] = interaction[2]
            row['Ligand'] = results[ligand]['prettyname']
            data.append(row)

    else:

        interactions = ResidueFragmentInteraction.objects.filter(
            structure_ligand_pair__structure__pdb_code__index=slug, structure_ligand_pair__annotated=True).order_by('rotamer__residue__sequence_number')
        data = []
        for interaction in interactions:
            row = {}
            row['Sequence Number'] = interaction.rotamer.residue.sequence_number
            row['Amino Acid'] = interaction.rotamer.residue.amino_acid
            if interaction.rotamer.residue.display_generic_number:
                row['Generic Number'] = interaction.rotamer.residue.display_generic_number.label
                row['Segment'] = interaction.rotamer.residue.protein_segment.slug
            else:
                row['Generic Number'] = 'N/A'
                row['Segment'] = '-'

            row['Interaction'] = interaction.interaction_type.name
            row['Interaction Slug'] = interaction.interaction_type.slug
            row['Ligand'] = interaction.structure_ligand_pair.ligand.name

            data.append(row)


    headers = ['Ligand','Amino Acid','Sequence Number','Generic Number','Segment','Interaction','Interaction Slug']

    #EXCEL SOLUTION
    output = BytesIO()
    workbook = xlsxwriter.Workbook(output)
    worksheet = workbook.add_worksheet()

    col = 0
    for h in headers:
        worksheet.write(0, col, h)
        col += 1
    row = 1
    for d in data:
        col = 0
        for h in headers:
            worksheet.write(row, col, d[h])
            col += 1
        row += 1
    workbook.close()
    output.seek(0)
    xlsx_data = output.read()

    response = HttpResponse(xlsx_data,content_type='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet')
    response['Content-Disposition'] = 'attachment; filename=Interaction_data_%s.xlsx' % slug
    return response

def ajax(request, slug, **response_kwargs):

    residuelist = Residue.objects.filter(protein_conformation__protein__entry_name=slug).prefetch_related('protein_segment','display_generic_number','generic_number')
    lookup = {}

    for r in residuelist:
        if r.generic_number:
            lookup[r.generic_number.label] = r.sequence_number

    interactions = ResidueFragmentInteraction.objects.filter(
        structure_ligand_pair__structure__protein_conformation__protein__parent__entry_name=slug, structure_ligand_pair__annotated=True).exclude(interaction_type__type ='hidden').order_by('rotamer__residue__sequence_number')
    # return HttpResponse("Hello, world. You're at the polls index. "+slug)
    jsondata = {}
    for interaction in interactions:
        if interaction.rotamer.residue.generic_number:
            sequence_number = interaction.rotamer.residue.sequence_number
            if interaction.rotamer.residue.generic_number.label in lookup:
                sequence_number = lookup[interaction.rotamer.residue.generic_number.label]

            #label = interaction.rotamer.residue.generic_number.label
            aa = interaction.rotamer.residue.amino_acid
            interactiontype = interaction.interaction_type.name
            if sequence_number not in jsondata:
                jsondata[sequence_number] = []
            jsondata[sequence_number].append([aa, interactiontype])

    jsondata = json.dumps(jsondata)
    response_kwargs['content_type'] = 'application/json'
    return HttpResponse(jsondata, **response_kwargs)


def ajaxLigand(request, slug, ligand, **response_kwargs):
    residuelist = Residue.objects.filter(protein_conformation__protein__entry_name=slug).prefetch_related('protein_segment','display_generic_number','generic_number')
    lookup = {}

    for r in residuelist:
        if r.generic_number:
            lookup[r.generic_number.label] = r.sequence_number


    interactions = ResidueFragmentInteraction.objects.filter(
        structure_ligand_pair__structure__protein_conformation__protein__parent__entry_name=slug, structure_ligand_pair__ligand__name=ligand).exclude(interaction_type__type ='hidden').order_by('rotamer__residue__sequence_number')
    # return HttpResponse("Hello, world. You're at the polls index. "+slug)
    jsondata = {}
    for interaction in interactions:
        sequence_number = interaction.rotamer.residue.sequence_number
        sequence_number = lookup[interaction.rotamer.residue.generic_number.label]
        aa = interaction.rotamer.residue.amino_acid
        interactiontype = interaction.interaction_type.name
        if sequence_number not in jsondata:
            jsondata[sequence_number] = []
        jsondata[sequence_number].append([aa, interactiontype])

    jsondata = json.dumps(jsondata)
    response_kwargs['content_type'] = 'application/json'
    return HttpResponse(jsondata, **response_kwargs)


def pdbfragment(request):
    #pdbname = request.GET.get('pdb')
    #ligand = request.GET.get('ligand')
    fragment = request.GET.get('fragment')

    result = ResidueFragmentInteraction.objects.filter(id=fragment).get()
    response = HttpResponse(result.rotamer.pdbdata.pdb +
                            result.fragment.pdbdata.pdb, content_type='text/plain')
    return response


def pdb(request):
    pdbname = request.GET.get('pdb')
    session = request.GET.get('session')
    # response = HttpResponse(mimetype='application/force-download')
    # #response['Content-Disposition'] = 'attachment; filename=%s' % smart_str(file_name)
    # response['Content-Disposition'] = 'attachment; filename=%s' % smart_str(pdbname+'_'+ligand+'.pdb')
    # mypath = module_dir+'/temp/results/'+pdbname+'/interaction/'+pdbname+'_'+ligand+'.pdb'
    # response['X-Sendfile'] = smart_str(mypath)
    if session:
        session = request.session.session_key
        pdbdata = open('/tmp/interactions/' + session +
                       '/pdbs/' + pdbname + '.pdb', 'r').read()
        response = HttpResponse(pdbdata, content_type='text/plain')
    else:
        web_resource = WebResource.objects.get(slug='pdb')
        web_link, created = WebLink.objects.get_or_create(
            web_resource=web_resource, index=pdbname)

        structure = Structure.objects.filter(pdb_code=web_link)
        if structure.exists():
            structure = Structure.objects.get(pdb_code=web_link)
        else:
            quit()  # quit!

        if structure.pdb_data is None:
            quit()

        response = HttpResponse(structure.pdb_data.pdb,
                                content_type='text/plain')
    return response

def complexpdb(request):
    pdbname = request.GET.get('pdb')
    session = request.GET.get('session')
    if session:
        session = request.session.session_key
        pdbdata = open('/tmp/interactions/complex/' + session +
                       '/pdbs/' + pdbname + '.pdb', 'r').read()
        response = HttpResponse(pdbdata, content_type='text/plain')
    else:
        web_resource = WebResource.objects.get(slug='pdb')
        web_link, created = WebLink.objects.get_or_create(
            web_resource=web_resource, index=pdbname)

        structure = Structure.objects.filter(pdb_code=web_link)
        if structure.exists():
            structure = Structure.objects.get(pdb_code=web_link)
        else:
            quit()  # quit!

        if structure.pdb_data is None:
            quit()

        response = HttpResponse(structure.pdb_data.pdb,
                                content_type='text/plain')
    return response
