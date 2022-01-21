#from django.db import connection
#from django.db import IntegrityError
from django.db.models import Q
from ligand.models import Ligand, LigandProperties, LigandType, LigandRole, TestLigand
from common.models import WebLink, WebResource
from common.tools import get_or_create_url_cache, fetch_from_web_api

import time
import os
import datamol as dm

def get_or_create_smallmolecule(name, ids = {}):
    """ This function tries to obtain a small molecule Ligand object.

        If the ligand already exists it will return the corresponding object. If
        not, it will processes the provided data and create the ligand object.
        When a ligand does not exist, or cannot be created, with the given input
        it will return None

        Parameters
        ----------
        name : str
            Name of the ligand
        ids : dict, optional
            Dictionary where the key is the id_type and the value is the id
            <id_type> is e.g. pubchem, gtp, chembl, etc. (slug of the web_resource)
            <id> is the Identifier itself

        Returns
        -------
        obj
            Ligand object or None

        Raises
        ------
        ...Error
            If no ...
    """

    ligand = None

    # Check and filter IDs
    for type_id in list(ids.keys()):
        if ids[type_id] == None or ids[type_id] == "None" or ids[type_id] == "":
            del ids[type_id]
        elif ids[type_id].isnumeric():
            ids[type_id] = int(ids[type_id])

    #print("Adding",name, ids)

    # No IDs are provided, so the only way forward is a name match
    if len(ids) == 0 and name.length > 0:
        results = TestLigand.objects.filter(name=name)
        if results.count() == 1:
            ligand = results.first()
            # DEBUGGING
            print("Semi-ambiguous matching using name")
        # DEBUGGING
        elif results.count() > 1:
            print("Ambiguous name (", name,") as it has ", results.count(), " corresponding entries")
        else:
            print("No ligand found with name",name)
    else:
        # PubChem ID, Guide to Pharmacology ID, ChEMBL ID
        for type in ["pubchem", "gtoplig", "chembl_ligand"]:
            if ligand == None and type in ids:
                ligand = get_ligand_by_id(type, ids[type])

        # InChiKey
        if ligand == None and "inchikey" in ids:
            ligand = get_ligand_by_inchikey(ids["inchikey"])

        # SMILES and calculated InChIKeys
        if ligand == None and "smiles" in ids:
            result = TestLigand.objects.filter(properties__smiles = ids["smiles"])
            if result.count() > 0:
                ligand = result.first()
                # DEBUGGING
                if result.count() > 1:
                    print("Ambiguous SMILES (", smiles,") as it has ", results.count(), " corresponding entries")
            elif not "inchikey" in ids:
                # calculate inchikey from given SMILES and repeat inchikey check
                input_mol = dm.to_mol(ids["smiles"], sanitize=False)
                ligand = get_ligand_by_inchikey(dm.to_inchikey(input_mol))
            elif ligand == None:
                # calculate cleaned InchiKey and repeate check
                clean_inchi = get_cleaned_inchikey(ids["smiles"])
                if clean_inchi != None:
                    ligand = get_ligand_by_inchikey(clean_inchi)

        # If the ligand has not been found => create ligand using the provided IDs
        # TODO merge this into one function
        if ligand == None:
            for type in ["smiles", "pubchem", "gtoplig", "chembl_ligand"]:
                if type in ids:
                    ligand = create_ligand_from_id(name, type, ids[type])
                    if ligand != None:
                        # Store provided InChIKey in the properties
                        if "inchikey" in ids:
                            ligand.properties.inchikey = ids["inchikey"]
                        break

        # Add missing IDs via weblinks to the ligand object
        if ligand != None:
            # Create list of existing weblinks
            current_ids = []
            for wl in ligand.properties.web_links.all():
                current_ids.append(wl.index)

            for type_id in ids:
                # If the ID does not yet exist => add
                if ids[type_id] not in current_ids and type not in ["smiles", "inchikey"]:
                    wr = WebResource.objects.get(slug=type)
                    wl, created = WebLink.objects.get_or_create(index=ids[type_id], web_resource=wr)
                    ligand.properties.web_links.add(wl)
    return ligand


def get_ligand_by_id(type, id):
    result = TestLigand.objects.filter(properties__web_links__index=id, properties__web_links__web_resource__slug=type)
    if result.count() > 0:
        if result.count() > 1:
            print("Multiple entries for the same ID - This should never happen - error", type, id)
        return result.first()
    else:
        return None

def create_ligand_from_id(name, type, id):
    # create new ligand
    ligand = TestLigand()
    ligand.name = name
    ligand.canonical = True # TODO: discuss if we should keep this
    ligand.ambiguous = False

    # create corresponding ligand properties object
    lp = LigandProperties()
    lp.ligand_type = LigandType.objects.get(slug='small-molecule')
    # Obtain external data (if needed)
    try:
        if type == "smiles":
            lp.smiles = id
        elif type == "pubchem":
            # check cache
            cache_dir = ['pubchem', 'cid', 'property']
            url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/$index/property/IsomericSMILES,InChIKey/json'
            pubchem = fetch_from_web_api(url, id, cache_dir)

            if pubchem != False and pubchem['PropertyTable']['Properties'][0]:
                lp.smiles = pubchem['PropertyTable']['Properties'][0]['IsomericSMILES']
                lp.inchikey = pubchem['PropertyTable']['Properties'][0]['InChIKey']
        elif type == "gtoplig":
            # check cache
            cache_dir = ['guidetopharmacology', 'ligand_structure']
            url = 'http://www.guidetopharmacology.org/services/ligands/$index/structure'
            gtop = fetch_from_web_api(url, id, cache_dir)

            if gtop:
                lp.smiles = gtop["smiles"]
                lp.inchikey = gtop["inchiKey"]
        elif type == "chembl_ligand":
            # check cache
            cache_dir = ['chembl', 'ligand_entry']
            url = 'https://www.ebi.ac.uk/chembl/api/data/molecule/$index/?format=json'
            chembl = fetch_from_web_api(url, id, cache_dir)

            # NOTE: canonical SMILES from ChEMBL are isomeric smiles
            if chembl and "molecule_structures" in chembl:
                lp.smiles = chembl["molecule_structures"]["canonical_smiles"]
                lp.inchikey = chembl["molecule_structures"]["standard_inchi_key"]
    except:
        skip=True

    # perform RDkit calculations
    if lp.smiles != None and lp.smiles != "":
        input_mol = dm.to_mol(lp.smiles, sanitize=True)
        if input_mol:
            # Check if InChIKey has been set
            if lp.inchikey == None:
                lp.inchikey = dm.to_inchikey(input_mol)

            # Check if cleaned InChIKey has been set
            if lp.clean_inchikey == None:
                lp.clean_inchikey = get_cleaned_inchikey(lp.smiles)

            # Calculate RDkit properties
            lp.mw = dm.descriptors.mw(input_mol)
            lp.rotatable_bonds = dm.descriptors.n_rotatable_bonds(input_mol)
            lp.hacc = dm.descriptors.n_hba(input_mol)
            lp.hdon = dm.descriptors.n_hbd(input_mol)
            lp.logp = dm.descriptors.clogp(input_mol)
        else:
            print("Issues with molecule", name)
            # TODO - try again if we have other IDs

        # Before storing - check one more time if we already have this ligand
        if lp.inchikey != None:
            checkligand = get_ligand_by_inchikey(lp.inchikey)
            if checkligand != None:
                return checkligand
        if lp.clean_inchikey != None and lp.clean_inchikey != lp.inchikey:
            checkligand = get_ligand_by_inchikey(lp.clean_inchikey)
            if checkligand != None:
                return checkligand

        # Store results and return
        lp.save()
        ligand.properties = lp
        ligand.save()
        return ligand
    else:
        # Ligand could not be created
        return None

def get_ligand_by_inchikey(inchikey):
    result = TestLigand.objects.filter(Q(properties__inchikey = inchikey) | Q(properties__clean_inchikey = inchikey))
    if result.count() > 0:
        if result.count() > 1:
            print("Multiple entries for the same InChIkey - This should never happen - error", inchikey)
        return result.first()
    else:
        return None

def get_cleaned_inchikey(smiles):
    inchikey = None
    try:
        input_mol = dm.to_mol(smiles, sanitize=True)
        if input_mol:
            input_mol = dm.to_neutral(input_mol, verbose=False)
            input_mol = dm.fix_mol(input_mol, verbose=False)
            cleaned_smiles = dm.standardize_smiles(dm.to_smiles(input_mol), tautomer=True, verbose=False)
            cleaned_mol = dm.to_mol(cleaned_smiles, sanitize=False)
            inchikey = dm.to_inchikey(cleaned_mol)
    except:
        skip = True
    return inchikey
