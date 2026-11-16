from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import inchi as rdinchi
import datamol as dm


# Tautomer canonicalization is ~99% of normalization time and grows with size
# (0.26 s/mol at <=50 heavy atoms, 1.3-3 s above), while above 50 heavy atoms it
# almost never changes the InChIKey connectivity block (measured on ChEMBL peptides).
# Larger molecules skip it; standard InChI still unifies mobile-H tautomers. The gate
# depends only on the molecule, so every source gets the same key for it.
MAX_TAUTOMER_HEAVY_ATOMS = 50


def _normalize_base(mol, keep_stereo=False):
    mol = rdMolStandardize.Cleanup(mol)
    mol = rdMolStandardize.MetalDisconnector().Disconnect(mol)
    mol = rdMolStandardize.LargestFragmentChooser().choose(mol)
    mol = rdMolStandardize.Normalizer().normalize(mol)
    mol = rdMolStandardize.Uncharger().uncharge(mol)
    if mol.GetNumHeavyAtoms() <= MAX_TAUTOMER_HEAVY_ATOMS:
        enumerator = rdMolStandardize.TautomerEnumerator()
        if keep_stereo:
            # By default canonicalization deletes sp3/E-Z stereo on atoms/bonds touched by
            # tautomerism (e.g. amino-acid alpha-CH), giving enantiomers the same key.
            enumerator.SetRemoveSp3Stereo(False)
            enumerator.SetRemoveBondStereo(False)
        mol = enumerator.Canonicalize(mol)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    return mol


def normalize_for_parent(smiles):
    """Normalize SMILES for parent record: strips stereo and isotopes.

    Stereo and isotopes are stripped *before* standardization: stereo influences
    tautomer canonicalization, so stripping it afterwards made the parent key depend
    on whether the input SMILES carried stereo (e.g. corticosterone got two keys).

    Returns dict with smiles, inchikey, clean_inchikey, or None on failure.
    Result is stored on the parent row.
    """
    mol = dm.to_mol(smiles, sanitize=True)
    if mol is None:
        return None
    try:
        Chem.RemoveStereochemistry(mol)
        for atom in mol.GetAtoms():
            atom.SetIsotope(0)
        mol = _normalize_base(mol)
        Chem.RemoveStereochemistry(mol)
        smi = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
        ik = rdinchi.InchiToInchiKey(rdinchi.MolToInchi(mol))
        return {
            "smiles": smi,
            "inchikey": ik,
            "clean_inchikey": ik.split("-")[0] if ik else None,
        }
    except Exception:
        return None


def normalize_for_child(smiles):
    """Normalize SMILES for child matching: preserves stereo and isotopes.

    Returns dict with inchikey and mol, or None on failure.
    Result is used for matching only — NOT stored (child stores raw SMILES).
    """
    mol = dm.to_mol(smiles, sanitize=True)
    if mol is None:
        return None
    try:
        mol = _normalize_base(mol, keep_stereo=True)
        ik = rdinchi.InchiToInchiKey(rdinchi.MolToInchi(mol))
        return {
            "inchikey": ik,
            "mol": mol,
        }
    except Exception:
        return None


def infer_stereo_status(mol):
    """Classify stereo status from a normalized mol object.

    Returns one of: 'full', 'partial', 'unspecified', or None. Structure alone cannot
    tell a racemate from undefined stereo; 'racemic' is assigned from source metadata
    in build_ligand_functions (_is_racemic_entry).
    """
    if mol is None:
        return None
    try:
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        stereo_info = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        if not stereo_info:
            return "unspecified"
        specified = sum(1 for _, s in stereo_info if s in ("R", "S"))
        unspecified = len(stereo_info) - specified
        if unspecified == 0:
            return "full"
        if specified == 0:
            return "unspecified"
        return "partial"
    except Exception:
        return None
