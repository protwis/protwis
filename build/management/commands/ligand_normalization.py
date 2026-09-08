from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import inchi as rdinchi
import datamol as dm


def _normalize_base(mol):
    mol = rdMolStandardize.Cleanup(mol)
    mol = rdMolStandardize.MetalDisconnector().Disconnect(mol)
    mol = rdMolStandardize.LargestFragmentChooser().choose(mol)
    mol = rdMolStandardize.Normalizer().normalize(mol)
    mol = rdMolStandardize.Uncharger().uncharge(mol)
    mol = rdMolStandardize.TautomerEnumerator().Canonicalize(mol)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    return mol


def normalize_for_parent(smiles):
    """Normalize SMILES for parent record: strips stereo and isotopes.

    Returns dict with smiles, inchikey, clean_inchikey, or None on failure.
    Result is stored on the parent row.
    """
    mol = dm.to_mol(smiles, sanitize=True)
    if mol is None:
        return None
    try:
        mol = _normalize_base(mol)
        Chem.RemoveStereochemistry(mol)
        for atom in mol.GetAtoms():
            atom.SetIsotope(0)
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
        mol = _normalize_base(mol)
        ik = rdinchi.InchiToInchiKey(rdinchi.MolToInchi(mol))
        return {
            "inchikey": ik,
            "mol": mol,
        }
    except Exception:
        return None


def infer_stereo_status(mol):
    """Classify stereo status from a normalized mol object.

    Returns one of: 'full', 'partial', 'unspecified', 'racemic', or None.
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
