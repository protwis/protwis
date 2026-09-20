import hashlib
from collections import OrderedDict
from types import SimpleNamespace

from django.contrib.postgres.aggregates import ArrayAgg

from common.selection import SelectionItem, SimpleSelection
from phylogenetic_trees.views import Treeclass
from protein.models import ProteinCouplings, ProteinSegment


GENERIC_CONSERVED_SEGMENT_FAMILY = "GPCR"
TREE_SETTINGS = ["0", "0", "0", "0"]
TREE_METHOD = "neighbor_joining"
SEGMENT_SOURCE = "generic_conserved"
BRANCH_MODE = "regular"


def _clean_annotation_text(value):
    return str(value or "").replace("<sub>", "").replace("</sub>", "").replace("<i>", "").replace("</i>", "").strip()


def _top_class_family(protein):
    family = getattr(protein, "family", None)
    ligand = getattr(family, "parent", None) if family else None
    receptor_family = getattr(ligand, "parent", None) if ligand else None
    return getattr(receptor_family, "parent", None) if receptor_family else None


def _receptor_family_node(protein):
    family = getattr(protein, "family", None)
    receptor_family = getattr(family, "parent", None) if family else None
    return receptor_family or family


def resolve_family_network_nodes(proteins, family_entry):
    family_obj = None
    class_family = None
    family_name = str(family_entry.get("name") or "").strip().lower()

    for protein in proteins:
        receptor_family = _receptor_family_node(protein)
        if not family_obj and receptor_family and str(receptor_family.name or "").strip().lower() == family_name:
            family_obj = receptor_family
        if not class_family:
            class_family = _top_class_family(protein)
        if family_obj and class_family:
            break

    if not family_obj and proteins:
        family_obj = _receptor_family_node(proteins[0])
    if not class_family and proteins:
        class_family = _top_class_family(proteins[0])
    return family_obj, class_family


def family_source_hash(proteins):
    digest = hashlib.sha1()
    for protein_id in sorted(int(protein.id) for protein in proteins):
        digest.update(str(protein_id).encode("ascii"))
        digest.update(b":")
    return digest.hexdigest()


def _generic_conserved_segments():
    segments = list(
        ProteinSegment.objects
        .filter(proteinfamily=GENERIC_CONSERVED_SEGMENT_FAMILY, fully_aligned=True)
        .order_by("id")
    )
    if not segments:
        segments = list(
            ProteinSegment.objects
            .filter(proteinfamily=GENERIC_CONSERVED_SEGMENT_FAMILY, partial=False)
            .order_by("id")
        )
    for segment in segments:
        segment.only_aligned_residues = True
    return segments


def _selection_for_proteins(proteins):
    selection = SimpleSelection()
    selection.targets = [SelectionItem("protein", protein) for protein in proteins]
    selection.segments = [SelectionItem("segment", segment) for segment in _generic_conserved_segments()]
    selection.tree_settings = list(TREE_SETTINGS)
    return selection


def _degenerate_newick(proteins):
    labels = [str(protein.entry_name) for protein in proteins if getattr(protein, "entry_name", None)]
    if not labels:
        return ""
    if len(labels) == 1:
        return "{};".format(labels[0])
    if len(labels) == 2:
        return "({0}:0,{1}:0);".format(labels[0], labels[1])
    return ""


def _build_tree_newick(proteins):
    if len(proteins) < 3:
        return _degenerate_newick(proteins)

    fake_request = SimpleNamespace(session={"selection": _selection_for_proteins(proteins)})
    tree = Treeclass()
    phylogeny_input, branches, ttype, total, legend, box, additional_info, buttons, protein_conformations = tree.Prepare_file(fake_request)
    if phylogeny_input in {"too big", "More_prots"}:
        raise ValueError("Unable to build family tree for group {}".format(total))
    return str(tree.phylip or "").replace("\n", "")


def _protein_annotation(protein):
    family = getattr(protein, "family", None)
    receptor_family = getattr(family, "parent", None) if family else None
    ligand_type = getattr(receptor_family, "parent", None) if receptor_family else None
    gpcr_class = getattr(ligand_type, "parent", None) if ligand_type else None
    return [
        _clean_annotation_text(getattr(protein, "name", "")),
        _clean_annotation_text(getattr(receptor_family, "name", "") or getattr(family, "name", "")),
        _clean_annotation_text(getattr(ligand_type, "name", "")),
        _clean_annotation_text(getattr(gpcr_class, "name", "")),
        str(getattr(family, "slug", "") or ""),
        _clean_annotation_text(getattr(protein, "name", "")),
    ]


def _build_annotations(proteins):
    annotations = OrderedDict()
    for protein in sorted(proteins, key=lambda item: str(item.entry_name or "")):
        annotations[str(protein.entry_name)] = _protein_annotation(protein)
    return annotations


def _build_gprotein_coupling(proteins):
    protein_slugs = sorted({str(getattr(protein.family, "slug", "") or "") for protein in proteins if getattr(protein, "family", None)})
    if not protein_slugs:
        return {}

    coupling = (
        ProteinCouplings.objects
        .filter(protein__family__slug__in=protein_slugs, source="GuideToPharma")
        .values_list("protein__family__slug", "transduction")
        .annotate(arr=ArrayAgg("g_protein__name"))
    )
    payload = {}
    for family_slug, transduction, g_proteins in coupling:
        if family_slug not in payload:
            payload[family_slug] = {}
        payload[family_slug][transduction] = g_proteins
    return payload


def build_family_tree_payload(proteins, family_entry, *, build_version="v1"):
    ordered_proteins = list(sorted(proteins, key=lambda item: str(item.entry_name or "")))
    note = (
        "Neighbor-joining family tree built from generic conserved GPCR segments."
        if len(ordered_proteins) >= 3
        else "Degenerate family tree stored because fewer than 3 receptors are available for this group."
    )
    tree_newick = _build_tree_newick(ordered_proteins)
    return {
        "tree": tree_newick,
        "annotations": _build_annotations(ordered_proteins),
        "Gprot_coupling": _build_gprotein_coupling(ordered_proteins),
        "meta": {
            "family": str(family_entry.get("label") or family_entry.get("name") or ""),
            "group_key": str(family_entry.get("key") or ""),
            "n_points": len(ordered_proteins),
            "classes": list(family_entry.get("class_labels") or []),
            "chemotypes": list(family_entry.get("chemotypes") or []),
            "modality_groups": list(family_entry.get("modality_groups") or []),
            "tree_method": TREE_METHOD,
            "segment_source": SEGMENT_SOURCE,
            "bootstrap": 0,
            "branch_mode": BRANCH_MODE,
            "build_version": str(build_version),
            "is_degenerate": len(ordered_proteins) < 3,
            "note": note,
        },
    }
