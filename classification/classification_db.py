"""
DB-driven classification-tree data, replacing the Excel-file-based tree builder that
`Classification_tree` (classification/views.py) used to run. `ProteinFamilyClassification`
(populated from the same Excel file by `build_classification_annotations`, but at build time
rather than per-request) is a family-level, verified-lossless equivalent of the Excel rows.
"""
from collections import OrderedDict

from django.core.cache import cache, caches
from django.db.models import Prefetch

from protein.models import Gene, Protein, ProteinFamily, ProteinFamilyClassification

try:
    cache_alignment = caches["alignments"]
except Exception:
    cache_alignment = cache

CLASSIFICATION_ROWS_CACHE_KEY = "classification:tree_rows:db:v2"
CLASSIFICATION_ROWS_CACHE_TIMEOUT = 60 * 60 * 24  # 24h

D3_EMPTY_LABEL = "Other / unknown"


def get_classification_rows(use_cache=True):
    """
    One dict per (ProteinFamilyClassification row x matching human SWISSPROT protein) --
    the DB equivalent of one Excel row. A family can carry 2 rows (chemotype/modality
    order 1 and 2), so a receptor with two tags appears twice, once per tag -- this is what
    lets it show up correctly under both in the Modality/Chemotype tabs.
    """
    cached = cache_alignment.get(CLASSIFICATION_ROWS_CACHE_KEY) if use_cache else None
    if cached is not None:
        return cached

    class_names = {
        fam.slug: fam.name
        for fam in ProteinFamily.objects.filter(parent__slug='000')
    }

    pfc_rows = list(
        ProteinFamilyClassification.objects
        .select_related('protein_family__parent', 'chemotype', 'modality', 'sense')
    )
    family_ids = {r.protein_family_id for r in pfc_rows}

    proteins_by_family = {}
    protein_qs = (
        Protein.objects
        .filter(family_id__in=family_ids, species__common_name='Human', source__name='SWISSPROT')
        .only('id', 'entry_name', 'name', 'family_id', 'sequence')
        .prefetch_related(
            Prefetch(
                'genes',
                queryset=Gene.objects.only('name', 'position').order_by('position'),
                to_attr='primary_genes_self',
            )
        )
        .order_by('entry_name')
    )
    for protein in protein_qs:
        proteins_by_family.setdefault(protein.family_id, []).append(protein)

    rows = []
    for pfc in pfc_rows:
        fam = pfc.protein_family
        class_slug = fam.slug.split('_')[0]
        receptor_family = fam.parent.name if fam.parent else ""
        for protein in proteins_by_family.get(fam.id, ()):
            genes = getattr(protein, 'primary_genes_self', None) or []
            try:
                protein_label = protein.short()
            except Exception:
                protein_label = protein.name or ""
            rows.append({
                "family_id": fam.id,
                "family_slug": fam.slug,
                "receptor_family": receptor_family,
                "class_slug": class_slug,
                "class_family_name": class_names.get(class_slug, ""),
                "chemotype": (pfc.chemotype.name if pfc.chemotype else ""),
                "chemotype_order": pfc.chemotype_order,
                "modality": (pfc.modality.name if pfc.modality else ""),
                "modality_order": pfc.modality_order,
                "sense": (pfc.sense.name if pfc.sense else ""),
                "protein_id": protein.id,
                "entry_name": protein.entry_name,
                "uniprot": protein.entry_name.split('_', 1)[0].upper(),
                "gene": genes[0].name if genes else "",
                "protein_label": protein_label,
                "protein_name": protein.name or "",
                "sequence": protein.sequence or "",
            })

    if use_cache:
        cache_alignment.set(CLASSIFICATION_ROWS_CACHE_KEY, rows, CLASSIFICATION_ROWS_CACHE_TIMEOUT)
    return rows


def get_primary_classification_rows(use_cache=True):
    """
    Like get_classification_rows(), but collapsed to one row per protein -- for consumers
    (GPCRBrowser, ClassificationWheel) that want exactly one row per receptor rather than one
    per (protein x tag). A protein whose family carries both a chemotype/modality order-1 and
    order-2 tag is deduped to its order-1 row (same "primary row" preference already used by
    StructureSim._get_pf_classification_map_db).
    """
    by_protein = OrderedDict()
    for row in get_classification_rows(use_cache=use_cache):
        existing = by_protein.get(row["protein_id"])
        if existing is None:
            by_protein[row["protein_id"]] = row
            continue
        is_primary = row["chemotype_order"] == 1 or row["modality_order"] == 1
        existing_is_primary = existing["chemotype_order"] == 1 or existing["modality_order"] == 1
        if is_primary and not existing_is_primary:
            by_protein[row["protein_id"]] = row
    return list(by_protein.values())


def get_pf_classification_map(use_cache=True):
    """
    {protein_family_id: [{'Chemotype', 'Modality', 'Sense', 'chemotype_order', 'modality_order'}, ...]}
    A list per family (not a single dict) since a family can carry up to 2 rows.
    """
    rows = get_classification_rows(use_cache=use_cache)
    mapping = {}
    for row in rows:
        mapping.setdefault(row["family_id"], []).append({
            "Chemotype": row["chemotype"],
            "Modality": row["modality"],
            "Sense": row["sense"],
            "chemotype_order": row["chemotype_order"],
            "modality_order": row["modality_order"],
        })
    return mapping


def d3_leaf(name):
    return OrderedDict([("name", name), ("value", 0), ("color", "")])


def d3_node(name, children):
    return OrderedDict([("name", name), ("value", 0), ("color", ""), ("children", children)])


def d3_root(children):
    return OrderedDict([("name", ""), ("value", 3000), ("color", ""), ("children", children)])


def build_grouped_children(rows, level_getters, level_sort_keys=None):
    """
    Generic replacement for the six Excel-specific `_build_nested*`/`_nested*_to_tree` builders:
    group `rows` by an ordered list of `level_getters` (row -> str), with `row["uniprot"]` as the
    leaf under the last level. `level_sort_keys` optionally overrides the default
    case-insensitive-name sort at a given depth (0-indexed) -- e.g. to sort the top level by a
    curated class order instead of alphabetically.
    """
    max_depth = len(level_getters)
    nested = {}
    for row in rows:
        node = nested
        for get in level_getters:
            key = (get(row) or "").strip() or D3_EMPTY_LABEL
            node = node.setdefault(key, {})
        node.setdefault("__leaves__", set()).add(row["uniprot"])

    def emit(mapping, depth):
        out = []
        sort_key = (level_sort_keys or {}).get(depth) or (lambda x: str(x).lower())
        for key in sorted(mapping.keys(), key=sort_key):
            sub = mapping[key]
            if depth == max_depth - 1:
                leaves = sub.get("__leaves__", set())
                children = [d3_leaf(u) for u in sorted(leaves, key=lambda x: str(x).upper())]
            else:
                children = emit(sub, depth + 1)
            out.append(d3_node(key, children))
        return out

    return emit(nested, 0)


def build_leaf_label_lookup(rows):
    """Pure-DB replacement for Classification_tree._build_leaf_label_lookup."""
    lookup = {}
    for row in rows:
        key = row["uniprot"]
        entry = lookup.setdefault(key, {"UniProt": key, "Gene": "", "Protein": ""})
        if row["gene"] and not entry["Gene"]:
            entry["Gene"] = row["gene"]
        if row["protein_label"] and not entry["Protein"]:
            entry["Protein"] = row["protein_label"]
    for key, entry in lookup.items():
        entry["Protein"] = entry["Protein"] or key
        entry["Gene"] = entry["Gene"] or key
    return lookup
