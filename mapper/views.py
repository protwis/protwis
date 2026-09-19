from django.http import HttpResponse, JsonResponse
from django.db.models import Q
from django.core.cache import cache
from django.views.generic import TemplateView, View
from protein.models import Protein, ProteinFamily
from classification.models import ClusterCoord

import json
from copy import deepcopy
from collections import OrderedDict
import sys


try:
    import importlib.metadata
except ImportError:
    sys.modules['importlib.metadata'] = __import__('importlib_metadata')

import pandas as pd
import numpy as np
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans
import re

class DataMapperHome(TemplateView):

    @staticmethod
    def keep_by_names(data, names_to_keep):
        data_copy = deepcopy(data)
        if isinstance(data_copy, list):
            # Process each item in the list
            kept_items = [DataMapperHome.keep_by_names(item, names_to_keep) for item in data_copy]
            # Return only non-None items
            return [item for item in kept_items if item is not None]
        elif isinstance(data_copy, OrderedDict):
            name = data_copy.get('name')
            if name not in names_to_keep.keys():
                if 'children' in data_copy:
                    # Recursively process children
                    data_copy['children'] = DataMapperHome.keep_by_names(data_copy['children'], names_to_keep)
                    # Remove the 'children' key if it's empty after processing
                    if not data_copy['children']:
                        return None
                else:
                    return None
            else:
                # If the name is in the keep list, update the 'value' field
                if 'Inner' in names_to_keep[name]:
                    data_copy['value'] = names_to_keep[name]['Inner']
                # Process children if present
                if 'children' in data_copy:
                    data_copy['children'] = DataMapperHome.keep_by_names(data_copy['children'], names_to_keep)
                    if not data_copy['children']:
                        del data_copy['children']
            return data_copy
        return data_copy

    # Query filters for GenerateGPCRomeDataStructure's two data_type variants. Pulled out so the
    # (previously near-duplicated) query setup can be shared by one helper -- the circle-bucketing
    # logic that follows is genuinely different between Classic and Odorant and stays separate.
    _GPCROME_QUERY_FILTERS = {
        "Classic": {
            "protein_q": Q(family_id__slug__startswith='0'),
            "protein_exclude_q": Q(family_id__slug__startswith='007') | Q(family_id__slug__startswith='008'),
            "family_q": Q(slug__startswith='0'),
            "family_exclude_q": Q(slug='000') | Q(slug__startswith='007') | Q(slug__startswith='008'),
        },
        "Odorant": {
            "protein_q": Q(family_id__slug__startswith='007') | Q(family_id__slug__startswith='008'),
            "protein_exclude_q": None,
            "family_q": Q(slug__startswith='007') | Q(slug__startswith='008'),
            "family_exclude_q": None,
        },
    }

    # Class A receptor families are split across two circles purely to keep the wheel visually
    # balanced (~half the families per circle) -- not derived from any biological or slug-based
    # grouping. Re-check by eye against the rendered wheel any time Class A's family list changes
    # (e.g. after a slug/family reclassification).
    _CLASS_A_FAMILIES_IN_CIRCLE_1 = 43

    # Class O2 odorant family numbers (parsed out of "Olfactory family N" display names) are
    # bucketed into 3 circles purely for wheel layout balance. Assumes family numbers currently
    # run 1-14 with no gaps -- a family number outside these ranges is silently dropped from the
    # wheel rather than erroring (see _bucket_odorant_circles).
    _ODORANT_O2_CIRCLE_RANGES = (
        ("Circle_1", range(1, 5)),
        ("Circle_2", range(5, 10)),
        ("Circle_3", range(10, 15)),
    )

    _GPCROME_CLASS_RENAME_MAP = {
        "Class A (Rhodopsin)": "A",
        "Class B1 (Secretin)": "B1",
        "Class B2 (Adhesion)": "B2",
        "Class C (Glutamate)": "C",
        "Class F (Frizzled)": "F",
        "Class O1 (Olfactory/extra-nasal 1)": "O1",
        "Class O2 (Olfactory/extra-nasal 2)": "O2",
        "Class T2 (Taste 2)": "T2",
        "Class V (Vomeronasal)": "V",
        "Unclassified": "U",
    }

    @staticmethod
    def _load_gpcrome_query_data(filters):
        proteins_qs = Protein.objects.filter(
            species_id=1,
            parent_id__isnull=True,
            accession__isnull=False,
        ).filter(filters["protein_q"])
        if filters["protein_exclude_q"] is not None:
            proteins_qs = proteins_qs.exclude(filters["protein_exclude_q"])
        proteins_qs = proteins_qs.prefetch_related('genes').order_by('entry_name')

        families = ProteinFamily.objects.filter(filters["family_q"])
        if filters["family_exclude_q"] is not None:
            families = families.exclude(filters["family_exclude_q"])

        # Single pass over the protein queryset (plus its genes prefetch) builds everything the
        # tree-building step needs -- no second round-trip re-fetching rows already in hand.
        valid_names = set()
        name_to_entry = {}
        entry_to_gene = {}
        for protein in proteins_qs:
            valid_names.add(protein.name)
            name_to_entry[protein.name] = protein.entry_name
            entry_to_gene[protein.entry_name] = next(
                (g.name for g in protein.genes.all() if g.position == 0), None
            )

        return families, valid_names, name_to_entry, entry_to_gene

    @staticmethod
    def _build_gpcrome_family_tree(families, name_to_entry, entry_to_gene):
        datatree = {}
        slug_to_name = {fam.slug: fam.name for fam in families}

        for item in families:
            slug_parts = item.slug.split('_')
            current_level = datatree

            for i in range(len(slug_parts)):
                full_slug = '_'.join(slug_parts[:i + 1])
                name = item.name if full_slug == item.slug else None

                if i == len(slug_parts) - 1:
                    if len(slug_parts) == 4:
                        if item.name in name_to_entry:
                            entry_name = name_to_entry[item.name]
                            entry_code = entry_name.split('_')[0].upper()
                            gene_symbol = entry_to_gene.get(entry_name)  # Entrez name, if available
                            current_level[item.name] = {
                                "Data": "Empty",
                                "EntryName": entry_code,
                                "Entrez": gene_symbol if gene_symbol else "UNKNOWN",
                                "Color": "#FFFFFF"
                            }
                    else:
                        current_level.setdefault(name, {})
                else:
                    if name:
                        current_level = current_level.setdefault(name, {})
                    else:
                        current_level = current_level.setdefault(full_slug, {})  # Temp key

        def convert_keys(tree):
            new_tree = {}
            for key, value in tree.items():
                new_key = slug_to_name.get(key, key)
                new_tree[new_key] = convert_keys(value) if isinstance(value, dict) else value
            return new_tree

        return convert_keys(datatree)

    @staticmethod
    def _prune_gpcrome_tree(tree, valid_leaves):
        pruned = {}
        for key, value in tree.items():
            if isinstance(value, dict):
                if 'Data' in value:
                    if key in valid_leaves:
                        pruned[key] = value
                else:
                    pruned_subtree = DataMapperHome._prune_gpcrome_tree(value, valid_leaves)
                    if pruned_subtree:
                        pruned[key] = pruned_subtree
        return pruned if pruned else None

    @staticmethod
    def _bucket_classic_circles(GPCRomeStructureDict, class_rename_map):
        GPCRome_dict = {
            "Circle_1": {},
            "Circle_2": {},
            "Circle_3": {},
            "Circle_4": {},
            "Circle_5": {}
        }

        class_A_receptor_families = 0

        for Class, ligand_types in GPCRomeStructureDict.items():
            renamed_class = class_rename_map.get(Class, Class)

            if Class == "Class A (Rhodopsin)":
                sorted_receptor_families = []

                for Ligand_type, receptor_families in ligand_types.items():
                    if Ligand_type == "Orphan receptors":
                        continue

                    for Receptor_Family in receptor_families:
                        sorted_receptor_families.append(
                            (Ligand_type, Receptor_Family, receptor_families[Receptor_Family])
                        )

                sorted_receptor_families.sort(key=lambda x: x[1])  # sort by receptor family name

                for Ligand_type, Receptor_Family, receptors in sorted_receptor_families:
                    target_circle = (
                        "Circle_1"
                        if class_A_receptor_families < DataMapperHome._CLASS_A_FAMILIES_IN_CIRCLE_1
                        else "Circle_2"
                    )

                    GPCRome_dict.setdefault(target_circle, {}).setdefault(renamed_class, {}).setdefault(Ligand_type, {})[Receptor_Family] = receptors
                    class_A_receptor_families += 1

                # Add Class A orphans to Circle_2 if present. The underlying ProteinFamily tree
                # nests these under a receptor-family node named "Orphan receptors" (not the
                # older "Class A orphans" name this lookup used to expect -- that name no longer
                # exists in the data, which silently dropped all ~80 Class A orphan receptors).
                # The output key is still relabelled "Class A orphans" since datamapper.js's
                # buildContentItems looks for that literal string to sort orphans to the end of
                # Class A's family list.
                orphans = ligand_types.get("Orphan receptors", {})
                if "Orphan receptors" in orphans:
                    GPCRome_dict["Circle_2"].setdefault(renamed_class, {}).setdefault("Orphan receptors", {})["Class A orphans"] = orphans["Orphan receptors"]

            elif Class in ["Class B1 (Secretin)", "Class B2 (Adhesion)"]:
                GPCRome_dict["Circle_3"].setdefault(renamed_class, {}).update(ligand_types)

            elif Class in ["Class C (Glutamate)", "Class F (Frizzled)"]:
                GPCRome_dict["Circle_4"].setdefault(renamed_class, {}).update(ligand_types)

            elif Class in ["Class T2 (Taste 2)", "Class V (Vomeronasal)", "Unclassified"]:
                GPCRome_dict["Circle_5"].setdefault(renamed_class, {}).update(ligand_types)

        # Remove the ligand type layer
        for circle in GPCRome_dict:
            for class_name in list(GPCRome_dict[circle].keys()):
                new_structure = {}

                for ligand_type in list(GPCRome_dict[circle][class_name].keys()):
                    for receptor_family, receptors in GPCRome_dict[circle][class_name][ligand_type].items():
                        new_structure[receptor_family] = receptors

                # Replace the old structure with the flattened one
                GPCRome_dict[circle][class_name] = new_structure

        # Final return (renamed and sorted into circles)
        return {
            "Data": GPCRome_dict
        }

    @staticmethod
    def _bucket_odorant_circles(GPCRomeStructureDict, class_rename_map):
        GPCRome_dict = {
            "Circle_1": {},  # Family 1-4 from Class O2
            "Circle_2": {},  # Family 5-9 from Class O2
            "Circle_3": {},  # Family 10-14 from Class O2
            "Circle_4": {}   # All of Class O1
        }

        for Class, ligand_types in GPCRomeStructureDict.items():
            renamed_class = class_rename_map.get(Class, Class)

            if Class == "Class O2 (Olfactory/extra-nasal 2)":
                sorted_families = []

                for Ligand_type, receptor_families in ligand_types.items():
                    for Receptor_Family, receptors in receptor_families.items():
                        try:
                            family_number = int(Receptor_Family.replace("Olfactory family", "").strip())
                        except ValueError:
                            continue

                        sorted_families.append(
                            (family_number, Ligand_type, Receptor_Family, receptors)
                        )

                sorted_families.sort(key=lambda x: x[0])  # Sort numerically by family number

                for family_number, Ligand_type, Receptor_Family, receptors in sorted_families:
                    renamed_family = f"Family {family_number}"
                    circle = next(
                        (name for name, rng in DataMapperHome._ODORANT_O2_CIRCLE_RANGES if family_number in rng),
                        None,
                    )
                    if circle is None:
                        continue  # family number outside the known ranges -- see _ODORANT_O2_CIRCLE_RANGES

                    GPCRome_dict.setdefault(circle, {}).setdefault(renamed_class, {}).setdefault(Ligand_type, {})[renamed_family] = receptors

            elif Class == "Class O1 (Olfactory/extra-nasal 1)":
                renamed_ligand_types = {}

                for Ligand_type, receptor_families in ligand_types.items():
                    renamed_receptor_families = {}

                    for Receptor_Family, receptors in receptor_families.items():
                        try:
                            family_number = int(Receptor_Family.replace("Olfactory family", "").strip())
                            renamed_family = f"Family {family_number}"
                        except ValueError:
                            renamed_family = Receptor_Family  # fallback to original if parsing fails

                        renamed_receptor_families[renamed_family] = receptors

                    renamed_ligand_types[Ligand_type] = renamed_receptor_families

                GPCRome_dict["Circle_4"].setdefault(renamed_class, {}).update(renamed_ligand_types)

        # Flatten ligand_type layer (remove ligand type level)
        for circle in GPCRome_dict:
            for class_name in list(GPCRome_dict[circle].keys()):
                new_structure = {}
                for ligand_type in GPCRome_dict[circle][class_name]:
                    for receptor_family, receptors in GPCRome_dict[circle][class_name][ligand_type].items():
                        new_structure[receptor_family] = receptors
                GPCRome_dict[circle][class_name] = new_structure

        return {
            "Data": GPCRome_dict
        }

    @staticmethod
    def GenerateGPCRomeDataStructure(data_type: str = "Classic"):
        # Reference data (Protein/ProteinFamily) only changes on scheduled data releases, and this
        # skeleton is otherwise rebuilt uncached on every page render (sometimes twice per request,
        # e.g. Classic + Odorant back-to-back) -- a week-long cache is the same low-ceremony pattern
        # already used elsewhere in this codebase for slow-changing reference lookups. Bump the
        # "v1" suffix if the skeleton shape ever changes, as a manual invalidation lever.
        cache_key = f"gpcrome_data_structure_v2_{data_type}"
        cached = cache.get(cache_key)
        if cached is not None:
            return deepcopy(cached)  # callers mutate the returned dict in place -- never hand out the cached object itself

        filters = DataMapperHome._GPCROME_QUERY_FILTERS.get(data_type)
        if filters is None:
            raise ValueError(f"Unsupported data structure type: {data_type}")

        families, valid_names, name_to_entry, entry_to_gene = DataMapperHome._load_gpcrome_query_data(filters)

        datatree = DataMapperHome._build_gpcrome_family_tree(families, name_to_entry, entry_to_gene)
        GPCRomeStructureDict = DataMapperHome._prune_gpcrome_tree(datatree, valid_names)

        if not GPCRomeStructureDict:
            return None

        if data_type == "Classic":
            result = DataMapperHome._bucket_classic_circles(GPCRomeStructureDict, DataMapperHome._GPCROME_CLASS_RENAME_MAP)
        else:
            result = DataMapperHome._bucket_odorant_circles(GPCRomeStructureDict, DataMapperHome._GPCROME_CLASS_RENAME_MAP)

        cache.set(cache_key, result, 60 * 60 * 24 * 7)
        return deepcopy(result)

    @staticmethod
    def update_nested_GPCRome_data(structure_dict, raw_data):
        def normalize_key(raw_key):
            return raw_key.split("_")[0].upper()

        # Normalize and extract only non-null values
        normalized_raw_data = {
            normalize_key(key): {
                'Data': val.get("Value1"),
                'Color': val.get("Value2")  # Might be None
            }
            for key, val in raw_data.items()
            if isinstance(val, dict) and "Value1" in val
        }

        def recursive_update(d):
            if isinstance(d, dict):
                if "EntryName" in d and "Data" in d:
                    entry = d["EntryName"].upper()
                    if entry in normalized_raw_data:
                        d["Data"] = normalized_raw_data[entry]["Data"]
                        color = normalized_raw_data[entry].get("Color")
                        if color is not None:
                            d["Color"] = color

                for value in d.values():
                    recursive_update(value)
            elif isinstance(d, list):
                for item in d:
                    recursive_update(item)

        recursive_update(structure_dict)
        return structure_dict

    @staticmethod
    def build_gpcrome_receptor_normalization_maps(include_odorant=False):
        """Build lookup tables for receptor normalization.

        Args:
            include_odorant: If True, include odorant receptors (family slugs 007, 008)
                             in the picker / autocomplete / resolve set.
                             The GPCRome Wheel page leaves this False; all other Mapper
                             pages pass True.
        """
        all_proteins = Protein.objects.filter(
            species_id=1,
            parent_id__isnull=True,
            accession__isnull=False,
            family_id__slug__startswith='0',
        ).values_list('entry_name', flat=True).distinct()

        if include_odorant:
            # All human GPCRs including odorant families 007 / 008
            proteins_gpcrome_tree = set(all_proteins)
        else:
            proteins_gpcrome_tree = set(
                Protein.objects.filter(
                    species_id=1,
                    parent_id__isnull=True,
                    accession__isnull=False,
                    family_id__slug__startswith='0',
                ).exclude(
                    family_id__slug__startswith='007',
                ).exclude(
                    family_id__slug__startswith='008',
                ).values_list('entry_name', flat=True).distinct()
            )

        proteins = Protein.objects.prefetch_related('genes').filter(entry_name__in=all_proteins)
        entry_to_gene = {
            protein.entry_name: next((g.name for g in protein.genes.all() if g.position == 0), None)
            for protein in proteins
        }
        gene_to_entry = {v.upper(): k for k, v in entry_to_gene.items() if v}
        entry_name_upper_to_entry = {k.upper(): k for k in entry_to_gene.keys()}
        entry_name_no_species_to_entry = {
            k.split('_')[0].upper(): k for k in entry_to_gene.keys()
        }
        iuphar_name_to_entry = {}
        for p in Protein.objects.filter(entry_name__in=proteins_gpcrome_tree).only('entry_name', 'name'):
            if p.name:
                plain = re.sub(r'<[^>]+>', '', p.name).strip()
                if plain:
                    iuphar_name_to_entry[plain.upper()] = p.entry_name
        return {
            'proteins_gpcrome_tree': proteins_gpcrome_tree,
            'gene_to_entry': gene_to_entry,
            'entry_name_upper_to_entry': entry_name_upper_to_entry,
            'entry_name_no_species_to_entry': entry_name_no_species_to_entry,
            'entry_to_gene': entry_to_gene,
            'iuphar_name_to_entry': iuphar_name_to_entry,
        }

    @staticmethod
    def gpcrome_receptor_select2_options(maps=None):
        if maps is None:
            maps = DataMapperHome.build_gpcrome_receptor_normalization_maps()
        entry_to_gene = maps['entry_to_gene']
        tree = sorted(maps['proteins_gpcrome_tree'])
        name_by_entry = {
            p['entry_name']: p['name']
            for p in Protein.objects.filter(entry_name__in=tree).values('entry_name', 'name')
        }
        options = []
        for entry_name in tree:
            gene = entry_to_gene.get(entry_name) or ''
            raw_name = name_by_entry.get(entry_name) or ''
            plain_name = re.sub(r'<[^>]+>', '', raw_name).strip() if raw_name else ''
            base = entry_name.split('_')[0].upper() if '_' in entry_name else entry_name.upper()
            parts = [plain_name, gene, base]
            label = ' | '.join([x for x in parts if x])
            if not label:
                label = entry_name
            search_parts = [
                entry_name,
                entry_name.upper(),
                base,
                gene.upper() if gene else '',
                plain_name.upper() if plain_name else '',
            ]
            search_text = ' '.join([x for x in search_parts if x]).upper()
            options.append({
                'id': entry_name,
                'text': label,
                'search_text': search_text,
                'name_plain': plain_name,
                'name_html': raw_name if raw_name else plain_name,
                'gene': gene,
                'uniprot': base,
            })
        return options

    @staticmethod
    def gpcrome_receptor_client_resolve_map(maps=None):
        """Uppercased lookup keys -> canonical entry_name (same rules as normalize where possible)."""
        if maps is None:
            maps = DataMapperHome.build_gpcrome_receptor_normalization_maps()
        tree = maps['proteins_gpcrome_tree']
        out = {}
        for en in tree:
            out[en.upper()] = en
            if '_' in en:
                out[en.split('_')[0].upper()] = en
        for gene_upper, en in maps['gene_to_entry'].items():
            if en in tree:
                out[gene_upper] = en
        for name_upper, en in (maps.get('iuphar_name_to_entry') or {}).items():
            if en in tree:
                out[name_upper] = en
        return out

    @staticmethod
    def _gpcrome_strip_markup(s):
        if not s:
            return ''
        return re.sub(r'<[^>]+>', '', str(s)).strip()

    # Same lineage idea as Drugs browser (target__family__parent__...) — prefetch enough FK hops
    # that in-memory traversal matches Protein.get_protein_class / get_protein_family without new queries.
    _GPCROME_FAM_PARENT_CHAIN = 'family__' + '__'.join(['parent'] * 14)

    @staticmethod
    def _gpcrome_class_display_name_after_family_traversal(protein_family):
        """Traversal copy of Protein.get_protein_class; uses cached FK links from select_related()."""
        if protein_family is None:
            return ''
        tmp = protein_family
        while tmp.parent is not None and tmp.parent.parent is not None:
            tmp = tmp.parent
        return tmp.name if tmp.name else ''

    @staticmethod
    def _gpcrome_receptor_family_display_name_after_traversal(protein_family):
        """Receptor family = one hop up from Protein.family (matches family__parent__name)."""
        if protein_family is None or protein_family.parent is None:
            return ''
        return protein_family.parent.name or ''

    @staticmethod
    def _gpcrome_ligand_type_display_name_after_family_traversal(protein_family):
        """Ligand type / chemotype = two hops up from Protein.family (matches family__parent__parent__name)."""
        if (
            protein_family is None
            or protein_family.parent is None
            or protein_family.parent.parent is None
        ):
            return ''
        return protein_family.parent.parent.name or ''

    @staticmethod
    def gpcrome_collapse_nonhuman_only(nonhuman_only):
        """One representative entry per stem (first species encountered), matching the
        receptor input table's collapsed-species convention: a single pick per receptor
        labeled "STEM (no human ortholog)" instead of one row per species suffixed
        "(Mouse only)" / "(Rat only)". Returns (extra_entry_names, stem_by_entry) for
        use with gpcrome_receptor_picker_table_rows.
        """
        seen_stems = set()
        extra_entry_names = []
        stem_by_entry = {}
        for info in nonhuman_only:
            stem = info['stem']
            if stem in seen_stems:
                continue
            seen_stems.add(stem)
            extra_entry_names.append(info['entry'])
            stem_by_entry[info['entry']] = stem.upper()
        return extra_entry_names, stem_by_entry

    @staticmethod
    def gpcrome_receptor_picker_table_rows(maps=None, extra_entry_names=None, entry_label_overrides=None):
        """Plain JSON rows for GPCR picker (wheel receptor set). Single batched Protein query.

        extra_entry_names: optional entry_names to include in addition to the human
            gpcrome tree (e.g. non-human-ortholog-only receptors for Tree/Heatmap/List).
        entry_label_overrides: optional {entry_name: STEM_UPPER} — for those entries the
            GtoPdb name column is fully replaced with "STEM (no human ortholog)", matching
            the receptor input table's collapsed-species convention.
        """
        if maps is None:
            maps = DataMapperHome.build_gpcrome_receptor_normalization_maps()
        entry_to_gene = dict(maps['entry_to_gene'])
        tree = set(maps['proteins_gpcrome_tree'])
        if extra_entry_names:
            tree |= set(extra_entry_names)
        tree = sorted(tree)
        if not tree:
            return []

        missing_gene_entries = [e for e in tree if e not in entry_to_gene]
        if missing_gene_entries:
            for protein in Protein.objects.prefetch_related('genes').filter(
                entry_name__in=missing_gene_entries
            ):
                entry_to_gene[protein.entry_name] = next(
                    (g.name for g in protein.genes.all() if g.position == 0), None
                )

        entry_label_overrides = entry_label_overrides or {}
        NHO_TAG = '(no human ortholog)'

        qs = Protein.objects.filter(entry_name__in=tree).select_related(
            DataMapperHome._GPCROME_FAM_PARENT_CHAIN
        )
        by_entry = {p.entry_name: p for p in qs}

        rows = []
        for entry_name in tree:
            p = by_entry.get(entry_name)
            if not p:
                continue
            gene = entry_to_gene.get(entry_name) or ''
            raw_name = p.name or ''
            plain_name = DataMapperHome._gpcrome_strip_markup(raw_name)

            override_stem = entry_label_overrides.get(entry_name)
            if override_stem:
                plain_name = '{} {}'.format(override_stem, NHO_TAG)
                raw_name = '{} <em>{}</em>'.format(override_stem, NHO_TAG)

            rfam = p.family
            family = DataMapperHome._gpcrome_strip_markup(
                DataMapperHome._gpcrome_receptor_family_display_name_after_traversal(rfam)
            )
            ligand_type = DataMapperHome._gpcrome_strip_markup(
                DataMapperHome._gpcrome_ligand_type_display_name_after_family_traversal(rfam)
            )
            prot_class = DataMapperHome._gpcrome_strip_markup(
                DataMapperHome._gpcrome_class_display_name_after_family_traversal(rfam)
            )

            entry_short = p.entry_short()
            gpcrdb_link = 'https://gpcrdb.org/protein/{}/'.format(entry_name)
            uniprot_link = 'https://www.uniprot.org/uniprot/{}'.format(entry_short) if entry_short else ''
            rows.append({
                'id': entry_name,
                'name_html': raw_name if raw_name else plain_name,
                'name_plain': plain_name or entry_short,
                'gene': gene or '',
                'uniprot': entry_short,
                'uniprot_link': uniprot_link,
                'family': family or '',
                'ligandtype': ligand_type or '',
                'class': prot_class or '',
                'gpcrdb_link': gpcrdb_link,
            })
        return rows

    @staticmethod
    def gpcrome_receptor_info_for_list(maps=None, extra_entry_names=None, entry_label_overrides=None):
        """Returns dict of entry_name → {class, ligandtype, family, name_plain, gene, uniprot}
        for all human non-odorant GPCRs, for building list plot data in the browser.

        extra_entry_names: optional entry_names to include in addition to the human
            gpcrome tree (e.g. non-human-ortholog-only receptors). Without an entry here,
            mapperListBuildData() has no class/ligandtype/family to bucket the receptor
            under and silently drops it from the plot.
        entry_label_overrides: optional {entry_name: STEM_UPPER}, mirrors
            gpcrome_receptor_picker_table_rows — replaces name_plain with
            "STEM (no human ortholog)" for those entries.
        """
        if maps is None:
            maps = DataMapperHome.build_gpcrome_receptor_normalization_maps()
        gpcrome_set = set(maps.get('proteins_gpcrome_tree', set()))
        if extra_entry_names:
            gpcrome_set |= set(extra_entry_names)
        select2_opts = {o['id']: o for o in DataMapperHome.gpcrome_receptor_select2_options(maps=maps)}
        entry_label_overrides = entry_label_overrides or {}
        NHO_TAG = '(no human ortholog)'

        entry_to_gene = dict(maps['entry_to_gene'])
        entry_short_by_name = {}
        missing_gene_entries = [e for e in gpcrome_set if e not in entry_to_gene]
        if missing_gene_entries:
            for protein in Protein.objects.prefetch_related('genes').filter(
                entry_name__in=missing_gene_entries
            ):
                entry_to_gene[protein.entry_name] = next(
                    (g.name for g in protein.genes.all() if g.position == 0), None
                )
                entry_short_by_name[protein.entry_name] = protein.entry_short()

        qs = Protein.objects.filter(
            entry_name__in=list(gpcrome_set)
        ).values_list(
            'entry_name',
            'family__parent__parent__parent__name',  # Class
            'family__parent__parent__name',           # Ligand type
            'family__parent__name',                   # Receptor family
        )

        result = {}
        for entry_name, cls, ligandtype, family in qs:
            opt = select2_opts.get(entry_name, {})
            override_stem = entry_label_overrides.get(entry_name)
            name_plain = '{} {}'.format(override_stem, NHO_TAG) if override_stem else opt.get('name_plain', '')
            result[entry_name] = {
                'class':      cls        or 'Other',
                'ligandtype': ligandtype or 'Other',
                'family':     family     or 'Other',
                'name_plain': name_plain,
                'gene':       opt.get('gene') or entry_to_gene.get(entry_name) or '',
                'uniprot':    opt.get('uniprot') or entry_short_by_name.get(entry_name, ''),
            }
        return result

    @staticmethod
    def build_ortholog_species_map():
        """
        Groups all SWISSPROT GPCR proteins by family to build a per-receptor
        species availability map for the Tree page species feature.
        Returns dict with:
          by_stem:          human_stem -> list of ortholog dicts (all species incl. human)
          nonhuman_only:    list of receptor dicts that have no human equivalent
          entry_to_species: entry_name -> {common, latin, label, is_human}
        """
        proteins = list(
            Protein.objects.filter(
                source__name='SWISSPROT',
                accession__isnull=False,
                parent_id__isnull=True,
                family__slug__startswith='0',
            ).values(
                'entry_name', 'family_id', 'species_id',
                'species__common_name', 'species__latin_name',
            ).order_by('family_id', 'species_id')
        )

        family_map = {}
        for p in proteins:
            fid = p['family_id']
            if fid not in family_map:
                family_map[fid] = []
            family_map[fid].append(p)

        by_stem = {}
        nonhuman_only = []
        entry_to_species = {}

        for _family_id, members in family_map.items():
            human_members = [m for m in members if m['species_id'] == 1]
            member_infos = []
            for m in members:
                common = (m['species__common_name'] or '').strip() or (m['species__latin_name'] or '').strip()
                latin = (m['species__latin_name'] or '').strip()
                if common and latin and common != latin:
                    label = '{} ({})'.format(common, latin)
                elif latin:
                    label = latin
                else:
                    label = m['entry_name']
                is_human = (m['species_id'] == 1)
                stem = m['entry_name'].split('_')[0]
                info = {
                    'entry': m['entry_name'],
                    'stem': stem,
                    'common': common,
                    'latin': latin,
                    'label': label,
                    'is_human': is_human,
                }
                member_infos.append(info)
                entry_to_species[m['entry_name']] = {
                    'common': common,
                    'latin': latin,
                    'label': label,
                    'is_human': is_human,
                }

            if human_members:
                human_stem = human_members[0]['entry_name'].split('_')[0]
                if human_stem not in by_stem:
                    by_stem[human_stem] = member_infos
            else:
                for info in member_infos:
                    info['label_suffix'] = '({} only)'.format(
                        info['common'] or info['latin'] or 'Unknown'
                    )
                    nonhuman_only.append(info)

        return {
            'by_stem': by_stem,
            'nonhuman_only': nonhuman_only,
            'entry_to_species': entry_to_species,
        }

    # Classification-tree config. Fully slug-driven: any ProteinFamily that is a direct child of
    # the root family ('000') is picked up as a top-level GPCR class automatically -- no hardcoded
    # per-class ProteinFamily.objects.get(name=...) lookups needed, so a new/renamed class needs no
    # code change here (unlike the old generate_tree_plot, which crashed once "Other GPCRs" was
    # renamed to "Unclassified").
    _CLASSIFICATION_ROOT_FAMILY_SLUG = '000'
    _CLASSIFICATION_TREE_DEPTH = 4  # class / ligand type / receptor family / receptor
    # Class D1 (slug '005') has exactly one SWISSPROT protein and no classification annotations
    # yet -- excluded until there's real data behind it.
    _CLASSIFICATION_EXCLUDED_CLASS_SLUGS = frozenset({'005'})
    _CLASSIFICATION_TREE_CACHE_KEY = 'mapper_classification_tree_skeleton_v1'

    @staticmethod
    def _load_classification_tree_data():
        """(families, leaf_by_family_slug) for every top-level GPCR class. families is a list of
        (slug, name) for every ProteinFamily at or below a non-excluded class, up to the tree's
        depth. leaf_by_family_slug maps a depth-4 family slug to the entry-name stem of its
        representative human SWISSPROT protein (one leaf per receptor family)."""
        class_slugs = set(
            ProteinFamily.objects
                .filter(parent__slug=DataMapperHome._CLASSIFICATION_ROOT_FAMILY_SLUG)
                .exclude(slug__in=DataMapperHome._CLASSIFICATION_EXCLUDED_CLASS_SLUGS)
                .values_list('slug', flat=True)
        )

        families = [
            (slug, name)
            for slug, name in ProteinFamily.objects
                .exclude(slug=DataMapperHome._CLASSIFICATION_ROOT_FAMILY_SLUG)
                .values_list('slug', 'name')
            if slug.split('_')[0] in class_slugs and len(slug.split('_')) <= DataMapperHome._CLASSIFICATION_TREE_DEPTH
        ]

        leaf_by_family = {}
        for fam_slug, entry_name in (
            Protein.objects
                .filter(source__name='SWISSPROT')
                .order_by('family__slug', 'species_id')
                .values_list('family__slug', 'entry_name')
        ):
            if fam_slug.split('_')[0] in class_slugs:
                leaf_by_family.setdefault(fam_slug, entry_name.split('_')[0])
        return families, leaf_by_family

    @staticmethod
    def _classification_node_label(name, level):
        # Ligand-type and receptor-family names carry markup/suffixes not meant for display;
        # class names and leaf (protein stem) names are used verbatim.
        if level in (2, 3):
            return name.replace('receptors', '').replace('<sub>', ' ').replace('</sub>', '').strip()
        return name

    @staticmethod
    def _build_classification_tree(families, leaf_by_family):
        def make_node(label, value=0):
            return OrderedDict([('name', label), ('value', value), ('color', ''), ('children', [])])

        root = make_node('', value=3000)
        nodes = {}
        # Sorting by slug guarantees a parent is emitted before its children, and gives a
        # deterministic (if cosmetic) sibling order with no curated list to maintain.
        for slug, fam_name in sorted(families, key=lambda f: f[0]):
            parts = slug.split('_')
            level = len(parts)
            parent = root if level == 1 else nodes.get('_'.join(parts[:-1]))
            if parent is None:
                continue  # parent excluded/absent -> drop the orphan subtree
            if level == DataMapperHome._CLASSIFICATION_TREE_DEPTH:
                leaf_name = leaf_by_family.get(slug)
                if not leaf_name:
                    continue  # family with no SWISSPROT member -> no leaf
                parent['children'].append(make_node(leaf_name))
            else:
                node = make_node(DataMapperHome._classification_node_label(fam_name, level))
                nodes[slug] = node
                parent['children'].append(node)
        return root

    @staticmethod
    def _prune_classification_tree(node, level=0):
        """Drop interior nodes left with no receptor leaves (e.g. a receptor-family branch whose
        only member has no SWISSPROT protein). Depth-gated rather than "no children -> leaf" so a
        genuinely childless interior node is distinguished from a real receptor leaf."""
        if level >= DataMapperHome._CLASSIFICATION_TREE_DEPTH:
            return node
        kept = [
            pruned for pruned in (
                DataMapperHome._prune_classification_tree(child, level + 1)
                for child in node['children']
            )
            if pruned is not None
        ]
        if not kept:
            return None
        node['children'] = kept
        return node

    @staticmethod
    def _classification_tree_skeleton():
        cached = cache.get(DataMapperHome._CLASSIFICATION_TREE_CACHE_KEY)
        if cached is not None:
            return deepcopy(cached)  # callers filter/mutate -- never hand out the cached object
        families, leaf_by_family = DataMapperHome._load_classification_tree_data()
        tree = DataMapperHome._prune_classification_tree(
            DataMapperHome._build_classification_tree(families, leaf_by_family)
        ) or OrderedDict([('name', ''), ('value', 3000), ('color', ''), ('children', [])])
        cache.set(DataMapperHome._CLASSIFICATION_TREE_CACHE_KEY, tree, 60 * 60 * 24 * 7)
        return deepcopy(tree)

    @staticmethod
    def _classification_tree_max_depth(node, level=0):
        if not node.get('children'):
            return level
        return max(DataMapperHome._classification_tree_max_depth(c, level + 1) for c in node['children'])

    @staticmethod
    def _classification_tree_options(tree, depth, anchor="tree_plot"):
        # branch_length holds, per level, the longest label actually present in the tree being
        # returned -- used by the renderers purely as a text-width proxy for ring spacing, so it
        # only needs to track real content rather than a hand-maintained literal.
        branch_length = {lvl: '' for lvl in range(1, depth + 1)}

        def walk(node, level):
            if 1 <= level < depth:
                name = node.get('name', '')
                if len(name) > len(branch_length[level]):
                    branch_length[level] = name
            for child in (node.get('children') or []):
                walk(child, level + 1)

        walk(tree, 0)
        return {
            'depth': depth,
            'branch_length': branch_length,
            'branch_trunc': 0,
            'leaf_offset': 30,
            'anchor': anchor,
            'label_free': [],
            'fontSize': {
                'class': "15px",
                'ligandtype': "14px",
                'receptorfamily': "13px",
                'receptor': "12px"
            }
        }

    @staticmethod
    def GenerateClassificationTreeData(input_data): #ADD AN INPUT FILTER DICTIONARY
        """
        Full GPCR classification skeleton (class / ligand type / receptor family / receptor) as a
        d3 list-tree (OrderedDict with name/value/color/children at every node), optionally
        filtered down to `input_data`'s receptors. `color` is left blank on every node -- neither
        consumer needs a backend-supplied color (Mapper_Tree.html's mapper_classification_tree.js
        already recomputes colors client-side; drugged_gpcrome.html's datamapper.js gets a small
        client-side fallback palette instead, see static/home/js/datamapper.js).
        """
        master_dict = DataMapperHome._classification_tree_skeleton()

        updated_data = {key.replace('_human', ''): value for key, value in input_data.items()}
        circles = {
            key.replace('_human', '').upper(): {k: v for k, v in value.items()}
            for key, value in input_data.items()
        }

        # Empty payloads: Mapper 2.0 loads the full skeleton client-side — do not prune the tree yet.
        if input_data:
            master_dict = DataMapperHome.keep_by_names(master_dict, updated_data)

        if isinstance(master_dict, dict) and master_dict.get('children') is not None and len(master_dict['children']) == 1:
            master_dict = master_dict['children'][0]

        depth = (
            DataMapperHome._classification_tree_max_depth(master_dict)
            if isinstance(master_dict, dict)
            else DataMapperHome._CLASSIFICATION_TREE_DEPTH
        )
        general_options = DataMapperHome._classification_tree_options(
            master_dict if isinstance(master_dict, dict) else {'children': []}, depth
        )

        entry_names_list = list(input_data.keys())  # or: input_data.keys()

        # Step 2: Trim protein queries
        whole_receptors = Protein.objects.filter(entry_name__in=entry_names_list).prefetch_related("family", "family__parent__parent__parent")

        protein_genes = Protein.objects.filter(entry_name__in=entry_names_list).prefetch_related('genes')

        entry_to_gene = {
            protein.entry_name: next((g.name for g in protein.genes.all() if g.position == 0), None)
            for protein in protein_genes
        }

        whole_rec_dict = {}
        entrez_label_dict = {}
        for rec in whole_receptors:
            rec_entryname = rec.entry_name
            rec_uniprot = rec.entry_short()
            rec_iuphar = rec.family.name.replace("receptor", '').replace("<i>", "").replace("</i>", "").strip()
            if (rec_iuphar[0].isupper()) or (rec_iuphar[0].isdigit()):
                whole_rec_dict[rec_uniprot] = [rec_iuphar]
            else:
                whole_rec_dict[rec_uniprot] = [rec_iuphar.capitalize()]
            # Entrez label — safely
            gene_name = entry_to_gene.get(rec_entryname)
            if gene_name:
                entrez_label_dict[rec_uniprot] = [gene_name]

        return master_dict, general_options, circles, whole_rec_dict, entrez_label_dict

    @staticmethod
    def recompute_position_layout(data):
        """
        Re-run t-SNE from user-submitted receptor "Position" values — the live
        recompute POSTed by ClusterRender.post whenever every row in the Position
        column is filled in. Class/Ligand type/Receptor family aren't returned —
        the client already holds that metadata (window.MAPPER_CLUSTER_ALL_POSITIONS)
        and merges it back in via getClusterPosMap() (mapper_cluster_page.js).
        """
        labels = [entry_name.replace('_human', '') for entry_name in data.keys()]
        positions = pd.Series(
            [row.get('Value2') for row in data.values()],
            index=labels,
            dtype=float,
        )
        positions = positions.fillna(positions.mean())

        distance_matrix = np.abs(positions.values.reshape(-1, 1) - positions.values.reshape(1, -1))

        n_points = distance_matrix.shape[0]
        perplexity = max(2, min(n_points / 5, 50))
        tsne = TSNE(n_components=2, metric='precomputed', init='random', random_state=42, perplexity=perplexity)
        coords = tsne.fit_transform(distance_matrix)
        clusters = KMeans(n_clusters=5, random_state=42).fit_predict(coords)

        result = pd.DataFrame(coords, columns=['x', 'y'])
        result['cluster'] = clusters
        result['label'] = positions.index
        return result.to_json(orient='records')

    # Global sequence-similarity t-SNE layout for all human GPCRs, sourced from
    # classification.ClusterCoord (built offline by build_clustercoord from ReceptorSimilarity)
    @staticmethod
    def generate_full_matrix():
        rows = (
            ClusterCoord.objects
            .filter(
                dataset_type=ClusterCoord.DATASET_SEQUENCE,
                plot_type=ClusterCoord.PLOT_TSNE,
                group_key='global',
                protein__species_id=1,
            )
            .select_related('protein')
            .order_by('protein__entry_name')
        )
        reduced_df = pd.DataFrame([
            {'label': r.protein.entry_name.replace('_human', ''), 'x': r.x, 'y': r.y}
            for r in rows.iterator()
        ])

        # add class/ligand_type/receptor_family clusters

        # Step 1: Fetch data
        proteins = Protein.objects.filter(
            parent_id__isnull=True, species_id=1
        ).values_list(
            'entry_name',
            "family__parent__parent__parent__name",  # To be renamed as 'Class'
            'family__parent__parent__name',  # To be renamed as 'Ligand type'
            'family__parent__name'  # To be renamed as 'Receptor family'
        )

        # Step 2: Convert to a DataFrame
        proteins_df = pd.DataFrame(list(proteins), columns=['entry_name', 'Class', 'Ligand type', 'Receptor family'])

        # Step 3: Remove '_human' suffix from 'entry_name'
        proteins_df['entry_name'] = proteins_df['entry_name'].str.replace('_human', '')

        # Step 4: Rename 'entry_name' to 'label'
        proteins_df = proteins_df.rename(columns={'entry_name': 'label'})

        # Step 5: Merge with reduced_df on 'label'
        merged_df = pd.merge(reduced_df, proteins_df, on='label', how='left')

        return merged_df

    @staticmethod
    def Label_conversion_info(data):
        # Get list of keys (UniProt-style entry names)
        Name_list = list(data.keys())

        # Fetch IUPHAR names
        names_dict = Protein.objects.filter(entry_name__in=Name_list).values('entry_name', 'name').order_by('entry_name')
        UniProt_to_IUPHAR_converter = {item['entry_name']: item['name'] for item in names_dict}
        IUPHAR_to_UniProt_converter = {item['name']: item['entry_name'] for item in names_dict}

        # Fetch genes with position == 0
        protein_genes = Protein.objects.prefetch_related('genes').filter(entry_name__in=Name_list)
        UniProt_to_Gene_converter = {
            protein.entry_name: next((g.name for g in protein.genes.all() if g.position == 0), None)
            for protein in protein_genes
        }

        # Create IUPHAR to Gene mapping using the two converters above
        IUPHAR_to_Gene_converter = {
            iuphar: UniProt_to_Gene_converter.get(entry)
            for iuphar, entry in IUPHAR_to_UniProt_converter.items()
            if UniProt_to_Gene_converter.get(entry) is not None
        }

        # Final converter dict
        Label_converter = {
            'UniProt_to_IUPHAR_converter': UniProt_to_IUPHAR_converter,
            'IUPHAR_to_UniProt_converter': IUPHAR_to_UniProt_converter,
            'UniProt_to_Gene_converter': UniProt_to_Gene_converter,
            'IUPHAR_to_Gene_converter': IUPHAR_to_Gene_converter
        }

        return Label_converter


class ClusterRender(View):
    def post(self, request, *args, **kwargs):
        Data_json = request.POST.get('Data')

        try:
            # Get data
            Data = json.loads(Data_json)
            # Calculate the plot
            output_seq = DataMapperHome.recompute_position_layout(Data)
            label_converter = DataMapperHome.Label_conversion_info(Data)

            return JsonResponse({
                'cluster_data_seq': json.loads(output_seq),
                'Label_converter': label_converter
            })

        except json.JSONDecodeError:
            return HttpResponse("Invalid JSON data")

class MapperLandingPageView(TemplateView):
    template_name = 'mapper/Mapper_landingPage.html'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        wheel_maps = DataMapperHome.build_gpcrome_receptor_normalization_maps()
        full_maps = DataMapperHome.build_gpcrome_receptor_normalization_maps(include_odorant=True)
        context['wheel_max_receptors'] = len(wheel_maps['proteins_gpcrome_tree'])
        context['full_max_receptors'] = len(full_maps['proteins_gpcrome_tree'])
        return context


class MapperGPCRomeWheelView(TemplateView):
    template_name = 'mapper/Mapper_GPCRomeWheel.html'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        base = DataMapperHome.GenerateGPCRomeDataStructure(data_type="Classic")
        structure = deepcopy(base["Data"])
        merged = DataMapperHome.update_nested_GPCRome_data(structure, {})
        context['GPCRomeData'] = json.dumps(merged)
        context['PlotType'] = 'Numeric'
        maps = DataMapperHome.build_gpcrome_receptor_normalization_maps()
        context['receptor_select2_json'] = json.dumps(
            DataMapperHome.gpcrome_receptor_select2_options(maps=maps)
        )
        context['gpcrome_resolve_json'] = json.dumps(
            DataMapperHome.gpcrome_receptor_client_resolve_map(maps=maps)
        )
        context['gpcrome_picker_rows_json'] = json.dumps(
            DataMapperHome.gpcrome_receptor_picker_table_rows(maps=maps)
        )
        return context


class MapperTreeView(TemplateView):
    template_name = 'mapper/Mapper_Tree.html'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        master_dict, general_options, circles, receptors, genes = DataMapperHome.GenerateClassificationTreeData({})
        context['tree'] = json.dumps(master_dict)
        context['tree_options'] = json.dumps(general_options)
        context['circles'] = json.dumps(circles if circles else {})
        context['Receptor_dict'] = json.dumps(receptors if receptors else {})
        context['Entrez_dict'] = json.dumps(genes if genes else {})
        context['PlotType'] = 'Numeric'
        context['Data'] = json.dumps({})
        maps = DataMapperHome.build_gpcrome_receptor_normalization_maps(include_odorant=True)
        context['receptor_select2_json'] = json.dumps(
            DataMapperHome.gpcrome_receptor_select2_options(maps=maps)
        )
        context['gpcrome_resolve_json'] = json.dumps(
            DataMapperHome.gpcrome_receptor_client_resolve_map(maps=maps)
        )
        ortholog_map = DataMapperHome.build_ortholog_species_map()
        context['ortholog_species_json'] = json.dumps(ortholog_map)
        nonhuman_entries, nonhuman_stem_by_entry = DataMapperHome.gpcrome_collapse_nonhuman_only(
            ortholog_map['nonhuman_only']
        )
        context['gpcrome_picker_rows_json'] = json.dumps(
            DataMapperHome.gpcrome_receptor_picker_table_rows(
                maps=maps, extra_entry_names=nonhuman_entries, entry_label_overrides=nonhuman_stem_by_entry
            )
        )
        return context


class MapperHeatmapView(TemplateView):
    template_name = 'mapper/Mapper_Heatmap.html'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        maps = DataMapperHome.build_gpcrome_receptor_normalization_maps(include_odorant=True)
        context['receptor_select2_json'] = json.dumps(
            DataMapperHome.gpcrome_receptor_select2_options(maps=maps)
        )
        context['gpcrome_resolve_json'] = json.dumps(
            DataMapperHome.gpcrome_receptor_client_resolve_map(maps=maps)
        )
        ortholog_map = DataMapperHome.build_ortholog_species_map()
        context['ortholog_species_json'] = json.dumps(ortholog_map)
        nonhuman_entries, nonhuman_stem_by_entry = DataMapperHome.gpcrome_collapse_nonhuman_only(
            ortholog_map['nonhuman_only']
        )
        context['gpcrome_picker_rows_json'] = json.dumps(
            DataMapperHome.gpcrome_receptor_picker_table_rows(
                maps=maps, extra_entry_names=nonhuman_entries, entry_label_overrides=nonhuman_stem_by_entry
            )
        )
        return context


class MapperListView(TemplateView):
    template_name = 'mapper/Mapper_List.html'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        maps = DataMapperHome.build_gpcrome_receptor_normalization_maps(include_odorant=True)
        context['receptor_select2_json'] = json.dumps(
            DataMapperHome.gpcrome_receptor_select2_options(maps=maps)
        )
        context['gpcrome_resolve_json'] = json.dumps(
            DataMapperHome.gpcrome_receptor_client_resolve_map(maps=maps)
        )
        ortholog_map = DataMapperHome.build_ortholog_species_map()
        context['ortholog_species_json'] = json.dumps(ortholog_map)
        nonhuman_entries, nonhuman_stem_by_entry = DataMapperHome.gpcrome_collapse_nonhuman_only(
            ortholog_map['nonhuman_only']
        )
        context['gpcrome_picker_rows_json'] = json.dumps(
            DataMapperHome.gpcrome_receptor_picker_table_rows(
                maps=maps, extra_entry_names=nonhuman_entries, entry_label_overrides=nonhuman_stem_by_entry
            )
        )
        context['receptor_info_json'] = json.dumps(
            DataMapperHome.gpcrome_receptor_info_for_list(
                maps=maps, extra_entry_names=nonhuman_entries, entry_label_overrides=nonhuman_stem_by_entry
            )
        )
        return context


class MapperClusterView(TemplateView):
    template_name = 'mapper/Mapper_Cluster.html'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        maps = DataMapperHome.build_gpcrome_receptor_normalization_maps(include_odorant=True)
        context['receptor_select2_json'] = json.dumps(
            DataMapperHome.gpcrome_receptor_select2_options(maps=maps)
        )
        context['gpcrome_resolve_json'] = json.dumps(
            DataMapperHome.gpcrome_receptor_client_resolve_map(maps=maps)
        )
        context['gpcrome_picker_rows_json'] = json.dumps(
            DataMapperHome.gpcrome_receptor_picker_table_rows(maps=maps)
        )
        full_matrix = DataMapperHome.generate_full_matrix()
        context['cluster_all_positions_json'] = full_matrix.to_json(orient='records')
        return context

