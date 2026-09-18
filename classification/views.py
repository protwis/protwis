from django.core.cache import cache, caches
from django.db.models import Case, F, IntegerField, Max, Prefetch, Q, When
from django.db.models.functions import Greatest, Least, Upper
from django.http import Http404, JsonResponse
from django.urls import reverse
from django.utils.text import slugify
from django.views import View
from django.views.generic import TemplateView

from classification import classification_db
from classification.models import ClusterCoord, ReceptorSimilarity, StructureSimilarity, TreeNetwork
from common.models import WebLink
from mapper.views import DataMapperHome
from protein.models import Gene, Protein, ProteinFamily, ProteinState
from collections import OrderedDict, defaultdict
import math
import json
import re
from string import Template
from urllib.parse import urlencode

import numpy as np
from sklearn.manifold import TSNE


try:
    cache_alignment = caches["alignments"]
except Exception:
    cache_alignment = cache


class ClassificationVisualizationMixin:
    GLOBAL_GROUP_KEY = "global"
    CLASS_GROUP_PREFIX = "class:"
    MODALITY_GROUP_LABELS = [
        "Orphan receptors",
        "Polypeptide receptors",
        "Small molecule receptors",
    ]
    ORPHAN_MODALITY_KEYS = {"orphan receptors"}
    POLYPEPTIDE_MODALITY_KEYS = {"peptide receptors", "protein receptors", "polypeptide receptors"}
    ORPHAN_MODALITY_GROUP_LABEL = "Orphan receptors"
    ORPHAN_SEARCH_SPLIT_LABELS = {"orphan receptors", "orphans receptors"}
    UNCLASSIFIED_BROWSER_LABEL = "Unclassified / Other GPCRs"
    UNCLASSIFIED_CLASS_NAME_CANDIDATES = ("Unclassified", "Classless", "Other GPCRs", "Unclassified / Other GPCRs")
    CLASS_KEY_ALIASES = {
        "UNCLASSIFIED": "U",
        "CLASSLESS": "U",
        "OTHERGPCRS": "U",
        "UNCLASSIFIEDOTHERGPCRS": "U",
    }
    # Alphabetical by symbol -- the general ordering rule for every class listing/dropdown/table
    # in this app -- with "U" (Unclassified) always forced last, never sorted in alphabetically.
    CLASS_VISUALIZATION_CONFIG = OrderedDict([
        ("A", {"label": "Class A", "title": "Class A (Rhodopsin)", "slug": "001"}),
        ("B1", {"label": "Class B1", "title": "Class B1 (Secretin)", "slug": "002"}),
        ("B2", {"label": "Class B2", "title": "Class B2 (Adhesion)", "slug": "003"}),
        ("C", {"label": "Class C", "title": "Class C (Glutamate)", "slug": "004"}),
        ("F", {"label": "Class F", "title": "Class F (Frizzled)", "slug": "006"}),
        ("O1", {"label": "Class O1", "title": "Class O1 (Fish-like olfactory receptors)", "slug": "007"}),
        ("O2", {"label": "Class O2", "title": "Class O2 (Tetrapod-specific olfactory receptors)", "slug": "008"}),
        ("T2", {"label": "Class T2", "title": "Class T2 (Taste 2)", "slug": "009"}),
        ("V", {"label": "Class V", "title": "Class V (Vomeronasal)", "slug": "010"}),
        ("U", {"label": "Unclassified", "title": "Unclassified", "slug": "011"}),
    ])
    TREE_DISABLED_CLASS_KEYS = {"U"}
    # Chemotypes that don't make sense as a standalone Chemotype-tab dataset.
    TREE_EXCLUDED_CHEMOTYPES = {"odorant receptors", "ion receptors"}
    # Editorial exceptions to the generic Class-tab tree shape -- verified against the DB (Class C
    # genuinely has chemotype annotations, so this can't be derived, only curated).
    CLASS_TREE_OVERRIDES = {
        "A": {"colorMode": "chemotype"},  # kept multi-color on purpose -- flat class color looked worse here
        "B2": {"skip_chemotype_layer": True},  # 6 sparse chemotypes over 42 receptors -- too crowded as its own ring
        "C": {"skip_chemotype_layer": True},
        "O2": {"family_sort": "natural"},  # "Family 1, 2, ..., 14", not "1, 10, 11, ..., 2, ..."
    }
    # Classes whose orphan-tagged receptors are split out into a legend rather than shown in the
    # main tree. Orphan-tagged receptors exist in B2/C/Unclassified too, but generalizing the split
    # is a visual redesign (it also requires a tree.js change) -- out of scope here.
    ORPHAN_SPLIT_CLASS_KEYS = {"A"}
    ORPHAN_CHEMOTYPE_KEYS = {"orphan receptors"}

    @classmethod
    def normalize_visualization_class_key(cls, raw_value):
        value = str(raw_value or "").strip().upper()
        if value in cls.CLASS_VISUALIZATION_CONFIG:
            return value
        compact_value = re.sub(r"[^A-Z0-9]", "", value)
        return cls.CLASS_KEY_ALIASES.get(compact_value)

    @classmethod
    def get_visualization_class_config(cls, class_key):
        key = cls.normalize_visualization_class_key(class_key)
        if not key:
            return None
        config = dict(cls.CLASS_VISUALIZATION_CONFIG[key])
        config["key"] = key
        config["group_key"] = cls.class_group_key(key)
        config["has_tree"] = key not in cls.TREE_DISABLED_CLASS_KEYS
        return config

    @classmethod
    def list_visualization_classes(cls):
        return [cls.get_visualization_class_config(key) for key in cls.CLASS_VISUALIZATION_CONFIG.keys()]

    @classmethod
    def class_group_key(cls, class_key):
        return f"{cls.CLASS_GROUP_PREFIX}{class_key}"

    @classmethod
    def class_key_from_group_key(cls, group_key):
        value = str(group_key or "").strip()
        if not value.startswith(cls.CLASS_GROUP_PREFIX):
            return None
        return cls.normalize_visualization_class_key(value[len(cls.CLASS_GROUP_PREFIX):])

    @classmethod
    def tree_class_keys(cls):
        """Ordered class keys eligible for classification trees -- single source of truth for
        all three tabs (Class / Modality / Chemotype)."""
        return [key for key in cls.CLASS_VISUALIZATION_CONFIG if key not in cls.TREE_DISABLED_CLASS_KEYS]

    @classmethod
    def class_slug_to_key_map(cls):
        """ProteinFamily class slug ('001') -> visualization key ('A')."""
        return {cfg["slug"]: key for key, cfg in cls.CLASS_VISUALIZATION_CONFIG.items()}

    def _group_modality_label(self, modality):
        key = str(modality or "").strip().lower()
        if key in self.ORPHAN_MODALITY_KEYS:
            return "Orphan receptors"
        if key in self.POLYPEPTIDE_MODALITY_KEYS:
            return "Polypeptide receptors"
        return "Small molecule receptors"

    @classmethod
    def _class_title_to_key_map(cls):
        mapping = {}
        legacy_title_aliases = {
            "Class O1 (fish-like)": "O1",
            "Class O1 (fish-like odorant)": "O1",
            "Class O1 (fish-like olfactory receptor)": "O1",
            "Class O1 (fish-like olfactory receptors)": "O1",
            "Class O2 (tetrapod specific)": "O2",
            "Class O2 (tetrapod specific odorant)": "O2",
            "Class O2 (tetrapod-specific olfactory receptor)": "O2",
            "Class O2 (tetrapod-specific olfactory receptors)": "O2",
            "Unclassified": "U",
            "Classless": "U",
            "Other GPCRs": "U",
            "Unclassified / Other GPCRs": "U",
        }
        for key, cfg in cls.CLASS_VISUALIZATION_CONFIG.items():
            title = cfg["title"]
            mapping[title] = key
            mapping[title.lower()] = key
        for alias, key in legacy_title_aliases.items():
            mapping[alias] = key
            mapping[alias.lower()] = key
        return mapping

    @staticmethod
    def _natural_sort_key(value):
        chunks = re.split(r"(\d+)", str(value or "").strip().lower())
        key = []
        for chunk in chunks:
            if chunk == "":
                continue
            if chunk.isdigit():
                key.append((0, int(chunk)))
            else:
                key.append((1, chunk))
        return tuple(key)

    @classmethod
    def _browser_class_key(cls, class_label, class_key=None):
        if class_key:
            return class_key
        label = str(class_label or "").strip()
        if not label:
            return None
        return "U"

    @classmethod
    def _browser_class_label(cls, browser_class_key):
        if browser_class_key == "U":
            return cls.UNCLASSIFIED_BROWSER_LABEL
        config = cls.get_visualization_class_config(browser_class_key)
        if config:
            return config["title"]
        return str(browser_class_key or "").strip()

    def _build_receptor_family_catalog(self):
        slug_to_key = self.class_slug_to_key_map()
        name_to_key = {}

        normalized_rows = []
        orphan_family_classes = defaultdict(OrderedDict)
        for row in classification_db.get_classification_rows():
            class_label = row["class_family_name"]
            chemotype = row["chemotype"]
            family_name = row["receptor_family"]
            modality = row["modality"]
            if not class_label or not chemotype or not family_name:
                continue
            class_key = slug_to_key.get(row["class_slug"])
            name_to_key[class_label] = class_key
            modality_group = self._group_modality_label(modality)
            normalized_rows.append({
                "class_label": class_label,
                "chemotype": chemotype,
                "family_name": family_name,
                "modality_group": modality_group,
                "receptor": row["uniprot"],
            })
            if modality_group == self.ORPHAN_MODALITY_GROUP_LABEL:
                orphan_family_classes[family_name][class_label] = True

        family_rows = OrderedDict()
        hierarchy = OrderedDict()
        for row in normalized_rows:
            class_label = row["class_label"]
            chemotype = row["chemotype"]
            family_name = row["family_name"]
            modality_group = row["modality_group"]
            class_key = name_to_key.get(class_label)
            browser_class_key = self._browser_class_key(class_label, class_key)
            class_config = self.get_visualization_class_config(class_key) if class_key else None
            split_orphan_by_class = (
                modality_group == self.ORPHAN_MODALITY_GROUP_LABEL
                and len(orphan_family_classes.get(family_name, {})) > 1
            )
            family_storage_key = (
                "{}::{}".format(family_name, class_label)
                if split_orphan_by_class else family_name
            )
            family_label = family_name
            if split_orphan_by_class and str(family_name).strip().lower() in self.ORPHAN_SEARCH_SPLIT_LABELS:
                family_label = "{} ({})".format(
                    family_name,
                    class_config["label"] if class_config else class_label,
                )

            entry = family_rows.setdefault(family_storage_key, {
                "storage_key": family_storage_key,
                "name": family_name,
                "label": family_label,
                "browser_label": family_name,
                "modality_groups": OrderedDict(),
                "classes": OrderedDict(),
                "chemotypes": OrderedDict(),
                "receptors": OrderedDict(),
                "split_by_class": split_orphan_by_class,
                "split_class_label": class_label if split_orphan_by_class else None,
                "split_class_key": class_key if split_orphan_by_class else None,
                "class_slug_prefix": class_config["slug"] if split_orphan_by_class and class_config else None,
            })
            entry["classes"][class_label] = True
            entry["chemotypes"][chemotype] = True
            entry["modality_groups"][modality_group] = True
            entry["receptors"][row["receptor"]] = True

            if not browser_class_key:
                continue
            class_bucket = hierarchy.setdefault(browser_class_key, {
                "label": self._browser_class_label(browser_class_key),
                "chemotypes": OrderedDict(),
            })
            chemotype_bucket = class_bucket["chemotypes"].setdefault(chemotype, {
                "modality_group": modality_group,
                "families": OrderedDict(),
            })
            chemotype_bucket["modality_group"] = modality_group
            chemotype_bucket["families"][family_name] = family_storage_key

        key_counts = defaultdict(int)
        family_entries = []
        family_lookup = OrderedDict()
        family_by_storage_key = {}
        for family_storage_key in sorted(
            family_rows.keys(),
            key=lambda value: (
                self._natural_sort_key(family_rows[value]["name"]),
                self._natural_sort_key(family_rows[value].get("split_class_label") or ""),
            ),
        ):
            raw = family_rows[family_storage_key]
            base_key_parts = [raw["name"]]
            if raw.get("split_class_key"):
                base_key_parts.append(raw["split_class_key"])
            elif raw.get("split_class_label"):
                base_key_parts.append(raw["split_class_label"])
            base_key = slugify("-".join(base_key_parts)) or "family"
            key_counts[base_key] += 1
            family_key = base_key if key_counts[base_key] == 1 else f"{base_key}-{key_counts[base_key]}"
            class_labels = list(raw["classes"].keys())
            chemotype_labels = list(raw["chemotypes"].keys())
            receptor_labels = list(raw["receptors"].keys())
            modality_groups = list(raw["modality_groups"].keys())
            class_keys = [name_to_key.get(label) or "U" for label in class_labels]
            entry = {
                "key": family_key,
                "name": raw["name"],
                "label": raw["label"],
                "browser_label": raw["browser_label"],
                "url": reverse("classification-visualizations-family", kwargs={"family_key": family_key}),
                "class_labels": class_labels,
                "class_keys": class_keys,
                "split_by_class": raw["split_by_class"],
                "split_class_label": raw["split_class_label"],
                "split_class_key": raw["split_class_key"],
                "class_slug_prefix": raw["class_slug_prefix"],
                "chemotypes": chemotype_labels,
                "modality_groups": modality_groups,
                "receptor_labels": receptor_labels,
                "receptor_count": len(receptor_labels),
            }
            family_entries.append(entry)
            family_lookup[family_key] = entry
            family_by_storage_key[family_storage_key] = entry

        browser_nodes = []
        browser_order = {key: idx for idx, key in enumerate(self.CLASS_VISUALIZATION_CONFIG.keys())}
        for browser_class_key in sorted(hierarchy.keys(), key=lambda value: browser_order.get(value, 999)):
            class_bucket = hierarchy[browser_class_key]
            grouped_chemotypes = []
            aggregated_families = OrderedDict()

            for chemotype in sorted(class_bucket["chemotypes"].keys(), key=self._natural_sort_key):
                families = []
                chemotype_node = class_bucket["chemotypes"][chemotype]
                modality_group = chemotype_node.get("modality_group") or "Small molecule receptors"
                modality_theme = slugify(modality_group)
                for family_name in sorted(chemotype_node["families"].keys(), key=self._natural_sort_key):
                    family_storage_key = chemotype_node["families"].get(family_name)
                    family_entry = family_by_storage_key.get(family_storage_key)
                    if not family_entry:
                        continue
                    family_node = {
                        "key": family_entry["key"],
                        "label": family_entry.get("browser_label") or family_name,
                        "url": family_entry["url"],
                        "receptor_count": family_entry["receptor_count"],
                    }
                    families.append(family_node)
                    aggregated_families[family_node["key"]] = family_node
                if families:
                    grouped_chemotypes.append({
                        "label": chemotype,
                        "modality_group": modality_group,
                        "modality_theme": modality_theme,
                        "families": families,
                    })

            if grouped_chemotypes:
                browser_nodes.append({
                    "label": class_bucket["label"],
                    "class_key": browser_class_key,
                    "interactive_chemotypes": browser_class_key == "A",
                    "chemotypes": grouped_chemotypes,
                    "families": sorted(
                        aggregated_families.values(),
                        key=lambda family: self._natural_sort_key(family.get("label") or "")
                    ),
                })

        return {
            "entries": family_entries,
            "lookup": family_lookup,
            "browser_nodes": browser_nodes,
        }

    def get_receptor_family_entry(self, family_key):
        catalog = self._build_receptor_family_catalog()
        return catalog["lookup"].get(str(family_key or "").strip()), catalog

    def _visualization_class_name_candidates(self, class_label):
        label = str(class_label or "").strip()
        if not label:
            return []
        if label.lower() in {"unclassified", "classless", "other gpcrs", "unclassified / other gpcrs"}:
            return list(self.UNCLASSIFIED_CLASS_NAME_CANDIDATES)
        return [label]

    def _visualization_class_family_q(self, class_key, family_path_prefix):
        config = self.get_visualization_class_config(class_key)
        if not config:
            return Q()

        prefix = str(family_path_prefix or "")
        q_obj = Q(**{prefix + "slug": config["slug"]})
        for candidate in self._visualization_class_name_candidates(config["title"]):
            q_obj |= Q(**{prefix + "name__iexact": candidate})
        return q_obj

    def _visualization_class_protein_q(self, class_key):
        return self._visualization_class_family_q(class_key, "family__parent__parent__parent__")

    def _visualization_class_clustercoord_q(self, class_key):
        return self._visualization_class_family_q(class_key, "protein__family__parent__parent__parent__")

    def _visualization_receptor_label_q(self, receptor_labels):
        q_obj = Q()
        for label in receptor_labels:
            clean_label = str(label or "").strip()
            if not clean_label:
                continue
            q_obj |= Q(accession__iexact=clean_label)
            q_obj |= Q(entry_name__istartswith="{}_".format(clean_label))
        return q_obj

    def _visualization_family_queryset_class_q(self, class_slug_prefix=None, class_label=None):
        q_obj = Q()
        if class_slug_prefix:
            q_obj |= Q(family_slug__startswith=class_slug_prefix)
        if class_label:
            for candidate in self._visualization_class_name_candidates(class_label):
                q_obj |= Q(family__parent__parent__parent__name__iexact=candidate)
        return q_obj

    def _get_visualization_family_queryset(self, family_name, class_slug_prefix=None, class_label=None, exact_family_name_only=False):
        qs = (
            Protein.objects
            .annotate(
                family_slug=F("family__slug"),
                entry_name_upper=Upper("entry_name"),
                accession_upper=Upper("accession"),
            )
            .filter(
                parent_id__isnull=True,
                species__common_name__iexact="Human",
            )
            .exclude(accession=None)
            .select_related("family__parent__parent__parent")
            .prefetch_related(
                Prefetch(
                    "genes",
                    queryset=Gene.objects.only("name", "position").order_by("position"),
                    to_attr="primary_genes_self",
                )
            )
            .order_by("entry_name")
            .distinct()
        )
        if exact_family_name_only:
            qs = qs.filter(
                Q(family__parent__name__iexact=family_name)
                | Q(family__name__iexact=family_name)
            )
        else:
            qs = qs.filter(
                Q(family__parent__name__iexact=family_name)
                | Q(family__name__iexact=family_name)
            )
        class_scope_q = self._visualization_family_queryset_class_q(
            class_slug_prefix=class_slug_prefix,
            class_label=class_label,
        )
        if class_scope_q:
            qs = qs.filter(class_scope_q)
        return qs

    def get_visualization_family_proteins(self, family_entry):
        family_name = family_entry["name"]
        receptor_labels = [
            str(label).strip().upper()
            for label in family_entry.get("receptor_labels", [])
            if str(label).strip()
        ]
        split_by_class = bool(family_entry.get("split_by_class"))
        class_slug_prefix = family_entry.get("class_slug_prefix")
        if split_by_class:
            if receptor_labels:
                receptor_q = self._visualization_receptor_label_q(receptor_labels)
                qs = (
                    Protein.objects
                    .annotate(
                        family_slug=F("family__slug"),
                    )
                    .filter(
                        parent_id__isnull=True,
                        species__common_name__iexact="Human",
                    )
                    .exclude(accession=None)
                    .filter(receptor_q)
                    .select_related("family__parent__parent__parent")
                    .prefetch_related(
                        Prefetch(
                            "genes",
                            queryset=Gene.objects.only("name", "position").order_by("position"),
                            to_attr="primary_genes_self",
                        )
                    )
                    .order_by("entry_name")
                    .distinct()
                )
                class_scope_q = self._visualization_family_queryset_class_q(
                    class_slug_prefix=class_slug_prefix,
                    class_label=family_entry.get("split_class_label"),
                )
                if class_scope_q:
                    qs = qs.filter(class_scope_q)
            else:
                qs = self._get_visualization_family_queryset(
                    family_name=family_name,
                    class_slug_prefix=class_slug_prefix,
                    class_label=family_entry.get("split_class_label"),
                    exact_family_name_only=True,
                )
        else:
            qs = self._get_visualization_family_queryset(
                family_name=family_name,
                class_slug_prefix=class_slug_prefix,
                class_label=family_entry.get("split_class_label"),
                exact_family_name_only=False,
            )
            if receptor_labels:
                qs = qs.filter(self._visualization_receptor_label_q(receptor_labels))
        return list(qs)

    def get_requested_visualization_scope(self, request=None, raise_404=False):
        request = request or self.request
        raw_class_key = request.GET.get("class")
        class_key = self.normalize_visualization_class_key(raw_class_key)
        if not raw_class_key:
            return {
                "class_key": None,
                "group_key": self.GLOBAL_GROUP_KEY,
                "config": None,
            }
        if not class_key:
            if raise_404:
                raise Http404("Unknown class visualization")
            raise ValueError("Unknown class visualization")
        config = self.get_visualization_class_config(class_key)
        if not config:
            if raise_404:
                raise Http404("Unknown class visualization")
            raise ValueError("Unknown class visualization")
        return {
            "class_key": class_key,
            "group_key": config["group_key"],
            "config": config,
        }

    @staticmethod
    def normalize_tree_selection_type(raw_value):
        value = str(raw_value or "").strip()
        if value.lower() == "modality":
            return "Modality"
        if value.lower() == "chemotype":
            return "Chemotype"
        return ""

    @staticmethod
    def _format_similarity_display(value):
        if value is None:
            return ""
        try:
            rounded = round(float(value), 1)
        except Exception:
            return ""
        if rounded == int(rounded):
            return str(int(rounded))
        return "{:.1f}".format(rounded)

    def resolve_tree_visualization_selection(self, tree_type, selection):
        normalized_type = self.normalize_tree_selection_type(tree_type)
        raw_selection = str(selection or "").strip()
        if not normalized_type or not raw_selection:
            raise ValueError("Unknown classification tree selection")

        requested_key = raw_selection.lower()
        slug_to_key = self.class_slug_to_key_map()
        allowed = set(self.tree_class_keys())

        matched = []
        normalized_selection = raw_selection
        for row in classification_db.get_classification_rows():
            class_key = slug_to_key.get(row["class_slug"])
            if class_key not in allowed:
                continue
            if normalized_type == "Modality":
                raw_modality = (row["modality"] or "").strip()
                if not raw_modality:
                    continue
                grouped_modality = self._group_modality_label(raw_modality)
                if (
                    grouped_modality.lower() != requested_key
                    and raw_modality.lower() != requested_key
                ):
                    continue
                normalized_selection = grouped_modality
            else:
                chemotype = (row["chemotype"] or "").strip()
                if not chemotype or chemotype.lower() in self.TREE_EXCLUDED_CHEMOTYPES:
                    continue
                if chemotype.lower() != requested_key:
                    continue
                normalized_selection = chemotype
            matched.append((class_key, row))

        if not matched:
            raise ValueError("Unknown classification tree selection")

        class_keys = OrderedDict()
        receptor_labels = OrderedDict()
        entry_names = OrderedDict()
        for class_key, row in matched:
            class_keys[class_key] = True
            label = row["uniprot"]
            receptor_labels[label] = True
            entry_names[row["entry_name"]] = True

        return {
            "type": normalized_type,
            "selection": normalized_selection,
            "class_keys": list(class_keys.keys()),
            "receptor_labels": list(receptor_labels.keys()),
            "entry_names": list(entry_names.keys()),
            "receptor_count": len(receptor_labels),
        }

    def build_tree_selection_matrix_payload(self, selection_info):
        entry_names = list((selection_info or {}).get("entry_names") or [])
        if not entry_names:
            return {"entities": [], "matrix": [], "meta": {"n_points": 0}}

        proteins = list(
            Protein.objects
            .filter(entry_name__in=entry_names)
            .only("id", "entry_name", "name", "accession")
            .prefetch_related(
                Prefetch(
                    "genes",
                    queryset=Gene.objects.only("name", "position").order_by("position"),
                    to_attr="primary_genes_self",
                )
            )
            .order_by("entry_name")
        )
        if not proteins:
            return {"entities": [], "matrix": [], "meta": {"n_points": 0}}

        protein_ids = [protein.id for protein in proteins]
        pair_identity = {}
        pair_similarity = {}
        pair_qs = (
            ReceptorSimilarity.objects
            .filter(protein_ref_id__in=protein_ids, protein_target_id__in=protein_ids)
            .values("protein_ref_id", "protein_target_id", "identity", "similarity")
        )
        for rec in pair_qs.iterator():
            try:
                a = int(rec["protein_ref_id"])
                b = int(rec["protein_target_id"])
            except Exception:
                continue
            if a == b:
                continue
            key = (a, b) if a < b else (b, a)
            try:
                identity = float(rec["identity"])
            except Exception:
                identity = None
            try:
                similarity = float(rec["similarity"])
            except Exception:
                similarity = None
            if identity is not None:
                previous = pair_identity.get(key)
                if previous is None or identity > previous:
                    pair_identity[key] = identity
            if similarity is not None:
                previous = pair_similarity.get(key)
                if previous is None or similarity > previous:
                    pair_similarity[key] = similarity

        entities = []
        for protein in proteins:
            primary_genes = getattr(protein, "primary_genes_self", None) or []
            gene_label = primary_genes[0].name if primary_genes else (
                str(protein.entry_name or "").split("_", 1)[0].upper()
            )
            entities.append({
                "symbol": protein.entry_name,
                "name": protein.short(),
                "short_label": protein.entry_short(),
                "gene_label": gene_label,
            })

        matrix = []
        for i, protein in enumerate(proteins):
            row = []
            for j, other_protein in enumerate(proteins):
                if i == j:
                    identity = 100.0
                    similarity = 100.0
                else:
                    key = (min(protein.id, other_protein.id), max(protein.id, other_protein.id))
                    identity = pair_identity.get(key)
                    similarity = pair_similarity.get(key)
                row.append({
                    "source": protein.entry_name,
                    "target": other_protein.entry_name,
                    "identity": float(identity) if identity is not None else None,
                    "identity_display": self._format_similarity_display(identity),
                    "similarity": float(similarity) if similarity is not None else None,
                    "similarity_display": self._format_similarity_display(similarity),
                })
            matrix.append(row)

        return {
            "entities": entities,
            "matrix": matrix,
            "meta": {
                "n_points": len(entities),
                "selection": (selection_info or {}).get("selection", ""),
                "type": (selection_info or {}).get("type", ""),
            },
        }


class ClassificationVisualizationsLanding(ClassificationVisualizationMixin, TemplateView):
    template_name = "classification/ClassificationVisualizations.html"

    CHEMOTYPE_PLACEHOLDERS = [
        "Adhesion receptors",
        "Alicarboxylic acid receptors",
        "Aminergic receptors",
        "Amino acid receptors",
        "Lipid receptors",
        "Melatonin receptors",
        "Nucleotide receptors",
        "Orphan receptors",
        "Peptide receptors",
        "Protein receptors",
        "Retinal receptors",
        "Steroid receptors",
        "Tastant receptors",
    ]
    def _tree_only_url(self, tree_type, selection):
        return "{}?{}".format(
            reverse("classification-visualizations-tree"),
            urlencode({
                "type": tree_type,
                "selection": selection,
            }),
        )

    def _build_modality_chemotype_branches(self):
        try:
            rows = classification_db.get_classification_rows()
        except Exception:
            return []

        chemotype_to_modality = OrderedDict()
        for row in rows:
            chemotype = row["chemotype"]
            modality = row["modality"] or "Other / unknown"
            if not chemotype:
                continue
            if chemotype.strip().lower() in self.TREE_EXCLUDED_CHEMOTYPES:
                continue
            if chemotype not in chemotype_to_modality:
                chemotype_to_modality[chemotype] = modality

        grouped = OrderedDict((label, []) for label in self.MODALITY_GROUP_LABELS)
        for chemotype in sorted(chemotype_to_modality.keys(), key=lambda value: str(value).lower()):
            modality_group = self._group_modality_label(chemotype_to_modality[chemotype])
            grouped[modality_group].append({
                "label": chemotype,
                "url": self._tree_only_url("Chemotype", chemotype),
                "modality_theme": slugify(modality_group),
            })

        branches = []
        for modality, children in grouped.items():
            if not children:
                continue
            branches.append({
                "label": modality,
                "url": self._tree_only_url("Modality", modality),
                "theme_key": slugify(modality),
                "chemotypes": children,
            })
        return branches

    def get_context_data(self, **kwargs):
        ctx = super().get_context_data(**kwargs)
        class_buttons = []
        for config in self.list_visualization_classes():
            if config["key"] == "U":
                continue
            class_buttons.append({
                "key": config["key"],
                "label": config["key"],
                "title": config["title"],
                "url": reverse(
                    "classification-visualizations-class",
                    kwargs={"class_key": config["key"]},
                ),
            })
        ctx["class_buttons"] = class_buttons
        modality_branches = self._build_modality_chemotype_branches()
        ctx["modality_branches"] = modality_branches
        family_catalog = self._build_receptor_family_catalog()
        ctx["receptor_family_entries"] = family_catalog["entries"]
        ctx["receptor_family_browser"] = family_catalog["browser_nodes"]
        ctx["receptor_family_entries_json"] = json.dumps({
            item["label"]: {
                "key": item["key"],
                "label": item["label"],
                "url": item["url"],
                "receptor_count": item["receptor_count"],
            }
            for item in family_catalog["entries"]
        })
        superfamily_url = reverse("classification-visualizations-superfamily")
        ctx["superfamily_url"] = superfamily_url
        ctx["landing_payload_json"] = json.dumps({
            "superfamily": {
                "label": "GPCR superfamily",
                "note": "(all classes)",
                "url": superfamily_url,
            },
            "classButtons": class_buttons,
            "receptorFamilies": {
                "entries": family_catalog["entries"],
                "entriesByLabel": {
                    item["label"]: {
                        "key": item["key"],
                        "label": item["label"],
                        "url": item["url"],
                        "receptor_count": item["receptor_count"],
                    }
                    for item in family_catalog["entries"]
                },
                "browserNodes": family_catalog["browser_nodes"],
            },
            "modalityBranches": modality_branches,
        })
        return ctx


class ClassificationVisualizationDetail(ClassificationVisualizationMixin, TemplateView):
    template_name = "classification/ClassificationVisualizationDetail.html"

    def dispatch(self, request, *args, **kwargs):
        class_key = kwargs.get("class_key")
        self.class_config = self.get_visualization_class_config(class_key)
        if not self.class_config:
            raise Http404("Unknown class visualization")
        return super(ClassificationVisualizationDetail, self).dispatch(request, *args, **kwargs)

    def get_context_data(self, **kwargs):
        ctx = super().get_context_data(**kwargs)
        class_key = self.class_config["key"]
        cluster_query = urlencode({
            "class": class_key,
            "embed": "1",
        })
        ctx["class_key"] = class_key
        ctx["class_title"] = self.class_config["title"]
        ctx["class_label"] = self.class_config["label"]
        ctx["has_tree"] = self.class_config["has_tree"]
        ctx["cluster_url"] = "{}?{}".format(
            reverse("classification-structuresim"),
            cluster_query,
        )
        if self.class_config["has_tree"]:
            tree_view = Classification_tree()
            tree_ctx, tree_sets = tree_view._build_tree_datasets()
            if "error" not in tree_ctx:
                tree_ctx.update(Classification_tree.build_selection_context(
                    tree_sets,
                    requested_type="Class",
                    requested_selection=class_key,
                    locked=True,
                ))
            ctx.update(tree_ctx)
        else:
            ctx["tree_note"] = (
                "Due to the unclassified nature of these receptors, a classification tree is unavailable."
            )
        return ctx


class ClassificationTreeVisualizationDetail(ClassificationVisualizationMixin, TemplateView):
    template_name = "classification/ClassificationTreeVisualizationDetail.html"

    def get_context_data(self, **kwargs):
        ctx = super().get_context_data(**kwargs)
        raw_type = self.request.GET.get("type")
        raw_selection = self.request.GET.get("selection")
        try:
            selection_info = self.resolve_tree_visualization_selection(raw_type, raw_selection)
        except ValueError as e:
            raise Http404(str(e))

        class_keys = selection_info["class_keys"]
        has_cluster = len(class_keys) == 1 and selection_info["receptor_count"] >= 2
        cluster_url = ""
        if has_cluster:
            cluster_query = urlencode({
                "class": class_keys[0],
                "filter_type": selection_info["type"],
                "filter_selection": selection_info["selection"],
                "embed": "1",
            })
            cluster_url = "{}?{}".format(reverse("classification-structuresim"), cluster_query)

        cluster_note = ""
        if not has_cluster:
            if len(class_keys) > 1:
                cluster_note = (
                    "Cluster is unavailable because this selection spans multiple GPCR classes."
                )
            elif selection_info["receptor_count"] < 2:
                cluster_note = (
                    "Cluster is unavailable because this selection has fewer than two receptors."
                )
            else:
                cluster_note = (
                    "Cluster is unavailable because this selection could not be mapped to a single GPCR class."
                )

        tree_view = Classification_tree()
        tree_ctx, tree_sets = tree_view._build_tree_datasets()
        if "error" not in tree_ctx:
            tree_ctx.update(Classification_tree.build_selection_context(
                tree_sets,
                requested_type=selection_info["type"],
                requested_selection=selection_info["selection"],
                locked=True,
            ))
        ctx.update(tree_ctx)

        ctx["page_title"] = selection_info["selection"]
        ctx["tree_type"] = selection_info["type"]
        ctx["tree_selection"] = selection_info["selection"]
        ctx["has_cluster"] = has_cluster
        ctx["cluster_url"] = cluster_url
        ctx["cluster_note"] = cluster_note
        ctx["matrix_payload_json"] = json.dumps(self.build_tree_selection_matrix_payload(selection_info))
        return ctx



class GPCRSuperfamilyVisualizationDetail(TemplateView):
    template_name = "classification/ClassificationSuperfamilyDetail.html"

    def get_context_data(self, **kwargs):
        ctx = super().get_context_data(**kwargs)
        ctx["page_title"] = "GPCR superfamily"
        ctx.update(ClassificationWheel.build_wheel_context())
        ctx.update(NewClassClusterTree.build_cluster_tree_context(self.request))
        ctx.update(CrossClassSimilarity().get_context_data())
        return ctx

class Classification(ClassificationVisualizationMixin, TemplateView):
    template_name = "classification/Classification.html"

    # Static class information for the Classes table
    CLASSES_TABLE_DATA = [
        {"code": "A", "name": "Rhodopsin", "species": "Yes", "non_sensory_share": "Majority", "sensory_function": "Vision & light-sensing"},
        {"code": "B1", "name": "Secretin", "species": "Yes", "non_sensory_share": "All", "sensory_function": "-"},
        {"code": "B2", "name": "Adhesion", "species": "Yes", "non_sensory_share": "Majority", "sensory_function": "-"},
        {"code": "C", "name": "Glutamate", "species": "Yes", "non_sensory_share": "Majority", "sensory_function": "Taste (sweet/umami)"},
        {"code": "D1", "name": "Fungal pheromone 1", "species": "Fungi", "non_sensory_share": "-", "sensory_function": "Pheromone-sensing"},
        {"code": "D2", "name": "Fungal pheromone 2", "species": "Fungi", "non_sensory_share": "-", "sensory_function": "Pheromone-sensing"},
        {"code": "E", "name": "Slime mold cAMP", "species": "Slime molds, amoebas", "non_sensory_share": "-", "sensory_function": "Pheromone-sensing (cAMP, in chemotaxis)"},
        {"code": "F", "name": "Frizzled", "species": "Yes", "non_sensory_share": "All", "sensory_function": "-"},
        {"code": "OP1", "name": "Olfactory-polyfunctional 1", "species": "Yes", "non_sensory_share": "Minority", "sensory_function": "Olfaction"},
        {"code": "OP2", "name": "Olfactory-polyfunctional 2", "species": "Yes", "non_sensory_share": "Minority", "sensory_function": "Olfaction"},
        {"code": "T2", "name": "Taste 2", "species": "Yes", "non_sensory_share": "-", "sensory_function": "Taste (bitter)"},
        {"code": "V1", "name": "Vomeronasal 1", "species": "Amphibia, reptiles & non-primate mammals", "non_sensory_share": "-", "sensory_function": "Pheromone-sensing"},
        {"code": "V2", "name": "Vomeronasal 2", "species": "Amphibia, reptiles & non-primate mammals", "non_sensory_share": "-", "sensory_function": "Pheromone-sensing"},
        {"code": "U", "name": "Unclassified", "species": "Yes", "non_sensory_share": "Unknown", "sensory_function": "Unknown"},
    ]


    def get_context_data(self, **kwargs):
        ctx = super().get_context_data(**kwargs)

        ctx["gpcr_rows"] = json.dumps(_build_gpcr_browser_rows())

        slug_to_key = self.class_slug_to_key_map()
        rows = [
            dict(row, class_symbol=slug_to_key.get(row["class_slug"]))
            for row in classification_db.get_classification_rows()
        ]
        rows = [row for row in rows if row["class_symbol"]]

        # nice ordered list of classes
        class_order = list(self.CLASS_VISUALIZATION_CONFIG.keys())
        def sort_classes(s):
            return sorted(s, key=lambda x: (class_order.index(x)
                                            if x in class_order else 999, x))

        # ---------- 1) Ligand type table ----------
        # Using Chemotype as ligand_type and Modality as ligand_group
        lt_agg = {}  # chemotype -> {"group": modality, "classes": set([...])}
        for row in rows:
            chemotype = row["chemotype"]
            if not chemotype:
                continue
            modality = row["modality"] or "Other / unknown"
            entry = lt_agg.setdefault(chemotype, {"group": modality, "classes": set()})
            entry["classes"].add(row["class_symbol"])

        ligand_type_rows = []
        for chemotype in sorted(lt_agg.keys(), key=str.lower):
            entry = lt_agg[chemotype]
            cls_list = sort_classes(entry["classes"])
            ligand_type_rows.append({
                "ligand_type": chemotype,
                "ligand_group": entry["group"],
                "classes": ", ".join(cls_list),
            })

        # ---------- 2) Receptor families: non-sensory vs sensory vs orphan ----------
        # Using Sense directly:
        # - "non-sensory" -> non-sensory
        # - "unknown" -> orphan
        # - any other non-empty Sense -> sensory

        non_sens_triples = set()       # (class_symbol, family, chemotype)
        sensory_map = {}              # (family, chemotype) -> set(classes)
        orphan_map = {}               # (family, chemotype) -> set(classes)

        for row in rows:
            fam = row["receptor_family"]
            chemotype = row["chemotype"]
            if not fam or not chemotype:
                continue

            symbol = row["class_symbol"]
            sense = (row["sense"] or "").strip().lower()

            if sense == "non-sensory":
                non_sens_triples.add((symbol, fam, chemotype))
            elif sense == "unknown":
                orphan_map.setdefault((fam, chemotype), set()).add(symbol)
            elif sense:  # Any other non-empty sense value is considered sensory
                sensory_map.setdefault((fam, chemotype), set()).add(symbol)
            # If sense is empty, skip it

        # non-sensory: Class / Receptor family / Chemotype
        rf_non_rows = [
            {
                "class_symbol": cs,
                "receptor_family": fam,
                "ligand_type": chemotype,  # Using chemotype as ligand_type for template
            }
            for (cs, fam, chemotype) in sorted(
                non_sens_triples,
                key=lambda t: (class_order.index(t[0])
                               if t[0] in class_order else 999,
                               t[0].lower(), t[1].lower())
            )
        ]

        # sensory: Receptor family / Chemotype / Found in classes
        rf_sens_rows = []
        for (fam, chemotype), classes in sensory_map.items():
            cls_list = sort_classes(classes)
            rf_sens_rows.append({
                "receptor_family": fam,
                "ligand_type": chemotype,  # Using chemotype as ligand_type for template
                "classes": ", ".join(cls_list),
            })

        rf_sens_rows.sort(key=lambda r: (r["receptor_family"].lower(),
                                         r["ligand_type"].lower()))

        # orphan: Receptor family / Chemotype / Found in classes
        rf_orphan_rows = []
        for (fam, chemotype), classes in orphan_map.items():
            cls_list = sort_classes(classes)
            rf_orphan_rows.append({
                "receptor_family": fam,
                "ligand_type": chemotype,  # Using chemotype as ligand_type for template
                "classes": ", ".join(cls_list),
            })

        rf_orphan_rows.sort(key=lambda r: (r["receptor_family"].lower(),
                                           r["ligand_type"].lower()))

        # ---------- 3) Classes table ----------
        # Use static data for the classes table, pre-sorted by Code
        classes_rows = self.CLASSES_TABLE_DATA.copy()

        # Custom sort function for Code (A, B1, B2, C, D1, D2, E, F, OP1, OP2, T2, V1, V2, Cl)
        def sort_code_key(item):
            code = item["code"]
            # Extract letter prefix and number suffix
            match = re.match(r'([A-Za-z]+)(\d*)', code)
            if match:
                letter_part = match.group(1)
                num_part = match.group(2)
                # Handle special cases
                if letter_part.upper() == "CL":
                    return (999, 0)  # Unclassified goes last
                # Convert letter part to sortable value
                letter_order = {"A": 1, "B": 2, "C": 3, "D": 4, "E": 5, "F": 6,
                               "OP": 7, "O": 7, "T": 8, "V": 9}
                letter_key = letter_order.get(letter_part.upper(), 999)
                # Convert number part to integer (empty string = 0)
                num_key = int(num_part) if num_part else 0
                return (letter_key, num_key)
            # Fallback for unexpected formats
            return (999, 0)

        # Sort by Code
        classes_rows.sort(key=sort_code_key)

        # Serialize to JSON strings for template (using |safe filter)
        ctx["classes_rows"] = json.dumps(classes_rows)
        ctx["ligand_type_rows"] = json.dumps(ligand_type_rows)
        ctx["rf_non_rows"] = json.dumps(rf_non_rows)
        ctx["rf_sens_rows"] = json.dumps(rf_sens_rows)
        ctx["rf_orphan_rows"] = json.dumps(rf_orphan_rows)

        return ctx


# Row-getters shared by Classification_tree's Class/Modality/Chemotype tree builders --
# classification_db.build_grouped_children groups classification_db.get_classification_rows()
# dicts by an ordered list of these.
CHEMOTYPE_GETTER = lambda row: row["chemotype"]
FAMILY_GETTER = lambda row: row["receptor_family"]
CLASS_KEY_GETTER = lambda row: row["class_key"]


class Classification_tree(ClassificationVisualizationMixin):
    """
    Plain helper class (not a Django view) — provides DB-derived tree datasets (via
    classification/classification_db.py) to whichever Detail page inlines the
    Classification-tree content.
    """

    def _build_tree_datasets(self, **kwargs):
        """
        Build the DB-derived tree datasets (via classification/classification_db.py) independent
        of the current request, so a parent page (e.g. ClassificationVisualizationDetail) can
        inline this content with its own fixed type/selection.

        Returns (ctx, tree_sets) where tree_sets is the raw dict (not JSON-serialized)
        needed by build_selection_context to validate/resolve a selection.
        """
        ctx = dict(kwargs)

        try:
            rows = classification_db.get_classification_rows()
        except Exception as e:
            ctx["error"] = f"Error loading classification annotations from the database: {e}"
            return ctx, {}

        slug_to_key = self.class_slug_to_key_map()
        allowed_keys = self.tree_class_keys()
        allowed = set(allowed_keys)

        rows = [dict(row, class_key=slug_to_key.get(row["class_slug"])) for row in rows]
        rows = [row for row in rows if row["class_key"] in allowed]
        if not rows:
            ctx["error"] = "No classification annotations found in the database."
            return ctx, {}

        base_tree_options = {
            "branch_trunc": 0,
            "leaf_offset": 30,
            "anchor": "",
            "label_free": [],
            "centerBadgeR": 0,
            "centerBadgePadding": 0,
            "firstRingExtra": 45,
        }
        tree_sets = {
            "Class": {"options": [], "plots": {}},
            "Modality": {"options": [], "plots": {}},
            "Chemotype": {"options": [], "plots": {}},
        }

        self._add_class_plots(tree_sets, rows, allowed_keys, base_tree_options)
        self._add_modality_plots(tree_sets, rows, base_tree_options)
        self._add_chemotype_plots(tree_sets, rows, base_tree_options)

        ctx["tree_sets"] = json.dumps(tree_sets)
        ctx["tree_leaf_label_lookup"] = json.dumps(classification_db.build_leaf_label_lookup(rows))
        return ctx, tree_sets

    def _add_class_plots(self, tree_sets, rows, allowed_keys, base_tree_options):
        by_class = {}
        for row in rows:
            by_class.setdefault(row["class_key"], []).append(row)

        for class_key in allowed_keys:  # iteration order == CLASS_VISUALIZATION_CONFIG order
            class_rows = by_class.get(class_key)
            if not class_rows:
                continue  # class has no annotation rows yet (e.g. Class V until it's populated)

            config = self.get_visualization_class_config(class_key)
            override = self.CLASS_TREE_OVERRIDES.get(class_key, {})
            label = config["label"]
            class_family_name = class_rows[0]["class_family_name"]

            # Orphan split (curated opt-in -- see ORPHAN_SPLIT_CLASS_KEYS).
            orphan_leaf_labels = []
            if class_key in self.ORPHAN_SPLIT_CLASS_KEYS:
                is_orphan = lambda row: (row["chemotype"] or "").strip().lower() in self.ORPHAN_CHEMOTYPE_KEYS
                orphan_rows = [row for row in class_rows if is_orphan(row)]
                class_rows = [row for row in class_rows if not is_orphan(row)]
                orphan_leaf_labels = sorted(
                    {row["uniprot"] for row in orphan_rows}, key=lambda x: str(x).upper()
                )
            if not class_rows:
                continue

            color_mode = override.get("colorMode", "class")
            family_sort_key = self._natural_sort_key if override.get("family_sort") == "natural" else None
            if override.get("skip_chemotype_layer"):
                sort_keys = {0: family_sort_key} if family_sort_key else None
                children = classification_db.build_grouped_children(class_rows, [FAMILY_GETTER], sort_keys)
                tree_opts = dict(base_tree_options, colorMode=color_mode)
                meta = {"title": label, "liftClassLayer": True, "collapseLabels": ["Family"]}
            else:
                sort_keys = {1: family_sort_key} if family_sort_key else None
                children = classification_db.build_grouped_children(class_rows, [CHEMOTYPE_GETTER, FAMILY_GETTER], sort_keys)
                tree_opts = dict(base_tree_options, colorMode=color_mode)
                meta = {"title": label, "liftClassLayer": True, "collapseLabels": ["Chemotype", "Family"]}

            if orphan_leaf_labels:
                meta["orphanLeafLabels"] = orphan_leaf_labels

            tree_sets["Class"]["options"].append({"key": class_key, "label": label})
            tree_sets["Class"]["plots"][class_key] = {
                "tree": classification_db.d3_root([classification_db.d3_node(class_family_name, children)]),
                "tree_options": tree_opts,
                "meta": meta,
            }

    def _add_modality_plots(self, tree_sets, rows, base_tree_options):
        buckets = {label: [] for label in self.MODALITY_GROUP_LABELS}
        for row in rows:
            modality = (row["modality"] or "").strip()
            if not modality:
                continue
            buckets[self._group_modality_label(modality)].append(row)

        for key in self.MODALITY_GROUP_LABELS:
            bucket_rows = buckets[key]
            if not bucket_rows:
                continue
            tree_sets["Modality"]["options"].append({"key": key, "label": key})
            tree_sets["Modality"]["plots"][key] = {
                "tree": classification_db.d3_root(
                    classification_db.build_grouped_children(
                        bucket_rows, [CLASS_KEY_GETTER, FAMILY_GETTER], {0: self._class_order_sort_key}
                    )
                ),
                "tree_options": dict(base_tree_options, colorMode="class"),
                "meta": {"title": key, "liftClassLayer": False, "collapseLabels": []},
            }

    def _add_chemotype_plots(self, tree_sets, rows, base_tree_options):
        by_chem = {}
        for row in rows:
            chem = (row["chemotype"] or "").strip()
            if not chem or chem.lower() in self.TREE_EXCLUDED_CHEMOTYPES:
                continue
            by_chem.setdefault(chem, []).append(row)

        for chem in sorted(by_chem.keys(), key=lambda x: str(x).lower()):
            chem_rows = by_chem[chem]
            class_syms = {row["class_key"] for row in chem_rows}
            if len(class_syms) > 1:
                tree_opts = dict(base_tree_options, colorMode="class")
            else:
                tree_opts = dict(base_tree_options, colorMode="chemotype", forceChemotype=chem)
            tree_sets["Chemotype"]["options"].append({"key": chem, "label": chem})
            tree_sets["Chemotype"]["plots"][chem] = {
                "tree": classification_db.d3_root(
                    classification_db.build_grouped_children(
                        chem_rows, [CLASS_KEY_GETTER, FAMILY_GETTER], {0: self._class_order_sort_key}
                    )
                ),
                "tree_options": tree_opts,
                "meta": {"title": chem, "liftClassLayer": False, "collapseLabels": ["Class", "Family"]},
            }

    def _class_order_sort_key(self, class_key):
        order = self.tree_class_keys()
        return order.index(class_key) if class_key in order else len(order)

    @staticmethod
    def build_selection_context(tree_sets, requested_type="Class", requested_selection="",
                                 locked=False):
        if requested_type not in tree_sets:
            requested_type = "Class"
        available_options = tree_sets.get(requested_type, {}).get("options", [])
        available_keys = [str(opt.get("key")) for opt in available_options]
        if requested_selection not in available_keys:
            if requested_type == "Class" and "A" in available_keys:
                requested_selection = "A"
            else:
                requested_selection = available_keys[0] if available_keys else ""
        tree_locked = locked and bool(requested_selection)
        return {
            "tree_initial_type": requested_type,
            "tree_initial_selection": requested_selection,
            "tree_locked": tree_locked,
            "tree_locked_type": requested_type if tree_locked else "",
            "tree_locked_selection": requested_selection if tree_locked else "",
        }

def _strip_tags(s):
    if not s:
        return ""
    return re.sub(r"<[^>]+>", "", str(s)).strip()


def _build_gpcr_browser_rows():
    """Flat, per-receptor row list backing the "GPCR list" tab on the
    Classification overview page (classification/Classification.html)."""
    return [
        {
            "uniprot": row["uniprot"],
            "entry_name": row["entry_name"],
            "gene": row["gene"],
            "protein_name_html": row["protein_name"],
            "protein_name_text": _strip_tags(row["protein_name"]),
            "class": row["class_family_name"],
            "family": row["receptor_family"],
            "modality": row["modality"],
            "chemotype": row["chemotype"],
            "sense": row["sense"],
            "sequence": row["sequence"],
        }
        for row in sorted(
            classification_db.get_primary_classification_rows(),
            key=lambda r: str(r["uniprot"]).lower(),
        )
    ]



class ClassificationWheel:
    """
    Plain helper class (not a Django view) — provides the GPCRome wheel data
    to whichever Detail page inlines the wheel content.
    """

    @staticmethod
    def build_wheel_context():
        """
        Build the wheel-data context, used by GPCRSuperfamilyVisualizationDetail
        to inline the wheel content directly.
        """
        context = {}

        # --- Step 1: DB-derived metadata, keyed by the wheel's EntryName (uppercased UniProt stem) ---
        meta_lookup = {
            row["uniprot"]: {
                "Class": row["class_family_name"],
                "Chemotype": row["chemotype"],
                "Receptor family": row["receptor_family"],
                "Modality": row["modality"],
                "Sense": row["sense"],
            }
            for row in classification_db.get_primary_classification_rows()
        }

        # --- Step 2: Helper to inject metadata into wheel structure ---
        def enrich_wheel_with_metadata(wheelstructure):
            def recurse(node, current_class=None):
                if isinstance(node, dict):
                    for k, v in node.items():
                        if isinstance(v, dict):
                            # If we're inside a Circle_X, the keys here are actual classes (A, B1, etc.)
                            if k.startswith("Circle_"):
                                recurse(v, current_class=None)  # reset class at start of a circle
                            elif current_class is None and not "EntryName" in v:
                                # This k is the class code (A, B1, etc.)
                                recurse(v, current_class=k)
                            elif "EntryName" in v:
                                entry_code = str(v.get("EntryName", "")).strip().upper()
                                meta = meta_lookup.get(entry_code, {})
                                v.update(meta)

                                # Add the class from one level above (A, B1, etc.)
                                v["Class"] = current_class

                                if "Color" not in v:
                                    v["Color"] = "#FFFFFF"
                                if "Data" not in v:
                                    v["Data"] = ""
                            else:
                                recurse(v, current_class=current_class)
                elif isinstance(node, list):
                    for item in node:
                        recurse(item, current_class=current_class)

            recurse(wheelstructure.get("Data", {}))
            return wheelstructure


        # --- Step 3: Build the wheels ---
        odorant_wheel = DataMapperHome.GenerateGPCRomeDataStructure(data_type="Odorant")
        classic_wheel = DataMapperHome.GenerateGPCRomeDataStructure(data_type="Classic")

        # Inject metadata into both
        updated_odorant = enrich_wheel_with_metadata(odorant_wheel)
        updated_classic = enrich_wheel_with_metadata(classic_wheel)

        # --- Step 4: Pass to template ---
        context['GPCRomeData'] = json.dumps(updated_classic['Data'])
        context['GPCRomeOdorantData'] = json.dumps(updated_odorant['Data'])

        return context

class ClassSimilarityDataMixin:
    # Core fixed classes (original order)
    CLASS_ORDER = [
        "Class A (Rhodopsin)",
        "Class B1 (Secretin)",
        "Class B2 (Adhesion)",
        "Class C (Glutamate)",
        "Class F (Frizzled)",
        "Class O1 (fish-like)",
        "Class O2 (tetrapod specific)",
        "Class T2 (Taste 2)",
        "Unclassified",
    ]

    # Map display names -> top-level family slug codes
    CLASS_CODE_BY_NAME = {
        "Class A (Rhodopsin)": "001",
        "Class B1 (Secretin)": "002",
        "Class B2 (Adhesion)": "003",
        "Class C (Glutamate)": "004",
        "Class F (Frizzled)":  "006",
        "Class O1 (fish-like)": "007",
        "Class O2 (tetrapod specific)": "008",
        "Class T2 (Taste 2)": "009",
        "Unclassified": "011",
    }

    CLASS_SYMBOL_BY_NAME = {
        "Class A (Rhodopsin)": "A",
        "Class B1 (Secretin)": "B1",
        "Class B2 (Adhesion)": "B2",
        "Class C (Glutamate)": "C",
        "Class F (Frizzled)": "F",
        "Class O1 (fish-like)": "O1",
        "Class O2 (tetrapod specific)": "O2",
        "Class T2 (Taste 2)": "T2",
        "Unclassified": "U",
    }

    CLASS_COLOR_BY_SYMBOL = {
        "A": "#1f78b4",
        "B1": "#33a02c",
        "B2": "#6A3D9A",
        "C": "#d62728",
        "F": "#FF7F0E",
        "O1": "#17becf",
        "O2": "#bc80bd",
        "T2": "#F7B6D2",
        "U": "#9e9e9e",
    }

    CLASS_CLUSTER_PAYLOAD_CACHE_KEY = "classclustertree:payload:v10"
    CLASS_CLUSTER_PAYLOAD_CACHE_TIMEOUT = 60 * 60 * 24
    CLASS_CLUSTER_PAYLOAD_VERSION = 3
    CLASS_CLUSTER_DATASET_KEY = "superfamily"
    CLASS_CLUSTER_SELECTION_KEY = "global"

    def _resolve_family_ids(self, slug_codes):
        qs = ProteinFamily.objects.filter(slug__in=slug_codes).only('id', 'slug', 'name')
        return {f.slug: f.id for f in qs}

    @staticmethod
    def _resolve_unclassified_family_ids():
        """
        build_receptor_similarity excludes the single top-level "Unclassified" (slug '011')
        family and instead builds ReceptorSimilarity rows against many individual bottom-level
        unclassified families. Resolve that full descendant set here (rather than the top-level
        id alone) so "Unclassified" can be matched against the rows that actually exist.
        """
        return list(
            ProteinFamily.objects.filter(slug__startswith='011')
            .exclude(slug='011')
            .values_list('id', flat=True)
        )

    @staticmethod
    def _max_pair_similarity(similarities):
        values = []
        for score in similarities or []:
            try:
                values.append(float(score))
            except Exception:
                continue
        return max(values) if values else None

    @staticmethod
    def _format_similarity_display(similarity):
        if similarity is None:
            return "n/a"
        try:
            rounded = round(float(similarity), 1)
        except Exception:
            return "n/a"
        if not math.isfinite(rounded):
            return "n/a"
        if float(int(rounded)) == float(rounded):
            return "{}%".format(int(rounded))
        return "{:.1f}%".format(rounded)

    def _build_class_only_similarity_data(self):
        code_to_famid = self._resolve_family_ids(list(self.CLASS_CODE_BY_NAME.values()))
        unclassified_family_ids = self._resolve_unclassified_family_ids()
        classes = []
        for display_name in self.CLASS_ORDER:
            if display_name == "Unclassified":
                # Excluded from the Cluster tab's own t-SNE input (not just its display) so
                # its row doesn't skew where the other classes land in the projection.
                continue
            slug_code = self.CLASS_CODE_BY_NAME.get(display_name)
            symbol = self.CLASS_SYMBOL_BY_NAME.get(display_name, display_name)
            if slug_code == "011":
                # "Unclassified" has no ReceptorSimilarity rows against its own top-level
                # family id — aggregate across the unclassified leaf families instead.
                if not unclassified_family_ids:
                    continue
                family_ids = list(unclassified_family_ids)
                family_id = code_to_famid.get(slug_code) or 0
            else:
                family_id = code_to_famid.get(slug_code)
                if not family_id:
                    continue
                family_ids = [family_id]
            classes.append({
                "name": display_name,
                "symbol": symbol,
                "slug": slug_code,
                "family_id": family_id,
                "family_ids": family_ids,
                "color": self.CLASS_COLOR_BY_SYMBOL.get(symbol, "#808080"),
            })

        family_id_to_class_index = {}
        for idx, class_row in enumerate(classes):
            for fam_id in class_row["family_ids"]:
                family_id_to_class_index[fam_id] = idx

        allowed_class_ids = list(family_id_to_class_index.keys())
        pair_scores = defaultdict(list)
        pair_qs = (
            ReceptorSimilarity.objects
            .filter(ref_class_id__in=allowed_class_ids, target_class_id__in=allowed_class_ids)
            .values_list('ref_class_id', 'target_class_id', 'similarity')
        )
        for ref_class_id, target_class_id, similarity in pair_qs.iterator():
            idx_a = family_id_to_class_index.get(ref_class_id)
            idx_b = family_id_to_class_index.get(target_class_id)
            if idx_a is None or idx_b is None or idx_a == idx_b:
                continue
            pair_key = (min(idx_a, idx_b), max(idx_a, idx_b))
            pair_scores[pair_key].append(similarity)

        n_classes = len(classes)
        similarity_matrix = [[100.0 if i == j else None for j in range(n_classes)] for i in range(n_classes)]
        distance_matrix = np.full((n_classes, n_classes), np.nan, dtype=float)
        np.fill_diagonal(distance_matrix, 0.0)

        seen_distances = []
        for i in range(n_classes):
            for j in range(i + 1, n_classes):
                similarity = self._max_pair_similarity(pair_scores.get((i, j)))
                if similarity is None:
                    continue

                distance = max(0.0, 100.0 - similarity)
                similarity_matrix[i][j] = similarity
                similarity_matrix[j][i] = similarity
                distance_matrix[i, j] = distance
                distance_matrix[j, i] = distance
                seen_distances.append(distance)

        fill_distance = float(max(seen_distances)) if seen_distances else 100.0

        missing_pairs = 0
        for i in range(n_classes):
            for j in range(i + 1, n_classes):
                if np.isnan(distance_matrix[i, j]):
                    distance_matrix[i, j] = fill_distance
                    distance_matrix[j, i] = fill_distance
                    missing_pairs += 1

        matrix_rows = []
        for i, class_row in enumerate(classes):
            row_values = []
            for j, other_row in enumerate(classes):
                similarity = similarity_matrix[i][j]
                distance = float(distance_matrix[i, j])
                row_values.append({
                    "source": class_row["symbol"],
                    "target": other_row["symbol"],
                    "similarity": float(similarity) if similarity is not None else None,
                    "similarity_display": self._format_similarity_display(similarity),
                    "distance": distance,
                })
            matrix_rows.append(row_values)

        return {
            "classes": classes,
            "matrix": matrix_rows,
            "distance_matrix": distance_matrix,
            "missing_pairs": missing_pairs,
            "fill_distance": fill_distance,
        }

    @staticmethod
    def _manual_tsne_fallback(distance_matrix):
        n_items = int(distance_matrix.shape[0])
        if n_items <= 0:
            return np.zeros((0, 2), dtype=float)
        if n_items == 1:
            return np.array([[0.0, 0.0]], dtype=float)
        angles = np.linspace(0.0, 2.0 * math.pi, num=n_items, endpoint=False)
        return np.column_stack((np.cos(angles), np.sin(angles)))

    def _compute_tsne_coords(self, distance_matrix):
        try:
            n = int(distance_matrix.shape[0])
            if n < 2:
                perplexity = 1.0
            else:
                perplexity = min(40.0, max(1.0, (n - 1) / 3.0))
                if perplexity >= (n - 1):
                    perplexity = float(max(1, n - 2))
            tsne = TSNE(
                n_components=2,
                metric="precomputed",
                perplexity=perplexity,
                random_state=42,
                init="random",
                learning_rate="auto",
            )
            return tsne.fit_transform(distance_matrix)
        except Exception:
            return self._manual_tsne_fallback(distance_matrix)

    def _build_projection_points(self, classes, coords):
        points = []
        for idx, class_row in enumerate(classes):
            points.append({
                "id": idx,
                "symbol": class_row["symbol"],
                "label": class_row["name"],
                "color": class_row["color"],
                "x": float(coords[idx, 0]),
                "y": float(coords[idx, 1]),
            })
        return points

    def _build_scatter_methods(self, classes, distance_matrix):
        return OrderedDict([
            ("tsne", {
                "label": "t-SNE",
                "points": self._build_projection_points(classes, self._compute_tsne_coords(distance_matrix)),
            }),
        ])

    def _build_class_tip_annotations(self, classes):
        annotations = OrderedDict()
        for class_row in classes:
            symbol = str(class_row.get("symbol") or "").strip()
            if not symbol:
                continue
            annotations[symbol] = {
                "symbol": symbol,
                "label": str(class_row.get("name") or symbol),
                "color": str(class_row.get("color") or "#808080"),
                "slug": str(class_row.get("slug") or ""),
                "family_id": int(class_row.get("family_id") or 0),
            }
        return annotations

    def _build_class_tree_storage_meta(self):
        return {
            "payload_version": self.CLASS_CLUSTER_PAYLOAD_VERSION,
            "dataset_key": self.CLASS_CLUSTER_DATASET_KEY,
            "selection_key": self.CLASS_CLUSTER_SELECTION_KEY,
            "entity_type": "gpcr_class",
            "source_model": ReceptorSimilarity._meta.db_table,
            "distance_metric": "100_minus_similarity",
        }

    def _build_class_cluster_tree_payload(self):
        cached_payload = cache_alignment.get(self.CLASS_CLUSTER_PAYLOAD_CACHE_KEY)
        if cached_payload:
            return cached_payload

        class_data = self._build_class_only_similarity_data()
        classes = class_data["classes"]
        tip_annotations = self._build_class_tip_annotations(classes)
        distance_matrix = np.array(class_data["distance_matrix"], dtype=float)
        scatter_methods = self._build_scatter_methods(classes, distance_matrix)

        payload = {
            "dataset": self._build_class_tree_storage_meta(),
            "classes": classes,
            "tip_annotations": tip_annotations,
            "matrix": class_data["matrix"],
            "scatter": {
                "default_method": "tsne",
                "methods": scatter_methods,
            },
            "meta": {
                "n_classes": len(classes),
                "missing_pairs": class_data["missing_pairs"],
                "fill_distance": class_data["fill_distance"],
            },
        }
        cache_alignment.set(
            self.CLASS_CLUSTER_PAYLOAD_CACHE_KEY,
            payload,
            self.CLASS_CLUSTER_PAYLOAD_CACHE_TIMEOUT,
        )
        return payload


class CrossClassSimilarity(ClassSimilarityDataMixin):
    """
    Plain helper class (not a Django view) — provides the cross-class similarity
    matrix data to whichever Detail page inlines the matrix content.
    """

    # Five single-protein “Unclassified” items as separate groups (display order)
    SINGLE_PROTEIN_LABELS = ["GPR107", "GPR137", "TPRA1", "GPR143", "GPR157"]

    # Exact entry_name per label (case-insensitive)
    SINGLE_PROTEIN_ENTRYNAMES = {
        "GPR107": "gp107_human",
        "GPR137": "g137a_human",
        "TPRA1":  "tpra1_human",
        "GPR143": "gp143_human",
        "GPR157": "gp157_human",
    }

    # ---------- helpers ----------
    @staticmethod
    def clean_name(nm: str) -> str:
        if not nm:
            return "-"
        return (nm.replace("receptor", "")
                  .replace("-adrenoceptor", "")
                  .replace("<i>", "").replace("</i>", "")
                  .strip()) or "-"

    @staticmethod
    def primary_gene_of(p: 'Protein') -> str:
        if getattr(p, "primary_genes_self", None):
            return p.primary_genes_self[0].name
        if p.entry_name:
            return p.entry_name.split("_")[0].upper()
        return "-"

    def build_gtop_url(self, wl):
        try:
            return Template(wl.web_resource.url).substitute(index=wl.index)
        except Exception:
            return None

    def pack_hover(self, p: 'Protein') -> dict:
        """Return per-protein metadata for the tooltip."""
        wl = p.gtop_links_self[0] if getattr(p, "gtop_links_self", None) else None
        return {
            "display_name": self.clean_name(p.name),
            "gtopdb_link": self.build_gtop_url(wl) or "",
            "uniprot": p.entry_name or "",
            "gene": self.primary_gene_of(p),
            "gpcrdb_link": f"/protein/{p.entry_name}" if p.entry_name else "",
            "uniprot_link": (f"https://www.uniprot.org/uniprot/{getattr(p, 'accession', '')}"
                             if getattr(p, "accession", None) else ""),
        }

    # ---- resolvers
    def _fetch_by_entry_names(self, entry_names_lower):
        if not entry_names_lower:
            return {}
        qs = (Protein.objects
              .filter(entry_name__in=entry_names_lower)
              .only('id', 'entry_name', 'name', 'accession'))
        return {(p.entry_name or "").lower(): p for p in qs}

    def _resolve_single_proteins(self):
        """Resolve SINGLE_PROTEIN_LABELS by entry_name first, then primary gene."""
        wanted_lc = {
            lbl: (self.SINGLE_PROTEIN_ENTRYNAMES.get(lbl) or "").lower()
            for lbl in self.SINGLE_PROTEIN_LABELS
        }
        entry_to_label = {en: lbl for lbl, en in wanted_lc.items() if en}

        found = {}
        if entry_to_label:
            by_en = self._fetch_by_entry_names(list(entry_to_label.keys()))
            for en, prot in by_en.items():
                lbl = entry_to_label.get(en)
                if lbl:
                    found[lbl] = prot

        missing = [lbl for lbl in self.SINGLE_PROTEIN_LABELS if lbl not in found]
        if missing:
            wanted_genes = set(missing)
            gene_qs = (
                Protein.objects
                .prefetch_related(
                    Prefetch('genes',
                             queryset=Gene.objects.filter(position=0),
                             to_attr='primary_genes_self')
                )
                .only('id', 'entry_name', 'name', 'accession')
            )
            for p in gene_qs:
                g = self.primary_gene_of(p)
                if g in wanted_genes and g not in found:
                    found[g] = p
        return found

    # ---- main
    def get_context_data(self, **kwargs):
        context = dict(kwargs)

        # 1) Build display list (no extra/non-human groups).
        # The matrix keeps unclassified receptors as toggleable single-protein rows,
        # but does not show "Unclassified" as a standalone class.
        base_names = [name for name in self.CLASS_ORDER if name != "Unclassified"]
        single_names = [f"{lab} (Unclassified)" for lab in self.SINGLE_PROTEIN_LABELS]
        display_names = base_names + single_names

        # 2) Resolve base class families for fast class↔class aggregation
        code_to_famid = self._resolve_family_ids(list(self.CLASS_CODE_BY_NAME.values()))
        name_to_famid = {
            name: code_to_famid[self.CLASS_CODE_BY_NAME[name]]
            for name in self.CLASS_ORDER
            if self.CLASS_CODE_BY_NAME[name] in code_to_famid
        }
        allowed_class_ids = list(name_to_famid.values())

        # 3) Resolve the 5 unclassified singles
        resolved_singles = self._resolve_single_proteins()
        unclassified_name_to_protein = {}
        for lab in self.SINGLE_PROTEIN_LABELS:
            key = f"{lab} (Unclassified)"
            p = resolved_singles.get(lab)
            if p:
                unclassified_name_to_protein[key] = p

        # 4) Final groups (skip unresolved safely)
        groups = []
        for name in display_names:
            if name in name_to_famid:
                groups.append({"display": name, "kind": "class", "id": name_to_famid[name]})
            elif name in unclassified_name_to_protein:
                groups.append({"display": name, "kind": "protein", "id": unclassified_name_to_protein[name].id})
        n = len(groups)

        # ------------------------------ OPTIMIZED AGGREGATION ------------------------------
        # A) class↔class maxima (2 queries)
        base_pairs = (
            ReceptorSimilarity.objects
            .filter(ref_class_id__in=allowed_class_ids, target_class_id__in=allowed_class_ids)
            .annotate(
                pair_a=Least('ref_class_id', 'target_class_id'),
                pair_b=Greatest('ref_class_id', 'target_class_id'),
            )
            .values('pair_a', 'pair_b')
        )
        cc_max = {
            (row['pair_a'], row['pair_b']): (row['max_id'], row['max_sim'])
            for row in base_pairs.annotate(
                max_id=Max('identity'),
                max_sim=Max('similarity')
            )
        }

        # B) class↔protein maxima (1 query)
        protein_ids = [g['id'] for g in groups if g['kind'] == 'protein']
        cp_max = {}
        if allowed_class_ids and protein_ids:
            qs_cp = (
                ReceptorSimilarity.objects
                .filter(
                    (Q(ref_class_id__in=allowed_class_ids, protein_target_id__in=protein_ids)) |
                    (Q(target_class_id__in=allowed_class_ids, protein_ref_id__in=protein_ids))
                )
                .annotate(
                    canon_class_id=Case(
                        When(ref_class_id__in=allowed_class_ids, then='ref_class_id'),
                        default='target_class_id',
                        output_field=IntegerField()
                    ),
                    canon_protein_id=Case(
                        When(protein_target_id__in=protein_ids, then='protein_target_id'),
                        default='protein_ref_id',
                        output_field=IntegerField()
                    ),
                )
                .values('canon_class_id', 'canon_protein_id')
                .annotate(
                    max_id=Max('identity'),
                    max_sim=Max('similarity')
                )
            )
            cp_max = {
                (row['canon_class_id'], row['canon_protein_id']): (row['max_id'], row['max_sim'])
                for row in qs_cp
            }

        # C) protein↔protein maxima (1 query)
        pp_max = {}
        if len(protein_ids) >= 2:
            qs_pp = (
                ReceptorSimilarity.objects
                .filter(protein_ref_id__in=protein_ids, protein_target_id__in=protein_ids)
                .annotate(
                    pair_a=Least('protein_ref_id', 'protein_target_id'),
                    pair_b=Greatest('protein_ref_id', 'protein_target_id'),
                )
                .values('pair_a', 'pair_b')
                .annotate(
                    max_id=Max('identity'),
                    max_sim=Max('similarity')
                )
            )
            pp_max = {
                (row['pair_a'], row['pair_b']): (row['max_id'], row['max_sim'])
                for row in qs_pp
            }

        # 5) Build value-only matrix
        matrix = [[None for _ in range(n)] for _ in range(n)]

        def best_for(a, b, metric):
            if a['kind'] == 'class' and b['kind'] == 'class':
                key = (min(a['id'], b['id']), max(a['id'], b['id']))
                tup = cc_max.get(key)
            elif a['kind'] == 'class' and b['kind'] == 'protein':
                tup = cp_max.get((a['id'], b['id']))
            elif a['kind'] == 'protein' and b['kind'] == 'class':
                tup = cp_max.get((b['id'], a['id']))
            else:
                key = (min(a['id'], b['id']), max(a['id'], b['id']))
                tup = pp_max.get(key)
            if not tup:
                return None
            return tup[0] if metric == 'identity' else tup[1]

        # Collect tie fetch specs so we can pull them in big batches later
        tie_specs = []
        for i in range(n):
            for j in range(n):
                if i == j:
                    matrix[i][j] = None
                    continue
                a, b = groups[i], groups[j]
                metric = 'identity' if i < j else 'similarity'
                best = best_for(a, b, metric)
                matrix[i][j] = {"value": int(best) if best is not None else None,
                                "type": metric,
                                "items": []}
                if best is not None:
                    if a['kind'] == 'class' and b['kind'] == 'class':
                        tie_specs.append(('cc', (min(a['id'], b['id']), max(a['id'], b['id'])), metric, int(best)))
                    elif a['kind'] == 'class' and b['kind'] == 'protein':
                        tie_specs.append(('cp', (a['id'], b['id']), metric, int(best)))
                    elif a['kind'] == 'protein' and b['kind'] == 'class':
                        tie_specs.append(('cp', (b['id'], a['id']), metric, int(best)))
                    else:
                        tie_specs.append(('pp', (min(a['id'], b['id']), max(a['id'], b['id'])), metric, int(best)))

        # 6) Batch-fetch tie rows (kept separate per metric to avoid mixing)
        def bucket_specs(specs):
            buckets = defaultdict(list)
            for kind, ids, metric, best in specs:
                buckets[(kind, metric)].append((ids, best))
            return buckets

        buckets = bucket_specs(tie_specs)

        # Collectors: dict-of-dicts keyed by metric
        cc_ties = {'identity': defaultdict(list), 'similarity': defaultdict(list)}
        cp_ties = {'identity': defaultdict(list), 'similarity': defaultdict(list)}
        pp_ties = {'identity': defaultdict(list), 'similarity': defaultdict(list)}

        # weblink prefetch (GtoPdb)
        gtop_links_qs = WebLink.objects.select_related('web_resource').filter(web_resource__slug='gtop')

        def fetch_cc(metric):
            pairs = buckets.get(('cc', metric), [])
            if not pairs:
                return
            q = Q()
            for (a_id, b_id), best in pairs:
                cond = (Q(ref_class_id=a_id, target_class_id=b_id) |
                        Q(ref_class_id=b_id, target_class_id=a_id))
                cond &= Q(**{metric: best})
                q |= cond
            if not q.children:
                return
            rows = (
                ReceptorSimilarity.objects
                .filter(q)
                .select_related('protein_ref', 'protein_target')
                .only(
                    'identity', 'similarity',
                    'protein_ref__id', 'protein_ref__entry_name', 'protein_ref__name', 'protein_ref__accession',
                    'protein_target__id', 'protein_target__entry_name', 'protein_target__name', 'protein_target__accession',
                    'ref_class_id', 'target_class_id'
                )
                .prefetch_related(
                    Prefetch('protein_ref__genes',
                             queryset=Gene.objects.filter(position=0),
                             to_attr='primary_genes_self'),
                    Prefetch('protein_target__genes',
                             queryset=Gene.objects.filter(position=0),
                             to_attr='primary_genes_self'),
                    Prefetch('protein_ref__web_links',
                             queryset=gtop_links_qs,
                             to_attr='gtop_links_self'),
                    Prefetch('protein_target__web_links',
                             queryset=gtop_links_qs,
                             to_attr='gtop_links_self'),
                )
            )
            for r in rows:
                a = min(r.ref_class_id, r.target_class_id)
                b = max(r.ref_class_id, r.target_class_id)
                cc_ties[metric][(a, b)].append(r)

        def fetch_cp(metric):
            pairs = buckets.get(('cp', metric), [])
            if not pairs:
                return
            q = Q()
            for (cls_id, prot_id), best in pairs:
                cond = (
                    Q(ref_class_id=cls_id, protein_target_id=prot_id) |
                    Q(target_class_id=cls_id, protein_ref_id=prot_id)
                )
                cond &= Q(**{metric: best})
                q |= cond
            if not q.children:
                return
            rows = (
                ReceptorSimilarity.objects
                .filter(q)
                .select_related('protein_ref', 'protein_target')
                .only(
                    'identity', 'similarity',
                    'protein_ref__id', 'protein_ref__entry_name', 'protein_ref__name', 'protein_ref__accession',
                    'protein_target__id', 'protein_target__entry_name', 'protein_target__name', 'protein_target__accession',
                    'ref_class_id', 'target_class_id', 'protein_ref_id', 'protein_target_id'
                )
                .prefetch_related(
                    Prefetch('protein_ref__genes',
                             queryset=Gene.objects.filter(position=0),
                             to_attr='primary_genes_self'),
                    Prefetch('protein_target__genes',
                             queryset=Gene.objects.filter(position=0),
                             to_attr='primary_genes_self'),
                    Prefetch('protein_ref__web_links',
                             queryset=gtop_links_qs,
                             to_attr='gtop_links_self'),
                    Prefetch('protein_target__web_links',
                             queryset=gtop_links_qs,
                             to_attr='gtop_links_self'),
                )
            )
            for r in rows:
                if r.ref_class_id is not None and r.protein_target_id is not None:
                    key = (r.ref_class_id, r.protein_target_id)
                else:
                    key = (r.target_class_id, r.protein_ref_id)
                cp_ties[metric][key].append(r)

        def fetch_pp(metric):
            pairs = buckets.get(('pp', metric), [])
            if not pairs:
                return
            q = Q()
            for (a_id, b_id), best in pairs:
                cond = (
                    Q(protein_ref_id=a_id, protein_target_id=b_id) |
                    Q(protein_ref_id=b_id, protein_target_id=a_id)
                )
                cond &= Q(**{metric: best})
                q |= cond
            if not q.children:
                return
            rows = (
                ReceptorSimilarity.objects
                .filter(q)
                .select_related('protein_ref', 'protein_target')
                .only(
                    'identity', 'similarity',
                    'protein_ref__id', 'protein_ref__entry_name', 'protein_ref__name', 'protein_ref__accession',
                    'protein_target__id', 'protein_target__entry_name', 'protein_target__name', 'protein_target__accession',
                    'protein_ref_id', 'protein_target_id'
                )
                .prefetch_related(
                    Prefetch('protein_ref__genes',
                             queryset=Gene.objects.filter(position=0),
                             to_attr='primary_genes_self'),
                    Prefetch('protein_target__genes',
                             queryset=Gene.objects.filter(position=0),
                             to_attr='primary_genes_self'),
                    Prefetch('protein_ref__web_links',
                             queryset=gtop_links_qs,
                             to_attr='gtop_links_self'),
                    Prefetch('protein_target__web_links',
                             queryset=gtop_links_qs,
                             to_attr='gtop_links_self'),
                )
            )
            for r in rows:
                a = min(r.protein_ref_id, r.protein_target_id)
                b = max(r.protein_ref_id, r.protein_target_id)
                pp_ties[metric][(a, b)].append(r)

        # Execute the 6 batched tie fetches
        for m in ('identity', 'similarity'):
            fetch_cc(m)
            fetch_cp(m)
            fetch_pp(m)

        # 7) Fill items for tooltips (include identity & similarity per row)
        def pack_rows(rows):
            return [{
                "ref":     self.pack_hover(r.protein_ref),
                "target":  self.pack_hover(r.protein_target),
                "identity": r.identity,
                "similarity": r.similarity,
            } for r in rows]

        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                a, b = groups[i], groups[j]
                metric = matrix[i][j]["type"]            # 'identity' for upper, 'similarity' for lower
                val = matrix[i][j]["value"]
                if val is None:
                    continue

                if a['kind'] == 'class' and b['kind'] == 'class':
                    key = (min(a['id'], b['id']), max(a['id'], b['id']))
                    rows = cc_ties[metric].get(key, [])
                elif a['kind'] == 'class' and b['kind'] == 'protein':
                    rows = cp_ties[metric].get((a['id'], b['id']), [])
                elif a['kind'] == 'protein' and b['kind'] == 'class':
                    rows = cp_ties[metric].get((b['id'], a['id']), [])
                else:
                    key = (min(a['id'], b['id']), max(a['id'], b['id']))
                    rows = pp_ties[metric].get(key, [])

                matrix[i][j]["items"] = pack_rows(rows)

        # 8) Context
        context["classes"] = [g["display"] for g in groups]
        context["matrix_json"] = json.dumps(matrix)
        return context


class NewClassClusterTree(ClassSimilarityDataMixin, View):
    """
    JSON API only — serves the class-cluster scatter+matrix payload consumed by
    the Cluster tab inlined into GPCRSuperfamilyVisualizationDetail's template.
    """

    def get(self, request, *args, **kwargs):
        try:
            payload = self._build_class_cluster_tree_payload()
        except Exception as e:
            return JsonResponse({"error": str(e)}, status=400)
        return JsonResponse(payload, safe=True)

    @classmethod
    def build_cluster_tree_context(cls, request):
        """
        Build the small bit of context a Detail page needs to inline the Cluster
        tab: the JSON API URL its own JS fetches on load.
        """
        base = request.build_absolute_uri(reverse("classification-newclassclustertree"))
        return {
            "cluster_embed_url": "{}?{}".format(base, urlencode({"format": "json"})),
        }


class ReceptorFamilyVisualizationDetail(ClassificationVisualizationMixin, ClassSimilarityDataMixin, TemplateView):
    template_name = "classification/ClassificationFamilyDetail.html"

    def dispatch(self, request, *args, **kwargs):
        family_key = kwargs.get("family_key")
        family_entry, family_catalog = self.get_receptor_family_entry(family_key)
        if not family_entry:
            raise Http404("Unknown receptor family visualization")
        self.family_entry = family_entry
        self.family_catalog = family_catalog
        return super(ReceptorFamilyVisualizationDetail, self).dispatch(request, *args, **kwargs)

    def _normalize_class_display(self, raw_name):
        name = str(raw_name or "").strip()
        if not name:
            return {"key": None, "label": "Unclassified", "color": "#708090"}
        config_map = self.CLASS_VISUALIZATION_CONFIG
        reverse_map = {cfg["title"]: key for key, cfg in config_map.items()}
        key = reverse_map.get(name)
        if key:
            cfg = config_map[key]
            return {
                "key": key,
                "label": cfg["title"],
                "color": ClassSimilarityDataMixin.CLASS_COLOR_BY_SYMBOL.get(key, "#708090"),
            }
        return {"key": None, "label": name, "color": "#708090"}

    def _get_family_proteins(self):
        return self.get_visualization_family_proteins(self.family_entry)

    def _build_family_similarity_dataset(self):
        proteins = list(sorted(self._get_family_proteins(), key=lambda protein: str(protein.entry_name or "")))
        family_name = self.family_entry["name"]
        family_label = self.family_entry.get("label") or family_name
        if not proteins:
            return {"family": family_label, "proteins": [], "distance_matrix": np.zeros((0, 0), dtype=float), "fill_distance": 100.0, "missing_pairs": 0, "pair_map": {}, "identity_pair_map": {}}

        protein_ids = [p.id for p in proteins]
        pair_map = {}
        identity_pair_map = {}
        pair_qs = (
            ReceptorSimilarity.objects
            .filter(protein_ref_id__in=protein_ids, protein_target_id__in=protein_ids)
            .values("protein_ref_id", "protein_target_id", "identity", "similarity")
        )
        for rec in pair_qs.iterator():
            try:
                a = int(rec["protein_ref_id"])
                b = int(rec["protein_target_id"])
            except Exception:
                continue
            if a == b:
                continue
            key = (a, b) if a < b else (b, a)
            try:
                sim = float(rec["similarity"])
            except Exception:
                sim = None
            try:
                identity = float(rec["identity"])
            except Exception:
                identity = None
            if sim is not None:
                prev = pair_map.get(key)
                if prev is None or sim > prev:
                    pair_map[key] = sim
            if identity is not None:
                prev_identity = identity_pair_map.get(key)
                if prev_identity is None or identity > prev_identity:
                    identity_pair_map[key] = identity

        n_points = len(proteins)
        distance_matrix = np.full((n_points, n_points), np.nan, dtype=float)
        np.fill_diagonal(distance_matrix, 0.0)
        seen_distances = []
        for i in range(n_points):
            for j in range(i + 1, n_points):
                key = (min(protein_ids[i], protein_ids[j]), max(protein_ids[i], protein_ids[j]))
                similarity = pair_map.get(key)
                if similarity is None:
                    continue
                distance = max(0.0, 100.0 - similarity)
                distance_matrix[i, j] = distance
                distance_matrix[j, i] = distance
                seen_distances.append(distance)

        fill_distance = float(max(seen_distances)) if seen_distances else 100.0
        missing_pairs = 0
        for i in range(n_points):
            for j in range(i + 1, n_points):
                if np.isnan(distance_matrix[i, j]):
                    distance_matrix[i, j] = fill_distance
                    distance_matrix[j, i] = fill_distance
                    missing_pairs += 1

        return {
            "family": family_label,
            "proteins": proteins,
            "distance_matrix": distance_matrix,
            "fill_distance": fill_distance,
            "missing_pairs": missing_pairs,
            "pair_map": pair_map,
            "identity_pair_map": identity_pair_map,
        }

    def _build_family_cluster_payload(self, similarity_data):
        proteins = similarity_data["proteins"]
        family_label = similarity_data["family"]
        if not proteins:
            return {
                "family": family_label,
                "points": [],
                "meta": {
                    "n_points": 0,
                    "note": "No human receptors were found for this receptor family.",
                    "classes": self.family_entry.get("class_labels", []),
                    "chemotypes": self.family_entry.get("chemotypes", []),
                    "modality_groups": self.family_entry.get("modality_groups", []),
                },
            }

        n_points = len(proteins)
        distance_matrix = similarity_data["distance_matrix"]
        fill_distance = similarity_data["fill_distance"]
        missing_pairs = similarity_data["missing_pairs"]

        coords = (
            self._manual_tsne_fallback(distance_matrix)
            if n_points <= 2
            else self._compute_tsne_coords(distance_matrix)
        )

        points = []
        for idx, protein in enumerate(proteins):
            fam = getattr(protein.family, "parent", None)
            lig = getattr(fam, "parent", None) if fam else None
            cls = getattr(lig, "parent", None) if lig else None
            class_info = self._normalize_class_display(getattr(cls, "name", ""))
            points.append({
                "id": protein.id,
                "label": protein.entry_short(),
                "entry_name": protein.entry_name,
                "display_name": protein.short(),
                "protein_url": f"/protein/{protein.entry_name}",
                "class_label": class_info["label"],
                "class_key": class_info["key"],
                "family": family_label,
                "chemotypes": self.family_entry.get("chemotypes", []),
                "modality_groups": self.family_entry.get("modality_groups", []),
                "color": class_info["color"],
                "x": float(coords[idx, 0]),
                "y": float(coords[idx, 1]),
            })

        return {
            "family": family_label,
            "points": points,
            "meta": {
                "n_points": n_points,
                "missing_pairs": missing_pairs,
                "fill_distance": fill_distance,
                "note": "Sequence-similarity t-SNE computed on load for this receptor family.",
                "classes": self.family_entry.get("class_labels", []),
                "chemotypes": self.family_entry.get("chemotypes", []),
                "modality_groups": self.family_entry.get("modality_groups", []),
            },
        }

    def _build_family_tree_ui_payload(self, tree_payload, similarity_data):
        payload = json.loads(json.dumps(tree_payload or {}))
        proteins = similarity_data["proteins"]
        if not proteins:
            payload.setdefault("tree", "")
            payload.setdefault("annotations", {})
            payload.setdefault("Gprot_coupling", {})
            payload["entities"] = []
            payload["matrix"] = []
            payload["meta"] = dict(payload.get("meta") or {})
            payload["meta"].update({
                "n_points": 0,
                "note": "No human receptors were found for this receptor family.",
                "classes": self.family_entry.get("class_labels", []),
                "chemotypes": self.family_entry.get("chemotypes", []),
                "modality_groups": self.family_entry.get("modality_groups", []),
            })
            return payload

        protein_ids = [protein.id for protein in proteins]
        pair_map = similarity_data["pair_map"]
        identity_pair_map = similarity_data.get("identity_pair_map", {})
        distance_matrix = similarity_data["distance_matrix"]

        entities = []
        for idx, protein in enumerate(proteins):
            fam = getattr(protein.family, "parent", None)
            lig = getattr(fam, "parent", None) if fam else None
            cls = getattr(lig, "parent", None) if lig else None
            class_info = self._normalize_class_display(getattr(cls, "name", ""))
            primary_genes = getattr(protein, "primary_genes_self", None) or []
            gene_label = primary_genes[0].name if primary_genes else (
                str(protein.entry_name or "").split("_", 1)[0].upper()
            )
            entities.append({
                "symbol": protein.entry_name,
                "name": protein.short(),
                "short_label": protein.entry_short(),
                "gene_label": gene_label,
                "subtitle": class_info["label"],
                "slug": getattr(protein.family, "slug", "") or "",
                "color": class_info["color"],
                "protein_url": f"/protein/{protein.entry_name}",
            })

        matrix_rows = []
        for i, protein in enumerate(proteins):
            row = []
            for j, other_protein in enumerate(proteins):
                if i == j:
                    similarity = 100.0
                    identity = 100.0
                else:
                    pair_key = (min(protein_ids[i], protein_ids[j]), max(protein_ids[i], protein_ids[j]))
                    similarity = pair_map.get(pair_key)
                    identity = identity_pair_map.get(pair_key)
                row.append({
                    "source": protein.entry_name,
                    "target": other_protein.entry_name,
                    "identity": float(identity) if identity is not None else None,
                    "identity_display": self._format_similarity_display(identity),
                    "similarity": float(similarity) if similarity is not None else None,
                    "similarity_display": self._format_similarity_display(similarity),
                    "distance": float(distance_matrix[i, j]),
                })
            matrix_rows.append(row)

        payload.setdefault("meta", {})
        payload["meta"].update({
            "n_points": len(proteins),
            "classes": self.family_entry.get("class_labels", []),
            "chemotypes": self.family_entry.get("chemotypes", []),
            "modality_groups": self.family_entry.get("modality_groups", []),
        })
        payload["entities"] = entities
        payload["matrix"] = matrix_rows
        return payload

    def get_context_data(self, **kwargs):
        ctx = super(ReceptorFamilyVisualizationDetail, self).get_context_data(**kwargs)
        similarity_data = self._build_family_similarity_dataset()
        tree_row = (
            TreeNetwork.objects
            .filter(group_key=str(self.family_entry.get("key") or ""))
            .only("payload", "protein_count", "updated_at")
            .first()
        )
        payload = tree_row.payload if tree_row and tree_row.payload else {
            "tree": "",
            "annotations": {},
            "Gprot_coupling": {},
            "meta": {
                "n_points": 0,
                "classes": self.family_entry.get("class_labels", []),
                "chemotypes": self.family_entry.get("chemotypes", []),
                "modality_groups": self.family_entry.get("modality_groups", []),
                "note": "No persisted family tree payload is available yet. Run build_treenetwork during the data build.",
            },
        }
        tree_payload = self._build_family_tree_ui_payload(payload, similarity_data)
        cluster_payload = self._build_family_cluster_payload(similarity_data)
        ctx["page_title"] = self.family_entry["label"]
        ctx["family_tree_payload_json"] = json.dumps(tree_payload)
        ctx["family_cluster_payload_json"] = json.dumps(cluster_payload)
        ctx["family_tree_has_payload"] = bool(tree_payload.get("tree"))
        return ctx


# ----------------------------- Shared helpers ------------------------------

class OrphanSelect2Mixin:
    ORPHAN_LT = 'Orphan receptors'

    @staticmethod
    def _clean_gtop_name(nm):
        if not nm:
            return "-"
        s = (nm
             .replace("receptor", "")
             .replace("-adrenoceptor", "")
             .replace("<i>", "").replace("</i>", "")
             .strip())
        return s or "-"

    def get_orphans_select2(self):
        # pull the whole lineage to read "Class ..."
        orphans_qs = (
            Protein.objects
            .filter(
                parent_id__isnull=True,
                species_id=1,
                family__parent__parent__name__iexact=self.ORPHAN_LT,
            )
            .select_related('family__parent__parent__parent')  # <-- add this
            .prefetch_related(
                Prefetch('genes',
                         queryset=Gene.objects.filter(position=0),
                         to_attr='primary_genes_self')
            )
            .order_by('entry_name')
        )

        data = []
        for p in orphans_qs:
            gene = (
                p.primary_genes_self[0].name
                if getattr(p, 'primary_genes_self', None)
                else (p.entry_name.split('_')[0].upper() if p.entry_name else "-")
            )
            nm = self._clean_gtop_name(p.name)

            # read raw class name from lineage (e.g. "Class A orphans")
            fam = getattr(p.family, 'parent', None)
            lig = getattr(fam, 'parent', None) if fam else None
            cls = getattr(lig, 'parent', None) if lig else None
            raw_class = getattr(cls, 'name', '') or ''

            # normalize directly in Python
            if raw_class.lower().startswith("other gpcr"):
                raw_class = "Unclassified"
            elif raw_class.lower().startswith("class "):
                # keep only the letter, e.g. "Class A orphans" → "Class A"
                m = re.search(r"class\s*([a-z])", raw_class, re.I)
                raw_class = f"Class {m.group(1).upper()}" if m else raw_class

            data.append({
                "id": p.id,
                "text": nm,
                "name": nm,
                "entry_name": p.entry_name,
                "gene": gene,
                "class": raw_class,
            })
        return data


def _get_ref_class_info(ref_id):
    """
    Return (class_id, class_name) for the reference protein, inferred directly
    from ReceptorSimilarity (works regardless of row direction).
    """
    row = (
        ReceptorSimilarity.objects
        .filter(Q(protein_ref_id=ref_id) | Q(protein_target_id=ref_id))
        .values('protein_ref_id', 'protein_target_id', 'ref_class_id', 'target_class_id')
        .first()
    )
    if not row:
        return (None, None)

    if row['protein_ref_id'] == ref_id:
        class_id = row['ref_class_id']
    else:
        class_id = row['target_class_id']

    class_name = None
    if class_id:
        try:
            cls = ProteinFamily.objects.only('id', 'name').get(id=class_id)
            # Normalize like your Select2 (Unclassified / "Class X")
            name = (cls.name or '').strip()
            if name.lower().startswith('other gpcr'):
                name = 'Unclassified'
            else:
                m = re.search(r'class\s*([a-z])', name, re.I)
                if m:
                    name = f'Class {m.group(1).upper()}'
            class_name = name
        except ProteinFamily.DoesNotExist:
            pass

    return (class_id, class_name)

def _build_similarity_rows(ref_id):
    """
    Reuses the exact logic from SimilarityTopAPI to produce the table rows.
    Returns: list[dict] (the 'results' list you already send today)
    """
    from string import Template

    pairs_qs = (
        ReceptorSimilarity.objects
        .filter(Q(protein_ref_id=ref_id) | Q(protein_target_id=ref_id))
        .annotate(
            other_id=Case(
                When(protein_ref_id=ref_id, then=F('protein_target_id')),
                default=F('protein_ref_id'),
                output_field=IntegerField(),
            )
        )
        .order_by('-similarity')
        .values('other_id', 'similarity', 'identity')
    )
    pairs = list(pairs_qs)
    if not pairs:
        return []

    other_ids_all = [p['other_id'] for p in pairs]

    LT_ORPHAN = 'Orphan receptors'
    lt_map = dict(
        Protein.objects
               .filter(id__in=other_ids_all)
               .values_list('id', 'family__parent__parent__name')
    )
    def is_orphan(pid):
        return (lt_map.get(pid) or '').strip().lower() == LT_ORPHAN.lower()

    liganded_rows = [p for p in pairs if not is_orphan(p['other_id'])]

    if liganded_rows:
        cutoff = liganded_rows[9]['similarity'] if len(liganded_rows) >= 10 else liganded_rows[-1]['similarity']
    else:
        cutoff = pairs[min(9, len(pairs) - 1)]['similarity']

    kept = [p for p in pairs if p['similarity'] >= cutoff]
    kept_ids = [p['other_id'] for p in kept]

    gtop_links_qs = WebLink.objects.select_related('web_resource').filter(web_resource__slug='gtop')
    proteins_qs = (
        Protein.objects
        .filter(id__in=kept_ids)
        .select_related('family__parent__parent__parent')
        .prefetch_related(
            Prefetch('genes',
                     queryset=Gene.objects.filter(position=0),
                     to_attr='primary_genes_self'),
            Prefetch('web_links',
                     queryset=gtop_links_qs,
                     to_attr='gtop_links_self'),
        )
    )
    proteins = {p.id: p for p in proteins_qs}

    endo_qs = (
        Endogenous_GTP.objects
        .filter(receptor_id__in=kept_ids)
        .select_related('ligand', 'ligand__ligand_type')
    )
    endo_by_receptor = {}
    for e in endo_qs:
        if e.ligand:
            endo_by_receptor.setdefault(e.receptor_id, []).append(e)

    def clean_iuphar_name(nm):
        if not nm:
            return "-"
        s = nm.replace("receptor", "").replace("-adrenoceptor", "").replace("<i>", "").replace("</i>", "").strip()
        return s or "-"

    def build_gtop_url(wl):
        try:
            return Template(wl.web_resource.url).substitute(index=wl.index)
        except Exception:
            return None

    sim_map = {p['other_id']: p['similarity'] for p in kept}
    idn_map = {p['other_id']: p['identity']   for p in kept}

    results = []
    for row in kept:
        pid = row['other_id']
        p = proteins.get(pid)
        if not p:
            continue

        gene_name = p.primary_genes_self[0].name if getattr(p, 'primary_genes_self', None) else (
            (p.entry_name.split('_')[0].upper()) if p.entry_name else "-")

        entry_name  = p.entry_name or None
        gpcrdb_link = f"/protein/{entry_name}" if entry_name else "-"
        uniprot_link = f"https://www.uniprot.org/uniprot/{p.accession}" if p.accession else None

        wl_self = p.gtop_links_self[0] if getattr(p, 'gtop_links_self', None) else None
        iuphar_link = build_gtop_url(wl_self) if wl_self else None
        iuphar_name = clean_iuphar_name(p.name)

        family_name = getattr(getattr(p.family, "parent", None), "name", None)
        ligand_type = getattr(getattr(getattr(p.family, "parent", None), "parent", None), "name", None)
        clazz       = getattr(getattr(getattr(getattr(p.family, "parent", None), "parent", None), "parent", None), "name", None)

        lig_items = endo_by_receptor.get(pid, [])
        seen, endo_ligands, lig_types = set(), [], set()
        for e in lig_items:
            lig = e.ligand
            if not lig:
                continue
            if lig.id not in seen:
                seen.add(lig.id)
                endo_ligands.append({"id": lig.id, "name": lig.name})
            if lig.ligand_type:
                lig_types.add(lig.ligand_type.name)
        endo_type = "<br>".join(sorted(lig_types)) if lig_types else "-"

        results.append({
            "other_id": pid,
            "Gene": gene_name,
            "entry_name": entry_name,
            "gpcrdb_link": gpcrdb_link,
            "uniprot_link": uniprot_link,
            "iuphar_name": iuphar_name,
            "iuphar_link": iuphar_link,
            "family": family_name,
            "ligand_type": ligand_type,
            "class": clazz,
            "similarity": sim_map.get(pid, 0),
            "identity": idn_map.get(pid, 0),
            "endo_ligands": endo_ligands,
            "endo_type": endo_type,
        })

    return results

def _build_embedding_payload(ref_id, *, top_n=50, metric='identity', exclude_orphans=True):
    """
    Embedding selection logic (per your new rules), then the same t-SNE pipeline.
    Rules:
      - Unclassified  -> top N across all classes
      - Class C    -> ALL Class C
      - Other      -> top N within the same class as ref
    """
    # 0) Figure out the reference class once (using the RS table)
    ref_class_id, ref_class_name = _get_ref_class_info(ref_id)

    # If we couldn't deduce the class, fall back to the broad top N across all
    if not ref_class_id:
        base_qs = (
            ReceptorSimilarity.objects
            .filter(Q(protein_ref_id=ref_id) | Q(protein_target_id=ref_id))
            .order_by('-similarity')[:max(10, min(200, int(top_n)))]
        )
    else:
        # 1) Build class-based neighbor selection
        #    Use indexed filters directly on RS:
        #    - same-class rows regardless of direction
        same_class_q = (
            Q(protein_ref_id=ref_id, target_class_id=ref_class_id) |
            Q(protein_target_id=ref_id, ref_class_id=ref_class_id)
        )

        if ref_class_name == 'Unclassified':
            # Top N across all classes
            base_qs = (
                ReceptorSimilarity.objects
                .filter(Q(protein_ref_id=ref_id) | Q(protein_target_id=ref_id))
                .order_by('-similarity')[:max(10, min(200, int(top_n)))]
            )
        elif ref_class_name == 'Class C':
            # ALL Class C (no slice)
            base_qs = (
                ReceptorSimilarity.objects
                .filter(same_class_q)
                .order_by('-similarity')
            )
        else:
            # Same class only, Top N
            base_qs = (
                ReceptorSimilarity.objects
                .filter(same_class_q)
                .order_by('-similarity')[:max(10, min(200, int(top_n)))]
            )

    # 2) Convert to a uniform shape with other_id + values (works both directions)
    pairs = list(
        base_qs.annotate(
            other_id=Case(
                When(protein_ref_id=ref_id, then=F('protein_target_id')),
                default=F('protein_ref_id'),
                output_field=IntegerField(),
            )
        ).values('other_id', 'similarity', 'identity')
    )

    if not pairs:
        return {"points": [], "ref": {"id": ref_id, "label": ""}, "meta": {"note": "No neighbors for ref"}}

    # Note: the old 'exclude_orphans' flag is redundant here because the
    # class-based selection already governs inclusion. Kept for compatibility.

    # 3) Cap the "ranked" set only if we have a slice-less case above (Class C keeps all already)
    ranked = pairs  # already sliced where needed
    kept_ids = [ref_id] + [p["other_id"] for p in ranked]

    # 4) Fetch protein annotation for labels & legend
    proteins = (
        Protein.objects
        .filter(id__in=kept_ids)
        .select_related("family__parent__parent__parent")
    )
    pmap = {p.id: p for p in proteins}
    if ref_id not in pmap:
        return {"error": "Reference protein not found"}

    def _lab(en):
        return (en or "").replace("_human", "")

    labels, ordered_ids = [], []
    for pid in kept_ids:
        p = pmap.get(pid)
        if not p:
            continue
        lab = _lab(p.entry_name or "")
        if lab and lab not in labels:
            labels.append(lab)
            ordered_ids.append(pid)

    N = len(labels)
    if N < 3:
        return {"points": [], "ref": {"id": ref_id, "label": labels[0] if labels else ""}, "meta": {"note": "Too few points", "n_points": N}}

    id_to_idx = {pid: i for i, pid in enumerate(ordered_ids)}

    # 5) Build dense pair set among kept ids for distances
    subpairs = list(
        ReceptorSimilarity.objects
        .filter(protein_ref_id__in=kept_ids, protein_target_id__in=kept_ids)
        .values("protein_ref_id", "protein_target_id", "similarity", "identity")
    )
    pv = {}
    for r in subpairs:
        a, b = r["protein_ref_id"], r["protein_target_id"]
        if a == b:
            continue
        key = (a, b) if a < b else (b, a)
        if key not in pv:
            pv[key] = r

    D = np.full((N, N), np.nan, dtype=float)
    np.fill_diagonal(D, 0.0)

    metric = (metric or "identity").lower()
    if metric not in ("identity", "similarity"):
        metric = "identity"

    def to_dist(rec):
        val = rec["identity"] if metric == "identity" else rec["similarity"]
        try:
            v = float(val)
        except Exception:
            return np.nan
        return max(0.0, min(1.0, 1.0 - v / 100.0))

    for (a, b), rec in pv.items():
        if a in id_to_idx and b in id_to_idx:
            i, j = id_to_idx[a], id_to_idx[b]
            d = to_dist(rec)
            D[i, j] = d
            D[j, i] = d

    if np.isnan(D).any():
        col_med = np.nanmedian(D, axis=0)
        inds = np.where(np.isnan(D))
        D[inds] = np.take(col_med, inds[1])
        D = 0.5 * (D + D.T)
        np.fill_diagonal(D, 0.0)

    perplexity = max(1.0, min(40.0, (N - 1) / 3.0, N - 1 - 1e-9))
    tsne = TSNE(
        n_components=2,
        metric="precomputed",
        perplexity=perplexity,
        random_state=42,
        init="random",
        learning_rate="auto",
        square_distances=True
    )
    coords = tsne.fit_transform(D)
    used_method = "tsne"

    # 6) Values vs ref (for hover/gradient)
    ref_field = "identity" if metric == "identity" else "similarity"
    val_vs_ref = {}
    for r in ranked:
        oid = r["other_id"]
        if oid in id_to_idx:
            p = pmap.get(oid)
            if p:
                lab = _lab(p.entry_name or "")
                val = r.get(ref_field)
                if val is not None:
                    val_vs_ref[lab] = float(val)

    # 7) Build points w/ annotation
    points = []
    for i, pid in enumerate(ordered_ids):
        p = pmap[pid]
        fam = getattr(p.family, "parent", None)
        lig = getattr(fam, "parent", None) if fam else None
        cls = getattr(lig, "parent", None) if lig else None

        clazz = getattr(cls, "name", "") if cls else ""
        lig_t = getattr(lig, "name", "") if lig else ""
        fam_n = getattr(fam, "name", "") if fam else ""
        lab = labels[i]

        # Normalize class label like before
        if clazz.lower().startswith('other gpcr'):
            clazz = 'Unclassified'
        else:
            m = re.search(r'class\s*([a-z])', clazz, re.I)
            if m:
                clazz = f'Class {m.group(1).upper()}'

        fill = 100.0 if pid == ref_id else val_vs_ref.get(lab)

        points.append({
            "id": pid,
            "label": lab,
            "x": float(coords[i, 0]),
            "y": float(coords[i, 1]),
            "Class": clazz,
            "Ligand type": lig_t,
            "Receptor family": fam_n,
            "fill": float(fill) if fill is not None else None,
            "is_ref": (pid == ref_id),
        })

    return {
        "points": points,
        "ref": {"id": ref_id, "label": labels[0]},
        "meta": {"method": used_method, "metric": metric, "top_n": len(points), "n_points": N, "ref_class": ref_class_name}
    }

# ------------------------------ New merged page -----------------------------

class OrphanSimilarityExplorer(OrphanSelect2Mixin, TemplateView):
    """
    Single page that will host tabs: (1) Neighbor Table, (2) Cluster Embedding.
    The template (to be added) will keep a hidden wrapper, and after a selection
    it will call the bundle API once, show BusyLoad, then reveal the tabs.
    """
    template_name = 'class_similarity/OrphanSimilarityExplorer.html'

    def get_context_data(self, **kwargs):
        ctx = super().get_context_data(**kwargs)
        ctx['orphans_select2'] = json.dumps(self.get_orphans_select2())
        return ctx


# ------------------------------- New merged API -----------------------------

class SimilarityBundleAPI(View):
    """
    GET /class_similarity/api/bundle?ref=<protein_id>&top_n=50&metric=identity&exclude_orphans=true
    Returns BOTH:
      - table.results[]   (same structure as SimilarityTopAPI)
      - embedding.{points,ref,meta} (same structure as SimilarityEmbeddingAPI)
    Optional query flags:
      - only=table   -> return only table part
      - only=embed   -> return only embedding part
    """

    @staticmethod
    def _truthy(v):
        return str(v).lower() not in ("", "0", "false", "no", "off", "none")

    def get(self, request):
        # ---- params ----
        ref_raw = request.GET.get('ref')
        try:
            ref_id = int(ref_raw)
        except (TypeError, ValueError):
            return JsonResponse({"error": "Missing or invalid 'ref' parameter"}, status=400)

        top_n = request.GET.get("top_n")
        try:
            top_n = int(top_n) if top_n is not None else 50
        except Exception:
            top_n = 50
        top_n = max(10, min(200, top_n))

        metric = (request.GET.get("metric") or "identity").lower()
        if metric not in ("identity", "similarity"):
            metric = "identity"

        exclude_orphans = self._truthy(request.GET.get("exclude_orphans", "true"))

        only = (request.GET.get("only") or "").strip().lower()

        payload = {}

        # Build parts according to 'only'
        if only in ("", "table"):
            table_rows = _build_similarity_rows(ref_id)
            payload["table"] = {"results": table_rows}

            # If user asked only table, short-circuit
            if only == "table":
                return JsonResponse(payload)

        if only in ("", "embed"):
            embed = _build_embedding_payload(
                ref_id,
                top_n=top_n,
                metric=metric,
                exclude_orphans=exclude_orphans,
            )
            payload["embedding"] = embed

            if only == "embed":
                return JsonResponse(payload)

        return JsonResponse(payload)



# (Removed) legacy Excel export endpoint.

# ------------------------------ Structure similarity -----------------------------

class StructureSim(ClassificationVisualizationMixin, TemplateView):
    """
    Serve combined t-SNE embeddings for:
      - sequence similarity (from `classification.ReceptorSimilarity`)
      - structure distances (from `classification.StructureSimilarity`, filtered by state)

    All points are annotated from database-backed classification metadata
    (ProteinFamily + ProteinFamilyClassification).
    """

    template_name = 'classification/StructureSim.html'

    # -------------------------- DB-backed embedding + annotations --------------------------

    @staticmethod
    def _plot_method_key(plot_type):
        return plot_type

    @staticmethod
    def _neighbor_dataset_types():
        return {
            ClusterCoord.DATASET_SEQUENCE,
        }

    @staticmethod
    def _neighbor_similarity_model(dataset_type):
        return ReceptorSimilarity

    @staticmethod
    def _entry_stem(entry_name):
        if not entry_name:
            return ""
        return str(entry_name).split("_", 1)[0].upper()

    @staticmethod
    def _gene_name(protein):
        """
        Best-effort gene symbol for a protein.
        Uses prefetched `genes` if available to avoid N+1 queries.
        """
        if protein is None:
            return ""
        try:
            cache = getattr(protein, '_prefetched_objects_cache', {}) or {}
            pref = cache.get('genes')
            if pref:
                g0 = pref[0]
                return getattr(g0, 'name', '') or ''
        except Exception:
            pass
        try:
            g = protein.genes.all().first()
            return g.name if g else ""
        except Exception:
            return ""

    @staticmethod
    def _first_gene_map(protein_ids):
        """
        Return {protein_id: gene_name} for the *first* gene per protein (by Gene.position),
        fetched in a single query. This avoids N+1 when iterating with `.iterator()`.
        """
        if not protein_ids:
            return {}
        gene_map = {}
        qs = (
            Gene.objects
            .filter(proteins__id__in=protein_ids, species_id=1)
            .values_list('proteins__id', 'name', 'position')
            .order_by('proteins__id', 'position')
        )
        for prot_id, gene_name, _pos in qs.iterator():
            if prot_id not in gene_map:
                gene_map[prot_id] = gene_name or ""
        return gene_map

    def _build_similarity_context(self, protein_ids, similarity_model, top_n=10, restrict_targets_to_refs=True):
        """
        Return sequence similarity context for StructureSim:
          - top-neighbor rows per protein
        """
        protein_ids = [int(pid) for pid in protein_ids if pid]
        if not protein_ids:
            return {
                "top_neighbors": {},
            }

        ref_ids = set(protein_ids)
        rows_by_ref = defaultdict(list)
        related_ids = set(ref_ids)

        pair_filter = Q(protein_ref_id__in=ref_ids, protein_ref__species_id=1, protein_target__species_id=1)
        if restrict_targets_to_refs:
            pair_filter &= Q(protein_target_id__in=ref_ids)
        else:
            pair_filter = (
                Q(protein_ref_id__in=ref_ids) |
                Q(protein_target_id__in=ref_ids)
            ) & Q(protein_ref__species_id=1, protein_target__species_id=1)

        pair_qs = (
            similarity_model.objects
            .filter(pair_filter)
            .values('protein_ref_id', 'protein_target_id', 'similarity', 'identity')
        )

        for rec in pair_qs.iterator():
            a = rec.get('protein_ref_id')
            b = rec.get('protein_target_id')
            if not a or not b or a == b:
                continue
            related_ids.add(a)
            related_ids.add(b)
            try:
                sim = float(rec.get('similarity'))
            except Exception:
                continue
            try:
                identity = float(rec.get('identity'))
            except Exception:
                identity = 0.0

            rows_by_ref[a].append({
                'other_id': b,
                'similarity': sim,
                'identity': identity,
            })
            rows_by_ref[b].append({
                'other_id': a,
                'similarity': sim,
                'identity': identity,
            })

        if not rows_by_ref:
            return {
                "top_neighbors": {pid: [] for pid in protein_ids},
            }

        proteins = (
            Protein.objects
            .filter(id__in=related_ids)
            .select_related('family__parent__parent__parent')
        )
        protein_map = {p.id: p for p in proteins}
        gene_map = self._first_gene_map(related_ids)

        LT_ORPHAN = 'Orphan receptors'

        def is_orphan(pid):
            p = protein_map.get(pid)
            if not p:
                return False
            try:
                ligand_type = getattr(getattr(getattr(p, 'family', None), 'parent', None), 'parent', None)
                ligand_type_name = getattr(ligand_type, 'name', '') or ''
            except Exception:
                ligand_type_name = ''
            return ligand_type_name.strip().lower() == LT_ORPHAN.lower()

        def sort_key(row):
            try:
                sim = float(row.get('similarity') or 0)
            except Exception:
                sim = 0.0
            try:
                identity = float(row.get('identity') or 0)
            except Exception:
                identity = 0.0
            other_id = row.get('other_id')
            op = protein_map.get(other_id)
            label = (getattr(op, 'name', None) or getattr(op, 'entry_name', None) or '')
            return (-sim, -identity, str(label).lower())

        neighbor_map = {}
        for ref_id in protein_ids:
            rows = sorted(rows_by_ref.get(ref_id, []), key=sort_key)
            if not rows:
                neighbor_map[ref_id] = []
                continue

            liganded_rows = [row for row in rows if not is_orphan(row.get('other_id'))]
            if liganded_rows:
                cutoff_idx = min(top_n - 1, len(liganded_rows) - 1)
                cutoff = liganded_rows[cutoff_idx].get('similarity')
            else:
                cutoff_idx = min(top_n - 1, len(rows) - 1)
                cutoff = rows[cutoff_idx].get('similarity')

            kept = [row for row in rows if row.get('similarity') >= cutoff]
            out = []
            for row in kept:
                other_id = row.get('other_id')
                op = protein_map.get(other_id)
                if not op:
                    continue
                stem = self._entry_stem(getattr(op, 'entry_name', None))
                out.append({
                    'id': other_id,
                    'label': (getattr(op, 'name', None) or stem or ''),
                    'gene': gene_map.get(other_id, ''),
                    'uniprot': stem,
                    'receptor_family': (
                        getattr(getattr(op, 'family', None), 'parent', None).name
                        if getattr(getattr(op, 'family', None), 'parent', None)
                        else ''
                    ),
                    'class_name': (
                        getattr(getattr(getattr(getattr(op, 'family', None), 'parent', None), 'parent', None), 'parent', None).name
                        if getattr(getattr(getattr(getattr(op, 'family', None), 'parent', None), 'parent', None), 'parent', None)
                        else ''
                    ),
                    'similarity': row.get('similarity') or 0,
                    'identity': row.get('identity') or 0,
                })
            neighbor_map[ref_id] = out

        return {
            "top_neighbors": neighbor_map,
        }

    @staticmethod
    def _rep_structure_pdb_map(state_slug, protein_ids):
        """
        Return {protein_id: pdb_code} for the representative structure per receptor/state.

        IMPORTANT: We derive this from `classification.StructureSimilarity`, which already stores
        the (protein ↔ representative structure) linkage used to build the structure distance
        matrices for the clustering datasets. This avoids relying on other “representative” flags
        that may not match the StructureSimilarity build inputs.
        """
        if not protein_ids:
            return {}

        state_obj = ProteinState.objects.only("id").get(slug=state_slug)
        mapping = {}

        # Ref side
        qs_ref = (
            StructureSimilarity.objects
            .filter(state_id=state_obj.id, protein_ref_id__in=protein_ids)
            .values_list("protein_ref_id", "structure_ref__pdb_code__index")
            .distinct()
        )
        for prot_id, pdb in qs_ref.iterator():
            if prot_id not in mapping:
                mapping[prot_id] = (pdb or "")

        # Target side
        qs_tgt = (
            StructureSimilarity.objects
            .filter(state_id=state_obj.id, protein_target_id__in=protein_ids)
            .values_list("protein_target_id", "structure_target__pdb_code__index")
            .distinct()
        )
        for prot_id, pdb in qs_tgt.iterator():
            if prot_id not in mapping:
                mapping[prot_id] = (pdb or "")

        return mapping

    def _get_pf_classification_map_db(self):
        """
        Map ProteinFamily.id -> {'Chemotype','Modality','Sense'} using ProteinFamilyClassification.
        Delegates to the shared classification_db loader (also used by Classification_tree) so
        there is one reader of ProteinFamilyClassification, not two that can drift. A family can
        carry up to 2 rows (order 1/2); this collapses to the order-1 row since StructureSim only
        needs one chemotype/modality per family, falling back to whichever row is present.
        """
        full = classification_db.get_pf_classification_map()
        mapping = {}
        for family_id, records in full.items():
            primary = next(
                (r for r in records if r.get('chemotype_order') == 1 or r.get('modality_order') == 1),
                records[0],
            )
            mapping[family_id] = {
                'Chemotype': primary.get('Chemotype', ''),
                'Modality': primary.get('Modality', ''),
                'Sense': primary.get('Sense', ''),
            }
        return mapping

    def _protein_annotations_db(self, protein, pf_class_map):
        """
        Return annotation keys expected by `StructureSim.html` JavaScript:
        Class, Receptor family, Chemotype, Modality, Sense.
        """
        if protein is None:
            return {'Class': "", 'Receptor family': "", 'Chemotype': "", 'Modality': "", 'Sense': ""}

        pf = getattr(protein, 'family', None)

        receptor_family = ""
        clazz = ""
        try:
            receptor_family = pf.parent.name if (pf and pf.parent) else ""
        except Exception:
            receptor_family = ""
        try:
            clazz = (
                pf.parent.parent.parent.name
                if (pf and pf.parent and pf.parent.parent and pf.parent.parent.parent)
                else ""
            )
        except Exception:
            clazz = ""

        chemotype = ""
        modality = ""
        sense = ""

        cur = pf
        for _ in range(12):
            if not cur:
                break
            rec = pf_class_map.get(cur.id)
            if rec:
                chemotype = rec.get('Chemotype', "") or ""
                modality = rec.get('Modality', "") or ""
                sense = rec.get('Sense', "") or ""
                break
            cur = getattr(cur, 'parent', None)

        # Temporary label tweaks (to be removed after DB rebuild)
        try:
            if str(clazz).strip().upper() in {"OTHER GPCRS", "CLASSLESS", "UNCLASSIFIED"}:
                clazz = "Unclassified"
            if str(receptor_family).strip().upper() == "OTHER GPCR ORPHANS":
                receptor_family = "Orphan receptor"
        except Exception:
            pass

        return {
            'Class': clazz,
            'Receptor family': receptor_family,
            'Chemotype': chemotype,
            'Modality': modality,
            'Sense': sense,
        }

    def _build_sequence_dataset_db(self, pf_class_map, plot_type, group_key, restrict_neighbor_targets_to_refs=True):
        """
        Load persisted coordinates for the sequence dataset from ClusterCoord.
        """
        return self._build_clustercoord_sequence_dataset_db(
            pf_class_map=pf_class_map,
            plot_type=plot_type,
            dataset_type=ClusterCoord.DATASET_SEQUENCE,
            point_dataset="sequence",
            group_key=group_key,
            restrict_neighbor_targets_to_refs=restrict_neighbor_targets_to_refs,
            filter_entry_names=None,
        )

    def _build_clustercoord_sequence_dataset_db(
        self,
        pf_class_map,
        plot_type,
        dataset_type,
        point_dataset,
        group_key,
        restrict_neighbor_targets_to_refs,
        filter_entry_names=None,
    ):
        """
        Shared loader for sequence-like datasets stored in ClusterCoord.
        """
        rows = (
            ClusterCoord.objects
            .filter(
                dataset_type=dataset_type,
                plot_type=plot_type,
                group_key=group_key,
                protein__species_id=1,
            )
            .select_related(
                'protein',
                'protein__family',
                'protein__family__parent',
                'protein__family__parent__parent',
                'protein__family__parent__parent__parent',
            )
            .order_by('protein__entry_name')
        )
        fallback_class_key = self.class_key_from_group_key(group_key)
        if fallback_class_key == "U" and not rows.exists():
            rows = (
                ClusterCoord.objects
                .filter(
                    dataset_type=dataset_type,
                    plot_type=plot_type,
                    group_key=self.GLOBAL_GROUP_KEY,
                    protein__species_id=1,
                )
                .filter(self._visualization_class_clustercoord_q(fallback_class_key))
                .select_related(
                    'protein',
                    'protein__family',
                    'protein__family__parent',
                    'protein__family__parent__parent',
                    'protein__family__parent__parent__parent',
                )
                .order_by('protein__entry_name')
            )

        if filter_entry_names:
            rows = rows.filter(protein__entry_name__in=filter_entry_names)

        # Avoid N+1 for gene lookups: build a single mapping.
        protein_ids = list(rows.values_list('protein_id', flat=True))
        gene_map = self._first_gene_map(protein_ids)
        neighbor_context = {
            "top_neighbors": {},
        }
        if dataset_type in self._neighbor_dataset_types():
            similarity_model = self._neighbor_similarity_model(dataset_type)
            neighbor_context = self._build_similarity_context(
                protein_ids,
                similarity_model,
                restrict_targets_to_refs=restrict_neighbor_targets_to_refs,
            )

        points = []
        for r in rows.iterator():
            p = getattr(r, 'protein', None)
            stem = self._entry_stem(getattr(p, 'entry_name', None))
            gene = gene_map.get(getattr(p, 'id', None), "")
            gtop = (getattr(p, 'name', None) or stem or "")
            ann = self._protein_annotations_db(p, pf_class_map)
            points.append({
                "id": getattr(p, 'id', None),
                "label": gtop,
                "gene": gene,
                "uniprot": stem,
                "entry_name": getattr(p, 'entry_name', None) or "",
                "accession": getattr(p, 'accession', None) or "",
                "x": float(r.x),
                "y": float(r.y),
                "cluster": None,
                "dataset": point_dataset,
                "top_neighbors": neighbor_context["top_neighbors"].get(getattr(p, 'id', None), []),
                **ann,
            })

        method = self._plot_method_key(plot_type)
        return {
            "method": method,
            "points": points,
            "n": len(points),
        }

    def _build_structure_dataset_db(self, state, pf_class_map, plot_type, group_key, filter_entry_names=None):
        """
        Load persisted coordinates for the structure dataset (active/inactive)
        from ClusterCoord.
        """
        dataset_type = (
            ClusterCoord.DATASET_STRUCTURE_ACTIVE
            if state == 'active'
            else ClusterCoord.DATASET_STRUCTURE_INACTIVE
        )

        rows = (
            ClusterCoord.objects
            .filter(
                dataset_type=dataset_type,
                plot_type=plot_type,
                group_key=group_key,
                protein__species_id=1,
            )
            .select_related(
                'protein',
                'protein__family',
                'protein__family__parent',
                'protein__family__parent__parent',
                'protein__family__parent__parent__parent',
            )
            .order_by('protein__entry_name')
        )
        fallback_class_key = self.class_key_from_group_key(group_key)
        if fallback_class_key == "U" and not rows.exists():
            rows = (
                ClusterCoord.objects
                .filter(
                    dataset_type=dataset_type,
                    plot_type=plot_type,
                    group_key=self.GLOBAL_GROUP_KEY,
                    protein__species_id=1,
                )
                .filter(self._visualization_class_clustercoord_q(fallback_class_key))
                .select_related(
                    'protein',
                    'protein__family',
                    'protein__family__parent',
                    'protein__family__parent__parent',
                    'protein__family__parent__parent__parent',
                )
                .order_by('protein__entry_name')
            )

        if filter_entry_names:
            rows = rows.filter(protein__entry_name__in=filter_entry_names)

        # Avoid N+1 for gene lookups: build a single mapping.
        protein_ids = list(rows.values_list('protein_id', flat=True))
        gene_map = self._first_gene_map(protein_ids)
        pdb_map = self._rep_structure_pdb_map(state, protein_ids)

        points = []
        for r in rows.iterator():
            p = getattr(r, 'protein', None)
            stem = self._entry_stem(getattr(p, 'entry_name', None))
            gene = gene_map.get(getattr(p, 'id', None), "")
            gtop = (getattr(p, 'name', None) or stem or "")
            pdb = pdb_map.get(getattr(p, 'id', None), "")
            ann = self._protein_annotations_db(p, pf_class_map)
            points.append({
                "id": getattr(p, 'id', None),
                "label": gtop,
                "gene": gene,
                "uniprot": stem,
                "entry_name": getattr(p, 'entry_name', None) or "",
                "accession": getattr(p, 'accession', None) or "",
                "pdb": pdb,
                "x": float(r.x),
                "y": float(r.y),
                "cluster": None,
                "dataset": f"struct_{state}",
                **ann,
            })

        method = self._plot_method_key(plot_type)
        return {"method": method, "points": points, "n": len(points), "state": state}

    def _build_payload_db(self, plot_type, group_key=None, restrict_neighbor_targets_to_refs=True, filter_entry_names=None):
        group_key = group_key or self.GLOBAL_GROUP_KEY
        is_class_scoped = group_key != self.GLOBAL_GROUP_KEY
        filter_entry_names = list(filter_entry_names or [])
        pf_class_map = self._get_pf_classification_map_db()
        payload = {
            "sequence": self._build_clustercoord_sequence_dataset_db(
                pf_class_map=pf_class_map,
                plot_type=plot_type,
                dataset_type=ClusterCoord.DATASET_SEQUENCE,
                point_dataset="sequence",
                group_key=group_key,
                restrict_neighbor_targets_to_refs=restrict_neighbor_targets_to_refs,
                filter_entry_names=filter_entry_names,
            ),
            "structure": {
                "inactive": self._build_structure_dataset_db("inactive", pf_class_map, plot_type, group_key=group_key, filter_entry_names=filter_entry_names),
                "active": self._build_structure_dataset_db("active", pf_class_map, plot_type, group_key=group_key, filter_entry_names=filter_entry_names),
            },
        }
        # DB-only mode: coordinates must exist; otherwise instruct user to build them.
        missing = []
        if not (payload.get("sequence", {}).get("n") or 0):
            missing.append("sequence")
        if (
            (not is_class_scoped)
            and not (payload.get("structure", {}).get("inactive", {}).get("n") or 0)
        ):
            missing.append("structure_inactive")
        if (
            (not is_class_scoped)
            and not (payload.get("structure", {}).get("active", {}).get("n") or 0)
        ):
            missing.append("structure_active")
        if missing:
            raise ValueError(
                "Missing ClusterCoord datasets: %s (plot_type=%s, group_key=%s). Run: python manage.py build_clustercoord"
                % (", ".join(missing), plot_type, group_key)
            )
        return payload

    def _build_payload(self):
        """
        Build combined payload with:
          - sequence t-SNE
          - structure inactive / active t-SNE
        DB-backed payload (no Excel/CSV).
        """
        return self._build_payload_db(ClusterCoord.PLOT_TSNE)

    # ---------- TemplateView overrides ----------

    def get(self, request, *args, **kwargs):
        """
        - HTML by default
        - JSON when ?format=json or Accept: application/json
        """
        want_json = (
            request.GET.get('format') == 'json'
            or request.GET.get('data') == '1'
            or 'application/json' in request.headers.get('Accept', '')
        )
        try:
            scope = self.get_requested_visualization_scope(request=request, raise_404=False)
        except ValueError as e:
            if want_json:
                return JsonResponse({"error": str(e)}, status=400)
            raise Http404(str(e))
        group_key = scope["group_key"]
        tree_filter = None
        if request.GET.get("filter_type") or request.GET.get("filter_selection"):
            try:
                tree_filter = self.resolve_tree_visualization_selection(
                    request.GET.get("filter_type"),
                    request.GET.get("filter_selection"),
                )
            except ValueError as e:
                if want_json:
                    return JsonResponse({"error": str(e)}, status=400)
                raise Http404(str(e))
        filter_entry_names = tree_filter.get("entry_names", []) if tree_filter else []
        restrict_neighbor_targets_to_refs = scope["class_key"] is None or bool(filter_entry_names)
        if want_json:
            try:
                payload = self._build_payload_db(
                    ClusterCoord.PLOT_TSNE,
                    group_key=group_key,
                    restrict_neighbor_targets_to_refs=restrict_neighbor_targets_to_refs,
                    filter_entry_names=filter_entry_names,
                )
            except Exception as e:
                return JsonResponse({"error": str(e)}, status=400)
            return JsonResponse(payload, safe=True)

        return super(StructureSim, self).get(request, *args, **kwargs)

    def get_context_data(self, **kwargs):
        """
        Expose a single URL that JS can fetch combined tsne data from.
        """
        ctx = super(StructureSim, self).get_context_data(**kwargs)
        scope = self.get_requested_visualization_scope(raise_404=True)
        tree_filter = None
        if self.request.GET.get("filter_type") or self.request.GET.get("filter_selection"):
            tree_filter = self.resolve_tree_visualization_selection(
                self.request.GET.get("filter_type"),
                self.request.GET.get("filter_selection"),
            )
        base = self.request.build_absolute_uri(self.request.path)
        params = self.request.GET.copy()
        params.pop("format", None)
        params.pop("data", None)
        params["format"] = "json"
        ctx["embed_url"] = "{}?{}".format(base, params.urlencode())
        ctx["structuresim_scope_class"] = scope["class_key"] or ""
        ctx["structuresim_scope_title"] = (
            tree_filter["selection"]
            if tree_filter
            else (scope["config"]["title"] if scope["config"] else "All classes")
        )
        ctx["structuresim_embed_mode"] = str(self.request.GET.get("embed") or "").strip().lower() in {"1", "true", "yes"}
        return ctx
