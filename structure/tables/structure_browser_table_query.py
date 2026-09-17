import re

from django.db.models import (
    Count, Value, IntegerField, FloatField, ExpressionWrapper,
    Exists, OuterRef, Subquery, Prefetch,
)
from django.db.models.functions import Length, Coalesce, Cast

from structure.models import Structure, StructureExtraProteins, StructureAuxiliarySmallMolecule
from interaction.models import StructureLigandInteraction, ResidueFragmentInteraction
from protein.models import Gene, IdentifiedSites
from common.models import WebLink

# ── stabilising-agent classification (moved out of the view) ───────────────
FUSION_RE = (
    r".*thase.*|PGS|BRIL|.*Lysozyme|.*b562.*|TrxA|Flavodoxin|Rubredoxin|"
    r"Sialidase|.*Thioredoxin.*|Endolysin|.*cytochrome.*|.*DARPin.*"
)
ANTIBODY_RE = (
    r".*bod.*|.*Ab.*|.*scFv.*|.*Fab.*|.*activity.*|.*RAMP.*|.*GRK.*|"
    r"Unidentified peptide|.*CD4.*|.*IgG.*|.*NB.*|.*Fv.*"
)
fusion_pat = re.compile(FUSION_RE, re.I)
antibody_pat = re.compile(ANTIBODY_RE, re.I)


# ── subquery: distinct annotated-ligand count per structure ────────────────
def _num_annotated_ligands_subquery():
    return (
        StructureLigandInteraction.objects
        .filter(structure=OuterRef("pk"), annotated=True)
        .values("structure")
        .annotate(c=Count("ligand", distinct=True))
        .values("c")
    )


# ── subquery: residue-fragment-interaction count per structure ─────────────
def _num_interactions_subquery():
    return (
        ResidueFragmentInteraction.objects
        .filter(structure_ligand_pair__structure=OuterRef("pk"))
        .exclude(interaction_type__slug="acc")
        .values("structure_ligand_pair__structure")
        .annotate(c=Count("id"))
        .values("c")
    )


def StructureBrowserTableRows():
    """
    Yields one flat dict per Structure for table_provider.GpcrStructureBrowserTable.

    Only fields that need real per-row work are computed here: regex-classified
    stabilising agents, picking "the" arrestin/G-alpha extra-protein row, building
    the ligand/endogenous-ligand lists, filtered gene/IUPHAR lookups, and a couple
    of genuine subqueries (coverage, sodium site, interaction counts). Everything
    else (family, class, species, state, pdb, resolution, ...) is a plain FK-chain
    scalar and is read live via select_related in StructureDataJsonView instead of
    being duplicated here.
    """
    structures = (
        Structure.objects
        .filter(structure_type__origin='experiment')
        .prefetch_related(
            Prefetch(
                "protein_conformation__protein__parent__genes",
                queryset=Gene.objects.filter(position=0).select_related("entrez_weblink"),
                to_attr="filtered_genes",
            ),
            Prefetch(
                "extra_proteins",
                queryset=StructureExtraProteins.objects.select_related(
                    "wt_protein__family__parent"),
                to_attr="prefetched_extras",
            ),
            Prefetch(                          # ligands for Python loop
                "ligands",
                queryset=StructureLigandInteraction.objects.select_related(
                    "ligand", "ligand_role", "ligand__ligand_type"),
                to_attr="prefetched_ligands",
            ),
            Prefetch(
                "auxiliary_small_molecules",
                queryset=StructureAuxiliarySmallMolecule.objects.all(),
                to_attr="prefetched_aux_small_molecules",
            ),
            Prefetch(
                "protein_conformation__protein__parent__web_links",
                queryset=WebLink.objects
                    .select_related("web_resource")
                    .filter(web_resource__slug="gtop"),
                to_attr="prefetched_gtop_links",
            ),
            "stabilizing_agents",               # regex filter
            "protein_conformation__protein__parent__endogenous_gtp_set__ligand__ligand_type",
        )
        .annotate(
            coverage_pct=ExpressionWrapper(
                Cast(Count("protein_conformation__residue", distinct=True),
                     FloatField()) * 100.0 /
                Coalesce(
                    Length("protein_conformation__protein__parent__sequence"),
                    Value(1.0)
                ),
                output_field=IntegerField(),
            ),
            has_sodium_site=Exists(
                IdentifiedSites.objects.filter(
                    protein_conformation=OuterRef("protein_conformation_id"),
                    site__slug="sodium_pocket")
            ),
            has_ligand_interactions=Exists(
                ResidueFragmentInteraction.objects.filter(
                    structure_ligand_pair__structure=OuterRef("pk"))
            ),
            num_annotated_ligands=Coalesce(
                Subquery(_num_annotated_ligands_subquery(), output_field=IntegerField()), 0),
            num_interactions=Coalesce(
                Subquery(_num_interactions_subquery(), output_field=IntegerField()), 0),
        )
    )

    for s in structures:
        pp = s.protein_conformation.protein.parent

        # arrestin / Gα — pick the row, store only its id (display formatting
        # happens at read time off the FK)
        arrestin = next(
            (ep for ep in getattr(s, "prefetched_extras", [])
             if ep.category in {"G alpha", "Arrestin"}), None
        )

        # stabilising agents – tiny regex pass
        fusions = "<br>".join(a.name for a in s.stabilizing_agents.all()
                              if fusion_pat.match(a.name)) or "-"
        antibodies = "<br>".join(a.name for a in s.stabilizing_agents.all()
                                 if antibody_pat.match(a.name)) or "-"

        # auxiliary small molecules – dedupe by title for the Name list,
        # distinct sorted sets for Type/Function (same shape as fusions/antibodies)
        aux_names, aux_types, aux_functions, seen_aux = [], set(), set(), set()
        for am in getattr(s, "prefetched_aux_small_molecules", []):
            title = am.title or am.name
            if title and title not in seen_aux:
                seen_aux.add(title)
                aux_names.append(title)
            if am.type:
                aux_types.add(am.type)
            if am.function:
                aux_functions.add(am.function)

        # ligands – build names/types/roles here
        lig_list, lig_types, lig_roles = [], set(), set()
        for li in getattr(s, "prefetched_ligands", []):
            if not li.ligand:
                continue
            lig_list.append({"id": li.ligand.id, "name": li.ligand.name})
            if li.ligand.ligand_type:
                lig_types.add(li.ligand.ligand_type.name)
            if li.ligand_role:
                lig_roles.add(li.ligand_role.name)

        # endogenous ligands from parent protein
        endos = getattr(pp, "endogenous_gtp_set", []).all() if hasattr(pp, "endogenous_gtp_set") else []
        endo_list, seen = [], set()
        for e in endos:
            lig = getattr(e, "ligand", None)
            if not lig:
                continue
            key = getattr(lig, "id", None) or lig.name
            if key in seen:
                continue
            seen.add(key)
            endo_list.append({"id": lig.id, "name": lig.name})

        # gene (filtered M2M — not a plain FK chain, so it's precomputed)
        gene_name = (pp.filtered_genes[0].name
                     if getattr(pp, "filtered_genes", []) else None)
        gene_entrez_url = (str(pp.filtered_genes[0].entrez_weblink)
                            if getattr(pp, "filtered_genes", []) and
                            pp.filtered_genes[0].entrez_weblink else None)

        # IUPHAR / GToP link (filtered reverse FK — same reasoning)
        iuphar_link, iuphar_index = None, None
        wl = next(iter(getattr(pp, "prefetched_gtop_links", [])), None)
        if wl and wl.web_resource and wl.index:
            iuphar_index = wl.index
            iuphar_link = wl.web_resource.url.replace("$index", str(wl.index))

        yield {
            "structure_id": s.id,
            "arrestin_extra_protein_id": arrestin.id if arrestin else None,
            "coverage": int(s.coverage_pct) if s.coverage_pct is not None else None,
            "sodium_site": bool(s.has_sodium_site),

            "fusions": fusions,
            "antibodies": antibodies,

            "auxiliary_molecules": "<br>".join(aux_names) or "-",
            "auxiliary_molecule_type": "<br>".join(sorted(aux_types)) or "-",
            "auxiliary_molecule_function": "<br>".join(sorted(aux_functions)) or "-",

            "ligands": lig_list,
            "ligand_type": "<br>".join(map(str, sorted(lig_types))) or "-",
            "ligand_role": "<br>".join(map(str, sorted(lig_roles))) or "-",

            "endo_ligands": endo_list,
            "endo_type": "<br>".join(sorted({e.ligand.ligand_type.name for e in endos
                                              if e.ligand and e.ligand.ligand_type})) or "-",

            "gene_name": gene_name,
            "gene_entrez_url": gene_entrez_url,
            "iuphar_link": iuphar_link,
            "iuphar_index": iuphar_index,

            "has_ligand_interactions": bool(s.has_ligand_interactions),
            "num_annotated_ligands": s.num_annotated_ligands,
            "num_interactions": s.num_interactions,
        }
