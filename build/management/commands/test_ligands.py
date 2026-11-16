from build.management.commands.base_build import Command as BaseBuild
import random
from collections import Counter, defaultdict

from rdkit import Chem

from django.db.models import Count, F, Q
from django.db.models.functions import Lower

from ligand.models import Ligand, LigandID, AssayExperiment, Endogenous_GTP
from build.management.commands.build_ligand_functions import get_radioactive_info
from build.management.commands.build_ligand_functions import load_gtp_source_data
from build.management.commands.ligand_normalization import normalize_for_parent
from build.management.commands.build_ligand_functions import _PEPTIDE_BACKBONE_SMARTS


class Command(BaseBuild):
    def handle(self, *args, **options):
        ligs_children = Ligand.objects.filter(parent__isnull=False).prefetch_related('parent', 'ids')
        ligs_all = Ligand.objects.all().prefetch_related('ligand_type')
        c=1

        print("=========================================")
        print("Checking for empty names")
        empty_names = ligs_all.filter(name="")
        if len(empty_names)>0:
            for e in empty_names:
                print(e.gpcrdb_id)
        else:
            print("None found")

        iexact_dict = {}
        radioactive_dict = {}

        # Parent clean_inchikey is the first block of the normalize_for_parent key; child
        # clean_inchikey is the full normalize_for_child key (same pipeline, stereo/isotopes kept),
        # so their first blocks should match. Only mismatches are re-normalized with RDKit.
        norm_cache = {}
        def _norm(smiles):
            if not smiles:
                return None
            if smiles not in norm_cache:
                norm_cache[smiles] = normalize_for_parent(smiles)
            return norm_cache[smiles]

        def _norm_key(smiles):
            norm = _norm(smiles)
            return norm.get("clean_inchikey") if norm else None

        def _formula(smiles):
            # Element counts (incl. H) of the parent-normalized structure
            norm = _norm(smiles)
            mol = Chem.MolFromSmiles(norm["smiles"]) if norm and norm.get("smiles") else None
            if mol is None:
                return None
            counts = Counter()
            for atom in mol.GetAtoms():
                counts[atom.GetSymbol()] += 1
                counts["H"] += atom.GetTotalNumHs()
            return counts

        def _formula_str(counts):
            if not counts:
                return "?"
            order = ["C", "H"] + sorted(k for k in counts if k not in ("C", "H"))
            return "".join(f"{k}{counts[k] if counts[k] != 1 else ''}" for k in order if counts[k])

        def _formula_diff_class(child, parent):
            # Classify child-vs-parent formula differences that still mean "same ligand family":
            #   'disulfide'  - only H differs, by an even number (open vs closed disulfide bridges/rings)
            #   'iodinated'  - radio-iodinated child: I added in place of H or another halogen
            cf, pf = _formula(child.smiles), _formula(parent.smiles)
            if not cf or not pf:
                return None, cf, pf
            diff = {k: cf[k] - pf[k] for k in set(cf) | set(pf) if cf[k] != pf[k]}
            if set(diff) == {"H"} and diff["H"] % 2 == 0:
                return "disulfide", cf, pf
            if (child.radioactive and child.radioactive.endswith("I") and diff.get("I", 0) > 0
                    and all(k in ("H", "F", "Cl", "Br") and v < 0 for k, v in diff.items() if k != "I")
                    and sum(v for k, v in diff.items() if k != "I") == -diff["I"]):
                return "iodinated", cf, pf
            return None, cf, pf

        roots_by_ck = defaultdict(list)
        for pid, ck, gid, pname in (Ligand.objects.filter(parent__isnull=True, clean_inchikey__isnull=False)
                                    .values_list('id', 'clean_inchikey', 'gpcrdb_id', 'name')):
            roots_by_ck[ck].append((pid, gid, pname))

        def _roots_str(*keys):
            hits = {gid: pname for k in keys if k for _, gid, pname in roots_by_ck.get(k, [])}
            return ", ".join(f"{gid} {pname!r}" for gid, pname in sorted(hits.items())) or "none"

        key_mismatches = []  # (child, child key, key came from raw inchikey)
        n_no_key = 0
        for l in ligs_children:
            ### Collect children whose key's first block differs from the parent's clean_inchikey
            lp = l.parent.clean_inchikey
            if lp:
                if l.clean_inchikey:
                    ck, from_raw = l.clean_inchikey.split('-')[0], False
                elif l.inchikey:
                    ck, from_raw = l.inchikey.split('-')[0], True
                else:
                    ck = None
                if ck is None:
                    n_no_key += 1
                elif ck != lp:
                    key_mismatches.append((l, ck, from_raw))

            isotopes = get_radioactive_info(l.name, [], l.smiles)
            if len(isotopes)>0:
                radioactive_dict[l] = isotopes

        # Re-normalize only the mismatches: if the stereo-free parent normalization of the child
        # SMILES gives the parent key, the assignment is right but normalize_for_child and
        # normalize_for_parent disagree (stereo influencing tautomer canonicalization).
        misassigned, inconsistent, disulfide, iodinated = [], [], [], []
        for l, ck, from_raw in key_mismatches:
            nk = _norm_key(l.smiles)
            if nk and nk == l.parent.clean_inchikey:
                inconsistent.append((l, ck, nk, from_raw))
                continue
            diff_class, cf, pf = _formula_diff_class(l, l.parent)
            if diff_class == "disulfide":
                disulfide.append((l, cf, pf))
            elif diff_class == "iodinated":
                iodinated.append((l, cf, pf))
            else:
                misassigned.append((l, ck, nk, from_raw, cf, pf))

        print("=========================================")
        print("Checking for incorrectly assigned parents (child clean_inchikey first block vs parent clean_inchikey)")
        print("child gpcrdb_id | name | source || current parent gpcrdb_id | name | clean_inchikey || child key (raw=from inchikey) | "
              "child parent-normalized key if different | same seq as parent | child formula / parent formula || root parents with either key")
        for l, ck, nk, from_raw, cf, pf in misassigned:
            same_seq = "same seq" if l.sequence and l.sequence == l.parent.sequence else "-"
            print(f"{l.gpcrdb_id} | {l.name!r} | {l.source} || {l.parent.gpcrdb_id} | {l.parent.name!r} | {l.parent.clean_inchikey} || "
                  f"{ck}{' (raw)' if from_raw else ''} | {nk if nk != ck else '-'} | {same_seq} | {_formula_str(cf)} / {_formula_str(pf)} || {_roots_str(ck, nk)}")
            c+=1
        print(f"Total: {len(misassigned)}  |  children without clean_inchikey/inchikey (not checked): {n_no_key}")

        print("=========================================")
        print("Checking for normalizer stereo/tautomer inconsistencies (child key differs, but parent-normalized child SMILES matches parent)")
        print("child gpcrdb_id | name | source || parent gpcrdb_id | name || child key (raw=from inchikey) | parent key")
        for l, ck, nk, from_raw in inconsistent:
            print(f"{l.gpcrdb_id} | {l.name!r} | {l.source} || {l.parent.gpcrdb_id} | {l.parent.name!r} || "
                  f"{ck}{' (raw)' if from_raw else ''} | {nk}")
        print(f"Total: {len(inconsistent)}")

        print("=========================================")
        print("Children differing from parent only in H count (open vs closed disulfide bridges; same ligand, expected)")
        print("child gpcrdb_id | name | source || parent gpcrdb_id | name || child formula / parent formula")
        for l, cf, pf in disulfide:
            print(f"{l.gpcrdb_id} | {l.name!r} | {l.source} || {l.parent.gpcrdb_id} | {l.parent.name!r} || {_formula_str(cf)} / {_formula_str(pf)}")
        print(f"Total: {len(disulfide)}")

        print("=========================================")
        print("Iodinated radioligands under their non-iodinated parent (expected)")
        print("child gpcrdb_id | name | radioactive || parent gpcrdb_id | name || child formula / parent formula")
        for l, cf, pf in iodinated:
            print(f"{l.gpcrdb_id} | {l.name!r} | {l.radioactive} || {l.parent.gpcrdb_id} | {l.parent.name!r} || {_formula_str(cf)} / {_formula_str(pf)}")
        print(f"Total: {len(iodinated)}")

        #     iexact_dict[l] = Ligand.objects.filter(name__iexact=l.name, parent__isnull=False)

        # print("===================================================")
        # print("Checking iexact name matches with different parents")
        # for l, li in iexact_dict.items():
        #     if len(li)>1 and len(li.order_by("parent").distinct("parent"))>1:
        #         print('===========')
        #         for lig in li:
        #             print(lig, lig.gpcrdb_id, lig.parent, lig.parent.gpcrdb_id)

        groups = defaultdict(list)

        for l in ligs_children.select_related("parent"):
            key = l.name.lower()
            groups[key].append(l)

        print("===================================================")
        print("Checking iexact name matches with different parents")

        for name, ligs in groups.items():
            parents = {lig.parent_id for lig in ligs}

            if len(ligs) > 1 and len(parents) > 1:
                print("===========")

                for lig in ligs:
                    print(
                        lig,
                        lig.gpcrdb_id,
                        lig.parent,
                        lig.parent.gpcrdb_id
                    )

        print("=====================================================================================")
        print("Checking for radioactive isotope ligand names missing radioactive isotope")
        for l, isotopes in radioactive_dict.items():
            if len(isotopes)>1:
                print(f'WARNING: more than 1 isotope data ({isotopes}) for {l} {l.gpcrdb_id}')
            if l.radioactive!=isotopes[0]:
                print(f"Missing radioactive info for {l} GPCRdbID {l.gpcrdb_id}, source: {l.source}")
                print(isotopes)
                print(l.radioactive)
        # Names saying "radiolabel" without an isotope should carry radioactive="Yes"
        for l in Ligand.objects.filter(parent__isnull=False, name__iregex=r"radiolabel", radioactive__isnull=True):
            print(f"Missing radioactive flag for radiolabel-named {l} GPCRdbID {l.gpcrdb_id}, source: {l.source}")

        print("================================================================")
        print("Checking for peptides assigned as small molecules and nan values")
        for l in ligs_all:
            ### Check for incorrect ligand type
            if l.sequence and l.ligand_type.slug not in ['peptide', 'protein']:
                print("PEPTIDE", l.name, l.sequence, l.helm, l.ligand_type.slug)

            ### Check for nan values
            count = False
            if l.clean_inchikey=='nan':
                l.clean_inchikey = None
                count = True
            if l.inchikey=='nan':
                l.inchikey = None
                count = True
            if l.smiles=='nan':
                l.smiles = None
                count = True
            if l.sequence=='nan':
                l.sequence = None
                count = True
            if l.helm=='nan':
                l.helm = None
                count = True
            if count:
                c+=1
                print("----nan", l, l.gpcrdb_id, l.source)
            # l.save()

            ### Check for missing gpcrdb id
            if l.gpcrdb_id==None or not isinstance(l.gpcrdb_id, int):
                print("++++no GPCRdb id", l, l.id)

        print("=====================================================================================")
        print("Checking for LigandID objects assigned to multiple ligand objects")

        ligand_ids = (
            LigandID.objects.all().exclude(web_resource__slug="uniprot")
            .select_related('web_resource', 'ligand')
        )

        groups = defaultdict(list)

        for lid in ligand_ids:
            key = (lid.web_resource.slug, lid.index)
            groups[key].append(lid)

        for (slug, index), d_set in groups.items():

            if (
                len(d_set) == 2 and
                (
                    (d_set[0].ligand.radioactive is not None and d_set[1].ligand.radioactive is None)
                    or
                    (d_set[1].ligand.radioactive is not None and d_set[0].ligand.radioactive is None)
                )
            ):
                # Radioactive entries
                pass

            elif len(d_set) > 1:
                print('--------------')

                for ds in d_set:
                    print(
                        ds.web_resource.slug,
                        ds.index,
                        ds.ligand.name,
                        ds.ligand.gpcrdb_id,
                        ds.ligand.radioactive
                    )

        print("=====================================================================================")
        print("Checking for GtP ligand IDs (with GPCR data) missing from database")
        gtp = load_gtp_source_data()
        source_ids = set(gtp["ligand_ids"]) | set(gtp["peptide_ligand_ids"])
        db_ids = set(
            LigandID.objects
            .filter(web_resource__slug="gtoplig")
            .values_list("index", flat=True)
        )
        missing = source_ids - db_ids
        if missing:
            name_map = (
                gtp["gtp_complete_ligands"]
                .set_index("ligand_id")[["name", "type"]]
                .to_dict("index")
            )
            for mid in sorted(missing, key=lambda x: int(x) if x.isdigit() else 0):
                info = name_map.get(mid, {})
                print(f"  Missing gtoplig {mid}: {info.get('name', '?')} (type={info.get('type', '?')})")
        else:
            print("None missing")

        print("=====================================================================================")
        print("Checking number of children and parents")
        print("Parents: ", ligs_all.filter(parent__isnull=True).count())
        print("Children: ", ligs_children.count())

        print("=====================================================================================")
        print("Checking for parent entries with no children")
        childless_parents = (
            Ligand.objects
            .filter(parent__isnull=True)
            .annotate(child_count=Count('children'))
            .filter(child_count=0)
            .select_related('ligand_type')
            .order_by('ligand_type__slug', 'name')
        )
        total_childless = childless_parents.count()
        print(f"Total: {total_childless}")
        for l in childless_parents:
            type_slug = l.ligand_type.slug if l.ligand_type else 'unknown'
            print(f"  [{type_slug}] gpcrdb_id={l.gpcrdb_id}  name={l.name!r}  "
                  f"uniprot={l.uniprot}  source={l.source}")

        print("=====================================================================================")
        print("Checking for child entries with no parent (parent_id IS NULL, no children, non-small-molecule)")
        orphaned_pep_prot = childless_parents.filter(
            ligand_type__slug__in=['peptide', 'protein']
        )
        print(f"Peptide/protein orphans: {orphaned_pep_prot.count()} of {total_childless} childless roots")
        for l in orphaned_pep_prot:
            print(f"  gpcrdb_id={l.gpcrdb_id}  name={l.name!r}  "
                  f"uniprot={l.uniprot}  sequence={str(l.sequence or '')[:40]!r}  source={l.source}")

        print("=====================================================================================")
        print("Checking for duplicate child entries with same SMILES or clean_inchikey under the same parent")
        dup_smiles = (
            Ligand.objects
            .filter(parent__isnull=False, smiles__isnull=False)
            .values("parent_id", "smiles")
            .annotate(cnt=Count("id"))
            .filter(cnt__gt=1)
        )
        dup_ck = (
            Ligand.objects
            .filter(parent__isnull=False, clean_inchikey__isnull=False)
            .values("parent_id", "clean_inchikey")
            .annotate(cnt=Count("id"))
            .filter(cnt__gt=1)
        )
        dup_found = False
        for row in dup_smiles:
            siblings = Ligand.objects.filter(
                parent_id=row["parent_id"], smiles=row["smiles"]
            ).prefetch_related("ids__web_resource")
            sibling_list = list(siblings)
            # Variants of one ligand with different radiolabels (incl. unlabelled) share the
            # unlabelled structure but are distinct ligands, e.g. [125I]SP and [3H]SP
            if len({s.radioactive for s in sibling_list}) == len(sibling_list):
                continue
            dup_found = True
            print(f"  Duplicate SMILES under parent_id={row['parent_id']}:")
            for s in siblings:
                pubchem_ids = [lid.index for lid in s.ids.all() if lid.web_resource.slug == "pubchem"]
                print(f"    gpcrdb_id={s.gpcrdb_id}  pubchem={pubchem_ids}  clean_inchikey={s.clean_inchikey}")
        for row in dup_ck:
            siblings = Ligand.objects.filter(
                parent_id=row["parent_id"], clean_inchikey=row["clean_inchikey"]
            ).prefetch_related("ids__web_resource")
            sibling_list = list(siblings)
            # Variants of one ligand with different radiolabels (incl. unlabelled) share the
            # unlabelled structure but are distinct ligands, e.g. [125I]SP and [3H]SP
            if len({s.radioactive for s in sibling_list}) == len(sibling_list):
                continue
            dup_found = True
            print(f"  Duplicate clean_inchikey={row['clean_inchikey']} under parent_id={row['parent_id']}:")
            for s in siblings:
                pubchem_ids = [lid.index for lid in s.ids.all() if lid.web_resource.slug == "pubchem"]
                print(f"    gpcrdb_id={s.gpcrdb_id}  pubchem={pubchem_ids}")
        if not dup_found:
            print("  None found")

        parents = Ligand.objects.filter(parent__isnull=True)
        children = Ligand.objects.filter(parent__isnull=False)

        print("=====================================================================================")
        print("Checking for parents sharing a clean_inchikey (duplicate parents)")
        dup_parents = (parents.exclude(clean_inchikey__isnull=True)
                       .values("clean_inchikey").annotate(n=Count("id")).filter(n__gt=1))
        for row in dup_parents:
            print(f"  {row['clean_inchikey']}: " + ", ".join(
                f"{p.gpcrdb_id} {p.name!r}" for p in parents.filter(clean_inchikey=row["clean_inchikey"])))
        print(f"Total: {dup_parents.count()}")

        print("=====================================================================================")
        print("Checking parent key self-consistency (random sample: stored clean_inchikey vs normalize_for_parent(smiles))")
        rng = random.Random(0)
        parent_ids = list(parents.filter(smiles__isnull=False).values_list("id", flat=True))
        sample_ids = rng.sample(parent_ids, min(1500, len(parent_ids)))
        n_inconsistent = 0
        for p in parents.filter(id__in=sample_ids).only("gpcrdb_id", "name", "smiles", "clean_inchikey"):
            nk = _norm_key(p.smiles)
            if nk and nk != p.clean_inchikey:
                n_inconsistent += 1
                print(f"  {p.gpcrdb_id} {p.name!r} stored={p.clean_inchikey} normalized={nk}")
        print(f"Total: {n_inconsistent} of {len(sample_ids)} sampled")

        print("=====================================================================================")
        print("Checking structure-less parents whose name matches exactly one structured family")
        amb_parents = {p.pk: p for p in parents.filter(ambiguous_alias=True)}
        amb_names = defaultdict(set)
        for pk, p in amb_parents.items():
            if p.name and p.name.strip():
                amb_names[pk].add(p.name.strip().lower())
        for pid, name in children.filter(parent_id__in=list(amb_parents)).values_list("parent_id", "name"):
            if name and name.strip():
                amb_names[pid].add(name.strip().lower())
        all_amb_names = set().union(*amb_names.values()) if amb_names else set()
        structured_roots = defaultdict(set)
        for lname, lid, pid, amb, pamb in (Ligand.objects.annotate(lname=Lower("name"))
                                           .filter(lname__in=all_amb_names)
                                           .values_list("lname", "id", "parent_id", "ambiguous_alias", "parent__ambiguous_alias")):
            root, root_amb = (pid, pamb) if pid else (lid, amb)
            if not root_amb:
                structured_roots[lname.strip()].add(root)
        n_single = 0
        for pk, names in amb_names.items():
            targets = set().union(*(structured_roots.get(n, set()) for n in names)) - {pk}
            if len(targets) == 1:
                n_single += 1
                target = Ligand.objects.get(pk=targets.pop())
                print(f"  {amb_parents[pk].gpcrdb_id} {amb_parents[pk].name!r} -> {target.gpcrdb_id} {target.name!r}")
        print(f"Total: {n_single} (of {len(amb_parents)} structure-less parents)")

        print("=====================================================================================")
        print("Checking for children whose ligand type differs from their parent's")
        type_diff = children.exclude(ligand_type=F("parent__ligand_type")).select_related("ligand_type", "parent__ligand_type")
        for l in type_diff[:50]:
            print(f"  {l.gpcrdb_id} {l.name!r} {l.ligand_type.slug if l.ligand_type else None} "
                  f"vs parent {l.parent.gpcrdb_id} {l.parent.ligand_type.slug if l.parent.ligand_type else None}")
        print(f"Total: {type_diff.count()}")

        print("=====================================================================================")
        print("Checking for parents named after a radiolabel")
        radio_named = parents.filter(name__regex=r"\[<sup>\d+m?</sup>[A-Z][a-z]?[],]")
        with_unlabelled = radio_named.filter(children__radioactive__isnull=True).distinct()
        for p in with_unlabelled[:50]:
            print(f"  {p.gpcrdb_id} {p.name!r} has an unlabelled child")
        print(f"Total: {radio_named.count()} radiolabel-named parents, {with_unlabelled.count()} with an unlabelled child")

        print("=====================================================================================")
        print("Checking for bioactivity / endogenous rows pointing to a parent ligand")
        print(f"AssayExperiment: {AssayExperiment.objects.filter(ligand__parent__isnull=True).count()}")
        print(f"Endogenous_GTP: {Endogenous_GTP.objects.filter(ligand__parent__isnull=True).count()}")

        print("=====================================================================================")
        print("Checking for malformed LigandID values")
        id_patterns = {"chembl_ligand": r"^CHEMBL\d+$", "pubchem": r"^\d+$", "drugbank": r"^DB\d{5}$",
                       "gtoplig": r"^\d+$", "drug_central": r"^\d+$"}
        n_malformed = 0
        for slug, pattern in id_patterns.items():
            bad = LigandID.objects.filter(web_resource__slug=slug).exclude(index__regex=pattern)
            for lid in bad.select_related("ligand")[:20]:
                print(f"  {slug} {lid.index!r} on {lid.ligand.gpcrdb_id} {lid.ligand.name!r} ({lid.ligand.source})")
            n_malformed += bad.count()
        print(f"Total: {n_malformed}")

        print("=====================================================================================")
        print("Checking for non-radioactive children with SMILES but no clean_inchikey")
        no_key = children.filter(smiles__isnull=False, clean_inchikey__isnull=True, radioactive__isnull=True)
        for l in no_key[:50]:
            print(f"  {l.gpcrdb_id} {l.name!r} {l.source}")
        print(f"Total: {no_key.count()}")

        print("=====================================================================================")
        print("Informational: parents typed small-molecule with >=3 peptide backbone units (by child source)")
        pep_like = Counter()
        for smiles, srcs in ((p.smiles, {c.source for c in p.children.all()})
                             for p in parents.filter(ligand_type__slug="small-molecule", smiles__isnull=False)
                             .only("smiles").prefetch_related("children").iterator(chunk_size=2000)):
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None and len(mol.GetSubstructMatches(_PEPTIDE_BACKBONE_SMARTS)) >= 3:
                for s in srcs:
                    pep_like[s] += 1
        print(f"By source: {dict(pep_like)}")
