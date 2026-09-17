from build.management.commands.base_build import Command as BaseBuild
from django.db.models import Count, Q

from ligand.models import Ligand, LigandID
from build.management.commands.build_ligand_functions import get_radioactive_info
from build.management.commands.build_ligand_functions import load_gtp_source_data


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

        print("=========================================")
        print("Checking for incorrectly assigned parents")
        for l in ligs_children:
            ### Check for incorrectly assigned parents
            if l.inchikey and l.parent.clean_inchikey:
                li = l.inchikey.split('-')[0]
                lp = l.parent.clean_inchikey
                try:
                    supposed_to_be_parent = Ligand.objects.get(clean_inchikey=li, parent__isnull=True)
                except Ligand.DoesNotExist:
                    # print(f'No right parent found: {l.gpcrdb_id}')
                    continue
                if li!=lp:
                    # print('Ligand ids')
                    # for lid in l.ids.all():
                    #     print(lid)
                    # print('Ligand parent ids')
                    # for plid in supposed_to_be_parent.ids.all():
                    #     print(plid)
                    print(l.gpcrdb_id, l.name, l.parent.name, supposed_to_be_parent.name, supposed_to_be_parent.gpcrdb_id)
                    c+=1

            isotopes = get_radioactive_info(l.name, [])
            if len(isotopes)>0:
                radioactive_dict[l] = isotopes

        #     iexact_dict[l] = Ligand.objects.filter(name__iexact=l.name, parent__isnull=False)

        # print("===================================================")
        # print("Checking iexact name matches with different parents")
        # for l, li in iexact_dict.items():
        #     if len(li)>1 and len(li.order_by("parent").distinct("parent"))>1:
        #         print('===========')
        #         for lig in li:
        #             print(lig, lig.gpcrdb_id, lig.parent, lig.parent.gpcrdb_id)

        from collections import defaultdict

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
            if (
                len(sibling_list) == 2 and
                (
                    (sibling_list[0].radioactive is not None and sibling_list[1].radioactive is None)
                    or
                    (sibling_list[1].radioactive is not None and sibling_list[0].radioactive is None)
                )
            ):
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
            if (
                len(sibling_list) == 2 and
                (
                    (sibling_list[0].radioactive is not None and sibling_list[1].radioactive is None)
                    or
                    (sibling_list[1].radioactive is not None and sibling_list[0].radioactive is None)
                )
            ):
                continue
            dup_found = True
            print(f"  Duplicate clean_inchikey={row['clean_inchikey']} under parent_id={row['parent_id']}:")
            for s in siblings:
                pubchem_ids = [lid.index for lid in s.ids.all() if lid.web_resource.slug == "pubchem"]
                print(f"    gpcrdb_id={s.gpcrdb_id}  pubchem={pubchem_ids}")
        if not dup_found:
            print("  None found")
