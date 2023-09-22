from django.core.management.base import BaseCommand, CommandError
from django.core.management import call_command
from django.conf import settings
from django.db import connection
from django.db.models import Q, F
from contactnetwork.models import *
from protein.models import ProteinFamily

import time


class Command(BaseCommand):

    help = "Build  contact representatives"


    def handle(self, *args, **options):

        self.receptor_representatives()
        self.class_level_contacts()
        self.class_based_representative()

    def receptor_representatives(self):
        print('Script to decide contact representative for a conformation. Maximising highest frequence of common contacts, while minimizing uncommon (50%)')

        structures = Structure.objects.all().exclude(structure_type__slug__startswith='af-').prefetch_related(
            "pdb_code",
            "state",
            "protein_conformation__protein__parent__family")

        distinct_proteins = {}

        resolution_lookup = {}
        for s in structures:
            pdb = s.pdb_code.index
            resolution_lookup[pdb] = s.resolution
            state = s.state.slug
            slug = s.protein_conformation.protein.parent.family.slug
            name = s.protein_conformation.protein.parent.family.name

            key = '{}_{}'.format(name,state)

            if key not in distinct_proteins:
                distinct_proteins[key] = []

            distinct_proteins[key].append(pdb)


        conformation_representative = []
        conformation_non_representative = []
        for conformation, pdbs in distinct_proteins.items():
            print(conformation, "PDBS:",pdbs)
            number_of_pdbs = len(pdbs)
            if (number_of_pdbs==1):
                # Do not care when only one PDB for a conformation rep
                print("REPRESENTATIVE:",pdbs[0])
                conformation_representative.append(pdbs[0])
                s = Structure.objects.get(pdb_code__index=pdbs[0])
                s.contact_representative_score = 1
                s.save()
                continue
            interactions = list(Interaction.objects.filter(
                interacting_pair__referenced_structure__pdb_code__index__in=pdbs
            ).values(
                'interaction_type',
                'interacting_pair__referenced_structure__pdb_code__index',
                'interacting_pair__res1__generic_number__label',
                'interacting_pair__res2__generic_number__label',
            ).exclude(interacting_pair__res1__generic_number__isnull=True).filter(interacting_pair__res1__pk__lt=F('interacting_pair__res2__pk')).distinct())

            contacts_in_pdb = {}
            contacts_count = {}
            for i in interactions:
                pair = "{}_{}".format(i['interacting_pair__res1__generic_number__label'],i['interacting_pair__res2__generic_number__label'])
                pdb = i['interacting_pair__referenced_structure__pdb_code__index']
                #print(pair,pdb)

                if pdb not in contacts_in_pdb:
                    contacts_in_pdb[pdb] = set()
                if pair not in contacts_count:
                    contacts_count[pair] = 0

                # Make sure not to count interactions several times when there are multiple interactions for a pair (WdV & Hydrophobic)
                if pair not in contacts_in_pdb[pdb]:
                    contacts_count[pair] += 1
                    contacts_in_pdb[pdb].add(pair)

            # print(contacts_in_pdb)
            # print(contacts_count)

            common_contacts = set()
            equalcommon_contacts = set()
            uncommon_contacts = set()
            for pair, count in contacts_count.items():
                if (count>number_of_pdbs/2):
                    common_contacts.add(pair)
                elif (count==number_of_pdbs/2):
                    equalcommon_contacts.add(pair)
                else:
                    uncommon_contacts.add(pair)

            print("common contacts count",len(common_contacts),"50% contacts count",len(equalcommon_contacts),"uncommon contacts count",len(uncommon_contacts))

            pdbs_scored = []
            for pdb in pdbs:
                pdb_common_contacts = common_contacts.intersection(contacts_in_pdb[pdb])
                pdb_uncommon_contacts = uncommon_contacts.intersection(contacts_in_pdb[pdb])
                score = len(pdb_common_contacts) + len(uncommon_contacts) - len(pdb_uncommon_contacts)
                if len(common_contacts)+len(uncommon_contacts) > 0:
                    score /= (len(common_contacts)+len(uncommon_contacts))
                score = 0
                print(pdb,score,resolution_lookup[pdb],len(pdb_common_contacts),len(pdb_uncommon_contacts),len(contacts_in_pdb[pdb]))
                pdbs_scored.append([pdb,score,resolution_lookup[pdb],len(pdb_common_contacts),len(pdb_uncommon_contacts),len(contacts_in_pdb[pdb])])
                s = Structure.objects.get(pdb_code__index=pdb)
                s.contact_representative_score = score
                s.save()

            pdbs_scored = sorted(pdbs_scored, key=lambda x: (-x[1],x[2]))
            print("REPRESENTATIVE:",pdbs_scored[0])
            conformation_representative.append(pdbs_scored[0][0])
            conformation_non_representative += [row[0] for row in pdbs_scored[1:]]
            #print(interactions)

            # break
        print('CONFORMATION REPRESENTATIVE')
        print(conformation_representative)
        for pdb in conformation_representative:
            s = Structure.objects.get(pdb_code__index=pdb)
            s.contact_representative = True
            s.save()

        print('CONFORMATION NON REPRESENTATIVE')
        print(conformation_non_representative)
        for pdb in conformation_non_representative:
            s = Structure.objects.get(pdb_code__index=pdb)
            s.contact_representative = False
            s.save()


    def class_level_contacts(self):

        class_level_contacts = {}
        structures = Structure.objects.filter(refined=False).exclude(structure_type__slug__startswith='af-').prefetch_related(
            "pdb_code",
            "state",
            "protein_conformation__protein__parent__family",
            "protein_conformation__protein__family",
            "protein_conformation__protein__species")

        distinct_proteins = {}

        pdb_to_protein_lookup = {}
        pdb_to_class_lookup = {}
        pdb_to_state_lookup = {}
        for s in structures:
            pdb = s.pdb_code.index
            gpcr_class = s.protein_conformation.protein.family.slug.split("_")[0]
            state = s.state
            slug = s.protein_conformation.protein.parent.family.slug
            name = s.protein_conformation.protein.parent.family.name
            species = s.protein_conformation.protein.species.common_name
            protein = s.protein_conformation.protein.parent

            pdb_to_protein_lookup[pdb] = protein
            pdb_to_class_lookup[pdb] = gpcr_class
            pdb_to_state_lookup[pdb] = state.slug

            if species != "Human":
                # Only use humans for this
                continue

            key = '{}_{}'.format(name,state)

            if key not in distinct_proteins:
                distinct_proteins[key] = []

            distinct_proteins[key].append(pdb)


            if gpcr_class not in class_level_contacts:
                class_level_contacts[gpcr_class] = {}
            if state not in class_level_contacts[gpcr_class]:
                class_level_contacts[gpcr_class][state] = {"pdbs":set(), "proteins":set()}

            class_level_contacts[gpcr_class][state]["pdbs"].add(pdb)
            class_level_contacts[gpcr_class][state]["proteins"].add(protein)


        #Reset data
        ConsensusInteraction.objects.all().delete()
        pdb_contacts = {}
        for gpcr_class in class_level_contacts:
            protein_class = ProteinFamily.objects.get(slug=gpcr_class)
            for state in class_level_contacts[gpcr_class]:
                pdbs = class_level_contacts[gpcr_class][state]["pdbs"]
                proteins_number = len(class_level_contacts[gpcr_class][state]["proteins"])
                print(gpcr_class, state, "PDBS:",pdbs)

                interactions = list(Interaction.objects.filter(
                    interacting_pair__referenced_structure__pdb_code__index__in=pdbs
                ).prefetch_related(
                    'interacting_pair__referenced_structure__pdb_code',
                    'interacting_pair__res1__generic_number',
                    'interacting_pair__res2__generic_number',
                ).exclude(interacting_pair__res1__generic_number__isnull=True
                ).filter(interacting_pair__res1__pk__lt=F('interacting_pair__res2__pk')).distinct('interacting_pair__referenced_structure__pdb_code__index',
                                                                                                  'interacting_pair__res1__generic_number__label',
                                                                                                  'interacting_pair__res2__generic_number__label'))

                consensus = {}
                for i in interactions:
                    gn1 = i.interacting_pair.res1.generic_number
                    gn2 = i.interacting_pair.res2.generic_number
                    pair = "{}_{}".format(gn1.label,gn2.label)
                    s = i.interacting_pair.referenced_structure
                    protein = pdb_to_protein_lookup[s.pdb_code.index]
                    if pair not in consensus:
                        consensus[pair] = {"pdbs": set(), "proteins": set(), "gn1":gn1, "gn2":gn2}
                    if s.pdb_code.index not in pdb_contacts:
                        pdb_contacts[s.pdb_code.index] = set()

                    pdb_contacts[s.pdb_code.index].add(pair)

                    consensus[pair]["pdbs"].add(s)
                    consensus[pair]["proteins"].add(protein)
                print(len(consensus), "consensus pairs to calculate")
                for pair in consensus:
                    # print(pair, len(consensus[pair]["pdbs"]),len(consensus[pair]["proteins"]),len(consensus[pair]["proteins"])/ proteins_number)
                    ci = ConsensusInteraction()

                    ci.gn1 = consensus[pair]["gn1"]
                    ci.gn2 = consensus[pair]["gn2"]
                    ci.protein_class = protein_class
                    ci.state = state
                    ci.frequency = len(consensus[pair]["proteins"])/ proteins_number
                    ci.save()

                    ci.structures.add(*consensus[pair]["pdbs"])
                    ci.proteins.add(*consensus[pair]["proteins"])
                #print(consensus)


        # Find state-specific contacts
        frequent_contacts = {}

        states = ['active','inactive']
        contacts = ConsensusInteraction.objects.filter(frequency__gt=0.4, state__slug__in = states).prefetch_related('gn1','gn2','state','protein_class').all()
        state_pairs = {}
        for c in contacts:
            class_slug = c.protein_class.slug
            pair = "{}_{}_{}".format(class_slug,c.gn1.label,c.gn2.label)
            state = c.state.slug
            if class_slug not in state_pairs:
                state_pairs[class_slug] = {'active':set(), 'inactive':set()}
            freq = c.frequency
            if pair not in frequent_contacts:
                frequent_contacts[pair] = {'active':[0,0], 'inactive':[0,0], }
            frequent_contacts[pair][state][0] = freq
            frequent_contacts[pair][state][1] = c

        for label, fc in frequent_contacts.items():
            class_slug = label.split("_")[0]
            pair = "_".join(label.split("_")[1:])
            if fc['active'][0]-fc['inactive'][0] > 0.4:
                fc['active'][1].state_specific = True
                fc['active'][1].save()
                state_pairs[class_slug]['active'].add(pair)
            elif fc['inactive'][0]-fc['active'][0] > 0.4:
                fc['inactive'][1].state_specific = True
                fc['inactive'][1].save()
                state_pairs[class_slug]['inactive'].add(pair)

        #print(state_pairs)
        for pdb in pdb_to_protein_lookup:
            if pdb in pdb_contacts:
                pairs = pdb_contacts[pdb]
            else:
                # Need to fetch pairs

                interactions = list(Interaction.objects.filter(
                    interacting_pair__referenced_structure__pdb_code__index=pdb
                ).prefetch_related(
                    'interacting_pair__referenced_structure__pdb_code',
                    'interacting_pair__res1__generic_number',
                    'interacting_pair__res2__generic_number',
                ).exclude(interacting_pair__res1__generic_number__isnull=True
                ).filter(interacting_pair__res1__pk__lt=F('interacting_pair__res2__pk')).distinct('interacting_pair__referenced_structure__pdb_code__index',
                                                                                                  'interacting_pair__res1__generic_number__label',
                                                                                                  'interacting_pair__res2__generic_number__label'))
                for i in interactions:
                    gn1 = i.interacting_pair.res1.generic_number
                    gn2 = i.interacting_pair.res2.generic_number
                    pair = "{}_{}".format(gn1.label,gn2.label)
                    if pdb not in pdb_contacts:
                        pdb_contacts[pdb] = set()
                    pdb_contacts[pdb].add(pair)
                pairs = pdb_contacts[pdb]

            gpcr_class = pdb_to_class_lookup[pdb]
            pdb_state = pdb_to_state_lookup[pdb]
            if gpcr_class in state_pairs:
                class_pairs = state_pairs[gpcr_class]
                active_intersection = class_pairs['active'].intersection(pairs)
                inactive_intersection = class_pairs['inactive'].intersection(pairs)
                if len(class_pairs['active']):
                    active_intersection_fraction = len(active_intersection)/len(class_pairs['active'])
                else:
                    active_intersection_fraction = 0
                if len(class_pairs['inactive']):
                    inactive_intersection_fraction = len(inactive_intersection)/len(class_pairs['inactive'])
                else:
                    inactive_intersection_fraction = 0
                print(pdb,pdb_state, gpcr_class,'active',len(active_intersection),active_intersection_fraction,'inactive',len(inactive_intersection),inactive_intersection_fraction)
                s = Structure.objects.get(pdb_code__index=pdb)
                s.active_class_contacts_fraction = active_intersection_fraction
                s.inactive_class_contacts_fraction = inactive_intersection_fraction
                s.save()
            else:
                s = Structure.objects.get(pdb_code__index=pdb)
                if s.state.slug == "active":
                    s.active_class_contacts_fraction = 1
                    s.inactive_class_contacts_fraction = 0
                    s.save()
                elif s.state.slug == "inactive":
                    s.active_class_contacts_fraction = 0
                    s.inactive_class_contacts_fraction = 1
                    s.save()
                print("Class has no pairs, setting", pdb," to 100%")

    def class_based_representative(self):
        print('Script to decide contact representative for a conformation. Maximising highest frequence of common contacts, while minimizing uncommon (50%)')

        structures = Structure.objects.filter(refined=False).exclude(structure_type__slug__startswith='af-').prefetch_related(
            "pdb_code",
            "state",
            "protein_conformation__protein__parent__family")

        distinct_proteins = {}

        resolution_lookup = {}
        for s in structures:
            pdb = s.pdb_code.index
            resolution_lookup[pdb] = s.resolution
            state = s.state.slug
            slug = s.protein_conformation.protein.parent.family.slug
            name = s.protein_conformation.protein.parent.family.name

            key = '{}_{}'.format(name,state)

            if key not in distinct_proteins:
                distinct_proteins[key] = []

            distinct_proteins[key].append(s)


        conformation_representative = []
        conformation_non_representative = []
        for conformation, pdbs in distinct_proteins.items():
            state = conformation.split("_")[-1]
            print(state,conformation)
            number_of_pdbs = len(pdbs)
            if (number_of_pdbs==1):
                # Do not care when only one PDB for a conformation rep
                print("REPRESENTATIVE:",pdbs[0])
                conformation_representative.append(pdbs[0])
                continue

            scores = []
            for s in pdbs:
                fraction_active = s.active_class_contacts_fraction
                fraction_inactive = s.inactive_class_contacts_fraction
                diff_in_fraction = fraction_inactive-fraction_active
                abs_diff = abs(diff_in_fraction)

                scores.append([s,fraction_inactive,fraction_active,diff_in_fraction,abs_diff])

            if state=='active':
                scores_sorted = sorted(scores, key=lambda x: (x[3]))
            elif state=='inactive':
                scores_sorted = sorted(scores, key=lambda x: (-x[3]))
            elif state=='intermediate':
                scores_sorted = sorted(scores, key=lambda x: (x[4]))

            print("REPRESENTATIVE:",scores_sorted[0])
            conformation_representative.append(scores_sorted[0][0])
            conformation_non_representative += [row[0] for row in scores_sorted[1:]]
            #print(interactions)

            # break
        print('CONFORMATION REPRESENTATIVE')
        print(conformation_representative)
        for s in conformation_representative:
            s.class_contact_representative = True
            s.save()

        print('CONFORMATION NON REPRESENTATIVE')
        print(conformation_non_representative)
        for s in conformation_non_representative:
            s.class_contact_representative = False
            s.save()
