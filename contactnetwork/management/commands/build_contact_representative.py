from django.core.management.base import BaseCommand, CommandError
from django.core.management import call_command
from django.conf import settings
from django.db import connection
from django.db.models import Q, F
from contactnetwork.models import *

import time


class Command(BaseCommand):

    help = "Build  contact representatives"


    def handle(self, *args, **options):

        print('Script to decide contact representative for a conformation. Maximising highest frequence of common contacts, while minimizing uncommon (50%)')

        structures = Structure.objects.filter(refined=False).prefetch_related(
            "pdb_code",
            "state",
            "protein_conformation__protein__parent__family")

        distinct_proteins = {}

        for s in structures:
            pdb = s.pdb_code.index
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
                continue
            interactions = list(Interaction.objects.filter(
                interacting_pair__referenced_structure__pdb_code__index__in=pdbs
            ).values(
                'interaction_type',
                'interacting_pair__referenced_structure__pdb_code__index',
                'interacting_pair__res1__generic_number__label',
                'interacting_pair__res2__generic_number__label',
            ).filter(interacting_pair__res1__pk__lt=F('interacting_pair__res2__pk')).distinct())

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
            uncommon_contacts = set()
            for pair, count in contacts_count.items():
                if (count>=number_of_pdbs/2):
                    common_contacts.add(pair)
                else:
                    uncommon_contacts.add(pair)

            print("common contacts count",len(common_contacts),"uncommon contacts count",len(uncommon_contacts))

            pdbs_scored = []
            for pdb in pdbs:
                pdb_common_contacts = common_contacts.intersection(contacts_in_pdb[pdb])
                pdb_uncommon_contacts = uncommon_contacts.intersection(contacts_in_pdb[pdb])
                print(pdb,len(pdb_common_contacts),len(pdb_uncommon_contacts),len(contacts_in_pdb[pdb]))
                pdbs_scored.append([pdb,len(pdb_common_contacts),len(pdb_uncommon_contacts),len(contacts_in_pdb[pdb])])

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