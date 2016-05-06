from django.core.management.base import BaseCommand, CommandError
from django.conf import settings

from residue.models import ResiduePositionSet, ResidueGenericNumberEquivalent
from interaction.models import ResidueFragmentInteraction

import logging

class Command(BaseCommand):
    help = 'Reads source data and creates common database tables'

    logger = logging.getLogger(__name__)

    def handle(self, *args, **options):
        functions = [
            'create_residue_sets',
        ]

        # execute functions
        for f in functions:
            try:
                getattr(self, f)()
            except Exception as msg:
                print(msg)
                self.logger.error(msg)

    def create_residue_sets(self):
        self.logger.info('CREATING RESIDUE SETS')

        residue_sets = {
            'Signalling protein pocket': ['gpcrdba', ['3x50', '3x53', '3x54', '3x55', '34x50', '34x51', '34x53', '34x54', '5x64', '5x67', '5x68', '6x29', '6x36', '7x55', '8x48', '8x49']],
                        }
        for set_name in residue_sets.keys():
            residues = []
            for res in residue_sets[set_name][1]:
                try:
                    residues.append(ResidueGenericNumberEquivalent.objects.get(label=res, scheme__slug=residue_sets[set_name][0]))
                except Exception as e: 
                    print(e)
            if residues:
                try:
                    rs = ResiduePositionSet.objects.create(name=set_name)
                    for res in residues:
                        rs.residue_position.add(res)
                except Exception as msg:
                    self.logger.debug('Failed to create residue set {}. Error: {}'.format(set_name, msg))
            self.logger.info('SET {} CREATED'.format(set_name))

        gpcr_class = {
            'gpcrdba':'001',
            #'gpcrdbb': '002',
            'gpcrdbc': '004',
            'gpcrdbf': '005',
            }
        set_names = {
            'gpcrdba': "Class A binding pocket",
#            'gpcrdbb': "Class B binding pocket",
            'gpcrdbc': "Class C binding pocket",
            'gpcrdbf': "Class F binding pocket",
            }
        for c in gpcr_class:
            class_interactions = ResidueFragmentInteraction.objects.filter(
                structure_ligand_pair__structure__protein_conformation__protein__family__slug__startswith=gpcr_class[c], structure_ligand_pair__annotated=True).prefetch_related(
                'rotamer__residue__display_generic_number','interaction_type',
                'structure_ligand_pair__structure__protein_conformation__protein__parent',
                'structure_ligand_pair__ligand__properities')

            generic = {}

            for i in class_interactions:
                if i.rotamer.residue.generic_number:
                    gn = i.rotamer.residue.generic_number.label
                else:
                    continue
                if gn not in generic.keys():
                    generic[gn] = 1
                else:
                    generic[gn] += 1
            try:
                rs = ResiduePositionSet.objects.create(name=set_names[c])
            except Exception as msg:
                print(msg)
                self.logger.debug('Failed to create residue set {}. Error: {}'.format(set_names[c], msg))
                continue
            for g in sorted(generic):
                if generic[g] >= 2:
                    try:
                        rs.residue_position.add(ResidueGenericNumberEquivalent.objects.get(label=g, scheme__slug=c))
                    except Exception as msg:
                        print(g)
                        print(msg)
                        continue
            self.logger.info('SET {} CREATED'.format(set_names[c]))
        self.logger.info('COMPLETED CREATING RESIDUE SETS')
