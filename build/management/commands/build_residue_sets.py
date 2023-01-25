from django.core.management.base import BaseCommand, CommandError
from django.conf import settings

from residue.models import ResiduePositionSet, ResidueGenericNumberEquivalent
from interaction.models import ResidueFragmentInteraction
from common.tools import test_model_updates
import logging
import django.apps

class Command(BaseCommand):
    help = 'Reads source data and creates common database tables'

    logger = logging.getLogger(__name__)
    tracker = {}
    all_models = django.apps.apps.get_models()[6:]
    test_model_updates(all_models, tracker, initialize=True)

    def purge_residue_sets_(self):
        try:
            ResiduePositionSet.objects.all().delete()
        except Exception as msg:
            print(msg)
            self.logger.warning('Existing data cannot be deleted')

    def handle(self, *args, **options):
        functions = [
            'create_residue_sets',
        ]

        try:
            self.purge_residue_sets_()
            test_model_updates(self.all_models, self.tracker, initialize=True)
        except Exception as msg:
            print(msg)
            self.logger.error(msg)

        # execute functions
        for f in functions:
            try:
                getattr(self, f)()
            except Exception as msg:
                print(msg)
                self.logger.error(msg)
            test_model_updates(self.all_models, self.tracker, check=True)

    def create_residue_sets(self):
        self.logger.info('CREATING RESIDUE SETS')

        residue_sets = {
            'G-protein interface': [
                'gpcrdba',
                ['3x50', '3x53', '3x54', '3x55', '3x56', '34x50', '34x51', '34x52', '34x53', '34x54', '34x55', '34x56', '34x57', '5x61', '5x64', '5x65', '5x66',
                 '5x67', '5x68', '5x69', '5x71', '5x72', '5x74', '5x75', '6x25', '6x26', '6x28', '6x29', '6x32', '6x33', '6x36', '6x37', '6x40', '7x55', '7x56', '8x47', '8x48', '8x49', '8x51'],
                 'Signal protein interfaces',
                 'gpcr'
            ],
            'State (micro-)switches': [
                'gpcrdba',
                ['3x49', '7x43', '3x50', '5x47', '5x50', '5x58', '6x30', '6x34', '6x48', '6x50', '7x49', '7x50', '7x53', '3x40', '6x44'],
                 'Conformation/state-stabilising',
                 'gpcr'
            ],
            'Sodium ion pocket': [
                'gpcrdba',
                ['1x50', '1x53', '2x46', '2x47', '2x49', '2x50', '3x39', '3x43', '6x44', '6x48', '7x45', '7x46', '7x49', '7x50', '7x53'],
                 'Conformation/state-stabilising',
                 'gpcr'
            ],
            'Gprotein Barcode': [
                'cgn',
                ['G.hns1.02','G.hns1.03','G.S1.02','G.s2s3.01','G.S3.01','G.S3.03','G.H4.16','G.H4.17','G.h4s6.03','G.h4s6.20','G.H5.08','G.H5.11','G.H5.12',
                 'G.H5.13','G.H5.15','G.H5.16','G.H5.17','G.H5.19','G.H5.20','G.H5.21','G.H5.22','G.H5.23','G.H5.24','G.H5.25','G.H5.26'],
                 'Barcode',
                 'gprotein'
            ],
            'YM binding site': [
                'cgn',
                ['G.H1.02', 'G.H1.05', 'G.H1.06', 'G.H1.09', 'G.h1ha.01', 'G.h1ha.04', 'H.HA.03', 'H.HA.06', 'H.HA.07', 'H.HA.10', 'G.hfs2.03', 'G.hfs2.05', 'G.hfs2.06', 'G.S2.01', 'G.S2.02', 'G.S2.03', 'G.S2.04'],
                 'Binding site',
                 'gprotein'
            ],
            'Arrestin interface': [
                'can',
                ['N.s5s6.05', 'N.s5s6.09', 'C.s15s16.03', 'C.S16.01', 'C.s15s16.02', 'C.S15.14', 'C.S16.03', 'N.S5.12', 'N.s8s9.09', 'N.S9.01', 'N.s8s9.02',
                'N.s8s9.10', 'N.S8.03', 'N.s5s6.10', 'N.S6.03', 'N.s9s10.08',
                'N.s5s6.11', 'N.s5s6.08', 'N.s5s6.07', 'N.s5s6.04', 'C.s17s18.39',
                'C.s17s18.41', 'N.S6.04', 'N.S6.05', 'C.s17s18.38'],
                'Interface',
                'arrestin'
            ]
            }

        for set_name in residue_sets.keys():
            residues = []
            for res in residue_sets[set_name][1]:
                try:
                    residues.append(ResidueGenericNumberEquivalent.objects.get(label=res, scheme__slug=residue_sets[set_name][0]))
                except Exception as e:
                    print('ERROR: ', res)
                    print(e)
            if residues:
                try:
                    rs = ResiduePositionSet.objects.create(name=set_name,set_type=residue_sets[set_name][2], protein_group=residue_sets[set_name][3])
                    for res in residues:
                        rs.residue_position.add(res)
                except Exception as msg:
                    self.logger.debug('Failed to create residue set {}. Error: {}'.format(set_name, msg))
            self.logger.info('SET {} CREATED'.format(set_name))

        gpcr_class = {
            'gpcrdba':'001',
            'gpcrdbb': '002',
            'gpcrdbb2': '003',
            'gpcrdbc': '004',
            'gpcrdbf': '006',
            }
        set_names = {
            'gpcrdba': "Class A",
            'gpcrdbb': "Class B1",
            'gpcrdbb2': "Class B2",
            'gpcrdbc': "Class C",
            'gpcrdbf': "Class F",
            }
        for c in gpcr_class:
            class_interactions = ResidueFragmentInteraction.objects.filter(
                structure_ligand_pair__structure__protein_conformation__protein__family__slug__startswith=gpcr_class[c], structure_ligand_pair__annotated=True).prefetch_related(
                'rotamer__residue__display_generic_number','interaction_type',
                'structure_ligand_pair__structure__protein_conformation__protein__parent')

            generic = {}

            for i in class_interactions:
                if i.rotamer.residue.generic_number:
                    gn = i.rotamer.residue.display_generic_number.label
                else:
                    continue
                if gn not in generic.keys():
                    generic[gn] = 1
                else:
                    generic[gn] += 1
            try:
                rs = ResiduePositionSet.objects.create(name=set_names[c],set_type='Ligand binding pocket', protein_group='gpcr')
            except Exception as msg:
                print(msg)
                self.logger.debug('Failed to create residue set {}. Error: {}'.format(set_names[c], msg))
                continue
            for g in sorted(generic):
                if generic[g] >= 2:
                    bw, gpcrdb = g.split('x')
                    h, pos = bw.split('.')
                    try:
                        scheme = c
                        if scheme == "gpcrdbb2":
                            scheme = "gpcrdbb"
                        rs.residue_position.add(ResidueGenericNumberEquivalent.objects.get(label='{}x{}'.format(h,gpcrdb), scheme__slug=scheme))
                    except Exception as msg:
                        print(g)
                        print(msg)
                        continue
            self.logger.info('SET {} CREATED'.format(set_names[c]))
        self.logger.info('COMPLETED CREATING RESIDUE SETS')
