from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.utils.text import slugify
from django.utils.html import strip_tags
from django.db import IntegrityError
from common.tools import test_model_updates
from build.management.commands.base_build import Command as BaseBuild
from protein.models import (Protein, ProteinConformation, ProteinState, ProteinSequenceType, ProteinSegment, ProteinSource)
from residue.models import Residue
from structure.functions import ParseStructureCSV

import os
import logging
import yaml
import django.apps


class Command(BaseBuild):
    help = 'Reads source data and creates protein records for constructs'
    #Setting the variables for the test tracking of the model upadates
    tracker = {}
    all_models = django.apps.apps.get_models()[6:]
    test_model_updates(all_models, tracker, initialize=True)
    def add_arguments(self, parser):
        parser.add_argument('-p', '--proc',
            type=int,
            action='store',
            dest='proc',
            default=1,
            help='Number of processes to run')
        parser.add_argument('-c', '--construct',
            # action='append',
            dest='construct',
            help='Constructs to import (PDB ID)',
            nargs='+')
        parser.add_argument('-u', '--purge',
            action='store_true',
            dest='purge',
            default=False,
            help='Purge existing construct records')

    def handle(self, *args, **options):
        # delete any existing construct data
        if options['purge']:
            try:
                self.purge_constructs()
                self.tracker = {}
                test_model_updates(self.all_models, self.tracker, initialize=True)
            except Exception as msg:
                print(msg)
                self.logger.error(msg)
        self.parsed_structures = ParseStructureCSV()

        # where pdb ids specified?
        if options['construct']:
            self.parsed_structures.pdb_ids = [i for i in self.parsed_structures.pdb_ids if i in options['construct'] or i.lower() in options['construct']]
        try:
            self.logger.info('CREATING CONSTRUCTS')
            self.prepare_input(options['proc'], self.parsed_structures.pdb_ids)
            test_model_updates(self.all_models, self.tracker, check=True)
            self.logger.info('COMPLETED CREATING CONSTRUCTS')
        except Exception as msg:
            print(msg)
            self.logger.error(msg)

    def purge_constructs(self):
        try:
            pst = ProteinSequenceType.objects.get(slug='mod')
            Protein.objects.filter(sequence_type=pst).delete()
        except ProteinSequenceType.DoesNotExist:
            self.logger.warning('ProteinSequenceType mod not found: nothing to delete.')

    def main_func(self, positions, iteration,count,lock):
        # setting up processes
        if not positions[1]:
            constructs = self.parsed_structures.pdb_ids[positions[0]:]
        else:
            constructs = self.parsed_structures.pdb_ids[positions[0]:positions[1]]

        for sd in constructs:
            sd = self.parsed_structures.structures[sd]
            # is a protein specified?
            if 'protein' not in sd:
                self.logger.error('Protein not specified for construct, skipping')
                continue

            # fetch the parent protein
            try:
                ppc = ProteinConformation.objects.prefetch_related('protein__family', 'protein__species',
                    'protein__residue_numbering_scheme').get(protein__entry_name=sd['protein'].lower(),
                    state__slug=settings.DEFAULT_PROTEIN_STATE)
            except ProteinConformation.DoesNotExist:
                # abort if parent protein is not found
                print('Parent protein {} for construct {} not found, aborting!'.format(
                    sd['protein'], sd['name']))
                self.logger.error('Parent protein {} for construct {} not found, aborting!'.format(
                    sd['protein'], sd['name']))
                continue
            # sequence type
            try:
                sequence_type, created = ProteinSequenceType.objects.get_or_create(slug='mod',
                    defaults={'name': 'Modified'})
                if created:
                    self.logger.info('Created sequence type {}'.format(sequence_type))
            except IntegrityError:
                sequence_type = ProteinSequenceType.objects.get(slug='mod')

            # protein source
            try:
                protein_source, created = ProteinSource.objects.get_or_create(name='OTHER')
                if created:
                    self.logger.info('Created protein source {}'.format(protein_source))
            except IntegrityError:
                protein_source = ProteinSource.objects.get(name='OTHER')


            if not Protein.objects.filter(name=sd['name']).exists():
                # create a protein record
                p = Protein()
                p.parent = ppc.protein
                p.family = ppc.protein.family
                p.species = ppc.protein.species
                p.residue_numbering_scheme = ppc.protein.residue_numbering_scheme
                p.sequence_type= sequence_type
                p.source = protein_source
                p.entry_name = slugify(strip_tags(sd['name']))
                p.name = sd['name']
                p.sequence = ppc.protein.sequence

                # save protein (construct)
                try:
                    p.save()
                    self.logger.info('Created construct {} with parent protein {}'.format(p.name,
                        ppc.protein.entry_name))
                except:
                    self.logger.error('Failed creating construct {} with parent protein {}'.format(p.name,
                        ppc.protein.entry_name))
                    continue
            else:
                p = Protein.objects.get(name=sd['name'])


            if not ProteinConformation.objects.filter(protein=p).exists():
                # create protein conformation record
                pc = ProteinConformation()
                pc.protein = p
                pc.state = ProteinState.objects.get(slug=settings.DEFAULT_PROTEIN_STATE)
                try:
                    pc.save()
                    self.logger.info('Created conformation {} of protein {}'.format(pc.state.name, p.name))
                except:
                    print('Failed creating conformation {} of protein {}'.format(pc.state.name,p.entry_name))
                    self.logger.error('Failed creating conformation {} of protein {}'.format(pc.state.name,
                        p.entry_name))
