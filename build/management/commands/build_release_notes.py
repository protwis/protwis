from django.core.management.base import BaseCommand, CommandError
from django.conf import settings

from common.models import ReleaseNotes, ReleaseStatistics, ReleaseStatisticsType
from protein.models import Protein, Species
from ligand.models import Ligand
from structure.models import Structure, StructureModel, StructureComplexModel
from mutation.models import MutationExperiment
from interaction.models import StructureLigandInteraction, ProteinLigandInteraction, ResidueFragmentInteraction

import logging
import shlex
import os
import yaml

class Command(BaseCommand):
    help = 'Builds release notes and stats'

    logger = logging.getLogger(__name__)

    release_notes_data_dir = os.sep.join([settings.DATA_DIR, 'release_notes'])

    def handle(self, *args, **options):
        try:
            self.create_release_notes()
        except Exception as msg:
            print(msg)
            self.logger.error(msg)

    def create_release_notes(self):
        self.logger.info('CREATING RELEASE NOTES')

        ReleaseNotes.objects.all().delete()
        
        # what files should be parsed?
        filenames = os.listdir(self.release_notes_data_dir)

        for source_file in filenames:
            self.logger.info('Parsing file ' + source_file)
            source_file_path = os.sep.join([self.release_notes_data_dir, source_file])

            if source_file.endswith('.yaml'):
                with open(source_file_path, 'r') as f:
                    ds = yaml.load(f, Loader=yaml.FullLoader)

                    release_notes, created = ReleaseNotes.objects.get_or_create(date=ds['date'])
                    if created:
                        self.logger.info('Created release notes for {}'.format(ds['date']))
                
                with open(source_file_path[:-4]+'html', 'r') as h:
                    release_notes.html = h.read()
                    release_notes.save()

                    if created:
                        self.logger.info('Updated html for release notes {}'.format(ds['date']))

        # generate statistics
        latest_release_notes = ReleaseNotes.objects.all()[0]

        stats = [
            ['Proteins', Protein.objects.filter(sequence_type__slug='wt').count()],
            ['Human proteins', Protein.objects.filter(sequence_type__slug='wt', species__common_name="Human").count()],
            ['Species', Species.objects.all().count()],
            ['Exp. structures', Structure.objects.filter(refined=False).count()],
            ['Refined structures', Structure.objects.filter(refined=True).count()],
            ['GPCR structure models', StructureModel.objects.all().count()],
            ['GPCR-G protein structure models', StructureComplexModel.objects.all().count()],
            ['Ligands', Ligand.objects.filter(canonical=True).count()],
            ['Mutants', MutationExperiment.objects.all().count()],
            ['Ligand interactions', ResidueFragmentInteraction.objects.all().count()],
        ]

        for stat in stats:
            stat_type, created = ReleaseStatisticsType.objects.get_or_create(name=stat[0])
            if created:
                self.logger.info('Created release statistics type {}'.format(stat[0]))

            rstat, created = ReleaseStatistics.objects.update_or_create(release=latest_release_notes, statistics_type=stat_type,
                defaults={'value': stat[1]})
            if created:
                self.logger.info('Created release stat {} {}'.format(latest_release_notes, stat[1]))

        self.logger.info('COMPLETED CREATING RELEASE NOTES')
