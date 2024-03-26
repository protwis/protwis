from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import connection

from build.management.commands.base_build import Command as BaseBuild
from protein.models import (Protein, ProteinConformation, ProteinSequenceType, ProteinSegment,
    ProteinConformationTemplateStructure)
from structure.models import Structure
from common.alignment import Alignment

from cProfile import Profile

class Command(BaseBuild):
    help = 'Goes though all protein records and finds the best structural templates for the whole 7TM domain, and ' \
        + 'each helix individually'

    def add_arguments(self, parser):
        parser.add_argument('-p', '--proc',
            type=int,
            action='store',
            dest='proc',
            default=1,
            help='Number of processes to run')
        parser.add_argument('-r', '--profile',
            action='store_true',
            dest='profile',
            default=False,
            help='Profile the script with cProfile')

    # segments
    segments = ProteinSegment.objects.filter(slug__in=settings.REFERENCE_POSITIONS)

    # fetch representative (inactive) structures FIXME add active structure??
    structures = Structure.objects.filter(representative=True,
        protein_conformation__state__slug=settings.DEFAULT_PROTEIN_STATE).exclude(structure_type__slug__startswith='af-').prefetch_related(
        'protein_conformation__protein__parent__family')

    # fetch all protein conformations
    pconfs = ProteinConformation.objects.all().prefetch_related('protein__family', 'template_structure')

    def handle(self, *args, **options):
        # run with profiling
        if options['profile']:
            profiler = Profile()
            profiler.runcall(self._handle, *args, **options)
            profiler.print_stats()
        # run without profiling
        else:
            self._handle(*args, **options)

    def _handle(self, *args, **options):
        try:
            self.logger.info('ASSIGNING STRUCTURE TEMPLATES FOR PROTEINS')
            self.prepare_input(options['proc'], self.pconfs)
            self.logger.info('COMPLETED ASSIGNING STRUCTURE TEMPLATES FOR PROTEINS')
        except Exception as msg:
            print(msg)
            self.logger.error(msg)

    def main_func(self, positions, iteration):
        # pconfs
        if not positions[1]:
            pconfs = self.pconfs[positions[0]:]
        else:
            pconfs = self.pconfs[positions[0]:positions[1]]

        # find templates
        for pconf in pconfs:
            # filter structure sequence queryset to include only sequences from within the same class
            pconf_class = pconf.protein.family.slug[:3]
            class_sps = self.structures.filter(
                protein_conformation__protein__parent__family__slug__startswith=pconf_class)
            sps = []
            sps_str = []
            if class_sps.exists():
                for structure in class_sps:
                    sps.append(structure.protein_conformation.protein.parent) # use the wild-type sequence for main tpl
                    sps_str.append(structure.protein_conformation.protein) # use the structure sequence for segment tpl
            else:
                for structure in self.structures:
                    sps.append(structure.protein_conformation.protein.parent)
                    sps_str.append(structure.protein_conformation.protein)

            # overall
            template = self.find_segment_template(pconf, sps, self.segments)
            template_structure = self.fetch_template_structure(self.structures, template.protein.entry_name)
            pconf.template_structure = template_structure
            pconf.save()
            self.logger.info("Assigned {} as overall template for {}".format(template_structure, pconf))

            # for each segment
            for segment in self.segments:
                template = self.find_segment_template(pconf, sps_str, [segment])
                template_structure = self.fetch_template_structure(self.structures, template.protein.parent.entry_name)
                pcts, created = ProteinConformationTemplateStructure.objects.get_or_create(protein_conformation=pconf,
                    protein_segment=segment, defaults={'structure': template_structure})
                if pcts.structure != template_structure:
                    pcts.structure = template_structure
                    pcts.save()
                self.logger.info("Assigned {} as {} template for {}".format(template_structure, segment, pconf))

    def find_segment_template(self, pconf, sconfs, segments):
            a = Alignment()
            a.load_reference_protein(pconf.protein)
            a.load_proteins(sconfs)
            a.load_segments(segments)
            a.build_alignment()
            a.calculate_similarity()

            return a.proteins[1]

    def fetch_template_structure(self, structures, template_protein):
        for structure in structures:
            structure_protein = structure.protein_conformation.protein.parent.entry_name
            if structure_protein == template_protein:
                return structure
