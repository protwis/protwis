from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import connection 

from protein.models import (Protein, ProteinConformation, ProteinSequenceType, ProteinSegment,
    ProteinConformationTemplateStructure)
from structure.models import Structure
from common.alignment import Alignment

import logging
from cProfile import Profile
from optparse import make_option
from multiprocessing import Queue, Process

class Command(BaseCommand):
    help = 'Goes though all protein records and finds the best structural templates for the whole 7TM domain, and ' \
        + 'each helix individually'

    def add_arguments(self, parser):
        parser.add_argument('--njobs', action='store', dest='njobs', help='Number of jobs to run')
        parser.add_argument('--profile', action='store_true', dest='profile', default=False,
            help='Profile the script with cProfile')

    logger = logging.getLogger(__name__)

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
        # how many jobs to run?
        if 'njobs' in options and options['njobs']:
            njobs = int(options['njobs'])
        else:
            njobs = 1

        try:
            self.prepare_input(njobs)
        except Exception as msg:
            print(msg)
            self.logger.error(msg)

    def prepare_input(self, njobs):
        self.logger.info('ASSIGNING STRUCTURE TEMPLATES FOR PROTEINS')
        
        q = Queue()
        procs = list()
        num_pconfs = self.pconfs.count()
        
        # make sure not to use more jobs than proteins (chunk size will be 0, which is not good)
        if njobs > num_pconfs:
            njobs = num_pconfs

        chunk_size = int(num_pconfs / njobs)
        connection.close()
        for i in range(0, njobs):
            first = chunk_size * i
            if i == njobs - 1:
                last = False
            else:
                last = chunk_size * (i + 1)
    
            p = Process(target=self.find_protein_templates, args=([(first, last)]))
            procs.append(p)
            p.start()

        for p in procs:
            p.join()

        self.logger.info('COMPLETED ASSIGNING STRUCTURE TEMPLATES FOR PROTEINS')

    def find_protein_templates(self, positions):
        # pconfs
        if not positions[1]:
            pconfs = self.pconfs[positions[0]:]
        else:
            pconfs = self.pconfs[positions[0]:positions[1]]

        # segments
        segments = ProteinSegment.objects.filter(slug__in=settings.REFERENCE_POSITIONS)

        # fetch wild-type sequences of receptors with available structures
        structures = Structure.objects.order_by('protein_conformation__protein__parent', 'resolution').distinct(
            'protein_conformation__protein__parent').prefetch_related('protein_conformation__protein__parent__family')

        # find templates
        for pconf in pconfs:
            # filter structure sequence queryset to include only sequences from within the same class
            pconf_class = pconf.protein.family.slug[:3]
            class_sps = structures.filter(protein_conformation__protein__parent__family__slug__startswith=pconf_class)
            sps = []
            if class_sps.exists():
                for structure in class_sps:
                    sps.append(structure.protein_conformation.protein.parent) # use the wild-type sequence
            else:
                for structure in structures:
                    sps.append(structure.protein_conformation.protein.parent) # use the wild-type sequence

            # overall
            template = self.find_segment_template(pconf, sps, segments)
            template_structure = self.fetch_template_structure(structures, template.protein.entry_name)
            pconf.template_structure = template_structure
            pconf.save()
            self.logger.info("Assigned {} as overall template for {}".format(template_structure, pconf))

            # for each segment
            for segment in segments:
                template = self.find_segment_template(pconf, sps, [segment])
                template_structure = self.fetch_template_structure(structures, template.protein.entry_name)
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