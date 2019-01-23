from build.management.commands.base_build import Command as BaseBuild
from django.conf import settings

from protein.models import Protein, ProteinConformation, ProteinAnomaly, ProteinState, ProteinSegment
from residue.models import Residue
from structure.models import *

from collections import OrderedDict
import os
import logging
import sys
from datetime import datetime, date
import traceback

class Command(BaseBuild):  
    help = 'Build automated chimeric GPCR homology models'
    
    def add_arguments(self, parser):
        super(Command, self).add_arguments(parser=parser)
        parser.add_argument('--verbose', help='Print specific outliers', default=False, action='store_true')
        
    def handle(self, *args, **options):
        structures = Structure.objects.filter(refined=False)
        structures_with_issue = []
        segments_query_obj = ProteinSegment.objects.filter(proteinfamily="GPCR")
        for s in structures:
            resis = Residue.objects.filter(protein_conformation=s.protein_conformation)
            c = 0
            segments = OrderedDict((i,[]) for i in segments_query_obj)
            for r in resis:
                if c==0:
                    this_segment = r.protein_segment
                if this_segment==None:
                    continue
                if r.protein_segment==None:
                    if s not in structures_with_issue:
                        structures_with_issue.append(s)
                    continue
                if len(segments[r.protein_segment])>0:
                    if r.sequence_number-segments[r.protein_segment][-1]>100:
                        if options['verbose']:
                            print('Outlier: {} {} {} following {}'.format(s, r.protein_segment, r.sequence_number, segments[r.protein_segment][-1]))
                        if s not in structures_with_issue:
                            structures_with_issue.append(s)
                    else:
                        segments[r.protein_segment].append(r.sequence_number)
                else:
                    segments[r.protein_segment].append(r.sequence_number)
                c+=1
        print(structures_with_issue)


