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
import yaml
import pprint

class Command(BaseBuild):  
    help = 'Build automated chimeric GPCR homology models'
    
    def add_arguments(self, parser):
        super(Command, self).add_arguments(parser=parser)
        parser.add_argument('--verbose', help='Print specific outliers', default=False, action='store_true')
        
    def handle(self, *args, **options):
        structures = Structure.objects.filter(refined=False).prefetch_related('protein_conformation__protein__parent','pdb_code')
        structures_with_issue = []
        missing_helices = {}
        segments_query_obj = ProteinSegment.objects.filter(proteinfamily="GPCR")
        with open(os.sep.join([settings.DATA_DIR, 'structure_data','annotation','xtal_segends.yaml']), 'r') as anno_f:
            annotations = yaml.load(anno_f)
        for s in structures:
            resis = Residue.objects.filter(protein_conformation=s.protein_conformation).prefetch_related('display_generic_number')
            c = 0
            segments = OrderedDict((i,[]) for i in segments_query_obj)
            x50s = [i.display_generic_number.label for i in resis.filter(display_generic_number__label__in=['1.50x50', '2.50x50', '3.50x50', '4.50x50', '5.50x50', '6.50x50', '7.50x50'])]
            parent_x50_resis = Residue.objects.filter(protein_conformation__protein=s.protein_conformation.protein.parent, 
                                                      display_generic_number__label__in=['1.50x50', '2.50x50', '3.50x50', '4.50x50', '5.50x50', '6.50x50', '7.50x50']).prefetch_related('display_generic_number')
            missing_helices[s] = []
            missing_a_helix = False
            for x in ['1.50x50', '2.50x50', '3.50x50', '4.50x50', '5.50x50', '6.50x50', '7.50x50']:
                if x not in x50s:
                    TMb = annotations[s.protein_conformation.protein.parent.entry_name+'_'+s.pdb_code.index][x[0]+'b']
                    TMe = annotations[s.protein_conformation.protein.parent.entry_name+'_'+s.pdb_code.index][x[0]+'e']
                    if TMe=='-':
                        continue
                    if s.pdb_code.index=='4L6R' and x=='4.50x50' and TMb==274:
                        continue
                    if TMe>parent_x50_resis.get(display_generic_number__label=x).sequence_number:
                        missing_helices[s].append(x)
                        missing_a_helix = True
            if not missing_a_helix:
                del missing_helices[s]
            else:
                print(s.pdb_code.index,'missing',missing_helices[s])
            for r in resis:
                if c==0:
                    this_segment = r.protein_segment
                if this_segment==None:
                    continue
                if r.protein_segment==None:
                    if len(Fragment.objects.filter(residue=r))==0 and s not in structures_with_issue:
                        structures_with_issue.append(s)
                    continue
                if len(segments[r.protein_segment])>0:
                    if r.sequence_number-segments[r.protein_segment][-1]>400:
                        if options['verbose']:
                            print('Outlier: {} {} {} following {}'.format(s, r.protein_segment, r.sequence_number, segments[r.protein_segment][-1]))
                        if s not in structures_with_issue:
                            structures_with_issue.append(s)
                    else:
                        segments[r.protein_segment].append(r.sequence_number)
                else:
                    segments[r.protein_segment].append(r.sequence_number)
                c+=1
        if len(missing_helices)>0 or len(structures_with_issue)>0:
            print_out = 'Missing helices: '+missing_helices+'\n'+'Structures with issue: '+structures_with_issue
        else:
            print_out = 'No structure residue issues detected'
        self.logger.info('Check structure residues: '+print_out)

