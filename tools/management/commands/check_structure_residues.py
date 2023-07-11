from build.management.commands.base_build import Command as BaseBuild
from django.conf import settings
from django.db import connection
from django.db.models import Count

from protein.models import Protein, ProteinConformation, ProteinAnomaly, ProteinState, ProteinSegment
from residue.models import Residue, ResidueGenericNumber
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
        structures = Structure.objects.all().exclude(structure_type__slug__startswith='af-').prefetch_related('protein_conformation__protein__parent','pdb_code').annotate(dc=Count('distances'))
        structures_with_issue = []
        missing_helices = {}
        segments_query_obj = ProteinSegment.objects.filter(proteinfamily="GPCR")
        all_gns = set(ResidueGenericNumber.objects.filter(scheme__slug='gpcrdb').all().values_list('label',flat=True))
        # print(all_gns)
        missing_gns_count = OrderedDict()
        # raise
        with open(os.sep.join([settings.DATA_DIR, 'structure_data','annotation','xtal_segends.yaml']), 'r') as anno_f:
            annotations = yaml.load(anno_f, Loader=yaml.FullLoader)
        for s in structures:
            resis = Residue.objects.filter(protein_conformation=s.protein_conformation).prefetch_related('display_generic_number','protein_segment')
            wt_resis = Residue.objects.filter(protein_conformation__protein=s.protein_conformation.protein.parent)

            wt_gn = wt_resis.exclude(display_generic_number=None).count()
            s_gn = resis.exclude(display_generic_number=None).count()
            fraction = s_gn / wt_gn

            s_gns = set(resis.exclude(display_generic_number=None).all().values_list('generic_number__label',flat=True))
            missing_gns = all_gns-s_gns
            for gn in missing_gns:
                if gn not in missing_gns_count:
                    missing_gns_count[gn] = {'count':0, 'pdbs' : []}
                missing_gns_count[gn]['count'] += 1
                missing_gns_count[gn]['pdbs'].append(str(s))

            print(s,"distances",s.dc,"WT residues",wt_resis.count(),"PDB residues",resis.count(),"WT GNs",wt_gn,"PDB GNs",s_gn,"Fraction",fraction)
            # print('missing gns',missing_gns)
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


        cut_off = 10
        print("=== GN PDBS MISSING, PDB list === CUTOFF", cut_off)
        missing_gns_count = OrderedDict(sorted(missing_gns_count.items(), key=lambda t: t[1]['count']))
        #missing_gns_count = OrderedDict(sorted(missing_gns_count.items(), key=lambda x: x['count']))
        pdbs_with_low_occuring_missing = OrderedDict()
        for m, v in missing_gns_count.items():
            if cut_off>10: break
            print(m,v['count'],v['pdbs'])
            for pdb in v['pdbs']:
                if pdb not in pdbs_with_low_occuring_missing:
                    pdbs_with_low_occuring_missing[pdb] = {'count':0,'gns':[]}
                pdbs_with_low_occuring_missing[pdb]['count'] += 1
                pdbs_with_low_occuring_missing[pdb]['gns'].append(m)


        print("=== MAIN PDB CULPRITS ===")
        pdbs_with_low_occuring_missing = OrderedDict(sorted(pdbs_with_low_occuring_missing.items(), key=lambda t: -t[1]['count']))
        for pdb, v in pdbs_with_low_occuring_missing.items():
            print(pdb,v)
        #print(pdbs_with_low_occuring_missing)
        #print(missing_gns_count)


        if len(missing_helices)>0 or len(structures_with_issue)>0:
            print_out = 'Missing helices: '+missing_helices+'\n'+'Structures with issue: '+structures_with_issue
        else:
            print_out = 'No structure residue issues detected'
        self.logger.info('Check structure residues: '+print_out)
