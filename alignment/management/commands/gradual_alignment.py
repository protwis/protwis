from build.management.commands.base_build import Command as BaseBuild
from django.db import connection
from django.conf import settings

from common.alignment import Alignment
from protein.models import Protein, ProteinSegment

import os


class Command(BaseBuild):
    help = 'Upload GPCRdb structure models from zips to db'

    def add_arguments(self, parser):
        super(Command, self).add_arguments(parser=parser)

    def handle(self, *args, **options):
        self.build_gradual_alignment()
        self.merge_segments()

    def build_gradual_alignment(self):
        ps = Protein.objects.filter(family__slug__startswith='001', accession__isnull=False, sequence_type__slug='wt').filter(species__common_name='Human')

        segments = ['N-term','TM1','ICL1','TM2','ECL1','TM3','ICL2','TM4','ECL2','TM5','ICL3','TM6','ECL3','TM7','ICL4','H8','C-term']
        for i, seg in enumerate(segments):
            a = Alignment()
            ######
            if seg=='N-term':
                gn_nums = ','
                seg_header = ','
            else:
                gn_nums = ''
                seg_header = ''
            ######

            segment = ProteinSegment.objects.get(slug=seg)
            a.load_proteins(ps)
            a.load_segments([segment])
            a.build_alignment()

            i = str(i)
            if len(i)==1:
                i = '0'+i

            with open('./predefined_alignments/{}_{}.csv'.format(i,seg), 'w') as f:
                seg_header+=seg+','*(len(list(a.segments.values())[0])-1)+',\n'
                f.write(seg_header)
                gn_nums+=','.join(list(a.segments.values())[0])+',\n'
                f.write(gn_nums)
                for p in a.proteins:
                    if seg=='N-term':
                        line = '{},{}'.format(p.protein.entry_name,','.join([j[2] for j in p.alignment[seg]]))+',\n'
                    else:
                        line = ','.join([j[2] for j in p.alignment[seg]])+',\n'
                    f.write(line)

    def merge_segments(self):
        files = sorted(os.listdir('./predefined_alignments'))
        master_array = []
        for f in files:
            with open('./predefined_alignments/'+f, 'r') as fi:
                lines = fi.readlines()
            for i, line in enumerate(lines):
                if f.startswith('00'):
                    if i==0:
                        line = line[:-1]
                    if i==1:
                        line_split = line.split(',')
                        new_gns = ''
                        for l in line_split:
                            if 'x' not in l:
                                new_gns+=','
                            else:
                                new_gns+=l+','
                        line = new_gns[:-1]
                    master_array.append(line.replace('\n',''))
                else:
                    if i==0:
                        line = line[:-1]
                    if i==1:
                        line_split = line.split(',')
                        new_gns = ''
                        for l in line_split:
                            if 'x' not in l:
                                new_gns+=','
                            else:
                                new_gns+=l+','
                        line = new_gns[:-1]
                    master_array[i]+=line.replace('\n','')

        with open('./master_array.csv', 'w') as ma:
            for l in master_array:
                ma.write(l+'\n')

