from django.core.management.base import BaseCommand, CommandError
from django.core.management import call_command
from django.conf import settings
from django.db import connection
from django.template.loader import render_to_string
from protein.models import Protein
from residue.models import ResidueGenericNumber, ResidueGenericNumberEquivalent
from common import definitions
from common.selection import SelectionItem
from common.alignment_gpcr import Alignment
import xlsxwriter, xlrd


import logging, json, os

class Alignment_cmd(Alignment):

    def evaluate_sites(self, site_defs):
        # go through all proteins and match against site definitions
        for protein in self.proteins:
            for k, segment in enumerate(protein.alignment.values(), start = 1):
                num_matched = 0
                min_match = site_defs[k]['min_match']
                for position in segment:
                    # position example: ['6x49', '6.49x49', 'L', 'GPCRdb(A)', 282, 282]
                    if position[2] in site_defs[k]['amino_acids'][position[0]]:
                        num_matched += 1
                        if num_matched >= min_match:
                            break
                else:
                    # if the protein sequence does not match the definitions, store it in non_matching_proteins
                    self.non_matching_proteins.append(protein)
                    break

        # remove non-matching proteins from protein list
        self.proteins = [p for p in self.proteins if p not in self.non_matching_proteins]


class Command(BaseCommand):

    help = "Function evaluating sequence compatibility with wildtype protein."

    logger = logging.getLogger(__name__)

    def add_arguments(self, parser):
        parser.add_argument('-f', dest='site_file', type=str, help="Input file with site definition")
        parser.add_argument('codes', nargs='+')

    def handle(self, *args, **options):

        #Get the proteins
        targets = []
        for up_code in options['codes']:
            print(up_code)
            try:
                targets.append(Protein.objects.get(entry_name=up_code.strip().lower()))
            except Exception as msg:
                print("Cannot find a protein {}.\n{}".format(up_code, msg))
        print(targets)
        #Get the data from excel
        segments = []
        selection_type = 'segments'
        selection_subtype = 'site_residue'

        workbook = xlrd.open_workbook(options['site_file'])
        worksheets = workbook.sheet_names()
        for worksheet_name in worksheets:

            site_defs = {}
            worksheet = workbook.sheet_by_name(worksheet_name)
            for row in worksheet.get_rows():

        #for row in ws.rows:
                if len(row) < 5:
                    continue
                group_id = int(row[0].value)
                min_match = int(row[1].value)
                try:
                    position = ResidueGenericNumberEquivalent.objects.get(label=row[2].value, scheme__slug=row[3].value)
                except e:
                    print(e)
                    continue
                feature = row[4].value

                properties = {}
                properties['feature'] = feature
                properties['site_residue_group'] = group_id
                properties['amino_acids'] = ','.join(definitions.AMINO_ACID_GROUPS[feature]) if feature != 'custom' else row[5].value
                if group_id not in site_defs:
                    site_defs[group_id] = {'min_match': min_match,
                                           'positions': {},
                                           'amino_acids': {},
                                           }
                site_defs[group_id]['positions'][position.label] = properties['feature']
                site_defs[group_id]['amino_acids'][position.label] = properties['amino_acids']
                segments.append(SelectionItem(selection_subtype, position, properties))

        #The alignment
        a = Alignment_cmd()

        # group residue selection (for site search)
        a.use_residue_groups = True

        # load data from selection into the alignment
        a.load_proteins(targets)
        a.load_segments(segments)
        # build the alignment data matrix
        a.build_alignment()

        # evaluate sites
        a.evaluate_sites(site_defs)

        num_of_sequences = len(a.proteins)
        num_of_non_matching_sequences = len(a.non_matching_proteins)
        num_residue_columns = len(a.positions) + len(a.segments)

        print("Seqs: {}\tNot matching: {}".format(num_of_sequences, num_of_non_matching_sequences))
        open("test_output,csv", "w").write(render_to_string('alignment/alignment_csv.html', {'a': a, 'num_of_sequences': num_of_sequences,
        'num_residue_columns': num_residue_columns}))
