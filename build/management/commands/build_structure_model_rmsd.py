# -*- coding: utf-8 -*-
"""
Created on 2020-11-13 10:16:25.482439

@author: Gaspar Pandy
"""
from build.management.commands.base_build import Command as BaseBuild
from django.conf import settings

# from common.alignment import Alignment
from protein.models import Protein, ProteinSegment
# from residue.models import Residue
from structure.models import Structure, StructureModelRMSD
from structure.sequence_parser import SequenceParser
from common.tools import test_model_updates
from datetime import datetime
import csv, os, pprint
import django.apps


starttime = datetime.now()

class Command(BaseBuild):

    #Setting the variables for the test tracking of the model upadates
    tracker = {}
    all_models = django.apps.apps.get_models()[6:]
    test_model_updates(all_models, tracker, initialize=True)
    def add_arguments(self, parser):
        super(Command, self).add_arguments(parser=parser)
        parser.add_argument('--verbose',
            action='store_true',
            default=False,
            help='Verbose')
        parser.add_argument('-u', '--purge',
            action='store_true',
            dest='purge',
            default=False,
            help='Purge existing records')

    def handle(self, *args, **options):
        bsmr = BuildStructureModelRMSD()
        if options['purge']:
            bsmr.purge()
            self.tracker = {}
            test_model_updates(self.all_models, self.tracker, initialize=True)
        bsmr.parse_data_file()
        bsmr.run_build()
        test_model_updates(self.all_models, self.tracker, check=True)


class BuildStructureModelRMSD():
    def __init__(self):
        self.data = []

    def purge(self):
        StructureModelRMSD.objects.all().delete()

    def parse_data_file(self):
        self.data_file = os.sep.join([settings.DATA_DIR, 'structure_data', 'model_rmsd.csv'])
        with open(self.data_file, newline='') as csvfile:
            reader = csv.reader(csvfile, delimiter=',')
            self.data = [i for i in reader]

    def run_build(self):
        for i in self.data[1:]:
            i = [None if j=='-' else float(j) if '.' in j and len(j)==3 else j for j in i]
            pdb, _, version, overall_all, overall_backbone, TM_all, TM_backbone, H8, ICL1, ECL1, ICL2, ECL2, ECL3, notes = i
            target_structure = Structure.objects.get(pdb_code__index=pdb.upper())
            # main_template = Structure.objects.get(pdb_code__index=main_temp.upper())
            # a = Alignment()
            # a.load_reference_protein(target_structure.protein_conformation.protein.parent)
            # a.load_proteins([main_template.protein_conformation.protein.parent])
            # resis = Residue.objects.filter(protein_conformation__protein=target_structure.protein_conformation.protein.parent)
            # segments = resis.order_by('protein_segment__id').distinct('protein_segment__id').values_list('protein_segment',flat=True)
            # a.load_segments(ProteinSegment.objects.filter(id__in=segments))
            # a.build_alignment()
            # a.remove_non_generic_numbers_from_alignment()
            # a.calculate_similarity()
            # seq_sim = a.proteins[1].similarity
            # seq_id = a.proteins[1].identity
            main_template = None
            seq_id, seq_sim = None, None
            smr, created = StructureModelRMSD.objects.get_or_create(target_structure=target_structure,
                                                                    main_template=main_template,
                                                                    version='{}-{}-{}'.format(version[-4:], version[3:5], version[:2]),
                                                                    seq_id=seq_id, seq_sim=seq_sim,
                                                                    overall_all=overall_all, overall_backbone=overall_backbone, TM_all=TM_all,
                                                                    TM_backbone=TM_backbone, H8=H8,
                                                                    ICL1=ICL1, ECL1=ECL1, ICL2=ICL2, ECL2=ECL2, ECL3=ECL3,
                                                                    notes=notes)
