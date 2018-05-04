from build.management.commands.base_build import Command as BaseBuild
from django.db.models import Q
from django.conf import settings

from protein.models import Protein, ProteinConformation, ProteinAnomaly, ProteinState, ProteinSegment
from residue.models import Residue
from residue.functions import dgn, ggn
from structure.models import *
from structure.functions import HSExposureCB, PdbStateIdentifier
from common.alignment import AlignedReferenceTemplate, GProteinAlignment
from common.definitions import *
from common.models import WebLink
from signprot.models import SignprotComplex
import structure.structural_superposition as sp
import structure.assign_generic_numbers_gpcr as as_gn
import structure.homology_models_tests as tests

import Bio.PDB as PDB
from Bio import pairwise2
from modeller import *
from modeller.automodel import *
from collections import OrderedDict
import os
import subprocess
import shlex
import logging
import pprint
from io import StringIO
import sys
import re
import zipfile
import shutil
import math
from copy import deepcopy
from datetime import datetime, date
import yaml
import traceback


startTime = datetime.now()


class Command(BaseBuild):  
	help = 'Test scripts'
	
	def add_arguments(self, parser):
		super(Command, self).add_arguments(parser=parser)
		parser.add_argument('--gns', help="Specifiy generic numbers involved in calculation for TestStateIdentifier", default=False, nargs='+')
		parser.add_argument('--only_xtals', help="Only run TestStateIdentifier on xtals", default=False, action='store_true')
		parser.add_argument('--cutoffs', help="Set inactive and intermediate cutoffs", default=False, nargs='+')
		parser.add_argument('--segment', help="Set protein segment label for StructuralStatistics", default=False)


	def handle(self, *args, **options):
		# strs = Structure.objects.filter(refined=False).exclude(protein_conformation__protein__parent__family__parent__parent__parent__slug__in=['002','005'])
		# for s in strs:
		# 	print(s,',',s.state.slug)
		# return 0
		# t = TestStateIdentifierSets(options['only_xtals'])
		# t.run()
		# with open('./structure/cutoff_test_02.csv','w') as fw:
		# with open('./structure/cutoff_test_04.csv','r') as f:
		# 	lines = f.readlines()
		# 	big_l = []
		# 	for line in lines:
		# 		li = line.split(',')
		# 		li[4] = float(li[4])
		# 		li[-1] = float(li[-1])
		# 		big_l.append(li)
		# 	# big_l = sorted(big_l, key=lambda x: (x[4]))
		# 	c = 0
		# 	for l in big_l:
		# 		c+=1
		# 		if c<104:
		# 			print("['{}','{}','{}','{}',{},{}],".format(l[0],l[1],l[2],l[3],l[7],l[8]))
		# self.data = [['2x39','6x35','3x47','7x53',-0.5,8],
		# 		['2x39','6x35','3x44','7x53',-2,6],
		# 		['2x39','6x38','3x47','7x52',-0.5,6],
		# 		['2x39','6x38','3x44','7x52',-2,4],
		# 		['2x39','6x35','3x47','7x53',1,8],
		# 		['2x39','6x37','3x46','7x53',-0.5,6],
		# 		['2x39','6x35','3x45','7x53',-2,6],
		# 		['2x39','6x37','3x46','7x53',-2,6],
		# 		['2x39','6x38','3x47','7x52',1,6],
		# 		['2x39','6x38','3x46','7x52',1,6],
		# 		['2x39','6x37','3x47','7x53',-2,4],
		# 		['2x39','6x35','3x46','7x53',1,8],
		# 		['2x40','6x35','3x47','7x53',1,8],
		# 		['2x39','6x38','3x47','7x53',-0.5,8],
		# 		['2x42','6x35','3x44','7x53',-0.5,8],
		# 		['2x39','6x38','3x45','7x53',-2,6],
		# 		['2x42','6x35','3x45','7x53',1,8],
		# 		['2x42','6x35','3x44','7x53',1,8],
		# 		['2x39','6x35','3x47','7x52',-0.5,6],
		# 		['2x39','6x35','3x45','7x53',-0.5,6],
		# 		['2x39','6x35','3x44','7x53',-0.5,6],
		# 		['2x39','6x36','3x46','7x53',-0.5,6],
		# 		['2x40','6x35','3x45','7x53',-0.5,6],
		# 		['2x40','6x37','3x46','7x53',-0.5,6],
		# 		['2x41','6x37','3x45','7x53',-2,6],
		# 		['2x41','6x37','3x44','7x53',-2,6],
		# 		['2x39','6x35','3x45','7x52',-2,4],
		# 		['2x39','6x35','3x44','7x52',-2,4],
		# 		['2x42','6x37','3x44','7x53',-2,4],
		# 		['2x41','6x37','3x46','7x53',1,8],
		# 		['2x39','6x35','3x45','7x53',1,6],
		# 		['2x39','6x36','3x46','7x53',1,6],
		# 		['2x41','6x38','3x45','7x52',1,6],
		# 		['2x42','6x38','3x46','7x52',1,6],
		# 		['2x39','6x38','3x45','7x53',-0.5,6],
		# 		['2x40','6x35','3x44','7x53',-0.5,6],
		# 		['2x41','6x37','3x45','7x53',-0.5,6],
		# 		['2x41','6x37','3x44','7x53',-0.5,6],
		# 		['2x39','6x35','3x45','7x52',-0.5,4],
		# 		['2x40','6x38','3x45','7x52',-0.5,4],
		# 		['2x42','6x38','3x44','7x52',-0.5,4],
		# 		['2x42','6x38','3x44','7x53',-2,6],
		# 		['2x39','6x35','3x47','7x50',-2,4],
		# 		['2x39','6x35','3x46','7x51',-2,4],
		# 		['2x39','6x38','3x45','7x52',-2,4],
		# 		['2x41','6x38','3x45','7x51',-2,4],
		# 		['2x39','6x38','3x47','7x53',1,8],
		# 		['2x41','6x37','3x47','7x53',1,8],
		# 		['2x42','6x38','3x47','7x53',1,8],
		# 		['2x39','6x35','3x44','7x53',1,6],
		# 		['2x39','6x37','3x46','7x53',1,6],
		# 		['2x40','6x35','3x45','7x53',1,6],
		# 		['2x40','6x35','3x44','7x53',1,6],
		# 		['2x40','6x37','3x46','7x53',1,6],
		# 		['2x39','6x35','3x46','7x50',-0.5,6],
		# 		['2x39','6x38','3x44','7x53',-0.5,6],
		# 		['2x41','6x35','3x44','7x51',-0.5,6],
		# 		['2x42','6x37','3x46','7x53',-0.5,6],
		# 		['2x42','6x38','3x45','7x53',-0.5,6],
		# 		['2x39','6x38','3x44','7x52',-0.5,4],
		# 		['2x39','6x38','3x44','7x53',-2,6],
		# 		['2x39','6x36','3x44','7x53',-2,4],
		# 		['2x39','6x35','3x47','7x52',1,6],
		# 		['2x39','6x35','3x46','7x52',1,6],
		# 		['2x40','6x38','3x47','7x52',1,6],
		# 		['2x41','6x38','3x46','7x51',1,6],
		# 		['2x42','6x35','3x45','7x53',1,6],
		# 		['2x40','6x35','3x44','7x53',-0.5,8],
		# 		['2x41','6x37','3x47','7x53',-0.5,8],
		# 		['2x42','6x38','3x47','7x53',-0.5,8],
		# 		['2x39','6x35','3x47','7x53',-0.5,6],
		# 		['2x39','6x38','3x46','7x52',-0.5,6],
		# 		['2x40','6x38','3x44','7x52',-0.5,6],
		# 		['2x42','6x35','3x44','7x52',-0.5,6],
		# 		['2x42','6x38','3x44','7x53',-0.5,6],
		# 		['2x39','6x35','3x46','7x51',-0.5,4],
		# 		['2x39','6x35','3x44','7x52',-0.5,4],
		# 		['2x39','6x37','3x47','7x53',-0.5,4],
		# 		['2x39','6x35','3x44','7x53',-2,4],
		# 		['2x39','6x35','3x44','7x50',-2,4],
		# 		['2x39','6x38','3x46','7x51',-2,4],
		# 		['2x40','6x37','3x44','7x53',-2,4],
		# 		['2x42','6x35','3x44','7x51',-2,4],
		# 		['2x42','6x37','3x45','7x53',-2,4],
		# 		['2x42','6x38','3x44','7x52',-2,4],
		# 		['2x42','6x38','3x44','7x50',-2,4],
		# 		['2x40','6x35','3x44','7x53',1,8],
		# 		['2x40','6x38','3x47','7x53',1,8],
		# 		['2x41','6x38','3x44','7x52',1,8],
		# 		['2x39','6x35','3x46','7x50',1,6],
		# 		['2x41','6x37','3x46','7x52',1,6],
		# 		['2x41','6x37','3x45','7x53',1,6],
		# 		['2x42','6x37','3x46','7x53',1,6],
		# 		['2x42','6x38','3x45','7x53',1,6],
		# 		['2x39','6x35','3x44','7x52',1,4],
		# 		['2x42','6x35','3x44','7x50',-0.5,6],
		# 		['2x39','6x35','3x47','7x52',-0.5,4],
		# 		['2x39','6x37','3x46','7x53',-0.5,4],
		# 		['2x42','6x37','3x44','7x53',-0.5,4],
		# 		['2x39','6x35','3x45','7x53',-2,4],
		# 		['2x39','6x37','3x46','7x53',-2,4],
		# 		['2x42','6x38','3x46','7x51',-2,4],
		# 		['2x42','6x38','3x45','7x52',-2,4]]

		# self.data = [['2x39','6x35','3x44','7x53',-2,6.5],
		# 			['2x39','6x35','3x47','7x53',-0.5,8.0],
		# 			['2x39','6x35','3x47','7x53',0,8.0],
		# 			['2x39','6x35','3x44','7x53',-2,6.0],
		# 			['2x39','6x35','3x47','7x53',-0.5,8.5],
		# 			['2x39','6x35','3x47','7x53',0,8.5],
		# 			['2x39','6x35','3x47','7x53',-0.5,7.5],
		# 			['2x39','6x35','3x47','7x53',0,7.5],
		# 			['2x39','6x35','3x44','7x53',-2,5.5],
		# 			['2x39','6x38','3x47','7x52',-1,6.0],
		# 			['2x39','6x38','3x47','7x52',-0.5,6.0],
		# 			['2x39','6x38','3x47','7x52',0,6.0],
		# 			['2x39','6x38','3x47','7x52',0.5,6.0],
		# 			['2x39','6x35','3x44','7x52',-2.5,5.0],
		# 			['2x39','6x35','3x44','7x52',-2,5.0],
		# 			['2x39','6x35','3x44','7x53',-1.5,6.5],
		# 			['2x39','6x35','3x45','7x53',-2.5,6.0],
		# 			['2x39','6x35','3x45','7x53',-1.5,6.0],
		# 			['2x39','6x35','3x46','7x53',1,8.5],
		# 			['2x39','6x35','3x46','7x53',1,9.0],
		# 			['2x39','6x35','3x47','7x53',0.5,8.0],
		# 			['2x39','6x37','3x46','7x53',-1.5,5.5],
		# 			['2x39','6x37','3x46','7x53',-1.5,6.0],
		# 			['2x39','6x38','3x44','7x52',-2,4.0],
		# 			['2x39','6x38','3x44','7x52',-2,5.0],
		# 			['2x40','6x35','3x44','7x53',0,7.5],
		# 			['2x42','6x35','3x44','7x53',0,8.5],
		# 			['2x39','6x35','3x44','7x52',-1.5,5.0],
		# 			['2x39','6x35','3x44','7x53',-2,7.0],
		# 			['2x39','6x35','3x44','7x53',-1.5,6.0],
		# 			['2x39','6x35','3x45','7x53',-2.5,5.5],
		# 			['2x39','6x35','3x45','7x53',-2.5,6.5],
		# 			['2x39','6x35','3x45','7x53',-2,6.0],
		# 			['2x39','6x35','3x45','7x53',-1.5,5.5],
		# 			['2x39','6x35','3x45','7x53',-1.5,6.5],
		# 			['2x39','6x35','3x45','7x53',-1,6.0],
		# 			['2x39','6x35','3x46','7x53',1.5,8.5],
		# 			['2x39','6x35','3x46','7x53',1.5,9.0],
		# 			['2x39','6x35','3x47','7x53',0.5,8.5],
		# 			['2x39','6x35','3x47','7x53',1,8.0],
		# 			['2x39','6x36','3x44','7x53',-3,3.0],
		# 			['2x39','6x37','3x46','7x53',-2,5.5],
		# 			['2x39','6x37','3x46','7x53',-2,6.0],
		# 			['2x39','6x37','3x46','7x53',-1,5.5],
		# 			['2x39','6x37','3x46','7x53',-1,6.0],
		# 			['2x39','6x37','3x46','7x53',-0.5,5.5],
		# 			['2x39','6x37','3x46','7x53',-0.5,6.0],
		# 			['2x39','6x37','3x47','7x53',-3,4.0],
		# 			['2x39','6x37','3x47','7x53',-3,4.5],
		# 			['2x39','6x37','3x47','7x53',-2.5,4.0],
		# 			['2x39','6x37','3x47','7x53',-2.5,4.5],
		# 			['2x39','6x38','3x44','7x52',-2,4.5],
		# 			['2x39','6x38','3x44','7x52',-1.5,4.0],
		# 			['2x39','6x38','3x44','7x52',-1.5,5.0],
		# 			['2x39','6x38','3x44','7x52',-1,4.0],
		# 			['2x39','6x38','3x44','7x52',-1,5.0],
		# 			['2x39','6x38','3x45','7x52',-2.5,3.0],
		# 			['2x39','6x38','3x45','7x52',-2.5,3.5],
		# 			['2x39','6x38','3x46','7x52',0.5,6.0],
		# 			['2x39','6x38','3x46','7x52',0.5,6.5],
		# 			['2x39','6x38','3x47','7x52',-1,5.5],
		# 			['2x39','6x38','3x47','7x52',-0.5,5.5],
		# 			['2x39','6x38','3x47','7x52',0,5.5],
		# 			['2x39','6x38','3x47','7x52',0.5,5.5],
		# 			['2x40','6x35','3x45','7x53',-0.5,7.0],
		# 			['2x40','6x38','3x44','7x52',0,5.5],
		# 			['2x40','6x38','3x47','7x53',1.5,8.5],
		# 			['2x41','6x37','3x44','7x53',-2,6.5],
		# 			['2x41','6x38','3x44','7x52',1,7.5],
		# 			['2x41','6x38','3x44','7x52',1.5,7.5],
		# 			['2x41','6x38','3x44','7x52',2,7.5],
		# 			['2x42','6x35','3x44','7x53',-0.5,8.5],
		# 			['2x42','6x35','3x44','7x53',0.5,8.5],
		# 			['2x42','6x38','3x47','7x53',0,8.0]]
		# self.only_xtals = options['only_xtals']
		# self.processors = options['proc']
		# self.data = [5, 2.5, 0]
		# self.prepare_input(options['proc'], self.data)

		# data_dir = '../../data/protwis/gpcr/structure_data/structures/'
		# files = os.listdir(data_dir)
		# for f in files:
		# 	with open(data_dir+f, 'r') as yf:
		# 		y = yaml.load(yf)
		# 		s = Structure.objects.get(pdb_code__index=f.split('.')[0])
		# 		if s.protein_conformation.protein.parent.family.parent.parent.parent.slug=='002':
		# 			p = PdbStateIdentifier(s)
		# 			p.run()
		# 			print(f, s.state, s.distance, p.state, p.activation_value)
			

		# t = TestStateIdentifier(options['gns'], options['only_xtals'], float(options['cutoffs'][0]), float(options['cutoffs'][1]))
		# t.run()

		# c = ChangeDistanceValues()
		# c.run()

		# ss = StructuralStatistics()
		# for s in Structure.objects.filter(protein_conformation__protein__family__parent__parent__parent__slug='001'):
		# 	print(s.protein_conformation.protein.parent.entry_name, s, s.state, options['segment'], ss.check_segment(s, options['segment']))

		# d = ['98', '99', '100', '101', '102', '103', '104', '105', '106', '107', '108', '109', '110', '111', '112', '113', '114', '115', '116', '117', '118', '119', '120', '121', '122', '123', '124', '125', '126', '127', '128', '129', '130', '131', '132', '133', '134', '135', '136', '137', '138', '139', '140', '141', '142', '143', '144', '145', '146', '147', '148', '149', '150', '151', '152', '153', '154', '155', '156', '157', '158', '159', '160', '161', '162', '163', '164', '165', '166', '167', '168', '169', '170', '171', '172', '173', '174', '175', '176', '177', '178', '179', '180', '181', '182', '183', '184', '185', '186', '187', '188', '189', '190', '191', '192', '193', '194', '195', '196', '197', '198', '199', '200', '201', '202', '203', '204', '205', '206', '207', '208', '209', '210', '211', '212', '213', '214', '215', '216', '217', '218', '219', '220', '221', '222', '223', '224', '225', '226', '227', '228', '229', '230', '231', '232', '233', '234', '235', '236', '237', '238', '239', '240', '241', '242', '243', '244', '245', '246', '247', '248', '249', '250', '251', '252', '253', '254', '255', '256', '257', '258', '259', '260', '261', '262', '263', '264', '265', '266', '267', '268', '269', '270', '271', '272', '273', '274', '275', '276', '277', '278', '279', '280', '281', '282', '283', '284', '285', '286', '287', '288', '289', '290', '291', '292', '293', '294', '295', '296', '297', '298', '299', '300', '301', '302', '303', '304', '305', '306', '307', '308', '309', '310', '311', '312', '313', '314', '315', '316', '317', '318', '319', '320', '321', '322', '323', '324', '325', '326', '327', '328', '329', '330', '331', '332', '333', '334', '335', '336', '337', '338', '339', '340', '341', '342', '343', '344', '345', '346', '347', '348', '349', '350', '351', '352', '353', '354', '355', '356', '357', '358', '359', '360', '361', '362', '363', '364', '365', '366', '367', '368', '369', '370', '371', '372', '373', '374', '375', '376', '377', '378', '379', '380', '381', '382', '383', '384', '385', '386', '387', '388', '389', '390', '391', '392', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35', '36', '37', '38', '39', '40', '41', '42', '43', '44', '45', '46', '47', '207', '208', '209', '210', '211', '212', '213', '214', '215', '216', '217', '218', '219', '220', '221', '222', '223', '224', '225', '226', '227', '228', '229', '230', '231', '232', '233', '234', '235', '236', '237', '238', '239', '240', '241', '242', '243', '244', '245', '246', '247', '248', '249', '250', '264', '265', '266', '267', '268', '269', '270', '271', '272', '273', '274', '275', '276', '277', '278', '279', '280', '281', '282', '283', '284', '285', '286', '287', '288', '289', '290', '291', '292', '309', '310', '311', '312', '313', '314', '315', '316', '317', '318', '319', '320', '321', '331', '332', '333', '334', '335', '336', '337', '338', '339', '340', '341', '342', '343', '344', '345', '346', '347', '348', '349', '350', '351', '352', '353', '354', '355', '356', '357', '358', '359', '360', '361', '362', '363', '364', '365', '370', '371', '372', '373', '374', '375', '376', '377', '378', '379', '380', '381', '382', '383', '384', '385', '386', '387', '388', '389', '390', '391', '392', '393', '394']
		
		print(StructureComplexModel.objects.defer('pdb').filter(receptor_protein__entry_name='gp139_human')[0])
		return 0

		s = SeqCompare()
		# s.compare('MAGCCCLSAEEKESQRISAEIERQLRRDKKDARRELKLLLLGTGESGKSTFIKQMRIIHGSGYSDEDRKGFTKLVYQNIFTAMQAMIRAMDTLRIQYVCEQNKENAQIIREVEVDKVSMLSREQVEAIKQLWQDPGIQECYDRRREYQLSDSAKYYLTDIDRIATPSFVPTQQDVLRVRVPTTGIIEYPFDLENIIFRMVDVGGQRSERRKWIHCFESVTSIIFLVALSEYDQVLAECDNENRMEESKALFKTIITYPWFLNSSVILFLNKKDLLEEKIMYSHLISYFPEYTGPKQDVRAARDFILKLYQDQNPDKEKVIYSHFTCATDTDNIRFVFAAVKDTILQLNLREFNLV',
		# 		  'MAGCCCLSAEEKESQRISAEIERQLRRDKKDARRELKLLLLGTGESGKSTFIKQMRIIHGSGYSDEDRKGFTKLVYQNIFTAMQAMIRAMDTLKIQYVCEQNKENAQLIREVEVDKVSTLSRDQVEAIKQLWQDPGIQECYDRRREYQLSDSAKYYLTDIDRIAMPAFVPTQQDVLRVRVPTTGIIEYPFDLENIIFRMVDVGGQRSERRKWIHCFESVTSIIFLVALSEYDQVLAECDNENRMEESKALFKTIITYPWFLNSSVILFLNKKDLLEEKIMYSHLISYFPEYTGPKQDVKAARDFILKLYQDQNPDKEKVIYSHFTCATDTENIRFVFAAVKDTILQLNLREFNLV')
		
		if True and True:
			print('true')
		s.align('MG------CTL------------------SAEERAALERSKAIEKNLKEDGISAAKDVKLLLLGAGESGKSTIVKQMKIIHEDGFS---------------GEDVKQYKPVVYSNTIQSLAAIVRAMDTLG--IEYGDKERKADAKMVCDVVSRMEDT------EPFSAELLSAMMRLWGDSGIQECFNRSREYQLNDSAKYYLDSLDRIGAADYQPTEQDILRTRVKTTGIVETHFTFKNLHFRLFDVGGQRSERKKWIHCFEDVTAIIFCVALSGYDQVLHEDETTNRMHESLMLFDSICNNKFFIDTSIILFLNKKDLFGEKIK--KSPLTICFPEYTGP-------------NTYEDAAA-YIQAQFESKNR----SPNKE--------IYCHMTCATDTNNIQVVFDAVTDIIIANNLRGCGLY--------------------------',
				'MGCAMSAEERAALARSRQIERNLREDGLQAAKDIKLLLLGAGESGKSTIVKQMKIIHESGFTSEDFKQYRPVVFSNTVQSLVAILRAMPNLGIGFGTNERETDAKMVLDVIQRMEDTEPFSEELLTAMKRLWADPGVQMCFSRSNEYQLNDSAKYFLDDLDRLGSKDYQPTEQDILRTRVKTTGIVEVHFSFKNLNFKLFDVGGQRSERKKWIHCFEDVTAIIFCVAMSEYDQVLHEDETTNRMQESLKLFDSICNNKWFTDTSIILFLNKKDLFEEKIKKSPLTICFPEYAGAQEYGEAAAYIQAQFEAKNKSTTKEIYCHMTCATDTNNIQFVFDAVTDVIIANNLRGCGLY')

		print(datetime.now()-startTime)

	def main_func(self, positions, iteration, count, lock):
		processor_id = round(self.processors*positions[0]/len(self.data))+1
		i = 0
		while count.value<len(self.data):
			i += 1
			with lock:
				d = self.data[count.value]
				count.value +=1
			t = TestStateIdentifierSets(self.only_xtals, d)
			t.run()
			# t = TestStateIdentifier([d[0],d[1],d[2],d[3]],self.only_xtals,d[4],d[5])
			# t.run()


class SeqCompare(object):
	def __init__(self):
		pass

	def compare(self, seq1, seq2):
		if len(seq1)==len(seq2):
			for i, j, k in zip(enumerate(seq1), seq1, seq2):
				if j!=k:
					print('Mismatch at',i)
					print(seq1[i[0]-3:i[0]+1])
					print(seq2[i[0]-3:i[0]+1], 'change here')

	def align(self, seq1, seq2):
		p = pairwise2.align.globalms(seq1, seq2, 1, 1, -18, -18)
		for i in p[0]:
			print(i)


class StructuralStatistics(object):
	def __init__(self):
		pass

	def check_segment(self, structure, segment):
		parent_res = Residue.objects.filter(protein_conformation__protein=structure.protein_conformation.protein.parent, protein_segment__slug=segment)
		parent_gns = [i for i in parent_res if i.display_generic_number!=None]

		struct_res = Residue.objects.filter(protein_conformation=structure.protein_conformation, protein_segment__slug=segment)
		struct_gns = [i for i in struct_res if i.display_generic_number!=None]

		if len(parent_gns)>0:
			if len(struct_gns)==len(parent_gns):
				return 'Conserved_ordered'
			elif len(struct_gns)==0:
				if len(parent_res)-2<=len(struct_res)<=len(parent_res)+2:
					return 'Conserved_disordered'
				elif len(struct_res)>len(parent_res):
					return 'Conserved_disordered'
				elif len(struct_res)<=len(parent_res)/2:
					return 'Conserved_missing'
				elif len(struct_res)<len(parent_res):
					return 'Conserved_missing'
			else:
				if len(struct_res)<=len(parent_res)/2:
					return 'Conserved_missing'
		else:
			if len(struct_res)<=len(parent_res)/2:
				return 'Non-conserved_missing'
			else:
				return 'Non-conserved_present'


class ChangeDistanceValues(object):
	def __init__(self):
		self.path = os.sep.join([settings.DATA_DIR, 'structure_data', 'structures'])

	def run(self):
		files = os.listdir(self.path)
		for f in files:
			try:
				with open(os.sep.join([self.path, f]), 'r') as yf:
					y = yaml.load(yf)
					ps = PdbStateIdentifier(Structure.objects.get(pdb_code__index=y['pdb']))
					ps.run()
					print(f, y['distance'], y['state'], round(float(ps.activation_value), 2), ps.state)
					y['distance'] = round(float(ps.activation_value), 2)
				with open(os.sep.join([self.path, f]), 'w') as syf:
					yaml.dump(y, syf, indent=4, default_flow_style=False)
			except Exception as msg:
				print(f, msg)


class TestStateIdentifierSets(object):
	def __init__(self, only_xtals=False, iac=2):
		self.only_xtals = only_xtals
		self.iac = iac

	def run(self):
		tm2 = ['2x39','2x40','2x41','2x42']
		tm6 = ['6x32','6x33','6x34','6x35'] # for class A ['6x35','6x36','6x37','6x38']
		tm3 = ['3x47','3x46','3x45','3x44']
		tm7 = ['7x53','7x52','7x51','7x50']
		inact_cutoffs = [1, -0.5, -2] # for class A [1, -0.5, -2]
		inter_cutoffs = [15, 12, 9, 6] # for class A [8, 6, 4] 
		best = 1000
		best_params = []
		counter = 0
		# with open('./structure/cutoff_test_01.csv', 'w') as f:
			# for iac in inact_cutoffs:
		iac = self.iac
		for inc in inter_cutoffs:
			for t2 in tm2:
				for t6 in tm6:
					for t3 in tm3:
						for t7 in tm7:
							counter+=1
							t = TestStateIdentifier([t2, t6, t3, t7], self.only_xtals, iac, inc)
							t.run()
							if t.mismatch<best:
								best = t.mismatch
								best_params = [t2, t6, t3, t7, iac, inc]
							if t.mismatch==0:
								print(counter, t2,t6,t3,t7, t.mismatch, t.match, t.exceptions, iac, inc)
							# f.write('{},{},{},{},{},{},{},{},{}\n'.format(t2, t6, t3, t7, t.mismatch, t.match, t.exceptions, iac, inc))
		print(best_params, best)


class TestStateIdentifierBestSets(TestStateIdentifierSets):
	def run(self, data):
		cutoff_finetune = [-1,-0.5,0,0.5,1]
		for plus_iac in cutoff_finetune:
			for plus_inc in cutoff_finetune:
				t = TestStateIdentifier([data[0],data[1],data[2],data[3]], self.only_xtals, data[4]+plus_iac, data[5]+plus_inc)
				t.run()
				print('{}-{}-{}-{},{},{},{},{},{}'.format(data[0],data[1],data[2],data[3],t.mismatch,t.match,t.exceptions,data[4]+plus_iac,data[5]+plus_inc))


class TestStateIdentifier(object):
	def __init__(self, gns, only_xtals=False, inact_cutoff=-1, inter_cutoff=8):
		self.tm2_gn, self.tm6_gn, self.tm3_gn, self.tm7_gn = '2x41', '6x33', '3x44', '7x51'
		self.inact_cutoff = inact_cutoff
		self.inter_cutoff = inter_cutoff
		if gns:
			for value in gns:
				if value.startswith('2'):
					self.tm2_gn = value
				elif value.startswith('6'):
					self.tm6_gn = value
				elif value.startswith('3'):
					self.tm3_gn = value
				elif value.startswith('7'):
					self.tm7_gn = value
		self.only_xtals = only_xtals

	def run(self):
		strs = Structure.objects.filter(refined=False).exclude(protein_conformation__protein__parent__family__parent__parent__parent__slug__in=['001','004','005'])
		self.match, self.mismatch, self.exceptions = 0,0,0
		for s in strs:
			try:
				if self.only_xtals:
					psis = PdbStateIdentifier(s, self.tm2_gn, self.tm6_gn, self.tm3_gn, self.tm7_gn, self.inact_cutoff, self.inter_cutoff)
					psis.run()
					if psis.state!=s.state:
						print(s, s.state, s.distance, psis.state, psis.activation_value, 'mismatch')
						self.mismatch+=1
					else:
						print(s,",",s.state,',', psis.state, psis.activation_value, self.inact_cutoff, self.inter_cutoff)
						self.match+=1
				else:
					r_s = Structure.objects.get(pdb_code__index=s.pdb_code.index+'_refined')
					psi = PdbStateIdentifier(r_s, self.tm2_gn, self.tm6_gn, self.tm3_gn, self.tm7_gn, self.inact_cutoff, self.inter_cutoff)
					psi.run()
					psis = PdbStateIdentifier(s, self.tm2_gn, self.tm6_gn, self.tm3_gn, self.tm7_gn, self.inact_cutoff, self.inter_cutoff)
					psis.run()
					if psi.state!=psis.state:
						print(s, psis.state.slug, psis.activation_value, psi.state.slug, psi.activation_value)
						self.mismatch+=1
					else:
						self.match+=1
			except:
				print('Exception: ', s)
				self.exceptions+=1
		if not self.only_xtals:
			hommods = StructureModel.objects.all().exclude(protein__family__parent__parent__parent__slug__in=['001','004','005'])
			for h in hommods:
				try:
					psih = PdbStateIdentifier(h, self.tm2_gn, self.tm6_gn, self.tm3_gn, self.tm7_gn, self.inact_cutoff, self.inter_cutoff)
					psih.run()
					psiss = PdbStateIdentifier(h.main_template, self.tm2_gn, self.tm6_gn, self.tm3_gn, self.tm7_gn, self.inact_cutoff, self.inter_cutoff)
					psiss.run()
					if psih.state!=psiss.state:
						print(h, psiss.state.slug, psiss.activation_value, psih.state.slug, psih.activation_value)
						self.mismatch+=1
					else:
						self.match+=1
				except:
					print('Exception hommod:', h)
					self.exceptions+=1
		print('match:', self.match, 'mismatch:', self.mismatch, 'exceptions:', self.exceptions)
		print('{}-{}-{}-{},{},{},{},{},{}'.format(self.tm2_gn, self.tm6_gn, self.tm3_gn, self.tm7_gn, self.mismatch, self.match, self.exceptions, self.inact_cutoff, self.inter_cutoff))
		return 0