# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 14:31:44 CET 2020

@author: Gaspar Pandy
"""
from build.management.commands.base_build import Command as BaseBuild
from django.db.models import Q
from django.conf import settings

from protein.models import Protein, ProteinConformation, ProteinSequenceType, ProteinSource, ProteinState
from residue.models import Residue
from structure.models import Structure, PdbData, StructureType
from structure.sequence_parser import SequenceParser
from structure.functions import PdbChainSelector, PdbStateIdentifier
from construct.functions import *
from common.models import WebResource, WebLink, Publication
from ligand.models import Ligand, LigandProperities, LigandType

import Bio.PDB as PDB
from datetime import datetime
import urllib
import re
import os
import xmltodict
import yaml
import shlex
import subprocess
import csv


startTime = datetime.now()

class Command(BaseBuild):

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
        el = EndogenousLigands(options['verbose'])
        if options['purge']:
            el.purge()
        el.parse_data_file()
        el.build_endogenous_ligands()


class EndogenousLigands():
    def __init__(self, verbose):
        self.data_file = None
        self.verbose = verbose

    def purge(self):
        ligands = Ligand.objects.filter(endogenous=True)
        for l in ligands:
            # for w in l.properities.web_links.all():
            #     w.delete()
            l.properities.delete()
        ligands.delete()

    def parse_data_file(self):
        self.data_file = os.sep.join([settings.DATA_DIR, 'ligand_data', '191107_endogenous_ligands.csv'])
        with open(self.data_file, newline='') as csvfile:
            reader = csv.reader(csvfile, delimiter=',')
            self.data = [i for i in reader]

    def build_endogenous_ligands(self):
        if self.verbose:
            print('Running build_endogenous_ligands...')
        for i, val in enumerate(self.data):
            if i!=0:
                protein = val[2]
                ligand_type = val[14]
                seq = val[23]
                name = val[12]
                if val[16]=='':
                    pubchemid = None
                else:
                    pubchemid = int(float(val[16]))
                try:
                    db_prot = Protein.objects.get(entry_name=protein)
                except Protein.DoesNotExist:
                    print('Missing {} from Protein table'.format(protein))
                    continue
                if ligand_type=='Peptide':
                    if len(seq)>50:
                        ligand_type, created = LigandType.objects.get_or_create(slug='protein', name='protein')
                    else:
                        ligand_type = LigandType.objects.get(slug='peptide')
                elif ligand_type=='Constitutive':
                    ligand_type, created = LigandType.objects.get_or_create(slug='constitutive', name='constitutive')
                else:
                    if seq=='':
                        ligand_type = LigandType.objects.get(slug='small-molecule')
                    else:
                        print('Error: unexpected ligand type: ',protein, db_prot, ligand_type, seq, name, pubchemid)
                if seq=='':
                    sequence = None
                else:
                    sequence = seq
                if not pubchemid and len(WebLink.objects.filter(index=pubchemid))==0:
                    created = True
                else:
                    wl, created = WebLink.objects.get_or_create(index=pubchemid, web_resource=WebResource.objects.get(slug='pubchem'))
                if created or (len(wl.ligandproperities_set.all())==0 and not created):
                    lp = LigandProperities(ligand_type=ligand_type, sequence=sequence)
                    if name=='succinic acid':
                        sa = LigandProperities.objects.get(inchikey='KDYFGRWQOYBRFD-UHFFFAOYSA-N')
                        sa.ligand_type = LigandType.objects.get(slug='small-molecule')
                        sa.save()
                    lp.save()
                    if pubchemid:
                        lp.web_links.add(wl)
                        lp.save()
                else:
                    lp = wl.ligandproperities_set.all()[0]
                try:
                    l = Ligand.objects.get(name=name)
                    l.endogenous = True
                except:
                    l = Ligand.objects.create(name=name, endogenous=True, properities=lp)
                l.properities = lp
                l.name = name
                # if self.verbose:
                #     print(l, lp, l.name)
                if pubchemid:
                    l.load_from_pubchem('cid', pubchemid, lp.ligand_type, ligand_title=l.name)
                l.save()
                orthologs = Protein.objects.filter(family=db_prot.family, accession__isnull=False)
                for o in orthologs:
                    o.endogenous_ligands.add(l)
                # if self.verbose:
                #     print(db_prot, db_prot.endogenous_ligands.all())
        if self.verbose:
            print('Runtime: ',datetime.now() - startTime)
