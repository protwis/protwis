import csv
import logging
import math
import os
import shlex
import subprocess
import sys
from collections import OrderedDict, defaultdict
from itertools import islice
from optparse import make_option
from urllib.request import urlopen
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
from pprint import pprint
import requests
import xlrd
import xmltodict, json
import yaml
import django.apps
import Bio.PDB as PDB
from Bio import SeqIO, pairwise2
from Bio.pairwise2 import format_alignment
from common.models import Publication, WebLink, WebResource
from django.conf import settings
from django.utils.text import slugify
from django.core.management.base import BaseCommand, CommandError
from django.core.management.color import no_style
from django.db import IntegrityError, connection
from protein.models import (Gene, Protein, ProteinAlias, ProteinConformation,
                            ProteinFamily, Species, ProteinSegment,
                            ProteinSequenceType, ProteinSource, ProteinState)
from residue.models import (Residue, ResidueGenericNumber,
                            ResidueGenericNumberEquivalent,
                            ResidueNumberingScheme)
from signprot.models import SignprotBarcode, SignprotComplex, SignprotStructure, SignprotStructureExtraProteins
from structure.models import Structure, StructureStabilizingAgent, StructureType, StructureExtraProteins
from common.tools import test_model_updates


class Command(BaseCommand):
    help = 'Build G proteins'

    # source files
    gprotein_data_path = os.sep.join([settings.DATA_DIR, 'g_protein_data'])
    if not os.path.exists(os.sep.join([settings.DATA_DIR, 'g_protein_data', 'PDB_UNIPROT_ENSEMBLE_ALL.txt'])):
        with open(os.sep.join([settings.DATA_DIR, 'g_protein_data', 'PDB_UNIPROT_ENSEMBLE_ALL.txt']), 'w') as f:
            f.write(
                'PDB_ID\tPDB_Chain\tPosition\tResidue\tCGN\tEnsembl_Protein_ID\tUniprot_ACC\tUniprot_ID\tsortColumn\n')
    gprotein_data_file = os.sep.join([settings.DATA_DIR, 'g_protein_data', 'PDB_UNIPROT_ENSEMBLE_ALL.txt'])
    barcode_data_file = os.sep.join([settings.DATA_DIR, 'g_protein_data', 'barcode_data.csv'])
    pdbs_path = os.sep.join([settings.DATA_DIR, 'g_protein_data', 'pdbs'])
    lookup = os.sep.join([settings.DATA_DIR, 'g_protein_data', 'CGN_lookup.csv'])
    alignment_file = os.sep.join([settings.DATA_DIR, 'g_protein_data', 'CGN_referenceAlignment.fasta'])
    ortholog_file = os.sep.join([settings.DATA_DIR, 'g_protein_data', 'gprotein_orthologs.csv'])
    local_uniprot_dir = os.sep.join([settings.DATA_DIR, 'g_protein_data', 'uniprot'])
    local_uniprot_beta_dir = os.sep.join([settings.DATA_DIR, 'g_protein_data', 'uniprot_beta'])
    local_uniprot_gamma_dir = os.sep.join([settings.DATA_DIR, 'g_protein_data', 'uniprot_gamma'])
    remote_uniprot_dir = 'https://www.uniprot.org/uniprot/'
    #Setting the variables for the test tracking of the model upadates
    tracker = {}
    all_models = django.apps.apps.get_models()[6:]
    test_model_updates(all_models, tracker, initialize=True)

    logger = logging.getLogger(__name__)

    def add_arguments(self, parser):
        parser.add_argument('--filename',
                            action='append',
                            dest='filename',
                            help='Filename to import. Can be used multiple times')
        parser.add_argument('--wt',
                            default=False,
                            type=str,
                            help='Add wild type protein sequence to data')
        parser.add_argument('--xtal',
                            default=False,
                            type=str,
                            help='Add xtal to data')
        parser.add_argument('--build_datafile',
                            default=False,
                            action='store_true',
                            help='Build PDB_UNIPROT_ENSEMBLE_ALL file')

    def handle(self, *args, **options):
        self.options = options
        if options['filename']:
            filenames = options['filename']
        else:
            filenames = False
        if self.options['wt']:
            self.add_entry()
            test_model_updates(self.all_models, self.tracker, check=True)
        elif self.options['build_datafile']:
            self.build_table_from_fasta()
            test_model_updates(self.all_models, self.tracker, check=True)
        else:
            # Add G-proteins from CGN-db Common G-alpha Numbering <https://www.mrc-lmb.cam.ac.uk/CGN/>
            try:
                self.purge_signprot_complex_data()
                self.logger.info('PASS: purge_signprot_complex_data')
                self.purge_cgn_proteins()
                self.logger.info('PASS: purge_cgn_proteins')
                self.purge_other_subunit_proteins()
                self.logger.info('PASS: purge_other_subunit_proteins')

                self.ortholog_mapping = OrderedDict()
                with open(self.ortholog_file, 'r') as ortholog_file:
                    ortholog_data = csv.reader(ortholog_file, delimiter=',')
                    for i, row in enumerate(ortholog_data):
                        if i == 0:
                            header = list(row)
                            continue
                        for j, column in enumerate(row):
                            if j in [0, 1]:
                                continue
                            if '_' in column:
                                self.ortholog_mapping[column] = row[0]
                            else:
                                if column == '':
                                    continue
                                self.ortholog_mapping[column + '_' + header[j]] = row[0]
                self.logger.info('PASS: ortholog_mapping')

                self.cgn_create_proteins_and_families()
                self.logger.info('PASS: cgn_create_proteins_and_families')

                human_and_orths = self.cgn_add_proteins()
                self.logger.info('PASS: cgn_add_proteins')
                self.update_protein_conformation(human_and_orths)
                self.logger.info('PASS: update_protein_conformation')
                self.create_barcode()
                self.logger.info('PASS: create_barcode')
                self.add_other_subunits()
                self.logger.info('PASS: add_other_subunits')
                test_model_updates(self.all_models, self.tracker, check=True)
            except Exception as msg:
                exc_type, exc_obj, exc_tb = sys.exc_info()
                fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                print(exc_type, fname, exc_tb.tb_lineno)
                self.logger.error(msg)

    def add_other_subunits(self):
        beta_fam, created = ProteinFamily.objects.get_or_create(slug='100_002', name='Beta',
                                                                parent=ProteinFamily.objects.get(name='G-Protein'))
        gigsgt, created = ProteinFamily.objects.get_or_create(slug='100_002_001', name='G(I)/G(S)/G(T)',
                                                              parent=beta_fam)
        gamma_fam, created = ProteinFamily.objects.get_or_create(slug='100_003', name='Gamma',
                                                                 parent=ProteinFamily.objects.get(name='G-Protein'))
        gigsgo, created = ProteinFamily.objects.get_or_create(slug='100_003_001', name='G(I)/G(S)/G(O)',
                                                              parent=gamma_fam)

        # create proteins
        self.create_beta_gamma_proteins(self.local_uniprot_beta_dir, gigsgt)
        self.create_beta_gamma_proteins(self.local_uniprot_gamma_dir, gigsgo)
        # create residues
        self.create_beta_gamma_residues(gigsgt)
        self.create_beta_gamma_residues(gigsgo)

    def create_beta_gamma_residues(self, proteinfamily):
        bulk = []
        fam = Protein.objects.filter(family=proteinfamily)
        for prot in fam:
            for i, j in enumerate(prot.sequence):
                prot_conf = ProteinConformation.objects.get(protein=prot)
                r = Residue(sequence_number=i + 1, amino_acid=j, display_generic_number=None, generic_number=None,
                            protein_conformation=prot_conf, protein_segment=None)
                bulk.append(r)
                if len(bulk) % 10000 == 0:
                    self.logger.info('Inserted bulk {} (Index:{})'.format(len(bulk), i))
                    # print(len(bulk),"inserts!",index)
                    Residue.objects.bulk_create(bulk)
                    # print('inserted!')
                    bulk = []
        Residue.objects.bulk_create(bulk)

    def create_beta_gamma_proteins(self, uniprot_dir, proteinfamily):
        files = os.listdir(uniprot_dir)
        for f in files:
            acc = f.split('.')[0]
            up = self.parse_uniprot_file(acc)
            pst = ProteinSequenceType.objects.get(slug='wt')
            try:
                species, created = Species.objects.get_or_create(latin_name=up['species_latin_name'],
                                                                 defaults={
                                                                     'common_name': up['species_common_name'],
                                                                 })
                if created:
                    self.logger.info('Created species ' + species.latin_name)
            except IntegrityError:
                species = Species.objects.get(latin_name=up['species_latin_name'])
            source = ProteinSource.objects.get(name=up['source'])
            try:
                name = up['names'][0].split('Guanine nucleotide-binding protein ')[1]
            except:
                name = up['names'][0]
            prot, created = Protein.objects.get_or_create(entry_name=up['entry_name'], accession=acc, name=name,
                                                          sequence=up['sequence'], family=proteinfamily, parent=None,
                                                          residue_numbering_scheme=None, sequence_type=pst,
                                                          source=source, species=species)
            state = ProteinState.objects.get(slug='active')
            prot_conf, created = ProteinConformation.objects.get_or_create(protein=prot, state=state)

    def fetch_missing_uniprot_files(self):
        """
        This function is not being used at the moment 23-04-2020
        TODO: Check if this function is needed.
        :return:
        """
        BASE = 'https://www.uniprot.org'
        KB_ENDPOINT = '/uniprot/'
        uniprot_files = os.listdir(self.local_uniprot_dir)
        new_uniprot_files = os.listdir(os.sep.join([settings.DATA_DIR, 'g_protein_data', 'uniprot']))
        for record in SeqIO.parse(self.alignment_file, 'fasta'):
            sp, accession, name, ens = record.id.split('|')
            if accession + '.txt' not in uniprot_files and accession + '.txt' not in new_uniprot_files:
                g_prot, species = name.split('_')
                result = requests.get(BASE + KB_ENDPOINT + accession + '.txt')
                with open(os.sep.join([settings.DATA_DIR, 'g_protein_data', 'uniprot', accession + '.txt']), 'w') as f:
                    f.write(result.text)
            else:
                try:
                    os.rename(os.sep.join([self.local_uniprot_dir, accession + '.txt']),
                              os.sep.join([settings.DATA_DIR, 'g_protein_data', 'uniprot', accession + '.txt']))
                except:
                    print('Missing: {}'.format(accession))

    def update_alignment(self):
        with open(self.lookup, 'r') as csvfile:
            lookup_data = csv.reader(csvfile, delimiter=',', quotechar='"')
            lookup_dict = OrderedDict([('-1', 'NA.N-terminal insertion.-1'), ('-2', 'NA.N-terminal insertion.-2'),
                                       ('-3', 'NA.N-terminal insertion.-3')])
            for row in lookup_data:
                lookup_dict[row[0]] = row[1:]

        residue_data = pd.read_csv(self.gprotein_data_file, sep="\t", low_memory=False)

        for i, row in residue_data.iterrows():
            try:
                residue_data['CGN'][i] = lookup_dict[str(int(residue_data['sortColumn'][i]))][6].replace('(',
                                                                                                         '').replace(
                    ')', '')
            except:
                pass
        residue_data['sortColumn'] = residue_data['sortColumn'].astype(int)

        residue_data.to_csv(path_or_buf=self.gprotein_data_path + '/test.txt', sep='\t', na_rep='NA', index=False)

    def build_lookup_file(self, segment_name=None, segment_value=None):
        segment_lengths = OrderedDict([('HN',53),('hns1',3),('S1',7),('s1h1',6),('H1',12),('h1ha',20),('HA',29),
                                        ('hahb',78),('HB',14),('hbhc',55),('HC',12),('hchd',2),('HD',12),('hdhe',5),
                                        ('HE',13),('hehf',7),('HF',6),('hfs2',7),('S2',8),('s2s3',2),('S3',8),('s3h2',3),
                                        ('H2',10),('h2s4',5),('S4',7),('s4h3',15),('H3',18),('h3s5',3),('S5',7),('s5hg',1),
                                        ('HG',17),('hgh4',21),('H4',17),('h4s6',20),('S6',5),('s6h5',5),('H5',26)])
        if segment_name and segment_value:
            segment_lengths[segment_name] = segment_value
        alt = {'S1':'S1', 's1h1':'P-loop', 'hfs2':'SwI', 's3h2':'SwII', 's4h3':'SwIII', 's6h5':'TCAT'}
        with open(os.sep.join([self.gprotein_data_path, "CGN_lookup_test1.csv"]), 'w') as f:
            f.write('"","Domain","GalphaSSE","CGN_pos","GalphaSSE_alternative","SSE_pos","CGN","CGN_new"\n')
            c = 1
            for s, l in segment_lengths.items():
                if s in ['HA','hahb','HB','hbhc','HC','hchd','HD','hdhe','HE','hehf','HF']:
                    domain = 'H'
                else:
                    domain = 'G'
                if s in alt:
                    alt_name = alt[s]
                else:
                    alt_name = "NA"
                if len(s)==4:
                    alt_s = s[:2]+'-'+s[2:]
                else:
                    alt_s = s
                for i in range(1, l+1):
                    line = '"{}","{}","{}",{},{},{},"({}).{}.{}","({}).{}.{}"\n'.format(c, domain, s, c, alt_name, i, domain, alt_s, c, domain, s, i)
                    f.write(line)
                    c+=1

    def add_new_orthologs(self):
        residue_data = pd.read_csv(self.gprotein_data_file, sep="\t", low_memory=False)
        with open(self.lookup, 'r') as csvfile:
            lookup_data = csv.reader(csvfile, delimiter=',', quotechar='"')
            lookup_dict = OrderedDict([('-1', 'NA.N-terminal insertion.-1'), ('-2', 'NA.N-terminal insertion.-2'),
                                       ('-3', 'NA.N-terminal insertion.-3')])
            for row in lookup_data:
                lookup_dict[row[0]] = row[1:]
        fasta_dict = OrderedDict()
        for record in SeqIO.parse(self.alignment_file, 'fasta'):
            sp, accession, name, ens = record.id.split('|')
            g_prot, species = name.split('_')
            if g_prot not in fasta_dict:
                fasta_dict[g_prot] = OrderedDict([(name, [accession, ens, str(record.seq)])])
            else:
                fasta_dict[g_prot][name] = [accession, ens, str(record.seq)]
        with open(self.ortholog_file, 'r') as ortholog_file:
            ortholog_data = csv.reader(ortholog_file, delimiter=',')
            for i, row in enumerate(ortholog_data):
                if i == 0:
                    header = list(row)
                    continue
                for j, column in enumerate(row):
                    if j in [0, 1]:
                        continue
                    if column != '':
                        if '_' in column:
                            entry_name = column
                            BASE = 'https://www.uniprot.org'
                            KB_ENDPOINT = '/uniprot/'
                            result = requests.get(BASE + KB_ENDPOINT,
                                                  params={'query': 'mnemonic:{}'.format(entry_name), 'format': 'list'})
                            accession = result.text.replace('\n', '')
                        else:
                            entry_name = '{}_{}'.format(column, header[j])
                            accession = column
                        if entry_name not in fasta_dict[row[0]]:
                            result = requests.get('https://www.uniprot.org/uniprot/{}.xml'.format(accession))
                            uniprot_entry = result.text
                            try:
                                entry_dict = xmltodict.parse(uniprot_entry)
                            except:
                                self.logger.warning('Skipped: {}'.format(accession))
                                continue
                            try:
                                ensembl = \
                                [i for i in entry_dict['uniprot']['entry']['dbReference'] if i['@type'] == 'Ensembl'][
                                    0]['@id']
                            except:
                                self.logger.warning('Missing Ensembl: {}'.format(accession))
                            sequence = entry_dict['uniprot']['entry']['sequence']['#text'].replace('\n', '')
                            seqc = SeqCompare()
                            aligned_seq = seqc.align(fasta_dict[row[0]][row[0] + '_HUMAN'][2], sequence)
                            fasta_dict[row[0]][entry_name] = [accession, ensembl, aligned_seq]

        with open(os.sep.join([settings.DATA_DIR, 'g_protein_data', 'fasta_test.fa']), 'w') as f:
            for g, val in fasta_dict.items():
                for i, j in val.items():
                    f.write('>sp|{}|{}|{}\n{}\n'.format(j[0], i, j[1], j[2]))

    def add_entry(self, pdbs_path=None):
        if not self.options['wt']:
            raise AssertionError('Error: Missing wt name')
        residue_data = pd.read_csv(self.gprotein_data_file, sep="\t", low_memory=False)
        try:
            if residue_data['Uniprot_ID'][self.options['wt']] and not self.options['xtal']:
                return 0
        except:
            pass
        with open(self.lookup, 'r') as csvfile:
            lookup_data = csv.reader(csvfile, delimiter=',', quotechar='"')
            lookup_dict = OrderedDict([('-1', 'NA.N-terminal insertion.-1'), ('-2', 'NA.N-terminal insertion.-2'),
                                       ('-3', 'NA.N-terminal insertion.-3')])
            for row in lookup_data:
                lookup_dict[row[0]] = row[1:]
        sequence = ''
        for record in SeqIO.parse(self.alignment_file, 'fasta'):
            sp, accession, name, ens = record.id.split('|')
            if name == self.options['wt'].upper():
                sequence = record.seq
                break

        with open(self.gprotein_data_file, 'a') as f:
            count_sort, count_res = 0, 0
            if self.options['xtal']:
                m = PDB.PDBParser(os.sep.join('st', [pdbs_path, self.options['xtal'] + '.pdb']))[0]

            else:
                for i in sequence:
                    count_sort += 1
                    if i == '-':
                        continue
                    count_res += 1
                    try:
                        cgn = lookup_dict[str(count_sort)][6].replace('(', '').replace(')', '')
                    except:
                        print('Dict error: {}'.format(self.options['wt']))
                    line = 'NA\tNA\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(count_res, i, cgn, ens, accession, name,
                                                                         count_sort)
                    f.write(line)

    def build_table_from_fasta(self):
        os.chdir('/vagrant/shared/sites/protwis')
        for record in SeqIO.parse(self.alignment_file, 'fasta'):
            sp, accession, name, ens = record.id.split('|')
            if len(record.seq) != 456:
                continue
            command = "/env/bin/python3 manage.py build_g_proteins --wt " + str(name.lower())
            subprocess.call(shlex.split(command))

    def purge_cgn_residues(self):
        try:
            Residue.objects.filter(generic_number_id__scheme__slug="cgn").delete()
        except:
            self.logger.warning('Existing Residue data cannot be deleted')

    def purge_other_subunit_residues(self):
        try:
            Residue.objects.filter(protein_conformation__protein__family__parent__name="Beta").delete()
        except:
            self.logger.warning('Existing Residue data cannot be deleted')
        try:
            Residue.objects.filter(protein_conformation__protein__family__parent__name="Gamma").delete()
        except:
            self.logger.warning('Existing Residue data cannot be deleted')

    def purge_signprot_complex_data(self):
        try:
            SignprotComplex.objects.all().delete()
        except:
            self.logger.warning('SignprotComplex data cannot be deleted')

    def create_barcode(self):

        barcode_data = pd.read_csv(self.barcode_data_file, low_memory=False)

        for index, entry in enumerate(barcode_data.iterrows()):

            similarity = barcode_data[index:index + 1]['aln_seqSim'].values[0]
            identity = barcode_data[index:index + 1]['aln_seqIdn'].values[0]
            entry_name = barcode_data[index:index + 1]['subfamily'].values[0].lower() + '_human'
            CGN = barcode_data[index:index + 1]['CGN'].values[0]
            paralog = barcode_data[index:index + 1]['paralog'].values[0]

            try:
                p = Protein.objects.get(entry_name=entry_name)
            except Protein.DoesNotExist:
                self.logger.warning('Protein not found for entry_name {}'.format(entry_name))
                continue

            try:
                cgn = Residue.objects.get(protein_conformation__protein=p, display_generic_number__label=CGN)
            except:
                # self.logger.warning('No residue number (GAP - position) for', CGN, "in ", p.name, "")
                continue

            if cgn:
                try:
                    barcode, created = SignprotBarcode.objects.get_or_create(protein=p, residue=cgn,
                                                                             seq_similarity=similarity,
                                                                             seq_identity=identity,
                                                                             paralog_score=paralog)
                    if created:
                        self.logger.info('Created barcode for ' + CGN + ' for protein ' + p.name)
                except IntegrityError:
                    self.logger.error('Failed creating barcode for ' + CGN + ' for protein ' + p.name)

    def purge_cgn_proteins(self):
        try:
            Protein.objects.filter(residue_numbering_scheme__slug='cgn').delete()
            ProteinAlias.objects.filter(protein__family__slug__startswith='100').delete()
        except:
            self.logger.info('Protein to delete not found')

    def purge_other_subunit_proteins(self):
        try:
            Protein.objects.filter(residue_numbering_scheme=None).delete()
        except:
            self.logger.info('Protein to delete not found')

    def add_cgn_residues(self, gprotein_list):
        # Parsing pdb uniprot file for residues
        self.logger.info('Start parsing PDB_UNIPROT_ENSEMBLE_ALL')
        self.logger.info('Parsing file ' + self.gprotein_data_file)
        residue_data = pd.read_csv(self.gprotein_data_file, sep="\t", low_memory=False)
        residue_data = residue_data.loc[residue_data['Uniprot_ACC'].isin(gprotein_list)]
        cgn_scheme = ResidueNumberingScheme.objects.get(slug='cgn')

        # Temp files to speed things up
        temp = {}
        temp['proteins'] = {}
        temp['rgn'] = {}
        temp['segment'] = {}
        temp['equivalent'] = {}
        bulk = []

        self.logger.info('Insert residues: {} rows'.format(len(residue_data)))
        for index, row in residue_data.iterrows():

            if row['Uniprot_ACC'] in temp['proteins']:
                pr = temp['proteins'][row['Uniprot_ACC']][0]
                pc = temp['proteins'][row['Uniprot_ACC']][1]
            else:
                # fetch protein for protein conformation
                pr, c = Protein.objects.get_or_create(accession=row['Uniprot_ACC'])

                # fetch protein conformation
                pc, c = ProteinConformation.objects.get_or_create(protein_id=pr)
                temp['proteins'][row['Uniprot_ACC']] = [pr, pc]

            # fetch residue generic number
            rgnsp = []

            if (int(row['CGN'].split('.')[2]) < 10):
                rgnsp = row['CGN'].split('.')
                rgn_new = rgnsp[0] + '.' + rgnsp[1] + '.0' + rgnsp[2]

                if rgn_new in temp['rgn']:
                    rgn = temp['rgn'][rgn_new]
                else:
                    rgn, c = ResidueGenericNumber.objects.get_or_create(label=rgn_new)
                    temp['rgn'][rgn_new] = rgn

            else:

                if row['CGN'] in temp['rgn']:
                    rgn = temp['rgn'][row['CGN']]
                else:
                    rgn, c = ResidueGenericNumber.objects.get_or_create(label=row['CGN'])
                    temp['rgn'][row['CGN']] = rgn

            d, s, n = row['CGN'].split('.')
            seg = d+'.'+s

            # fetch protein segment id
            if seg in temp['segment']:
                ps = temp['segment'][seg]
            else:
                ps = ProteinSegment.objects.get(slug=seg, proteinfamily='Alpha')
                temp['segment'][seg] = ps

            try:
                bulk_r = Residue(sequence_number=row['Position'], protein_conformation=pc, amino_acid=row['Residue'],
                                 generic_number=rgn, display_generic_number=rgn, protein_segment=ps)
                # self.logger.info("Residues added to db")
                bulk.append(bulk_r)
            except:
                self.logger.error("Failed to add residues")
            if len(bulk) % 10000 == 0:
                self.logger.info('Inserted bulk {} (Index:{})'.format(len(bulk), index))
                # print(len(bulk),"inserts!",index)
                Residue.objects.bulk_create(bulk)
                # print('inserted!')
                bulk = []

            # Add also to the ResidueGenericNumberEquivalent table needed for single residue selection
            try:
                if rgn.label not in temp['equivalent']:
                    ResidueGenericNumberEquivalent.objects.get_or_create(label=rgn.label, default_generic_number=rgn,
                                                                         scheme=cgn_scheme)
                    temp['equivalent'][rgn.label] = 1
                # self.logger.info("Residues added to ResidueGenericNumberEquivalent")

            except:
                self.logger.error("Failed to add residues to ResidueGenericNumberEquivalent")
        self.logger.info('Inserted bulk {} (Index:{})'.format(len(bulk), index))
        Residue.objects.bulk_create(bulk)

    def update_protein_conformation(self, gprotein_list):
        # gprotein_list=['gnaz_human','gnat3_human', 'gnat2_human', 'gnat1_human', 'gnas2_human', 'gnaq_human', 'gnao_human', 'gnal_human', 'gnai3_human', 'gnai2_human','gnai1_human', 'gna15_human', 'gna14_human', 'gna12_human', 'gna11_human', 'gna13_human']
        state = ProteinState.objects.get(slug='active')

        # add new cgn protein conformations
        for g in gprotein_list:
            gp = Protein.objects.get(accession=g)

            try:
                pc, created = ProteinConformation.objects.get_or_create(protein=gp, state=state)
                self.logger.info('Created protein conformation')
            except:
                self.logger.error('Failed to create protein conformation')

        self.update_genericresiduenumber_and_proteinsegments(gprotein_list)

    def update_genericresiduenumber_and_proteinsegments(self, gprotein_list):

        # Parsing pdb uniprot file for generic residue numbers
        self.logger.info('Start parsing PDB_UNIPROT_ENSEMBLE_ALL')
        self.logger.info('Parsing file ' + self.gprotein_data_file)
        residue_data = pd.read_csv(self.gprotein_data_file, sep="\t", low_memory=False)

        residue_data = residue_data[residue_data.Uniprot_ID.notnull()]

        # residue_data = residue_data[residue_data['Uniprot_ID'].str.contains('_HUMAN')]

        residue_data = residue_data[residue_data['Uniprot_ACC'].isin(gprotein_list)]

        # filtering for human gproteins using list above
        residue_generic_numbers = residue_data['CGN']

        # Residue numbering scheme is the same for all added residue generic numbers (CGN)

        cgn_scheme = ResidueNumberingScheme.objects.get(slug='cgn')

        # purge line
        # ResidueGenericNumber.objects.filter(scheme_id=12).delete()

        for rgn in residue_generic_numbers.unique():
            d, s, n = rgn.split('.')
            seg = d+'.'+s
            ps = ProteinSegment.objects.get(slug=seg, proteinfamily='Alpha')

            rgnsp = []

            if (int(rgn.split('.')[2]) < 10):
                rgnsp = rgn.split('.')
                rgn_new = rgnsp[0] + '.' + rgnsp[1] + '.0' + rgnsp[2]
            else:
                rgn_new = rgn

            try:
                res_gen_num, created = ResidueGenericNumber.objects.get_or_create(label=rgn_new, scheme=cgn_scheme,
                                                                                  protein_segment=ps)
                self.logger.info('Created generic residue number')

            except:
                self.logger.error('Failed creating generic residue number')

        self.add_cgn_residues(gprotein_list)

    def cgn_add_proteins(self):
        self.logger.info('Start parsing PDB_UNIPROT_ENSEMBLE_ALL')
        self.logger.info('Parsing file ' + self.gprotein_data_file)

        # parsing file for accessions
        df = pd.read_csv(self.gprotein_data_file, sep="\t", low_memory=False)
        prot_type = 'purge'
        pfm = ProteinFamily()

        # Human proteins from CGN with families as keys: https://www.mrc-lmb.cam.ac.uk/CGN/about.html
        cgn_dict = {}
        cgn_dict['G-Protein'] = ['Gs', 'Gi/o', 'Gq/11', 'G12/13', 'GPa1 family']
        cgn_dict['100_001_001_001'] = ['GNAS2_HUMAN']
        cgn_dict['100_001_001_002'] = ['GNAL_HUMAN']
        cgn_dict['100_001_002_001'] = ['GNAI1_HUMAN']
        cgn_dict['100_001_002_002'] = ['GNAI2_HUMAN']
        cgn_dict['100_001_002_003'] = ['GNAI3_HUMAN']
        cgn_dict['100_001_002_004'] = ['GNAT1_HUMAN']
        cgn_dict['100_001_002_005'] = ['GNAT2_HUMAN']
        cgn_dict['100_001_002_006'] = ['GNAT3_HUMAN']
        cgn_dict['100_001_002_007'] = ['GNAZ_HUMAN']
        cgn_dict['100_001_002_008'] = ['GNAO_HUMAN']
        cgn_dict['100_001_003_001'] = ['GNAQ_HUMAN']
        cgn_dict['100_001_003_002'] = ['GNA11_HUMAN']
        cgn_dict['100_001_003_003'] = ['GNA14_HUMAN']
        cgn_dict['100_001_003_004'] = ['GNA15_HUMAN']
        cgn_dict['100_001_004_001'] = ['GNA12_HUMAN']
        cgn_dict['100_001_004_002'] = ['GNA13_HUMAN']
        cgn_dict['100_001_005_001'] = ['GPA1_YEAST']

        # list of all 16 proteins
        cgn_proteins_list = []
        for k in cgn_dict.keys():
            for p in cgn_dict[k]:
                if p.endswith('_HUMAN') or p.endswith('_YEAST'):
                    cgn_proteins_list.append(p)

        # print(cgn_proteins_list)
        # GNA13_HUMAN missing from cambridge file

        accessions = df.loc[df['Uniprot_ID'].isin(cgn_proteins_list)]
        accessions = accessions['Uniprot_ACC'].unique()

        # Create new residue numbering scheme
        self.create_cgn_rns()

        # purging one cgn entry
        # ResidueNumberingScheme.objects.filter(name='cgn').delete()

        rns = ResidueNumberingScheme.objects.get(slug='cgn')

        for a in accessions:
            up = self.parse_uniprot_file(a)

            # Fetch Protein Family for gproteins
            for k in cgn_dict.keys():
                name = str(up['entry_name']).upper()

                if name in cgn_dict[k]:
                    pfm = ProteinFamily.objects.get(slug=k)

            # Create new Protein
            self.cgn_creat_gproteins(pfm, rns, a, up)

        ###################ORTHOLOGS###############
        orthologs_pairs = []
        orthologs = []

        # Orthologs for human gproteins
        allprots = list(df.Uniprot_ID.unique())
        allprots = list(set(allprots) - set(cgn_proteins_list))

        # for gp in cgn_proteins_list:
        for p in allprots:
            # if str(p).startswith(gp.split('_')[0]):
            if str(p) in self.ortholog_mapping:
                # orthologs_pairs.append((str(p), self.ortholog_mapping[str(p)]+'_HUMAN'))
                orthologs.append(str(p))

        accessions_orth = df.loc[df['Uniprot_ID'].isin(orthologs)]
        accessions_orth = accessions_orth['Uniprot_ACC'].unique()
        for a in accessions_orth:
            up = self.parse_uniprot_file(a)
            # Fetch Protein Family for gproteins
            for k in cgn_dict.keys():
                name = str(up['entry_name']).upper()
                name = name.split('_')[0] + '_' + 'HUMAN'

                if name in cgn_dict[k]:
                    pfm = ProteinFamily.objects.get(slug=k)
                else:
                    try:
                        if self.ortholog_mapping[str(up['entry_name']).upper()] + '_HUMAN' in cgn_dict[k]:
                            pfm = ProteinFamily.objects.get(slug=k)
                    except:
                        pass

            # Create new Protein
            self.cgn_creat_gproteins(pfm, rns, a, up)

        # human gproteins
        orthologs_lower = [x.lower() for x in orthologs]
        # print(orthologs_lower)

        # orthologs to human gproteins
        cgn_proteins_list_lower = [x.lower() for x in cgn_proteins_list]

        # all gproteins
        gprotein_list = cgn_proteins_list_lower + orthologs_lower
        accessions_all = list(accessions_orth) + list(accessions)
        return list(accessions_all)

    def cgn_creat_gproteins(self, family, residue_numbering_scheme, accession, uniprot):
        # get/create protein source
        try:
            source, created = ProteinSource.objects.get_or_create(name=uniprot['source'],
                                                                  defaults={'name': uniprot['source']})
            if created:
                self.logger.info('Created protein source ' + source.name)
        except IntegrityError:
            source = ProteinSource.objects.get(name=uniprot['source'])

        # get/create species
        try:
            species, created = Species.objects.get_or_create(latin_name=uniprot['species_latin_name'],
                                                             defaults={
                                                                 'common_name': uniprot['species_common_name'],
                                                             })
            if created:
                self.logger.info('Created species ' + species.latin_name)
        except IntegrityError:
            species = Species.objects.get(latin_name=uniprot['species_latin_name'])

        # get/create protein sequence type
        # Wild-type for all sequences from source file
        try:
            sequence_type, created = ProteinSequenceType.objects.get_or_create(slug='wt',
                                                                               defaults={
                                                                                   'slug': 'wt',
                                                                                   'name': 'Wild-type',
                                                                               })
            if created:
                self.logger.info('Created protein sequence type Wild-type')
        except:
            self.logger.error('Failed creating protein sequence type Wild-type')

        # create protein
        p = Protein()
        p.family = family
        p.species = species
        p.source = source
        p.residue_numbering_scheme = residue_numbering_scheme
        p.sequence_type = sequence_type

        if accession:
            p.accession = accession
        p.entry_name = uniprot['entry_name'].lower()
        try:
            p.name = uniprot['names'][0].split('Guanine nucleotide-binding protein ')[1]
        except:
            p.name = uniprot['names'][0]
        p.sequence = uniprot['sequence']

        try:
            p.save()
            self.logger.info('Created protein {}'.format(p.entry_name))
        except Exception as msg:
            self.logger.error('Failed creating protein {} ({})'.format(p.entry_name, msg))
            p = Protein.objects.get(entry_name=p.entry_name)

        # protein aliases
        for i, alias in enumerate(uniprot['names']):
            pcgn = Protein.objects.get(entry_name=uniprot['entry_name'].lower())
            a = ProteinAlias()
            a.protein = pcgn
            a.name = alias
            a.position = i

            try:
                a.save()
                self.logger.info('Created protein alias ' + a.name + ' for protein ' + p.name)
            except:
                self.logger.error('Failed creating protein alias ' + a.name + ' for protein ' + p.name)

        # genes
        for i, gene in enumerate(uniprot['genes']):
            g = False
            try:
                g, created = Gene.objects.get_or_create(name=gene, species=species, position=i)
                if created:
                    self.logger.info('Created gene ' + g.name + ' for protein ' + p.name)
            except IntegrityError:
                g = Gene.objects.get(name=gene, species=species, position=i)

            if g:
                pcgn = Protein.objects.get(entry_name=uniprot['entry_name'].lower())
                g.proteins.add(pcgn)

        # structures
        # for i, structure in enumerate(uniprot['structures']):
        #     created = False
        #     try:
        #         res = structure[1]
        #         if res == '-':
        #             res = 0
        #         print(structure[0], structure)
        #         if len(SignprotStructure.objects.filter(PDB_code=structure[0]))==0:
        #             wl, wl_created = WebLink.objects.get_or_create(web_resource=WebResource.objects.get(slug='pdb'), index=structure[0])
        #             structure, created = SignprotStructure.objects.get_or_create(pdb_code=wl, resolution=res,
        #                                                                          protein=p, id=self.signprot_struct_ids())
        #             if created:
        #                 self.logger.info('Created structure ' + structure.PDB_code + ' for protein ' + p.name)
        #     except IntegrityError:
        #         self.logger.error('Failed creating structure ' + structure.PDB_code + ' for protein ' + p.name)
        #     if created:
        #         if g:
        #             pcgn = Protein.objects.get(entry_name=uniprot['entry_name'].lower())
        #             structure.protein = p
        #             structure.save()

    def signprot_struct_ids(self):
        structs = Structure.objects.exclude(structure_type__slug__startswith='af-').count()
        s_structs = SignprotStructure.objects.count()
        offset = 1000
        if s_structs == None:
            return structs+1+offset
        else:
            return structs+s_structs+1+offset

    def cgn_parent_protein_family(self):

        pf_cgn, created_pf = ProteinFamily.objects.get_or_create(slug='100', defaults={
            'name': 'G-Protein'}, parent=ProteinFamily.objects.get(slug='000'))

        pff_cgn = ProteinFamily.objects.get(slug='100', name='G-Protein')

        # Changed name "No Ligands" to "Gprotein"
        pf1_cgn = ProteinFamily.objects.get_or_create(slug='100_001', name='Alpha', parent=pff_cgn)

    def create_cgn_rns(self):
        rns_cgn, created = ResidueNumberingScheme.objects.get_or_create(slug='cgn', short_name='CGN', defaults={
            'name': 'Common G-alpha numbering scheme'})

    def cgn_create_proteins_and_families(self):

        # Creating single entries in "protein_family' table
        ProteinFamily.objects.filter(slug__startswith="100").delete()
        self.cgn_parent_protein_family()

        # Human proteins from CGN: https://www.mrc-lmb.cam.ac.uk/CGN/about.html
        cgn_dict = {}

        levels = ['2', '3']
        keys = ['Alpha', 'Gs', 'Gi/o', 'Gq/11', 'G12/13', '001']
        slug1 = '100'
        slug3 = ''
        i = 1

        cgn_dict['Alpha'] = ['001']
        cgn_dict['001'] = ['Gs', 'Gi/o', 'Gq/11', 'G12/13', 'GPa1 family']

        # Protein families to be added
        # Key of dictionary is level in hierarchy
        cgn_dict['1'] = ['Alpha']
        cgn_dict['2'] = ['001']
        cgn_dict['3'] = ['Gs', 'Gi/o', 'Gq/11', 'G12/13', 'GPa1 family']

        # Protein lines not to be added to Protein families
        cgn_dict['4'] = ['GNAS2', 'GNAL', 'GNAI1', 'GNAI2', 'GNAI3', 'GNAT1', 'GNAT2', 'GNAT3', 'GNAZ', 'GNAO', 'GNAQ',
                         'GNA11', 'GNA14', 'GNA15', 'GNA12', 'GNA13', 'GPa1']

        grouped_subtypes = OrderedDict()
        for j in cgn_dict['3']:
            grouped_subtypes[j] = []
            for k in cgn_dict['4']:
                if j == 'Gs' and k in ['GNAS2', 'GNAL']:
                    grouped_subtypes[j].append(k)
                elif j == 'Gi/o' and k in ['GNAI1', 'GNAI2', 'GNAI3', 'GNAT1', 'GNAT2', 'GNAT3', 'GNAZ', 'GNAO']:
                    grouped_subtypes[j].append(k)
                elif j == 'Gq/11' and k in ['GNAQ', 'GNA11', 'GNA14', 'GNA15']:
                    grouped_subtypes[j].append(k)
                elif j == 'G12/13' and k in ['GNA12', 'GNA13']:
                    grouped_subtypes[j].append(k)
                elif j =='GPa1 family' and k=='GPa1':
                    grouped_subtypes[j].append(k)

        for entry in cgn_dict['001']:

            name = entry

            slug2 = '_001'
            slug3 = '_00' + str(i)

            slug = slug1 + slug2 + slug3

            slug3 = ''
            i = i + 1

            pff_cgn = ProteinFamily.objects.get(slug='100_001')

            new_pf, created = ProteinFamily.objects.get_or_create(slug=slug, name=entry, parent=pff_cgn)
            j = 1
            for en in grouped_subtypes[entry]:
                ort_slug = slug + '_00' + str(j)
                new_ort_fam, created = ProteinFamily.objects.get_or_create(slug=ort_slug, name=en, parent=new_pf)
                j += 1

        # function to create necessary arguments to add protein entry
        self.cgn_add_proteins()

    def parse_uniprot_file(self, accession):
        filename = accession + '.txt'
        local_file_path = os.sep.join([self.local_uniprot_dir, filename])
        remote_file_path = self.remote_uniprot_dir + filename

        up = {}
        up['genes'] = []
        up['names'] = []
        up['structures'] = []

        read_sequence = False
        remote = False

        # record whether organism has been read
        os_read = False

        # should local file be written?
        local_file = False

        try:
            if os.path.isfile(local_file_path):
                uf = open(local_file_path, 'r')
                self.logger.info('Reading local file ' + local_file_path)
            else:
                uf = urlopen(remote_file_path)
                remote = True
                self.logger.info('Reading remote file ' + remote_file_path)
                local_file = open(local_file_path, 'w')

            for raw_line in uf:
                # line format
                if remote:
                    line = raw_line.decode('UTF-8')
                else:
                    line = raw_line

                # write to local file if appropriate
                if local_file:
                    local_file.write(line)

                # end of file
                if line.startswith('//'):
                    break

                # entry name and review status
                if line.startswith('ID'):
                    split_id_line = line.split()
                    up['entry_name'] = split_id_line[1].lower()
                    review_status = split_id_line[2].strip(';')
                    if review_status == 'Unreviewed':
                        up['source'] = 'TREMBL'
                    elif review_status == 'Reviewed':
                        up['source'] = 'SWISSPROT'

                # species
                elif line.startswith('OS') and not os_read:
                    species_full = line[2:].strip().strip('.')
                    species_split = species_full.split('(')
                    up['species_latin_name'] = species_split[0].strip()
                    if len(species_split) > 1:
                        up['species_common_name'] = species_split[1].strip().strip(')')
                    else:
                        up['species_common_name'] = up['species_latin_name']
                    os_read = True

                # names
                elif line.startswith('DE'):
                    split_de_line = line.split('=')
                    if len(split_de_line) > 1:
                        split_segment = split_de_line[1].split('{')
                        up['names'].append(split_segment[0].strip().strip(';'))

                # genes
                elif line.startswith('GN'):
                    split_gn_line = line.split(';')
                    for segment in split_gn_line:
                        if '=' in segment:
                            split_segment = segment.split('=')
                            split_segment = split_segment[1].split(',')
                            for gene_name in split_segment:
                                split_gene_name = gene_name.split('{')
                                up['genes'].append(split_gene_name[0].strip())

                # structures
                elif line.startswith('DR') and 'PDB' in line and not 'sum' in line:
                    split_gn_line = line.split(';')
                    up['structures'].append([split_gn_line[1].lstrip(), split_gn_line[3].lstrip().split(" A")[0]])

                # sequence
                elif line.startswith('SQ'):
                    split_sq_line = line.split()
                    seq_len = int(split_sq_line[2])
                    read_sequence = True
                    up['sequence'] = ''
                elif read_sequence == True:
                    up['sequence'] += line.strip().replace(' ', '')

            # close the Uniprot file
            uf.close()
        except:
            return False

        # close the local file if appropriate
        if local_file:
            local_file.close()

        return up


class SeqCompare(object):
    def __init__(self):
        pass

    def align(self, seq1, seq2):
        for p in pairwise2.align.globalms(seq1, seq2, 8, 5, -5, -5):
            return format_alignment(*p).split('\n')[2]
