from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import connection
from django.db import IntegrityError

from common.models import WebResource, WebLink

from protein.models import (Protein, ProteinConformation, ProteinState, ProteinFamily, ProteinAlias, ProteinSequenceType, Species, Gene, ProteinSource, ProteinSegment)

from residue.models import (ResidueNumberingScheme, ResidueGenericNumber, Residue, ResidueGenericNumberEquivalent)

from signprot.models import SignprotStructure
import pandas as pd

from optparse import make_option

import requests
from xml.etree.ElementTree import fromstring

import math, os, csv
import numpy as np
import logging
import requests

from urllib.request import urlopen

class Command(BaseCommand):
    help = 'Build Arrestin proteins'

    # source file directory
    arrestin_data_file = os.sep.join([settings.DATA_DIR, 'arrestin_data', 'TableS2.xlsx'])

    local_uniprot_dir = os.sep.join([settings.DATA_DIR, 'protein_data', 'uniprot'])
    remote_uniprot_dir = 'http://www.uniprot.org/uniprot/'
    logger = logging.getLogger(__name__)

    def add_arguments(self, parser):
        parser.add_argument('--filename', action='append', dest='filename',
            help='Filename to import. Can be used multiple times')

    def handle(self, *args, **options):
        if options['filename']:
            filenames = options['filename']
        else:
            filenames = False

        ## create protein families and infrastructure
        try:
            self.purge_can_residues()
            self.purge_can_proteins()

            self.can_create_proteins_and_families()

        except Exception as msg:
            print(msg)
            self.logger.error(msg)

        ## add residues from can db
        try:
            human_and_orths = self.can_add_proteins()

            self.update_protein_conformation(human_and_orths)
        except Exception as msg:
            self.logger.error(msg)

    def purge_can_residues(self):
        try:
            Residue.objects.filter(generic_number_id__scheme__slug="can").delete()
        except:
            self.logger.warning('Existing Residue data cannot be deleted')

    def purge_can_proteins(self):
        try:
            Protein.objects.filter(residue_numbering_scheme__slug='can').delete()
        except:
            self.logger.info('Protein to delete not found')

    def add_can_residues(self, arrestin_list):

        # Parsing pdb uniprot file for residues
        self.logger.info('Start parsing ARRESTIN RESIDUES')
        self.logger.info('Parsing file ' + self.arrestin_data_file)
        residue_data =  pd.read_excel(self.arrestin_data_file, sheetname=1)
        # residue_data = residue_data.loc[residue_data['Uniprot_ACC'].isin(arrestin_list)]

        ## Example data to populate a table for further infrastructure work
        residue_data = residue_data[residue_data['pdb_id']=='1ayr']

        can_scheme = ResidueNumberingScheme.objects.get(slug='can')

        for index, row in residue_data.iterrows():
            # fetch protein for protein conformation
            pr, c = Protein.objects.get_or_create(accession='P10523')

            # fetch protein conformation
            pc, c = ProteinConformation.objects.get_or_create(protein_id=pr)

            # fetch residue generic number
            rgnsp = []
            if(int(row['CAN'].split('.')[2])<10):
                rgnsp = row['CAN'].split('.')
                rgn_new = rgnsp[0] + '.' + rgnsp[1] + '.0' + rgnsp[2]
                rgn, c = ResidueGenericNumber.objects.get_or_create(label=rgn_new)

            else:
                rgn, c = ResidueGenericNumber.objects.get_or_create(label=row['CAN'])

            ps, c = ProteinSegment.objects.get_or_create(slug=row['CAN'].split('.')[1], proteinfamily='Arrestin')
            try:
                Residue.objects.get_or_create(sequence_number=row['pdbPos'], protein_conformation=pc, amino_acid=row['res_id'][0], generic_number=rgn, display_generic_number=rgn, protein_segment=ps)

            except:
                print("failed to add residue")
                self.logger.error("Failed to add residues")

             # Add also to the ResidueGenericNumberEquivalent table needed for single residue selection
            try:
                ResidueGenericNumberEquivalent.objects.get_or_create(label=rgn.label,default_generic_number=rgn, scheme=can_scheme) ## Update scheme_id

            except:
                print("failed to add residue generic number")
                self.logger.error("Failed to add residues to ResidueGenericNumberEquivalent")

    def update_protein_conformation(self, arrestin_list):

        state = ProteinState.objects.get(slug='active')

        # add new can protein conformations
        for p in arrestin_list:
            arrestin = Protein.objects.get(accession=p)

            try:
                pc, created = ProteinConformation.objects.get_or_create(protein=arrestin, state=state, template_structure=None)
                self.logger.info('Created protein conformation')
            except:
                self.logger.error('Failed to create protein conformation')

        self.update_genericresiduenumber_and_proteinsegments(arrestin_list)

    def update_genericresiduenumber_and_proteinsegments(self, arrestin_list):

        # Parsing pdb uniprot file for generic residue numbers
        self.logger.info('Start parsing ARRESTIN RESIDUE FILE')
        self.logger.info('Parsing file ' + self.arrestin_data_file)
        residue_data =  pd.read_excel(self.arrestin_data_file, sheetname=1)

        # residue_data = residue_data[residue_data.Uniprot_ID.notnull()]
        # residue_data = residue_data[residue_data['Uniprot_ACC'].isin(arrestin_list)]

        # filtering for human arrestins using list above
        residue_generic_numbers = residue_data['CAN'].unique()

        can_scheme = ResidueNumberingScheme.objects.get(slug='can')

        for rgn in residue_generic_numbers:
            ps, c = ProteinSegment.objects.get_or_create(slug=rgn.split('.')[1], proteinfamily='Arrestin')

            rgnsp = []

            if(int(rgn.split('.')[2])<10):
                rgnsp = rgn.split('.')
                rgn_new = rgnsp[0] + '.' + rgnsp[1] + '.0' + rgnsp[2]
            else:
                rgn_new = rgn

            try:
                res_gen_num, created =  ResidueGenericNumber.objects.get_or_create(label=rgn_new, scheme=can_scheme, protein_segment=ps)
                self.logger.info('Created generic residue number')

            except:
                self.logger.error('Failed creating generic residue number')

        self.add_can_residues(arrestin_list)

    def get_uniprot_accession_id(self, response_xml):
        root = fromstring(response_xml)
        return next(
            el for el in root.getchildren()[0].getchildren()
            if el.attrib['dbSource'] == 'UniProt'
        ).attrib['dbAccessionId']

    def map_pdb_to_uniprot(self, pdb_id):
        pdb_mapping_url = 'http://www.rcsb.org/pdb/rest/das/pdb_uniprot_mapping/alignment'
        pdb_mapping_response = requests.get(
            pdb_mapping_url, params={'query': pdb_id}
        ).text
        uniprot_id = self.get_uniprot_accession_id(pdb_mapping_response)
        return uniprot_id

    def can_add_proteins(self):

        self.logger.info('Start adding ARRESTIN proteins')
        self.logger.info('Parsing file ' + self.arrestin_data_file)

        # parsing file for accessions
        df = pd.read_excel(self.arrestin_data_file, sheetname=1)
        pfm = ProteinFamily()

        # PDB to Uniprot_ID
        accessions = []
        for pdb_id in df.pdb_id.unique():
            accessions.append(self.map_pdb_to_uniprot(pdb_id))

        # 4 arrestin subtypes, two of which are primarily expressed in the retina and bind only to visual opsins (arrestin 1 and arrestin 4), while the other two (β-arrestin 1 and β-arrestin 2) interact with the remaining ~800 GPCRs
        translation = {'Beta':'200_000_001', 'Visual':'200_000_002'}

        can_dict = {}
        can_dict['Arrestin']=['Beta','Visual']
        can_dict['200_000_001']=['arrb1_human','arrb2_human']
        can_dict['200_000_002']=['arrc_human','arrs_human']

        # Create new residue numbering scheme
        self.create_can_rns()

        accessions = ['P49407','P10523','P32121','P36575']
        rns = ResidueNumberingScheme.objects.get(slug='can')

        for accession in list(set(accessions)):
            up = self.parse_uniprot_file(accession)

            #Fetch Protein Family for arrestins
            for k in can_dict.keys():
                entry_name = str(up['entry_name']).lower()

                if entry_name in can_dict[k]:
                    pfm = ProteinFamily.objects.get(slug=k)

            # Create new Protein
            self.can_create_arrestins(pfm, rns, accession, up)

        ################## ORTHOLOGS ##############
        ## Add orthologs of arrestins

        orthologs_pairs =[]
        orthologs =[]
        #
        # # Orthologs for human arrestins
        # allprots = list(df.Uniprot_ID.unique())
        # allprots = list(set(allprots) - set(can_proteins_list))
        #
        # for gp in can_proteins_list:
        #     for p in allprots:
        #         if str(p).startswith(gp.split('_')[0]):
        #             orthologs_pairs.append((str(p), gp))
        #             orthologs.append(str(p))
        #
        # accessions_orth = df.loc[df['Uniprot_ID'].isin(orthologs)]
        # accessions_orth = accessions_orth['Uniprot_ACC'].unique()

        # accessions_orth = []
        # for accession in accessions_orth:
        #     up = self.parse_uniprot_file(a)
        #
        #     # Fetch Protein Family for arrestins
        #     for k in can_dict.keys():
        #         name = str(up['entry_name']).upper()
        #         name = name.split('_')[0]+'_'+'HUMAN'
        #
        #         if name in can_dict[k]:
        #             pfm = ProteinFamily.objects.get(slug=k)
        #
        #     # Create new Protein
        #     self.can_create_arrestins(pfm, rns, accession, up)
        #
        # # human arrestins
        # orthologs_lower = [x.lower() for x in orthologs]
        #
        # # orthologs to human arrestins
        # can_proteins_list_lower = [x.lower() for x in can_proteins_list]
        #
        # # all arrestins
        # accessions_all = list(accessions_orth) + list(accessions)

        return list(accessions)

    def can_create_arrestins(self, family, residue_numbering_scheme, accession, uniprot):
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
        p.name = uniprot['names'][0]
        p.sequence = uniprot['sequence']

        try:
            p.save()
            self.logger.info('Created protein {}'.format(p.entry_name))
        except:
            self.logger.error('Failed creating protein {}'.format(p.entry_name))

        # protein aliases
        for i, alias in enumerate(uniprot['names']):
            pcan = Protein.objects.get(entry_name=uniprot['entry_name'].lower())
            a = ProteinAlias()
            a.protein = pcan
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
                pcan = Protein.objects.get(entry_name=uniprot['entry_name'].lower())
                g.proteins.add(pcan)

        # structures
        for i, structure in enumerate(uniprot['structures']):
            try:
                res = structure[1]
                if res == '-':
                    res = 0

                structure, created = SignprotStructure.objects.get_or_create(PDB_code=structure[0], resolution=res)
                if created:
                    self.logger.info('Created structure ' + structure.PDB_code + ' for protein ' + p.name)
            except IntegrityError:
                self.logger.error('Failed creating structure ' + structure.PDB_code + ' for protein ' + p.name)

            if g:
                pcan = Protein.objects.get(entry_name=uniprot['entry_name'].lower())
                structure.origin.add(pcan)
                structure.save()

    def can_parent_protein_family(self):
        ## New protein family entry

        pf_can, created_pf = ProteinFamily.objects.get_or_create(slug='200', defaults={
            'name': 'Arrestins'})

        pff_can = ProteinFamily.objects.get(slug='200', name='Arrestins')
        pf1_can = ProteinFamily.objects.get_or_create(slug='200_000', name='Arrestin', parent=pff_can)

    def create_can_rns(self):
        ## New numbering scheme entry_name

        rns_can, created= ResidueNumberingScheme.objects.get_or_create(slug='can', short_name='CAN', defaults={
            'name': 'Common arrestin numbering scheme'})

    def can_create_proteins_and_families(self):

        # Purge and create arrestin in "protein_family'
        ProteinFamily.objects.filter(slug__startswith="200").delete()
        self.can_parent_protein_family()

        can_dict = {}
        can_dict['Arrestin'] = ['000']
        can_dict['000'] = ['Beta','Visual']

        # Protein families to be added
        # Key of dictionary is level in hierarchy
        can_dict['1'] = ['Arrestins']
        can_dict['2'] = ['000']
        can_dict['3'] = ['Beta','Visual']

        # Protein lines not to be added to Protein families
        can_dict['4'] = ['ARRB2','ARRB1','ARRS','ARRC']

        for i,entry in enumerate(can_dict['000']):

            # slug for the different levels
            slug = '200' + '_000' + '_00' + str(i+1)

            pff_can = ProteinFamily.objects.get(slug='200_000')

            new_pf, created = ProteinFamily.objects.get_or_create(slug=slug, name=entry, parent=pff_can)

        ## function to create necessary arguments to add protein entry
        self.can_add_proteins()

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
                    up['structures'].append([split_gn_line[1].lstrip(),split_gn_line[3].lstrip().split(" A")[0]])

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
