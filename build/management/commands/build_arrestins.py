import logging
import os
import sys
from urllib.request import urlopen
from xml.etree.ElementTree import fromstring

import pandas as pd
import requests
import xlrd
import django.apps
from django.conf import settings
from django.core.management.base import BaseCommand
from django.core.management.color import no_style
from django.db import IntegrityError, connection
from common.tools import test_model_updates, parse_uniprot_file, build_entrezgeneid_lookup_dict, select_entrez_id
from common.models import WebLink, WebResource
from protein.models import (Gene, Protein, ProteinAlias, ProteinConformation,
                            ProteinFamily, ProteinSegment,
                            ProteinSequenceType, ProteinSource, ProteinState, Species)
from residue.models import (Residue, ResidueGenericNumber,
                            ResidueGenericNumberEquivalent,
                            ResidueNumberingScheme)
from signprot.models import SignprotStructure
from structure.models import Structure


class Command(BaseCommand):
    help = 'Build Arrestin proteins'

    # source files
    arrestin_data_file = os.sep.join([settings.DATA_DIR, 'arrestin_data', 'ortholog_alignment.xlsx'])
    local_uniprot_dir = os.sep.join([settings.DATA_DIR, 'protein_data', 'uniprot'])
    entrezid_source_file = os.sep.join([settings.DATA_DIR, 'gene_data', 'entrez_id_lookup.txt'])
    
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

    def handle(self, *args, **options):
        self.entrez_lookup = build_entrezgeneid_lookup_dict(self.entrezid_source_file, self.logger)
        self.options = options
        if options['filename']:
            filenames = options['filename']
        else:
            filenames = False

        try:
            self.purge_can_residues()
            self.logger.info('PASS: purge_can_residues')
            self.purge_can_proteins()
            self.logger.info('PASS: purge_can_proteins')

            # add proteins
            self.can_create_families()
            self.logger.info('PASS: can_create_families')
            self.can_add_proteins()
            self.logger.info('PASS: can_add_proteins')

            # add residues
            self.add_can_residues()
            self.logger.info('PASS: add_can_residues')

            #Perform model update check
            test_model_updates(self.all_models, self.tracker, check=True)
        except Exception as msg:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print(exc_type, exc_obj, fname, exc_tb.tb_lineno)
            self.logger.error(msg)

    def purge_can_residues(self):
        """Purge residues."""
        try:
            Residue.objects.filter(generic_number_id__scheme__slug="can").delete()
        except Exception as msg:
            self.logger.warning('Existing Residue data cannot be deleted', msg)

    def purge_can_proteins(self):
        """Purge proteins."""
        try:
            Protein.objects.filter(residue_numbering_scheme__slug='can').delete()
        except Exception as msg:
            self.logger.warning('Protein to delete not found' + str(msg))

    def add_can_residues(self):
        """Add CAN residues from source file provided by Andrija Sente."""
        # Parsing pdb uniprot file for residues
        self.logger.info('Start parsing ARRESTIN RESIDUES')
        self.logger.info('Parsing file ' + self.arrestin_data_file)
        residue_data = pd.read_excel(self.arrestin_data_file)

        can_scheme = ResidueNumberingScheme.objects.get(slug='can')

        can_dict = residue_data[residue_data.ID == 'CAN_id'].iloc[:, 3:].to_dict('list')

        # Loop over data table, but skip "CAN_posand" and "CAN_id" from current input file
        for index, row in residue_data[2:].iterrows():

            try:
                # for now only allow for ortholog with uniprot entries:
                if not row['AccessionID'].startswith('ENS'):
                    # fetch protein for protein conformation
                    pr, c = Protein.objects.get_or_create(accession=row['AccessionID'])

                    # fetch protein conformation
                    pc, c = ProteinConformation.objects.get_or_create(protein_id=pr)
                else:
                    continue
            except:
                print('error making/getting protein', row['AccessionID'])
                continue

            # loop over residue generic number
            sequence_number = 1
            for aln_pos in can_dict:

                canId = can_dict[aln_pos][0]

                # Add '0' infront of single digit positions
                if (int(canId.split('.')[2]) < 10):
                    rgnsp = canId.split('.')
                    canId = rgnsp[0] + '.' + rgnsp[1] + '.0' + rgnsp[2]

                ps, c = ProteinSegment.objects.get_or_create(slug=canId.split('.')[1], proteinfamily='Arrestin')

                rgn, c = ResidueGenericNumber.objects.get_or_create(label=canId, scheme=can_scheme, protein_segment=ps)

                # only add AA information if not gap
                if not row[aln_pos] == '-':

                    try:
                        Residue.objects.get_or_create(sequence_number=sequence_number, protein_conformation=pc,
                                                      amino_acid=row[aln_pos], generic_number=rgn, display_generic_number=rgn,
                                                      protein_segment=ps)
                        sequence_number += 1
                    except Exception as msg:
                        print("failed to add residue", msg)
                        self.logger.error("Failed to add residues", msg)

                    # Add also to the ResidueGenericNumberEquivalent table needed for single residue selection
                    try:
                        ResidueGenericNumberEquivalent.objects.get_or_create(label=rgn.label, default_generic_number=rgn,
                                                                             scheme=can_scheme)  # Update scheme_id
                    except Exception as msg:
                        print("failed to add residue generic number", msg)
                        self.logger.error("Failed to add residues to ResidueGenericNumberEquivalent")

    def get_uniprot_accession_id(self, response_xml):
        # TODO: This function seems to be legacy, perhaps it can be deleted?
        """Get Uniprot accession ID."""
        root = fromstring(response_xml)
        return next(
            # el for el in root.getchildren()[0].getchildren()
            el for el in root[0]
            if el.attrib['dbSource'] == 'UniProt'
        ).attrib['dbAccessionId']

    def map_pdb_to_uniprot(self, pdb_id):
        # TODO: This function seems to be legacy, perhaps it can be deleted?
        # Due to the new RCSB API this doesn't even work after December 2020
        """Get uniprot ID from PDB ID."""
        pdb_mapping_url = 'https://www.rcsb.org/pdb/rest/das/pdb_uniprot_mapping/alignment'
        pdb_mapping_response = requests.get(
            pdb_mapping_url, params={'query': pdb_id}
        ).text
        uniprot_id = self.get_uniprot_accession_id(pdb_mapping_response)
        return uniprot_id

    def can_add_proteins(self):
        """Add arrestin proteins."""
        self.logger.info('Start adding ARRESTIN proteins')
        self.logger.info('Parsing file ' + self.arrestin_data_file)

        # Import ortholog alignment as pandas dataframe
        residue_data = pd.read_excel(self.arrestin_data_file)

        # Create new residue numbering scheme
        self.create_can_rns()

        rns = ResidueNumberingScheme.objects.get(slug='can')
        state = ProteinState.objects.get(slug='active')

        arrestins = residue_data[2:].Ortholog.unique()

        for arrestin in arrestins:

            pfm = ProteinFamily.objects.get(name=arrestin)

            for accession in residue_data[residue_data.Ortholog == arrestin].AccessionID.unique():
                # only allow uniprot accession:
                if not accession.startswith('ENS'):
                    up = parse_uniprot_file(accession= accession, logger=self.logger, local_uniprot_dir=self.local_uniprot_dir)
                    #     if len(up['genes']) == 0:
                    #         print('There is no GN field in the uniprot!', accession)
                    #         self.logger.error('There is no GN field in the uniprot! {}'.format(accession))
                    #         continue
                    if up == False or not 'source' in up:
                        print('No source found, probably deprecated!', accession)
                        self.logger.error('No source found, probably deprecated! {}'.format(accession))
                        continue

                    # Create new Protein
                    self.can_create_arrestins(pfm, rns, accession, up)

                    # add new can protein conformations
                    try:
                        arrestin = Protein.objects.get(accession=accession)

                        pc, created = ProteinConformation.objects.get_or_create(protein=arrestin, state=state)
                        self.logger.info('Created protein conformation')
                    except Exception as msg:
                        self.logger.error('Failed to create protein conformation', msg)

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
                entrez_geneid = select_entrez_id(i, uniprot, self.entrez_lookup)
                
                if entrez_geneid is not None:
                    resource = WebResource.objects.get(slug='entrez_gene')
                    entrez_link, created = WebLink.objects.get_or_create(web_resource=resource, index=entrez_geneid)
                else:
                    entrez_link = None

                g, created = Gene.objects.get_or_create(name=gene, species=species, position=i, 
                                                        entrez_id=entrez_geneid, entrez_weblink=entrez_link)                    
                
                if created:
                    self.logger.info('Created gene ' + g.name + ' for protein ' + p.name)
            except IntegrityError:
                g = Gene.objects.get(name=gene, species=species, position=i)

            if g:
                pcgn = Protein.objects.get(entry_name=uniprot['entry_name'].lower())
                g.proteins.add(pcgn)

        # structures
        # for i, structure in enumerate(uniprot['structures']):
        #     # try:
        #     res = structure[1]
        #     if res == '-':
        #         res = 0

        #     structure, created = SignprotStructure.objects.get_or_create(PDB_code=structure[0], resolution=res, protein = p, id=self.signprot_struct_ids())
        #     if created:
        #         self.logger.info('Created structure ' + structure.PDB_code + ' for protein ' + p.name)

    def signprot_struct_ids(self):
        structs = Structure.objects.filter(structure_type__origin='experiment').count()
        s_structs = SignprotStructure.objects.count()
        offset = 1000
        if s_structs == None:
            return structs + 1 + offset
        else:
            return structs + s_structs + 1 + offset

    def create_can_rns(self):
        """Add new numbering scheme entry_name."""
        #        rns_can, created= ResidueNumberingScheme.objects.get_or_create(slug='can', short_name='CAN', defaults={
        #            'name': 'Common arrestin numbering scheme'})

        try:
            rns_can, created = ResidueNumberingScheme.objects.get_or_create(slug='can', short_name='CAN',
                                                                            defaults={'name': 'Common arrestin numbering scheme'})
            if created:
                self.logger.info('Created Arrestin Numbering ' + rns_can.slug)
        except IntegrityError:
            rns_can = ResidueNumberingScheme.objects.get(slug='can')
            self.logger.info('Integrity Error on creating can numbering')

    def can_create_families(self):
        """Purge and create arrestin in protein_family."""
        ProteinFamily.objects.filter(slug__startswith="200").delete()

        # 4 arrestin subtypes, two of which are primarily expressed in the retina and bind only to visual opsins (arrestin 1 and arrestin 4), while the other two (β-arrestin 1 and β-arrestin 2) interact with the remaining ~800 GPCRs

        can_dict = {}
        can_dict['Arrestin'] = ['Beta', 'Visual']
        can_dict['Beta'] = ['ARRB1', 'ARRB2']
        can_dict['Visual'] = ['ARRC', 'ARRS']

        pff_can, created_pf = ProteinFamily.objects.get_or_create(slug='200', defaults={
            'name': 'Arrestins'})
        pf1_can = ProteinFamily.objects.get_or_create(slug='200_000', name='Arrestin', parent=pff_can)

        for i, family in enumerate(can_dict['Arrestin']):

            # slug for the different levels
            fam_slug = '200_000_00' + str(i + 1)

            pff_can = ProteinFamily.objects.get(slug='200_000')
            new_pf, created = ProteinFamily.objects.get_or_create(slug=fam_slug, name=family, parent=pff_can)

            for i, protein in enumerate(can_dict[family]):
                prot_slug = fam_slug + '_00' + str(i + 1)

                pff_fam = ProteinFamily.objects.get(slug=fam_slug)
                new_pf, created = ProteinFamily.objects.get_or_create(slug=prot_slug, name=protein, parent=pff_fam)

