from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import connection
from django.db import IntegrityError

from build.management.commands.base_build import Command as BaseBuild
from protein.models import (Protein, ProteinConformation, ProteinState, ProteinFamily, ProteinAlias,
        ProteinSequenceType, Species, Gene, ProteinSource, ProteinSegment)
from common.models import WebResource, WebLink
from residue.models import ResidueNumberingScheme
from common.tools import test_model_updates

import django.apps
import shlex
import os
from urllib.request import urlopen

import pandas as pd
import numpy  as np
import math
import yaml


class Command(BaseBuild):
    help = 'Reads source data and creates protein families, proteins, and associated tables'

    protein_source_file = os.sep.join([settings.DATA_DIR, 'protein_data', 'proteins_and_families.txt'])
    local_uniprot_dir = os.sep.join([settings.DATA_DIR, 'protein_data', 'uniprot'])
    with open(os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation', 'sequences.yaml']), 'r') as f:
        excel_sequences = yaml.load(f, Loader=yaml.FullLoader)
    remote_uniprot_dir = 'http://www.uniprot.org/uniprot/'
    #Setting the variables for the test tracking of the model upadates
    tracker = {}
    all_models = django.apps.apps.get_models()[6:]
    test_model_updates(all_models, tracker, initialize=True)

    def handle(self, *args, **options):
        # use a smaller protein file if in test mode
        if options['test']:
            self.protein_source_file = os.sep.join([settings.DATA_DIR, 'protein_data',
                'proteins_and_families_test.txt'])

        # create parent protein family, 000
        try:
            self.create_parent_protein_family()
            test_model_updates(self.all_models, self.tracker, check=True)
        except Exception as msg:
            # print(msg)
            self.logger.error(msg)

        # create proteins and families
        try:
            self.create_proteins_and_families()
            test_model_updates(self.all_models, self.tracker, check=True)
        except Exception as msg:
            print(msg)
            self.logger.error(msg)


    def create_parent_protein_family(self):
        pf = ProteinFamily.objects.get_or_create(slug='000', defaults={
            'name': 'Parent family'})

    def create_proteins_and_families(self):
        self.logger.info('CREATING PROTEINS')
        self.logger.info('Parsing file ' + self.protein_source_file)

        # get/create protein sequence type
        # Wild-type for all sequences from source file, isoforms handled separately
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

        with open(self.protein_source_file, "r", encoding='UTF-8') as protein_file:
            # family hierarchy is determined by indent
            spaces_per_indent_level = 4
            last_indent = 0
            level_family_counter = [0]
            parent_family = [0]
            residue_numbering_scheme = False

            for row in protein_file:
                # determine the level of indent
                indent = len(row) - len(row.lstrip(' '))
                indent = int(indent / spaces_per_indent_level)

                # has the indent changed
                if indent != last_indent:
                    # did the level increase or decrease?
                    if indent > last_indent:
                        parent_family.append(0)
                        level_family_counter.append(0)
                    elif indent < last_indent:
                        for j in range(last_indent-indent):
                            parent_family.pop()
                            level_family_counter.pop()

                    last_indent = indent

                # process the line
                # is this a family or protein line?

                ###########
                # family/family type
                ###########
                create_protein = False
                if row.strip().startswith('"'): # protein row
                    split_row = row.strip().split('","')
                    family_name = split_row[4]
                    create_protein = True
                else: # family row
                    # check for residue numbering scheme
                    split_row = row.strip().split('|')
                    if len(split_row) > 1:
                        try:
                            rns = split_row[1].strip()
                            residue_numbering_scheme = ResidueNumberingScheme.objects.get(slug=rns)
                        except:
                            # abort if residue numbering scheme is not found in db
                            raise Exception('Residue numbering scheme ' + rns + ' not found, aborting')
                    else:
                        if not residue_numbering_scheme:
                            # abort if no residue numbering scheme is specified in the protein source file
                            raise Exception('No residue numbering scheme specified in source data, aborting')

                    family_name = split_row[0].strip()

                    # create the protein family
                    created_family = self.create_protein_family(family_name, indent, parent_family,
                        level_family_counter)
                    if created_family:
                        parent_family = created_family['parent_family']
                        level_family_counter = created_family['level_family_counter']
                    else:
                        continue

                ###########
                # protein
                ###########
                if create_protein:
                    protein_name = split_row[4]

                    # accession codes for human, mouse and rat receptors (from IU-PHAR)
                    # accessions = [split_row[15], split_row[31], split_row[23]]
                    accessions = [split_row[15]]

                    # create a family for this protein
                    created_family = self.create_protein_family(protein_name, indent, parent_family,
                        level_family_counter)
                    if created_family:
                        pf = created_family['pf']
                        parent_family = created_family['parent_family']
                        level_family_counter = created_family['level_family_counter']
                    else:
                        continue

                    for protein_accession in accessions:
                        # skip protein if there is no accession code
                        if not protein_accession:
                            self.logger.error('No accession code for protein ' + protein_name + ', skipping')
                            continue

                        # skip protein if accession code already exists
                        if Protein.objects.filter(accession=protein_accession).count() > 0:
                            self.logger.error('Protein accession ' + protein_accession + ' already exists, skipping')
                            continue

                        # parse uniprot file for this protein
                        self.logger.info('Parsing uniprot file for protein ' + protein_name)
                        up = self.parse_uniprot_file(protein_accession)
                        if not up:
                            self.logger.error('Failed parsing uniprot file for protein ' + protein_name + ', skipping')
                            continue

                        self.create_protein(protein_name, pf, sequence_type, residue_numbering_scheme,
                            protein_accession, up)

        self.logger.info('COMPLETED CREATING PROTEINS')

    def create_protein(self, name, family, sequence_type, residue_numbering_scheme, accession, uniprot):
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

        # create protein
        p = Protein()
        p.family = family
        p.species = species
        p.source = source
        p.residue_numbering_scheme = residue_numbering_scheme
        p.sequence_type = sequence_type
        if accession:
            p.accession = accession
        p.entry_name = uniprot['entry_name']
        p.name = name
        p.sequence = uniprot['sequence']

        try:
            p.save()
            self.logger.info('Created protein {}'.format(p.entry_name))
        except Exception as e:
            self.logger.error('Failed creating protein {} {}'.format(p.entry_name, str(e)))
            print('WARNING:', p.family, p.species, p.source, p.residue_numbering_scheme, p.sequence_type, p.accession, p.entry_name, p.name, p.sequence, e)

        # protein conformations
        try:
            ps, created = ProteinState.objects.get_or_create(slug=settings.DEFAULT_PROTEIN_STATE,
                defaults={'name': settings.DEFAULT_PROTEIN_STATE.title()})
        except IntegrityError:
            ps = ProteinState.objects.get(slug=settings.DEFAULT_PROTEIN_STATE)

        pc = ProteinConformation.objects.create(protein=p, state=ps)

        # protein uniprot links
        if accession:
            for ac in uniprot['accessions'][1:]:
                resource = WebResource.objects.get(slug='uniprot')

                # create a link
                link, created = WebLink.objects.get_or_create(web_resource=resource, index=ac)

                # add the link to this protein
                p.web_links.add(link)

        # protein aliases
        for i, alias in enumerate(uniprot['names']):
            a = ProteinAlias()
            a.protein = p
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
                g.proteins.add(p)

    def create_protein_family(self, family_name, indent, parent_family, level_family_counter):
        # find the parent family
        if indent == 0:
            try:
                ppf = ProteinFamily.objects.get(parent__isnull=True)
            except ProteinFamily.DoesNotExist:
                raise Exception('Family 000 not found, aborting')
        else:
            parent_family_id = parent_family[indent-1]
            try:
                ppf = ProteinFamily.objects.get(pk=parent_family_id)
            except ProteinFamily.DoesNotExist:
                self.logger.error('Parent family of ' + family_name + ' not found, skipping')
                return False

        # does this family already exists in db?
        try:
            pf = ProteinFamily.objects.get(name=family_name, parent=ppf)
        except ProteinFamily.DoesNotExist:
            # increment the family counter for the current indent level
            level_family_counter[indent] += 1

            # protein family slug
            family_slug = []
            for level in level_family_counter:
                family_slug.append(str(level).zfill(3))
            family_slug = '_'.join(family_slug)

            # create the protein family
            pf = ProteinFamily()
            pf.parent = ppf
            pf.slug = family_slug
            pf.name = family_name
            try:
                pf.save()
                parent_family[indent] = pf.id

                self.logger.info('Created protein family ' + family_name)
            except Exception as msg:
                self.logger.error('Failed creating protein family' + family_name,msg)

        return {
            'pf': pf,
            'parent_family': parent_family,
            'level_family_counter': level_family_counter,
        }

    def parse_uniprot_file(self, accession):
        filename = accession + '.txt'
        local_file_path = os.sep.join([self.local_uniprot_dir, filename])
        remote_file_path = self.remote_uniprot_dir + filename

        up = {}
        up['genes'] = []
        up['names'] = []
        up['accessions'] = []


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

                # accessions
                elif line.startswith('AC'):
                    sline = line.split()
                    for ac in sline:
                        up['accessions'].append(ac.strip(';'))

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
            try:
                up['sequence'] = self.excel_sequences[up['entry_name']]['Sequence']
            except:
                pass
        except:
            return False

        # close the local file if appropriate
        if local_file:
            local_file.close()

        return up
