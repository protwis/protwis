from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import connection

from protein.models import Protein
from protein.models import ProteinFamily
from protein.models import ProteinAlias
from protein.models import ProteinSegment
from protein.models import Species
from protein.models import Gene
from protein.models import ProteinSource

import logging
import shlex


class Command(BaseCommand):
    help = 'Reads source data and creates protein families, proteins, and associated tables'

    logger = logging.getLogger(__name__)

    protein_source_file = settings.DATA_DIR + '/protein_data/proteins_and_families.txt'
    segment_source_file = settings.DATA_DIR + '/protein_data/segments.txt'

    def handle(self, *args, **options):
        # delete any existing protein data
        try:
            self.truncate_protein_tables()
        except Exception as msg:
            # print(msg)
            self.logger.error(msg)

        # create parent protein family, 000
        try:
            self.create_parent_protein_family()
        except Exception as msg:
            # print(msg)
            self.logger.error(msg)

        # create protein segments
        try:
            self.create_protein_segments()
        except Exception as msg:
            print(msg)
            self.logger.error(msg)

        # create proteins and families
        try:
            self.create_proteins_and_families()
        except Exception as msg:
            print(msg)
            self.logger.error(msg)

    def truncate_protein_tables(self):
        cursor = connection.cursor()
        
        tables_to_truncate = [
            'protein_alias',
            'protein_family',
            'protein_links',
            'protein_reosurce', # ಠ_ಠ
            'protein_segment',
            'protein_set',
            'protein_proteinset_protein',
            'protein_source',
            'gene',
            'protein_gene_proteins',
            'species',
            'protein',
        ]

        for table in tables_to_truncate:
            cursor.execute("TRUNCATE TABLE " + table + " CASCADE")

    def create_parent_protein_family(self):
        pf = ProteinFamily()
        pf.slug = '000'
        pf.name = 'Parent family'
        pf.save()

    def create_protein_segments(self):
        self.logger.info('Parsing file ' + self.segment_source_file)
        self.logger.info('CREATING PROTEIN SEGMENTS')

        with open(self.segment_source_file, "r", encoding='UTF-8') as segment_file:
            for i, row in enumerate(segment_file):
                split_row = shlex.split(row)

                # create segment
                s = ProteinSegment()
                s.slug = split_row[0]
                s.category = split_row[1]
                s.name = split_row[2]
                s.position = i

                try:
                    s.save()
                    self.logger.info('Created protein segment ' + s.name)
                except:
                    self.logger.error('Failed creating protein segment ' + s.name)
                    continue

        self.logger.info('COMPLETED CREATING PROTEIN SEGMENTS')


    def create_proteins_and_families(self):
        self.logger.info('Parsing file ' + self.protein_source_file)
        self.logger.info('CREATING PROTEINS')

        with open(self.protein_source_file, "r", encoding='UTF-8') as protein_file:
            # family hierarchy is determined by indent
            spaces_per_indent_level = 4
            last_indent = 0
            level_family_counter = [0]
            parent_family = [0]

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
                # family
                ###########
                create_protein = False
                if row.strip().startswith('"'): # protein row
                    split_row = row.strip().split('","')
                    family_name = split_row[4]
                    create_protein = True
                else: # family row
                    family_name = row.strip()

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
                    accessions = [split_row[15], split_row[31], split_row[23]]

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
                            self.logger.error('Protein with accession ' + protein_accession + ' already exists, skipping')
                            continue

                        # parse uniprot file for this protein
                        self.logger.info('Parsing uniprot file for protein ' + protein_name)
                        uniprot_file = settings.DATA_DIR + '/uniprot/txt/' + protein_accession + '.txt'
                        up = self.parse_uniprot_file(uniprot_file)
                        if not up:
                            self.logger.error('Failed parsing uniprot file for protein ' + protein_name + ', skipping')
                            continue

                        # get/create protein source
                        try:
                            source = ProteinSource.objects.get(name=up['source'])
                        except ProteinSource.DoesNotExist:
                            source = ProteinSource()
                            source.name = up['source']
                            
                            try:
                                source.save()
                                self.logger.info('Created protein source ' + source.name)
                            except:
                                self.logger.error('Failed creating protein source ' + source.name)

                        # get/create species
                        try:
                            species = Species.objects.get(latin_name=up['species_latin_name'])
                        except Species.DoesNotExist:
                            species = Species()
                            species.latin_name = up['species_latin_name']
                            species.common_name = up['species_common_name']

                            try:
                                species.save()
                                self.logger.info('Created species ' + species.latin_name)
                            except:
                                self.logger.error('Failed creating species ' + species.latin_name)

                        # create protein
                        p = Protein()
                        p.family = pf
                        p.species = species
                        p.source = source
                        p.accession = protein_accession
                        p.entry_name = up['entry_name']
                        p.name = protein_name
                        p.sequence = up['sequence']

                        try:
                            p.save()
                            self.logger.info('Created protein ' + p.name)
                        except:
                            self.logger.error('Failed creating protein ' + p.name)

                        # protein aliases
                        for i, alias in enumerate(up['names']):
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
                        for i, gene in enumerate(up['genes']):
                            g = Gene()
                            g.species = species
                            g.name = gene
                            g.position = i

                            try:
                                g.save()
                                g.proteins.add(p)
                                self.logger.info('Created gene ' + g.name + ' for protein ' + p.name)
                            except:
                                self.logger.error('Failed creating gene ' + g.name + ' for protein ' + p.name)

        self.logger.info('COMPLETED CREATING PROTEINS')

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
            pf = ProteinFamily.objects.get(name=family_name)
        except ProteinFamily.DoesNotExist:
            # increment the famliy counter for the current indent level
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
            except:
                self.logger.error('Failed creating protein family' + family_name)

        return {
            'pf': pf,
            'parent_family': parent_family,
            'level_family_counter': level_family_counter,
        }

    def parse_uniprot_file(self, file_path):
        up = {}
        up['genes'] = []
        up['names'] = []
        read_sequence = False
        try:
            with open(file_path) as uf:
                for line in uf:
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
                    elif line.startswith('OS'):
                        species_full = line[2:].strip().strip('.')
                        species_split = species_full.split('(')
                        up['species_latin_name'] = species_split[0]
                        if len(species_split) > 1:
                            up['species_common_name'] = species_split[1].strip(')')
                        else:
                            up['species_common_name'] = up['species_latin_name']

                    # names
                    elif line.startswith('DE'):
                        split_de_line = line.split('=')
                        if len(split_de_line) > 1:
                            split_segment = split_de_line[1].split('{')
                            up['names'].append(split_segment[0].strip())

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
                        seq_len = split_sq_line[2]
                        read_sequence = True
                        up['sequence'] = ''
                    elif line.startswith('//'):
                        read_sequence = False
                        if len(sequence) != seq_len:
                            self.logger.error('Sequence length does not match specification')
                    elif read_sequence == True:
                        up['sequence'] += line.strip().replace(' ', '')
        except:
            return False

        return up
