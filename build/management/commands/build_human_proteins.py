from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import connection
from django.db import IntegrityError

from build.management.commands.base_build import Command as BaseBuild
from protein.models import (Protein, ProteinConformation, ProteinState, ProteinFamily, ProteinAlias,
        ProteinSequenceType, Species, Gene, ProteinSource, ProteinSegment)

from residue.models import (ResidueNumberingScheme, ResidueGenericNumber, Residue)

import shlex
import os
from urllib.request import urlopen

import pandas as pd
import numpy  as np
import math


class Command(BaseBuild):
    help = 'Reads source data and creates protein families, proteins, and associated tables'

    protein_source_file = os.sep.join([settings.DATA_DIR, 'protein_data', 'proteins_and_families.txt'])
    local_uniprot_dir = os.sep.join([settings.DATA_DIR, 'protein_data', 'uniprot'])
    remote_uniprot_dir = 'http://www.uniprot.org/uniprot/'

    gprotein_data_file = os.sep.join([settings.DATA_DIR, 'g_protein_data', 'PDB_UNIPROT_ENSEMBLE_ALL.txt'])

    def handle(self, *args, **options):
        # use a smaller protein file if in test mode
        if options['test']:
            self.protein_source_file = os.sep.join([settings.DATA_DIR, 'protein_data',
                'proteins_and_families_test.txt'])

        # create parent protein family, 000
        try:
            self.create_parent_protein_family()
        except Exception as msg:
            # print(msg)
            self.logger.error(msg)

        # create proteins and families
        try:
            self.create_proteins_and_families()
        except Exception as msg:
            print(msg)
            self.logger.error(msg)

        #add gproteins from cgn db
        try:
            self.cgn_create_proteins_and_families()

            #delete added g-proteins
            #self.purge_cgn_proteins()
        except Exception as msg:
            print(msg)
            self.logger.error(msg)

        #add residues from cgn db
        try:
            self.update_protein_conformation()

        except Exception as msg:
            print(msg)
            self.logger.error(msg)


    def purge_cgn_proteins(self):
        try:
            Protein.objects.filter(prot_type='purge').delete()
        except:
            self.logger.info('Protein to delete not found')

    def add_cgn_residues(self, residue_generic_numbers_list):

        #gproteins list (lower case)
        gprotein_list=['gnaz_human','gnat3_human', 'gnat2_human', 'gnat1_human', 'gnas2_human', 'gnaq_human', 'gnao_human', 'gnal_human', 'gnai3_human', 'gnai2_human','gnai1_human', 'gna15_human', 'gna14_human', 'gna12_human', 'gna11_human']
        i=0

        for gp in gprotein_list:
            gprotein_list[i] = gp.upper()
            i=i+1

        #Parsing pdb uniprot file for residues
        self.logger.info('Start parsing PDB_UNIPROT_ENSEMBLE_ALL')
        self.logger.info('Parsing file ' + self.gprotein_data_file)
        residue_data =  pd.read_table(self.gprotein_data_file, sep="\t", low_memory=False)
        residue_data = residue_data.loc[residue_data['Uniprot_ID'].isin(gprotein_list)]


        for index, row in residue_data.iterrows():
            #fetch protein for protein conformation
            pr, c= Protein.objects.get_or_create(entry_name=row['Uniprot_ID'].lower())

            #fetch protein conformation
            pc, c= ProteinConformation.objects.get_or_create(protein_id=pr)

            #fetch residue generic number
            rgn, c= ResidueGenericNumber.objects.get_or_create(label=row['CGN'])

            #fetch protein segment id
            ps, c= ProteinSegment.objects.get_or_create(slug=row['CGN'].split(".")[1])

            try:
                Residue.objects.get_or_create(sequence_number=row['Position'], protein_conformation=pc, amino_acid=row['Residue'], generic_number=rgn, display_generic_number=rgn, protein_segment=ps)
                self.logger.info("Residues added to db")

            except:
                self.logger.error("Failed to add residues")
            


    def update_protein_conformation(self):
        gprotein_list=['gnaz_human','gnat3_human', 'gnat2_human', 'gnat1_human', 'gnas2_human', 'gnaq_human', 'gnao_human', 'gnal_human', 'gnai3_human', 'gnai2_human','gnai1_human', 'gna15_human', 'gna14_human', 'gna12_human', 'gna11_human']
        state = ProteinState.objects.get(slug='active')

        #add new cgn protein conformations
        for g in gprotein_list:
            gp = Protein.objects.get(entry_name=g)

            try:
                pc, created= ProteinConformation.objects.get_or_create(protein=gp, state=state, template_structure=None)
                self.logger.info('Created protein conformation')
            except:
                self.logger.error('Failed to create protein conformation')

        self.update_genericresiduenumber_and_proteinsegments()


    def update_genericresiduenumber_and_proteinsegments(self):

        gprotein_list=['gnaz_human','gnat3_human', 'gnat2_human', 'gnat1_human', 'gnas2_human', 'gnaq_human', 'gnao_human', 'gnal_human', 'gnai3_human', 'gnai2_human','gnai1_human', 'gna15_human', 'gna14_human', 'gna12_human', 'gna11_human']
        i=0

        for gp in gprotein_list:
            gprotein_list[i] = gp.upper()
            i=i+1

        #Parsing pdb uniprot file for generic residue numbers
        self.logger.info('Start parsing PDB_UNIPROT_ENSEMBLE_ALL')
        self.logger.info('Parsing file ' + self.gprotein_data_file)
        residue_data =  pd.read_table(self.gprotein_data_file, sep="\t", low_memory=False)

        #filtering for human gproteins using list above
        residue_data = residue_data.loc[residue_data['Uniprot_ID'].isin(gprotein_list)]
        residue_generic_numbers= residue_data['CGN']

        #add protein segment entries:

        segments =[]
        cgns = residue_data['CGN'].unique()

        for s in cgns:
            segments.append(s.split(".")[1])


        #Commit protein segments in db

        #purge line
        #ProteinSegment.objects.filter(slug=np.unique(segments)).delete()

        for s in np.unique(segments):

            if(s.startswith('s') or s.startswith('S')):
                category = 'sheet'
            else:
                category = 'helix'

            try:
                ProteinSegment.objects.get_or_create(slug=s, name=s, category=category, fully_aligned=True)
                self.logger.info('Created protein segment')

            except:
                self.logger.error('Failed to create protein segment')

        #Residue numbering scheme is the same for all added residue generic numbers (CGN)

        cgn_scheme = ResidueNumberingScheme.objects.get(slug='cgn')

        #purge line
        ResidueGenericNumber.objects.filter(scheme_id=12).delete()

        for rgn in residue_generic_numbers.unique():

            ps, c= ProteinSegment.objects.get_or_create(slug=rgn.split('.')[1])

            try:
                res_gen_num, created= ResidueGenericNumber.objects.get_or_create(label=rgn, scheme=cgn_scheme, protein_segment=ps)
                self.logger.info('Created generic residue number')

            except:
                self.logger.error('Failed creating generic residue number')


        self.add_cgn_residues(residue_generic_numbers)


    def cgn_add_proteins(self):

        self.logger.info('Start parsing PDB_UNIPROT_ENSEMBLE_ALL')
        self.logger.info('Parsing file ' + self.gprotein_data_file)

        #parsing file for accessions
        df =  pd.read_table(self.gprotein_data_file, sep="\t", low_memory=False)
        prot_type = 'purge'
        pfm = ProteinFamily()

        #Human proteins from CGN with families as keys: http://www.mrc-lmb.cam.ac.uk/CGN/about.html
        cgn_dict = {}
        cgn_dict['G-Protein']=['Gs', 'Gi', 'Gq', 'G12']
        cgn_dict['100_000_001']=['GNAS2_HUMAN', 'GNAL_HUMAN']
        cgn_dict['100_000_003']=['GNAQ_HUMAN', 'GNA11_HUMAN', 'GNA14_HUMAN', 'GNA15_HUMAN']
        cgn_dict['100_000_002']=['GNAI2_HUMAN', 'GNAI1_HUMAN', 'GNAI3_HUMAN', 'GNAT2_HUMAN', 'GNAT1_HUMAN', 'GNAT3_HUMAN', 'GNAZ_HUMAN', 'GNAO_HUMAN' ]
        cgn_dict['100_000_004']=['GNA12_HUMAN', 'GNA13_HUMAN']

        #list of all 16 proteins
        cgn_proteins_list=[]
        for k in cgn_dict.keys():
            for p in cgn_dict[k]:
                if p.endswith('_HUMAN'):
                    cgn_proteins_list.append(p)

        #print(cgn_proteins_list)

        #GNA13_HUMAN missing from cambridge file
        accessions= df.loc[df['Uniprot_ID'].isin(cgn_proteins_list)]
        accessions= accessions['Uniprot_ACC'].unique()

        #Create new residue numbering scheme
        self.create_cgn_rns()

        #purging one cgn entry
        #ResidueNumberingScheme.objects.filter(name='cgn').delete()

        rns = ResidueNumberingScheme.objects.get(slug='cgn')

        for a in accessions:
            up = self.parse_uniprot_file(a)

            #Fetch Protein Family for gproteins
            for k in cgn_dict.keys():
                name=up['entry_name'].upper()

                if name in cgn_dict[k]:
                    pfm = ProteinFamily.objects.get(slug=k)

            #Create new Protein
            self.cgn_creat_gproteins(pfm, rns, a, up, 'purge')



    def cgn_creat_gproteins(self, family, residue_numbering_scheme, accession, uniprot, prot_type):

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
        p.prot_type = 'purge'
        p.sequence = uniprot['sequence']

        try:
            p.save()
            self.logger.info('Created protein {}'.format(p.entry_name))
        except:
            self.logger.error('Failed creating protein {}'.format(p.entry_name))


    def cgn_parent_protein_family(self):

        pf_cgn, created_pf = ProteinFamily.objects.get_or_create(slug='100', defaults={
            'name': 'G-Protein'})

        pff_cgn = ProteinFamily.objects.get(slug='100', name='G-Protein')
        pf1_cgn = ProteinFamily.objects.get_or_create(slug='100_000', name='No Ligands', parent=pff_cgn)

    def create_cgn_rns(self):
        rns_cgn, created= ResidueNumberingScheme.objects.get_or_create(slug='cgn', short_name='CGN', defaults={
            'name': 'Common G-alpha numbering scheme'})

    
    def cgn_create_proteins_and_families(self):

        #Creating single entries in "protein_family' table
        self.cgn_parent_protein_family()

        #Human proteins from CGN: http://www.mrc-lmb.cam.ac.uk/CGN/about.html
        cgn_dict = {}

        levels = ['2', '3']
        keys = ['Gprotein', 'Gs', 'Gi', 'Gq', 'G12', '000']
        slug1='100'
        slug3= ''
        i=1

        cgn_dict['Gprotein']=['000']
        cgn_dict['000']=['Gs', 'Gi', 'Gq', 'G12']
        
        #Protein families to be added
        #Key of dictionary is level in hierarchy
        cgn_dict['1']=['Gprotein']
        cgn_dict['2']=['000']
        cgn_dict['3']=['Gs', 'Gi', 'Gq', 'G12']

        #Protein lines not to be added to Protein families
        cgn_dict['4']=['GNAS2', 'GNAL', 'GNAI1', 'GNAI1', 'GNAI3', 'GNAT2', 'GNAT1', 'GNAT3', 'GNAO', 'GNAZ', 'GNAQ', 'GNA11', 'GNA14', 'GNA15', 'GNA12', 'GNA13']

        for entry in cgn_dict['000']:

            name = entry

            slug2= '_000'
            slug3= '_00' + str(i)

            slug = slug1 + slug2 + slug3

            slug3 = ''
            i = i+1

            pff_cgn = ProteinFamily.objects.get(slug='100_000')

            new_pf, created = ProteinFamily.objects.get_or_create(slug=slug, name=entry, parent=pff_cgn)

        #function to create necessary arguments to add protein entry
        self.cgn_add_proteins()

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
        except:
            self.logger.error('Failed creating protein {}'.format(p.entry_name))

        # protein conformations
        try:
            ps, created = ProteinState.objects.get_or_create(slug=settings.DEFAULT_PROTEIN_STATE,
                defaults={'name': settings.DEFAULT_PROTEIN_STATE.title()})
        except IntegrityError:
            ps = ProteinState.objects.get(slug=settings.DEFAULT_PROTEIN_STATE)

        pc = ProteinConformation.objects.create(protein=p, state=ps)

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
            except:
                self.logger.error('Failed creating protein family' + family_name)

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