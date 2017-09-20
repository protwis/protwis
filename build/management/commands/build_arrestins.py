from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import connection
from django.db import IntegrityError


from common.models import WebResource, WebLink

from protein.models import (Protein, ProteinGProtein,ProteinGProteinPair, ProteinConformation, ProteinState, ProteinFamily, ProteinAlias,
        ProteinSequenceType, Species, Gene, ProteinSource, ProteinSegment)

from residue.models import (ResidueNumberingScheme, ResidueGenericNumber, Residue, ResidueGenericNumberEquivalent)

from signprot.models import SignprotStructure, SignprotBarcode
import pandas as pd

from optparse import make_option

import pandas as pd
import math
import numpy as np
import logging
import csv
import os

from urllib.request import urlopen

class Command(BaseCommand):
    help = 'Build Arrestin proteins'

    # source file directory
    arrestin_data_file = os.sep.join([settings.DATA_DIR, 'arrestin_data', 'ArrestinLookup.txt'])

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

        #add protein
        try:
            self.purge_can_residues()
            self.purge_can_proteins()

            self.create_arrestins(filenames)
            self.can_create_proteins_and_families()

        except Exception as msg:
            print(msg)
            self.logger.error(msg)

        # add residues from can db
        try:
            human_and_orths = self.cgn_add_proteins()

            self.update_protein_conformation(human_and_orths)
        except Exception as msg:
            self.logger.error(msg)

    def purge_can_residues(self):
        try:
            Residue.objects.filter(generic_number_id__scheme__slug="can").delete()
        except:
            self.logger.warning('Existing Residue data cannot be deleted')


    def create_arrestins(self, filenames=False):
        self.logger.info('CREATING ARRESTINS')

        # 4 arrestin subtypes, two of which are primarily expressed in the retina and bind only to visual opsins (arrestin 1 and arrestin 4), while the other two (β-arrestin 1 and β-arrestin 2) interact with the remaining ~800 GPCRs
        translation = {'Beta':'200_000_001', 'Visual':'200_000_002', }

        # read source files
        if not filenames:
            filenames = [fn for fn in os.listdir(self.arrestin_data_file) if fn.endswith('ArrestinLookup.csv')]

        for filename in filenames:
            filepath = os.sep.join([self.arrestin_data_file, filename])

            self.logger.info('Reading filename' + filename)

            with open(filepath, 'r') as f:
                reader = csv.reader(f)
                for row in reader:

                    ## TO BE UPDATED! PANDAS?
                    entry_name = row[4]
                    primary = row[8]
                    secondary = row[9]

                    # fetch protein
                    try:
                        p = Protein.objects.get(entry_name=entry_name)
                    except Protein.DoesNotExist:
                        self.logger.warning('Protein not found for entry_name {}'.format(entry_name))
                        continue

                    primary = primary.replace("G protein (identity unknown)","None") #replace none
                    primary = primary.split(", ")

                    secondary = secondary.replace("G protein (identity unknown)","None") #replace none
                    secondary = secondary.split(", ")

                    if primary=='None' and secondary=='None':
                        print('no data for ',entry_name)
                        continue

                    for gp in primary:
                        if gp in ['None','_-arrestin','Arrestin','G protein independent mechanism']: #skip bad ones
                            continue
                        g = ProteinGProtein.objects.get_or_create(name=gp, slug=translation[gp])[0]
                        gpair = ProteinGProteinPair(protein=p, g_protein=g, transduction='primary')
                        gpair.save()

                    for gp in secondary:
                        if gp in ['None','_-arrestin','Arrestin','G protein independent mechanism', '']: #skip bad ones
                            continue
                        if gp in primary: #sip those that were already primary
                             continue
                        g = ProteinGProtein.objects.get_or_create(name=gp, slug=translation[gp])[0]
                        gpair = ProteinGProteinPair(protein=p, g_protein=g, transduction='secondary')
                        gpair.save()


        self.logger.info('COMPLETED CREATING G PROTEINS')

    def purge_can_proteins(self):
        try:
            Protein.objects.filter(residue_numbering_scheme_id=13).delete()
        except:
            self.logger.info('Protein to delete not found')

    def add_cgn_residues(self, gprotein_list):

        #Parsing pdb uniprot file for residues
        self.logger.info('Start parsing PDB_UNIPROT_ENSEMBLE_ALL')
        self.logger.info('Parsing file ' + self.gprotein_data_file)
        residue_data =  pd.read_table(self.gprotein_data_file, sep="\t")
        residue_data = residue_data.loc[residue_data['Uniprot_ACC'].isin(gprotein_list)]

        for index, row in residue_data.iterrows():
            #fetch protein for protein conformation
            pr, c= Protein.objects.get_or_create(accession=row['Uniprot_ACC'])

            #fetch protein conformation
            pc, c= ProteinConformation.objects.get_or_create(protein_id=pr)

            #fetch residue generic number
            rgnsp=[]
            if(int(row['CGN'].split('.')[2])<10):
                rgnsp = row['CGN'].split('.')
                rgn_new = rgnsp[0]+'.'+rgnsp[1]+'.0'+rgnsp[2]
                rgn, c= ResidueGenericNumber.objects.get_or_create(label=rgn_new)

            else:
                rgn, c= ResidueGenericNumber.objects.get_or_create(label=row['CGN'])

            #fetch protein segment id
            ps, c= ProteinSegment.objects.get_or_create(slug=row['CGN'].split(".")[1])

            try:
                Residue.objects.get_or_create(sequence_number=row['Position'], protein_conformation=pc, amino_acid=row['Residue'], generic_number=rgn, display_generic_number=rgn, protein_segment=ps)
                # self.logger.info("Residues added to db")

            except:
                self.logger.error("Failed to add residues")


             # Add also to the ResidueGenericNumberEquivalent table needed for single residue selection
            try:
                ResidueGenericNumberEquivalent.objects.get_or_create(label=rgn.label,default_generic_number=rgn, scheme_id=12)
                # self.logger.info("Residues added to ResidueGenericNumberEquivalent")

            except:
                self.logger.error("Failed to add residues to ResidueGenericNumberEquivalent")

    def update_protein_conformation(self, gprotein_list):
        #gprotein_list=['gnaz_human','gnat3_human', 'gnat2_human', 'gnat1_human', 'gnas2_human', 'gnaq_human', 'gnao_human', 'gnal_human', 'gnai3_human', 'gnai2_human','gnai1_human', 'gna15_human', 'gna14_human', 'gna12_human', 'gna11_human', 'gna13_human']
        state = ProteinState.objects.get(slug='active')

        #add new cgn protein conformations
        for g in gprotein_list:
            gp = Protein.objects.get(accession=g)

            try:
                pc, created= ProteinConformation.objects.get_or_create(protein=gp, state=state, template_structure=None)
                self.logger.info('Created protein conformation')
            except:
                self.logger.error('Failed to create protein conformation')

        self.update_genericresiduenumber_and_proteinsegments(gprotein_list)

    def update_genericresiduenumber_and_proteinsegments(self, gprotein_list):

        #Parsing pdb uniprot file for generic residue numbers
        self.logger.info('Start parsing PDB_UNIPROT_ENSEMBLE_ALL')
        self.logger.info('Parsing file ' + self.gprotein_data_file)
        residue_data =  pd.read_table(self.gprotein_data_file, sep="\t")


        residue_data = residue_data[residue_data.Uniprot_ID.notnull()]

        #residue_data = residue_data[residue_data['Uniprot_ID'].str.contains('_HUMAN')]

        residue_data = residue_data[residue_data['Uniprot_ACC'].isin(gprotein_list)]

        #filtering for human gproteins using list above
        residue_generic_numbers= residue_data['CGN']

        #add protein segment entries:

        segments =[]
        cgns = residue_data['CGN'].unique()

        for s in cgns:
            segments.append(s.split(".")[1])

        #Commit protein segments in db

        #purge line
        #ProteinSegment.objects.filter(slug__in=np.unique(segments)).delete()

        for s in np.unique(segments):

            if s.startswith('S') and len(s) == 2:
                category = 'sheet'
            elif s.startswith('H') and len(s) == 2:
                category = 'helix'
            else:
                category = 'loop'

            try:
                ProteinSegment.objects.get_or_create(slug=s, name=s, category=category, fully_aligned=True)
                self.logger.info('Created protein segment')

            except:
                self.logger.error('Failed to create protein segment')

        #Residue numbering scheme is the same for all added residue generic numbers (CGN)

        cgn_scheme = ResidueNumberingScheme.objects.get(slug='cgn')

        #purge line
        #ResidueGenericNumber.objects.filter(scheme_id=12).delete()

        for rgn in residue_generic_numbers.unique():
            ps, c= ProteinSegment.objects.get_or_create(slug=rgn.split('.')[1])

            rgnsp=[]

            if(int(rgn.split('.')[2])<10):
                rgnsp = rgn.split('.')
                rgn_new = rgnsp[0]+'.'+rgnsp[1]+'.0'+rgnsp[2]
            else:
                rgn_new = rgn

            try:
                res_gen_num, created= ResidueGenericNumber.objects.get_or_create(label=rgn_new, scheme=cgn_scheme, protein_segment=ps)
                self.logger.info('Created generic residue number')

            except:
                self.logger.error('Failed creating generic residue number')


        self.add_cgn_residues(gprotein_list)

    def cgn_add_proteins(self):

        self.logger.info('Start parsing PDB_UNIPROT_ENSEMBLE_ALL')
        self.logger.info('Parsing file ' + self.gprotein_data_file)

        #parsing file for accessions
        df =  pd.read_table(self.gprotein_data_file, sep="\t")
        prot_type = 'purge'
        pfm = ProteinFamily()

        #Human proteins from CGN with families as keys: http://www.mrc-lmb.cam.ac.uk/CGN/about.html
        cgn_dict = {}
        cgn_dict['G-Protein']=['Gs', 'Gi/o', 'Gq/11', 'G12/13']
        cgn_dict['100_000_001']=['GNAS2_HUMAN', 'GNAL_HUMAN']
        cgn_dict['100_000_002']=['GNAI2_HUMAN', 'GNAI1_HUMAN', 'GNAI3_HUMAN', 'GNAT2_HUMAN', 'GNAT1_HUMAN', 'GNAT3_HUMAN', 'GNAZ_HUMAN', 'GNAO_HUMAN' ]
        cgn_dict['100_000_003']=['GNAQ_HUMAN', 'GNA11_HUMAN', 'GNA14_HUMAN', 'GNA15_HUMAN']
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
                name=str(up['entry_name']).upper()

                if name in cgn_dict[k]:
                    pfm = ProteinFamily.objects.get(slug=k)

            #Create new Protein
            self.cgn_creat_gproteins(pfm, rns, a, up)

        ###################ORTHOLOGS###############
        orthologs_pairs =[]
        orthologs =[]

        #Orthologs for human gproteins
        allprots = list(df.Uniprot_ID.unique())
        allprots = list(set(allprots) - set(cgn_proteins_list))

        for gp in cgn_proteins_list:
            for p in allprots:
                if str(p).startswith(gp.split('_')[0]):
                    orthologs_pairs.append((str(p), gp))
                    orthologs.append(str(p))

        accessions_orth= df.loc[df['Uniprot_ID'].isin(orthologs)]
        accessions_orth= accessions_orth['Uniprot_ACC'].unique()


        for a in accessions_orth:
            up = self.parse_uniprot_file(a)

            #Fetch Protein Family for gproteins
            for k in cgn_dict.keys():
                name=str(up['entry_name']).upper()
                name = name.split('_')[0]+'_'+'HUMAN'

                if name in cgn_dict[k]:
                    pfm = ProteinFamily.objects.get(slug=k)

            #Create new Protein
            self.cgn_creat_gproteins(pfm, rns, a, up)

        #human gproteins
        orthologs_lower = [x.lower() for x in orthologs]
        #print(orthologs_lower)

        #orthologs to human gproteins
        cgn_proteins_list_lower = [x.lower() for x in cgn_proteins_list]

        #all gproteins
        gprotein_list = cgn_proteins_list_lower + orthologs_lower
        accessions_all = list(accessions_orth) + list(accessions)

        return list(accessions_all)

    def can_creat_gproteins(self, family, residue_numbering_scheme, accession, uniprot):

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
            'name': 'Arrestin'})

        pff_can = ProteinFamily.objects.get(slug='200', name='Arrestin')

        #Changed name "No Ligands" to "Gprotein"
        pf1_can = ProteinFamily.objects.get_or_create(slug='200_000', name='Arrestin', parent=pff_can)

    def create_can_rns(self):
        ## New numbering scheme entry_name

        rns_can, created= ResidueNumberingScheme.objects.get_or_create(slug='can', short_name='CAN', defaults={
            'name': 'Common arrestin numbering scheme'})

    def can_create_proteins_and_families(self):

        #Creating single entries in "protein_family' table
        ProteinFamily.objects.filter(slug__startswith="200").delete()
        self.can_parent_protein_family()

        can_dict = {}

        levels = ['2', '3']
        keys = ['Beta','Visual']
        slug1='200'
        slug3= ''
        i=1

        can_dict['Arrestin']=['000']
        can_dict['000']=['Gs', 'Gi/o', 'Gq/11', 'G12/13']

        # Protein families to be added
        # Key of dictionary is level in hierarchy
        can_dict['1']=['Arrestin']
        can_dict['2']=['000']
        can_dict['3']=['Beta','Visual']

        # Protein lines not to be added to Protein families
        can_dict['4']=['ARRB2','ARRB1','ARRS','ARRC']

        for entry in can_dict['000']:

            name = entry

            slug2= '_000'
            slug3= '_00' + str(i)

            slug = slug1 + slug2 + slug3

            slug3 = ''
            i = i+1

            pff_can = ProteinFamily.objects.get(slug='200_000')

            new_pf, created = ProteinFamily.objects.get_or_create(slug=slug, name=entry, parent=pff_can)

        #function to create necessary arguments to add protein entry
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
