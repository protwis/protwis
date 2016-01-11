from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import connection
from django.utils.text import slugify
from django.db import IntegrityError

from build.management.commands.base_build import Command as BaseBuild
from protein.models import (Protein, ProteinConformation, ProteinState, ProteinAnomaly, ProteinAnomalyType,
    ProteinSegment)
from residue.models import ResidueGenericNumber, ResidueNumberingScheme, Residue
from common.models import WebLink, WebResource, Publication
from structure.models import (Structure, StructureType, StructureSegment, StructureStabilizingAgent,PdbData,
    Rotamer, StructureSegmentModeling, StructureCoordinates, StructureCoordinatesDescription, StructureEngineering,
    StructureEngineeringDescription)
from ligand.models import Ligand, LigandType, LigandRole, LigandProperities
from interaction.models import *
from interaction.views import runcalculation,parsecalculation

import logging
import os
import re
import yaml
from collections import OrderedDict
import json
from urllib.request import urlopen
from Bio.PDB import parse_pdb_header


## FOR VIGNIR ORDERED DICT YAML IMPORT/DUMP
_mapping_tag = yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG
def dict_constructor(loader, node):
    return OrderedDict(loader.construct_pairs(node))

def represent_ordereddict(dumper, data):
    value = []

    for item_key, item_value in data.items():
        node_key = dumper.represent_data(item_key)
        node_value = dumper.represent_data(item_value)

        value.append((node_key, node_value))

    return yaml.nodes.MappingNode(u'tag:yaml.org,2002:map', value)

yaml.add_representer(OrderedDict, represent_ordereddict)
yaml.add_constructor(_mapping_tag, dict_constructor)

class Command(BaseBuild):
    help = 'Reads source data and creates pdb structure records'
    
    def add_arguments(self, parser):
        parser.add_argument('-p', '--proc',
            type=int,
            action='store',
            dest='proc',
            default=1,
            help='Number of processes to run')
        parser.add_argument('-f', '--filename',
            action='append',
            dest='filename',
            help='Filename to import. Can be used multiple times')
        parser.add_argument('-u', '--purge',
            action='store_true',
            dest='purge',
            default=False,
            help='Purge existing records')

    # source file directory
    structure_data_dir = os.sep.join([settings.DATA_DIR, 'structure_data', 'structures'])
    pdb_data_dir = os.sep.join([settings.DATA_DIR, 'structure_data', 'pdbs'])

    # read source files
    filenames = os.listdir(structure_data_dir)

    def handle(self, *args, **options):
        # delete any existing structure data
        if options['purge']:
            try:
                self.purge_structures()
            except Exception as msg:
                print(msg)
                self.logger.error(msg)

        # where filenames specified?
        if options['filename']:
            self.filenames = options['filename']

        try:
            self.logger.info('CREATING STRUCTURES')

            # run the function twice (once for representative structures, once for non-representative)
            iterations = 2
            for i in range(1,iterations+1):
                self.prepare_input(options['proc'], self.filenames, i)

            self.logger.info('COMPLETED CREATING STRUCTURES')
        except Exception as msg:
            print(msg)
            self.logger.error(msg)
    
    def purge_structures(self):
        Structure.objects.all().delete()

    def create_rotamers(self, structure):
        AA = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D',
     'CYS':'C', 'GLN':'Q', 'GLU':'E', 'GLY':'G',
     'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K',
     'MET':'M', 'PHE':'F', 'PRO':'P', 'SER':'S',
     'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V'}
        pdb = structure.pdb_data.pdb
        protein_conformation=structure.protein_conformation
        preferred_chain = structure.preferred_chain
        temp = ''
        check = 0
        errors = 0
        for line in pdb.splitlines():
            if line.startswith('ATOM'): #If it is a residue
                residue_number = line[22:26]
                chain = line[21]
                if preferred_chain and chain!=preferred_chain: #If perferred is defined and is not the same as the current line, then skip
                    continue
                preferred_chain = chain #If reached this point either the preferred_chain is specified or needs to be the first one and is thus specified now.
                if check==0 or residue_number==check: #If this is either the begining or the same as previous line add to current rotamer
                    temp += line + "\n"
                else: #if this is a new residue
                    try:
                        residue=Residue.objects.get(protein_conformation=protein_conformation, sequence_number=check)
                        if not residue_name==residue.three_letter():
                            #print('Residue not found in '+structure.pdb_code.index+', skipping! '+residue_name+' vs '+residue.three_letter()+ ' at position '+residue_number)
                            #errors += 1
                            self.logger.error("Changing WT residue amino_acid for sequence_number "+str(residue.sequence_number)+" from "+residue.amino_acid+" to "+AA[residue_name.upper()])
                            residue.amino_acid = AA[residue_name.upper()]
                            residue.save()
                        else:
                            rotamer_data, created = PdbData.objects.get_or_create(pdb=temp)
                            rotamer, created = Rotamer.objects.get_or_create(residue=residue, structure=structure, pdbdata=rotamer_data)
                    except Residue.DoesNotExist:
                        #print('Residue not found in '+structure.pdb_code.index+' ' +residue_name+' at position '+residue_number)
                        residue = None
                        errors += 1
                        #self.logger.error("No residue found for sequence_number "+str(check))
                    
                    temp = line + "\n"
                check = residue_number
                residue_name = line[17:20].title() #use title to get GLY to Gly so it matches
        try:
            residue=Residue.objects.get(protein_conformation=protein_conformation, sequence_number=check)
            if not residue_name==residue.three_letter():
                #print('Residue not found in '+structure.pdb_code.index+', skipping! '+residue_name+' vs '+residue.three_letter()+ ' at position '+residue_number)
                #prerrors += 1
                self.logger.error("Changing WT residue amino_acid for sequence_number "+str(residue.sequence_number)+" from "+residue.amino_acid+" to "+AA[residue_name.upper()])
                residue.amino_acid = AA[residue_name.upper()]
                residue.save()
            else:
                rotamer_data, created = PdbData.objects.get_or_create(pdb=temp)
                rotamer, created = Rotamer.objects.get_or_create(residue=residue, structure=structure, pdbdata=rotamer_data)
        except Residue.DoesNotExist:
            #print('Residue not found in '+structure.pdb_code.index+' ' +residue_name+' at position '+residue_number)
            residue = None
            errors += 1
            #self.logger.error("No residue found for sequence_number "+str(check))
        if errors:
            self.logger.error(structure.pdb_code.index + " had " + str(errors) + " residues that did not match in the database")
        return None

    def main_func(self, positions, iteration):
        # filenames
        if not positions[1]:
            filenames = self.filenames[positions[0]:]
        else:
            filenames = self.filenames[positions[0]:positions[1]]

        for source_file in filenames:
            source_file_path = os.sep.join([self.structure_data_dir, source_file])
            if os.path.isfile(source_file_path) and source_file[0] != '.':
                self.logger.info('Reading file {}'.format(source_file_path))
                # read the yaml file
                with open(source_file_path, 'r') as f:
                    sd = yaml.load(f)
                    
                    # is this a representative structure (will be used to guide structure-based alignments)?
                    representative = False
                    if 'representative' in sd and sd['representative']:
                        representative = True

                    # only process representative structures on first iteration
                    if not representative and iteration == 1:
                        continue

                    # skip representative structures on second iteration
                    if representative and iteration == 2:
                        continue

                    # is there a construct?
                    if 'construct' not in sd:
                        self.logger.error('No construct specified, skipping!')
                        continue

                    # does the construct exists?
                    try:
                        con = Protein.objects.get(entry_name=sd['construct'])
                    except Protein.DoesNotExist:
                        self.logger.error('Construct {} does not exists, skipping!'.format(sd['construct']))
                        continue

                    # create a structure record
                    try:
                        s = Structure.objects.get(protein_conformation__protein=con)
                    except Structure.DoesNotExist:
                        s = Structure()
                        s.representative = representative

                    # protein state
                    if 'state' not in sd:
                        self.logger.warning('State not defined, using default state {}'.format(
                            settings.DEFAULT_PROTEIN_STATE))
                        state = settings.DEFAULT_STATE.title()
                    else:
                        state = sd['state']
                    state_slug = slugify(state)
                    try:
                        ps, created = ProteinState.objects.get_or_create(slug=state_slug, defaults={'name': state})
                        if created:
                            self.logger.info('Created protein state {}'.format(ps.name))
                    except IntegrityError:
                        ps = ProteinState.objects.get(slug=state_slug)
                    s.state = ps

                    # protein conformation
                    try:
                        s.protein_conformation = ProteinConformation.objects.get(protein=con)
                    except ProteinConformation.DoesNotExist:
                        self.logger.error('Protein conformation for construct {} does not exists'.format(con))
                        continue
                    if s.protein_conformation.state is not state:
                        ProteinConformation.objects.filter(protein=con).update(state=ps)

                    # get the PDB file and save to DB
                    sd['pdb'] = sd['pdb'].upper()
                    if not os.path.exists(self.pdb_data_dir):
                        os.makedirs(self.pdb_data_dir)
                    
                    pdb_path = os.sep.join([self.pdb_data_dir, sd['pdb'] + '.pdb'])
                    if not os.path.isfile(pdb_path):
                        self.logger.info('Fetching PDB file {}'.format(sd['pdb']))
                        url = 'http://www.rcsb.org/pdb/files/%s.pdb' % sd['pdb']
                        pdbdata_raw = urlopen(url).read().decode('utf-8')
                        with open(pdb_path, 'w') as f:
                            f.write(pdbdata_raw)
                    else:
                        with open(pdb_path, 'r') as pdb_file:
                            pdbdata_raw = pdb_file.read()
                    
                    pdbdata, created = PdbData.objects.get_or_create(pdb=pdbdata_raw)
                    s.pdb_data = pdbdata

                    # UPDATE HETSYN with its PDB reference instead + GRAB PUB DATE, PMID, DOI AND RESOLUTION
                    hetsyn = {}
                    hetsyn_reverse = {}
                    for line in pdbdata_raw.splitlines():
                        if line.startswith('HETSYN'): 
                            m = re.match("HETSYN[\s]+([\w]{3})[\s]+(.+)",line) ### need to fix bad PDB formatting where col4 and col5 are put together for some reason -- usually seen when the id is +1000
                            if (m):
                                hetsyn[m.group(2).strip()] = m.group(1).upper()
                                hetsyn_reverse[m.group(1)] = m.group(2).strip().upper()
                        if line.startswith('HETNAM'): 
                            m = re.match("HETNAM[\s]+([\w]{3})[\s]+(.+)",line) ### need to fix bad PDB formatting where col4 and col5 are put together for some reason -- usually seen when the id is +1000
                            if (m):
                                hetsyn[m.group(2).strip()] = m.group(1).upper()
                                hetsyn_reverse[m.group(1)] = m.group(2).strip().upper()
                        if line.startswith('REVDAT   1'):
                            sd['publication_date'] = line[13:22]
                        if line.startswith('JRNL        PMID'):
                            sd['pubmed_id'] = line[19:].strip()
                        if line.startswith('JRNL        DOI'):
                            sd['doi_id'] = line[19:].strip()

                    if len(hetsyn) == 0:
                        self.logger.info("PDB file contained NO hetsyn")

                    with open(pdb_path,'r') as header:
                        header_dict = parse_pdb_header(header)
                    sd['publication_date'] = header_dict['release_date']
                    sd['resolution'] = str(header_dict['resolution']).strip()
                    sd['structure_method'] = header_dict['structure_method']

                    # structure type
                    if 'structure_method' in sd and sd['structure_method']:
                        structure_type = sd['structure_method'].capitalize()
                        structure_type_slug = slugify(sd['structure_method'])
                        
                        try:
                            st, created = StructureType.objects.get_or_create(slug=structure_type_slug,
                                defaults={'name': structure_type})
                            if created:
                                self.logger.info('Created structure type {}'.format(st))
                        except IntegrityError:
                            st = StructureType.objects.get(slug=structure_type_slug)
                        s.structure_type = st
                    else:
                        self.logger.warning('No structure type specified in PDB file {}'.format(sd['pdb']))

                    matched = 0
                    if 'ligand' in sd and sd['ligand']:
                        if isinstance(sd['ligand'], list):
                            ligands = sd['ligand']
                        else:
                            ligands = [sd['ligand']]
                        for ligand in ligands:
                            if 'name' in ligand:
                                if ligand['name'].upper() in hetsyn:
                                    self.logger.info('Ligand {} matched to PDB records'.format(ligand['name']))
                                    matched = 1
                                    ligand['name'] = hetsyn[ligand['name'].upper()]
                                elif ligand['name'].upper() in hetsyn_reverse:
                                    matched = 1

                    if matched==0 and len(hetsyn)>0:
                        self.logger.info('No ligand names found in HET in structure {}'.format(sd['pdb']))

                    # REMOVE? can be used to dump structure files with updated ligands
                    # yaml.dump(sd, open(source_file_path, 'w'), indent=4)

                    # pdb code
                    if 'pdb' in sd:
                        try:
                            web_resource = WebResource.objects.get(slug='pdb')
                        except:
                            # abort if pdb resource is not found
                            raise Exception('PDB resource not found, aborting!')
                        s.pdb_code, created = WebLink.objects.get_or_create(index=sd['pdb'],
                            web_resource=web_resource)
                    else:
                        self.logger.error('PDB code not specified for structure {}, skipping!'.format(sd['pdb']))
                        continue

                    # insert into plain text fields
                    if 'preferred_chain' in sd:
                        s.preferred_chain = sd['preferred_chain']
                    else:
                        self.logger.warning('Preferred chain not specified for structure {}'.format(sd['pdb']))
                    if 'resolution' in sd:
                        s.resolution = float(sd['resolution'])
                    else:
                        self.logger.warning('Resolution not specified for structure {}'.format(sd['pdb']))
                    if 'publication_date' in sd:
                        s.publication_date = sd['publication_date']
                    else:
                        self.logger.warning('Publication date not specified for structure {}'.format(sd['pdb']))

                    # publication
                    try:                     
                        if 'doi_id' in sd:
                            try:
                                s.publication = Publication.objects.get(web_link__index=sd['doi_id'])
                            except Publication.DoesNotExist as e:
                                p = Publication()
                                try:
                                    p.web_link = WebLink.objects.get(index=sd['doi_id'], web_resource__slug='doi')
                                except WebLink.DoesNotExist:
                                    wl = WebLink.objects.create(index=sd['doi_id'],
                                        web_resource = WebResource.objects.get(slug='doi'))
                                    p.web_link = wl
                                p.update_from_doi(doi=sd['doi_id'])
                                p.save()
                                s.publication = p
                        elif 'pubmed_id' in sd:
                            try:
                                s.publication = Publication.objects.get(web_link__index=sd['pubmed_id'])
                            except Publication.DoesNotExist as e:
                                p = Publication()
                                try:
                                    p.web_link = WebLink.objects.get(index=sd['pubmed_id'],
                                        web_resource__slug='pubmed')
                                except WebLink.DoesNotExist:
                                    wl = WebLink.objects.create(index=sd['pubmed_id'],
                                        web_resource = WebResource.objects.get(slug='pubmed'))
                                    p.web_link = wl
                                p.update_from_pubmed_data(index=sd['pubmed_id'])
                                p.save()
                                s.publication = p
                    except:
                        self.logger.error('Error saving publication'.format(ps.name))

                    # save structure before adding M2M relations
                    s.save()

                    # endogenous ligand(s)
                    default_ligand_type = 'Small molecule'
                    if representative and 'endogenous_ligand' in sd and sd['endogenous_ligand']:
                        if isinstance(sd['endogenous_ligand'], list):
                            endogenous_ligands = sd['endogenous_ligand']
                        else:
                            endogenous_ligands = [sd['endogenous_ligand']]
                        for endogenous_ligand in endogenous_ligands:
                            if endogenous_ligand['type']:
                                lt, created = LigandType.objects.get_or_create(slug=slugify(endogenous_ligand['type']),
                                    defaults={'name': endogenous_ligand['type']})
                            else:
                                lt, created = LigandType.objects.get_or_create(slug=slugify(default_ligand_type),
                                    defaults={'name': default_ligand_type})
                            ligand = Ligand()

                            if 'iupharId' not in endogenous_ligand:
                                endogenous_ligand['iupharId'] = 0

                            ligand = ligand.load_by_gtop_id(endogenous_ligand['name'], endogenous_ligand['iupharId'],
                                lt)
                            try:
                                s.protein_conformation.protein.parent.endogenous_ligands.add(ligand)
                            except IntegrityError:
                                self.logger.info('Endogenous ligand for protein {}, already added. Skipping.'.format(
                                    s.protein_conformation.protein.parent))

                    # ligands
                    if 'ligand' in sd and sd['ligand']:
                        if isinstance(sd['ligand'], list):
                            ligands = sd['ligand']
                        else:
                            ligands = [sd['ligand']]
                        for ligand in ligands:
                            l = False
                            if ligand['name'] and ligand['name'] != 'None': # some inserted as none.

                                # use annoted ligand type or default type
                                if ligand['type']:
                                    lt, created = LigandType.objects.get_or_create(slug=slugify(ligand['type']),
                                        defaults={'name': ligand['type']})
                                else:
                                    lt, created = LigandType.objects.get_or_create(
                                        slug=slugify(default_ligand_type), defaults={'name': default_ligand_type})

                                # set pdb reference for structure-ligand interaction
                                pdb_reference = ligand['name']

                                # use pubchem_id
                                if 'pubchemId' in ligand and ligand['pubchemId'] and ligand['pubchemId'] != 'None':
                                    # create ligand
                                    l = Ligand()


                                    # update ligand by pubchem id
                                    ligand_title = False
                                    if 'title' in ligand and ligand['title']:
                                        ligand_title = ligand['title']
                                    l = l.load_from_pubchem('cid', ligand['pubchemId'], lt, ligand_title)


                                # if no pubchem id is specified, use name
                                else:
                                    # use ligand title, if specified
                                    if 'title' in ligand and ligand['title']:
                                        ligand['name'] = ligand['title']

                                    # create empty properties
                                    lp = LigandProperities.objects.create()
                                    
                                    # create the ligand
                                    try:
                                        l, created = Ligand.objects.get_or_create(name=ligand['name'], canonical=True,
                                            defaults={'properities': lp, 'ambigious_alias': False})
                                        if created:
                                            self.logger.info('Created ligand {}'.format(ligand['name']))
                                        else:
                                            continue
                                    except IntegrityError:
                                        l = Ligand.objects.get(name=ligand['name'], canonical=True)

                                    # save ligand
                                    l.save()
                            else:
                                continue

                            # structure-ligand interaction
                            if l and ligand['role']:
                                role_slug = slugify(ligand['role'])
                                try:
                                    lr, created = LigandRole.objects.get_or_create(slug=role_slug,
                                    defaults={'name': ligand['role']})
                                    if created:
                                        self.logger.info('Created ligand role {}'.format(ligand['role']))
                                except IntegrityError:
                                    lr = LigandRole.objects.get(slug=role_slug)

                                i, created = StructureLigandInteraction.objects.get_or_create(structure=s,
                                    ligand=l, ligand_role=lr, annotated=True,
                                    defaults={'pdb_reference': pdb_reference})
                    
                    # structure segments
                    if 'segments' in sd and sd['segments']:
                        for segment, positions in sd['segments'].items():
                            # fetch (create if needed) sequence segment
                            try:
                                protein_segment = ProteinSegment.objects.get(slug=segment)
                            except ProteinSegment.DoesNotExist:
                                self.logger.error('Segment {} not found'.format(segment))
                                continue

                            struct_seg, created = StructureSegment.objects.update_or_create(structure=s,
                                protein_segment=protein_segment, defaults={'start': positions[0], 'end': positions[1]})
                    # all representive structures should have defined segments
                    elif representative:
                        self.logger.warning('Segments not defined for representative structure {}'.format(sd['pdb']))

                    # structure segments for modeling
                    if 'segments_in_structure' in sd and sd['segments_in_structure']:
                        for segment, positions in sd['segments_in_structure'].items():
                            # fetch (create if needed) sequence segment
                            try:
                                protein_segment = ProteinSegment.objects.get(slug=segment)
                            except ProteinSegment.DoesNotExist:
                                self.logger.error('Segment {} not found'.format(segment))
                                continue

                            struct_seg_mod, created = StructureSegmentModeling.objects.update_or_create(structure=s,
                                protein_segment=protein_segment, defaults={'start': positions[0], 'end': positions[1]})

                    # structure coordinates
                    if 'coordinates' in sd and sd['coordinates']:
                        for segment, coordinates in sd['coordinates'].items():
                            # fetch (create if needed) sequence segment
                            try:
                                protein_segment = ProteinSegment.objects.get(slug=segment)
                            except ProteinSegment.DoesNotExist:
                                self.logger.error('Segment {} not found'.format(segment))
                                continue

                            # fetch (create if needed) coordinates description
                            try:
                                description, created = StructureCoordinatesDescription.objects.get_or_create(
                                    text=coordinates)
                                if created:
                                    self.logger.info('Created structure coordinate description {}'.format(coordinates))
                            except IntegrityError:
                                description = StructureCoordinatesDescription.objects.get(text=coordinates)

                            sc = StructureCoordinates()
                            sc.structure = s
                            sc.protein_segment = protein_segment
                            sc.description = description
                            sc.save()

                    # structure engineering
                    if 'engineering' in sd and sd['engineering']:
                        for segment, engineering in sd['engineering'].items():
                            # fetch (create if needed) sequence segment
                            try:
                                protein_segment = ProteinSegment.objects.get(slug=segment)
                            except ProteinSegment.DoesNotExist:
                                self.logger.error('Segment {} not found'.format(segment))
                                continue

                            # fetch (create if needed) engineering description
                            try:
                                description, created = StructureEngineeringDescription.objects.get_or_create(
                                    text=engineering)
                                if created:
                                    self.logger.info('Created structure coordinate description {}'.format(engineering))
                            except IntegrityError:
                                description = StructureEngineeringDescription.objects.get(text=engineering)

                            se = StructureEngineering()
                            se.structure = s
                            se.protein_segment = protein_segment
                            se.description = description
                            se.save()

                    # protein anomalies
                    scheme = s.protein_conformation.protein.residue_numbering_scheme
                    if 'bulges' in sd and sd['bulges']:
                        pa_slug = 'bulge'
                        try:
                            pab, created = ProteinAnomalyType.objects.get_or_create(slug=pa_slug, defaults={
                                'name': 'Bulge'})
                            if created:
                                self.logger.info('Created protein anomaly type {}'.format(pab))
                        except IntegrityError:
                            pab = ProteinAnomalyType.objects.get(slug=pa_slug)
                        
                        for segment, bulges in sd['bulges'].items():
                            for bulge in bulges:
                                try:
                                    gn, created = ResidueGenericNumber.objects.get_or_create(label=bulge,
                                        scheme=scheme, defaults={'protein_segment': ProteinSegment.objects.get(
                                        slug=segment)})
                                    if created:
                                        self.logger.info('Created generic number {}'.format(gn))
                                except IntegrityError:
                                    gn =  ResidueGenericNumber.objects.get(label=bulge, scheme=scheme)

                                try:
                                    pa, created = ProteinAnomaly.objects.get_or_create(anomaly_type=pab,
                                        generic_number=gn)
                                    if created:
                                        self.logger.info('Created protein anomaly {}'.format(pa))
                                except IntegrityError:
                                    pa, created = ProteinAnomaly.objects.get(anomaly_type=pab, generic_number=gn)

                                s.protein_anomalies.add(pa)
                    if 'constrictions' in sd and sd['constrictions']:
                        pa_slug = 'constriction'
                        try:
                            pac, created = ProteinAnomalyType.objects.get_or_create(slug=pa_slug, defaults={
                                'name': 'Constriction'})
                            if created:
                                self.logger.info('Created protein anomaly type {}'.format(pac))
                        except IntegrityError:
                            pac = ProteinAnomalyType.objects.get(slug=pa_slug)
                        
                        for segment, constrictions in sd['constrictions'].items():
                            for constriction in constrictions:
                                try:
                                    gn, created = ResidueGenericNumber.objects.get_or_create(label=constriction,
                                        scheme=scheme, defaults={'protein_segment': ProteinSegment.objects.get(
                                        slug=segment)})
                                    if created:
                                        self.logger.info('Created generic number {}'.format(gn))
                                except IntegrityError:
                                    gn =  ResidueGenericNumber.objects.get(label=constriction, scheme=scheme)

                                try:
                                    pa, created = ProteinAnomaly.objects.get_or_create(anomaly_type=pac,
                                        generic_number=gn)
                                    if created:
                                        self.logger.info('Created protein anomaly {}'.format(pa))
                                except IntegrityError:
                                    pa, created = ProteinAnomaly.objects.get(anomaly_type=pac, generic_number=gn)

                                s.protein_anomalies.add(pa)
                    
                    # stabilizing agents, FIXME - redesign this!
                    # fusion proteins moved to constructs, use this for G-proteins and other agents?
                    aux_proteins = []
                    if 'signaling_protein' in sd and sd['signaling_protein'] and sd['signaling_protein'] != 'None':
                        aux_proteins.append('signaling_protein')
                    if 'auxiliary_protein' in sd and sd['auxiliary_protein'] and sd['auxiliary_protein'] != 'None':
                        aux_proteins.append('auxiliary_protein')
                    for index in aux_proteins:
                        if isinstance(sd[index], list):
                            aps = sd[index]
                        else:
                            aps = [sd[index]]
                        for aux_protein in aps:
                            aux_protein_slug = slugify(aux_protein)[:50]
                            try:
                                sa, created = StructureStabilizingAgent.objects.get_or_create(
                                    slug=aux_protein_slug, defaults={'name': aux_protein})
                            except IntegrityError:
                                sa = StructureStabilizingAgent.objects.get(slug=aux_protein_slug)
                            s.stabilizing_agents.add(sa)

                    # save structure
                    s.save()

                    self.create_rotamers(s)
                    self.logger.info('Calculate interactions')

                    runcalculation(sd['pdb'])
                    try:
                        parsecalculation(sd['pdb'],False)
                    except:
                        self.logger.error('Error parsing interactions output for {}'.format(sd['pdb']))
