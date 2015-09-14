from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import connection
from django.utils.text import slugify

from protein.models import (Protein, ProteinConformation, ProteinState, ProteinAnomaly, ProteinAnomalyType,
    ProteinSegment)
from residue.models import ResidueGenericNumber, ResidueNumberingScheme, Residue
from common.models import WebLink, WebResource, Publication
from structure.models import (Structure, StructureType, StructureSegment, StructureStabilizingAgent,PdbData,
    Rotamer)
from ligand.models import Ligand, LigandType, LigandRole, LigandProperities
from interaction.models import *
from interaction.views import runcalculation,parsecalculation

from optparse import make_option
from datetime import datetime
import logging, os, re
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

class Command(BaseCommand):
    help = 'Reads source data and creates pdb structure records'
    
    def add_arguments(self, parser):
        parser.add_argument('--filename', action='append', dest='filename',
            help='Filename to import. Can be used multiple times')
        parser.add_argument('--purge', action='store_true', dest='purge', default=False,
            help='Purge existing construct records')

    logger = logging.getLogger(__name__)

    # source file directory
    structure_data_dir = os.sep.join([settings.DATA_DIR, 'structure_data', 'structures'])

    def handle(self, *args, **options):
        # delete any existing structure data
        if options['purge']:
            try:
                self.purge_structures()
            except Exception as msg:
                print(msg)
                self.logger.error(msg)
        # import the structure data
        try:
            self.create_structures(options['filename'])
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

    def create_structures(self, filenames):
        self.logger.info('CREATING PDB STRUCTURES')
        
        # what files should be parsed?
        if not filenames:
            filenames = os.listdir(self.structure_data_dir)

        for source_file in filenames:
            source_file_path = os.sep.join([self.structure_data_dir, source_file])
            if os.path.isfile(source_file_path) and source_file[0] != '.':
                self.logger.info('Reading file {}'.format(source_file_path))
                # read the yaml file
                with open(source_file_path, 'r') as f:
                    sd = yaml.load(f)

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

                    # is this a representative structure (will be used to guide structure-based alignments)?
                    representative = False
                    if 'representative' in sd and sd['representative']:
                        representative = True
                        s.representative = True

                    # get the PDB file and save to DB
                    if not os.path.exists(self.structure_data_dir+'/../pdbs/'):
                        os.makedirs(self.structure_data_dir+'/../pdbs/')
                    pdb_path = self.structure_data_dir+'/../pdbs/'+sd['pdb']+'.pdb'
                    if not os.path.isfile(pdb_path):
                        self.logger.info('Fetching PDB file {}'.format(sd['pdb']))
                        url = 'http://www.rcsb.org/pdb/files/%s.pdb' % sd['pdb']
                        pdbdata_raw = urlopen(url).read().decode('utf-8')
                        f=open(pdb_path,'w')
                        f.write(pdbdata_raw)
                        f.close();
                    else:
                        pdbdata_raw = open(pdb_path, 'r').read()
                    
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

                    header = open(pdb_path,'r')
                    header_dict = parse_pdb_header(header)
                    header.close()
                    sd['publication_date'] = header_dict['release_date']
                    sd['resolution'] = str(header_dict['resolution']).strip()

                    matched = 0
                    if 'ligand' in sd:
                        if isinstance(sd['ligand'], list):
                            ligands = sd['ligand']
                        else:
                            ligands = [sd['ligand']]
                        for ligand in ligands:
                            if ligand['name'].upper() in hetsyn:
                                self.logger.info('match')
                                matched = 1
                                ligand['name'] = hetsyn[ligand['name'].upper()]
                            elif ligand['name'].upper() in hetsyn_reverse:
                                matched = 1

                    if matched==0 and len(hetsyn)>0:
                        self.logger.info('No ligand names found in HET in structure')
                        self.logger.info(hetsyn)
                        self.logger.info(ligands)

                    # endogenous ligand(s)
                    if 'endogenous_ligand' in sd and sd['endogenous_ligand']:
                        # name, id
                        pass

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

                    # protein state
                    if 'state' not in sd:
                        self.logger.error('State not defined, using default state {}'.format(
                            settings.DEFAULT_PROTEIN_STATE))
                        state = settings.DEFAULT_STATE.title()
                    else:
                        state = sd['state']
                    ps, created = ProteinState.objects.get_or_create(slug=slugify(state), defaults={'name': state})
                    if created:
                        self.logger.info('Created protein state {}'.format(ps.name))
                    s.state = ps

                    # protein conformation
                    try:
                        s.protein_conformation, created = ProteinConformation.objects.get_or_create(protein=con,
                            state=ps)
                        if created:
                            self.logger.info('Created conformation {} for protein {}'.format(ps.name, con.name))
                    except Protein.DoesNotExist:
                        self.logger.error('Construct {} for structure {} does not exists, skipping!'.format(
                            sd['construct'], sd['pdb']))
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

                    # structure type
                    st, created = StructureType.objects.get_or_create(slug='xray', defaults={'name': 'X-ray'})
                    s.structure_type = st

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
                                    p.web_link = WebLink.objects.get(index=sd['pubmed_id'], web_resource__slug='pubmed')
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
                    
                    # ligands
                    if 'ligand' in sd:
                        if isinstance(sd['ligand'], list):
                            ligands = sd['ligand']
                        else:
                            ligands = [sd['ligand']]
                        for ligand in ligands:
                            l = False
                            if ligand['name'] and ligand['name']!='None': #some inserted as none.

                                if ligand['name'] in hetsyn_reverse: ligand['name'] = hetsyn_reverse[ligand['name']] #USE HETSYN NAME, not 3letter pdb reference
                                pdb_reference = None
                                if ligand['name'] in hetsyn: pdb_reference = hetsyn[ligand['name']]
                                if Ligand.objects.filter(name=ligand['name'], canonical=True).exists(): #if this name is canonical and it has a ligand record already
                                    try:
                                        l = Ligand.objects.get(name=ligand['name'], canonical=True)
                                    except: 
                                        try:
                                            l = Ligand.objects.filter(name=ligand['name'], canonical=True, properities__inchikey__isnull=False)[0]
                                        except:
                                            self.logger.error('Skipping '+ligand['name']+" for "+sd['pdb'] +' Something wrong with getting ligand from DB')
                                            continue
                                elif Ligand.objects.filter(name=ligand['name'], canonical=False, ambigious_alias=False).exists(): #if this matches an alias that only has "one" parent canonical name - eg distinct
                                    l = Ligand.objects.get(name=ligand['name'], canonical=False, ambigious_alias=False)
                                elif Ligand.objects.filter(name=ligand['name'], canonical=False, ambigious_alias=True).exists(): #if this matches an alias that only has several canonical parents, must investigate, start with empty.
                                    lp = LigandProperities()
                                    lp.save()
                                    l = Ligand()
                                    l.properities = lp
                                    l.name = ligand['name']
                                    l.canonical = False
                                    l.ambigious_alias = True
                                    l.save()
                                    l.load_by_name(ligand['name'])
                                else: #if niether a canonical or alias exists, create the records. Remember to check for canonical / alias status.
                                    self.logger.info('Inserting '+ligand['name']+" for "+sd['pdb'])
                                    lp = LigandProperities()
                                    lp.save()
                                    l = Ligand()
                                    l.properities = lp
                                    l.name = ligand['name']
                                    l.canonical = True
                                    l.ambigious_alias = False
                                    l.save()
                                    l.load_by_name(ligand['name'])
                            # save ligand
                                l.save()

                            # elif ligand['inchi']:
                            #     pass # FIXME write!
                            else:
                                continue


                            if l:
                                if ligand['role']:
                                    lr, created = LigandRole.objects.get_or_create(slug=slugify(ligand['role']),
                                        defaults={'name': ligand['role']})
                                    i, created = StructureLigandInteraction.objects.get_or_create(structure=s, ligand=l,
                                        ligand_role=lr, annotated=True, defaults={'pdb_reference':pdb_reference})
                    
                    # structure segments
                    if 'segments' in sd and sd['segments']:
                        for segment, positions in sd['segments'].items():
                            # fetch (create if needed) sequence segment
                            try:
                                protein_segment = ProteinSegment.objects.get(slug=segment)
                            except ProteinSegment.DoesNotExist:
                                try:
                                    parent_slug = segment.split('_')[0]
                                    parent_protein_segment = ProteinSegment.objects.get(slug=parent_slug)
                                    protein_segment = ProteinSegment.objects.create(slug=segment,
                                        name=parent_protein_segment.name, category=parent_protein_segment.category,
                                        partial=True)
                                    self.logger.info('Created protein segment {}'.format(segment))
                                except:
                                    self.logger.error('Protein segment {} could not be created'.format(segment))
                                    continue
                            ps = StructureSegment()
                            ps.structure = s
                            ps.protein_segment = ProteinSegment.objects.get(slug=segment)
                            ps.start = positions[0]
                            ps.end = positions[1]
                            ps.save()
                    elif representative:
                        self.logger.warning('Segments not defined for representative structure {}'.format(sd['pdb']))

                    # protein anomalies
                    if 'bulges' in sd and sd['bulges']:
                        scheme = s.protein_conformation.protein.residue_numbering_scheme
                        pab, created = ProteinAnomalyType.objects.get_or_create(slug='bulge', defaults={
                            'name': 'Bulge'})
                        for segment, bulges in sd['bulges'].items():
                            for bulge in bulges:
                                gn, created = ResidueGenericNumber.objects.get_or_create(label=bulge, scheme=scheme,
                                    defaults={'protein_segment': ProteinSegment.objects.get(slug=segment)})
                                pa, created = ProteinAnomaly.objects.get_or_create(anomaly_type=pab, generic_number=gn)
                                s.protein_anomalies.add(pa)
                    if 'constrictions' in sd and sd['constrictions']:
                        pac, created = ProteinAnomalyType.objects.get_or_create(slug='constriction', defaults={
                            'name': 'Constriction'})
                        for segment, constrictions in sd['constrictions'].items():
                            for constriction in constrictions:
                                gn, created = ResidueGenericNumber.objects.get_or_create(label=constriction, defaults={
                                    'protein_segment': ProteinSegment.objects.get(slug=segment),
                                    'scheme': ResidueNumberingScheme.objects.get(slug=scheme.slug)})
                                pa, created = ProteinAnomaly.objects.get_or_create(anomaly_type=pac, generic_number=gn)
                                s.protein_anomalies.add(pa)
                    
                    # stabilizing agents
                    # fusion proteins moved to constructs, use this for G-proteins and other agents?
                    # if 'fusion_protein' in sd:
                    #     if isinstance(sd['fusion_protein'], list):
                    #         fusion_proteins = sd['fusion_protein']
                    #     else:
                    #         fusion_proteins = [sd['fusion_protein']]
                    #     for fusion_protein in fusion_proteins:
                    #         sa, created = StructureStabilizingAgent.objects.get_or_create(slug=slugify(fusion_protein),
                    #             name=fusion_protein)
                    #         s.stabilizing_agents.add(sa)

                    # save structure
                    s.save()

                    self.create_rotamers(s)
                    self.logger.info('Calculate interactions')

                    runcalculation(sd['pdb'])
                    try:
                        parsecalculation(sd['pdb'],False)
                    except:
                        self.logger.error('Error parsing interactions output for {}'.format(sd['pdb']))

        self.logger.info('COMPLETED CREATING PDB STRUCTURES')