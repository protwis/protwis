from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import connection
from django.utils.text import slugify

from protein.models import (Protein, ProteinConformation, ProteinState, ProteinAnomaly, ProteinAnomalyType,
    ProteinSegment)
from residue.models import ResidueGenericNumber, ResidueNumberingScheme
from common.models import WebLink, WebResource, Publication
from structure.models import Structure, StructureType, StructureSegment, StructureStabilizingAgent
from ligand.models import Ligand, LigandType, LigandRole
from interaction.models import StructureLigandInteraction

from optparse import make_option
from datetime import datetime
import logging, os
import yaml
import json
from urllib.request import urlopen


class Command(BaseCommand):
    help = 'Reads source data and creates pdb structure records'
    
    def add_arguments(self, parser):
        parser.add_argument('--filename', action='append', dest='filename',
            help='Filename to import. Can be used multiple times')
        parser.add_argument('--purge', action='store_true', dest='purge', default=False,
            help='Purge existing construct records')

    logger = logging.getLogger(__name__)

    # source file directory
    structure_data_dir = os.sep.join([settings.DATA_DIR, 'structure_data', 'pdb_structures'])

    def handle(self, *args, **options):
        # delete any existing structure data
        if options['purge']:
            try:
                self.truncate_structure_tables()
            except Exception as msg:
                print(msg)
                self.logger.error(msg)
        # import the structure data
        try:
            self.create_structures(options['filename'])
        except Exception as msg:
            print(msg)
            self.logger.error(msg)
    
    def truncate_structure_tables(self):
        cursor = connection.cursor()
        
        tables_to_truncate = [
            #Following the changes in the models - SM
            'structure',                
        ]

        for table in tables_to_truncate:
            cursor.execute("TRUNCATE TABLE {!s} CASCADE".format(table))

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
                    if 'pubmed_id' in sd:
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
                            if ligand['name']:
                                l, created = Ligand.objects.get_or_create(name=ligand['name'])
                                if created:
                                    l.load_by_name(ligand['name'])
                            # elif ligand['inchi']:
                            #     pass # FIXME write!
                            else:
                                continue

                            # save ligand
                            l.save()

                            if l:
                                if ligand['role']:
                                    lr, created = LigandRole.objects.get_or_create(slug=slugify(ligand['role']),
                                        defaults={'name': ligand['role']})
                                    i = StructureLigandInteraction.objects.create(structure=s, ligand=l,
                                        ligand_role=lr)
                    
                    # structure segments
                    if 'segments' in sd and sd['segments']:
                        for segment, positions in sd['segments'].items():
                            ps = StructureSegment()
                            ps.structure = s
                            ps.protein_segment = ProteinSegment.objects.get(slug=segment)
                            ps.start = positions[0]
                            ps.end = positions[1]
                            ps.save()

                    # protein anomalies
                    if 'bulges' in sd and sd['bulges']:
                        scheme = s.protein_conformation.protein.residue_numbering_scheme
                        pab, created = ProteinAnomalyType.objects.get_or_create(slug='bulge', defaults={
                            'name': 'Bulge'})
                        for segment, bulges in sd['bulges'].items():
                            for bulge in bulges:
                                gn, created = ResidueGenericNumber.objects.get_or_create(label=bulge, defaults={
                                    'protein_segment': ProteinSegment.objects.get(slug=segment),
                                    'scheme': ResidueNumberingScheme.objects.get(slug=scheme.slug)})
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

        self.logger.info('COMPLETED CREATING PDB STRUCTURES')