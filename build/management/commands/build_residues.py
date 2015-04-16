from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import connection

from protein.models import Protein, ProteinSegment
from residue.models import Residue, ResidueGenericNumber, ResidueNumberingScheme

from optparse import make_option
import logging, os
import shlex


class Command(BaseCommand):
    

    logger = logging.getLogger(__name__)

    #avoiding pathing shenanigans
    dump_source_dir = os.sep.join([settings.DATA_DIR, 'residue_data', 'dump'])
    generic_numbers_source_dir = os.sep.join([settings.DATA_DIR, 'residue_data', 'generic_numbers'])
    help = 'Creates residues from protein records from the filenames specifeid as arguments. The files are looked up in the directory {}'.format(dump_source_dir)
    option_list = BaseCommand.option_list + (
        make_option('--purge_tables',
                    action='store_true',
                    dest='purge',
                    default=False,
                    help='Truncate all the associated tables before inserting new records'),
    #    make_option('--update-generic-numbers',
    #                action='store_true',
    #                dest='generic',
    #                default=False,
    #                help='Update the residue records with generic numbers extracted from the old gpcrdb'),
        )
    

    def handle(self, *args, **options):
        if options['purge']:
            try:
                self.truncate_residue_tables()
            except Exception as msg:
                print(msg)
                self.logger.error(msg)

        # create residue records for all proteins
        try:
            self.create_residues(args)
        except Exception as msg:
            print(msg)
            self.logger.error(msg)

    def truncate_residue_tables(self):
        cursor = connection.cursor()
        
        tables_to_truncate = [
            #Following the changes in the models - SM
            'generic_number',
            'residue_set',
            'residue',
                
        ]

        for table in tables_to_truncate:
            cursor.execute("TRUNCATE TABLE {!s} CASCADE".format(table))

    def create_residues(self, args):
        self.logger.info('CREATING RESIDUES')
        for arg in args:
            residue_data = {}
            if os.path.exists(os.sep.join([self.dump_source_dir, arg])):
                residue_data = self.parse_residue_data_file(os.sep.join([self.dump_source_dir, arg]))
                self.logger.info('USING DATA FROM OLD GPCRB')
                
                # fetch schemes and conversion tables
                schemes = {
                    'gpcrdb': {'type': False},
                    'gpcrdba': {
                        'type': 'structure',
                        'seq_based': 'bw',
                    },
                    'gpcrdbb': {
                        'type': 'structure',
                        'seq_based': 'woot',
                    },
                    'gpcrdbc': {
                        'type': 'structure',
                        'seq_based': 'pin',
                    },
                    'gpcrdbf': {
                        'type': 'structure',
                        'seq_based': 'wang',
                    },
                    'bw': {'type': 'sequence'},
                    'woot': {'type': 'sequence'},
                    'pin': {'type': 'sequence'},
                    'wang': {'type': 'sequence'},
                }
                for scheme_name, scheme in schemes.items():
                    schemes[scheme_name]['obj'] = ResidueNumberingScheme.objects.get(slug=scheme_name)
                    mapping_file = os.sep.join([self.generic_numbers_source_dir, 'mapping_' + scheme_name + '.txt'])
                    if os.path.isfile(mapping_file):
                        with open(mapping_file, "r", encoding='UTF-8') as scheme_table_file:
                            schemes[scheme_name]['table'] = {}
                            for row in scheme_table_file:
                                split_row = shlex.split(row)
                                schemes[scheme_name]['table'][split_row[0]] = split_row[1]
            else:
                print("Can't find {!s}".format(os.sep.join([self.dump_source_dir, arg])))
            proteins = Protein.objects.all()

            for protein in proteins:
                for i, aa in enumerate(protein.sequence):
                    r = Residue()
                    r.protein = protein
                    r.sequence_number = i+1
                    r.amino_acid = aa
                    generic_numbers = []
                
                    if protein.entry_name in residue_data.keys():
                        try:
                            r.save()
                            self.logger.info('Created residue {:n}{!s} for protein {!s}'.format(i, aa,
                                protein.entry_name))
                        except Exception as msg:
                            print(msg)
                            self.logger.error('Failed to create residue {:n}{!s} for protein {!s}'.format(i, aa,
                                protein.entry_name))
                  
                        for res_record in residue_data[protein.entry_name]:
                            if int(res_record[0]) == r.sequence_number and res_record[1] == r.three_letter():
                                # residue segment
                                dump_segment = res_record[6]
                                r.protein_segment = ProteinSegment.objects.get(slug=dump_segment)
                                
                                # separate bulge number (1241 - > 124 + 1)
                                bulge_prime = ''
                                dump_oliveira = str(res_record[2])
                                if len(dump_oliveira) == 4:
                                    bulge_prime = dump_oliveira[3]
                                    dump_oliveira = dump_oliveira[:3]
                                dump_gpcrdb = res_record[3][:4]
                                dump_seq_based = res_record[4]

                                # default gpcrdb number
                                def_gpcrdb = False
                                if dump_oliveira in schemes[settings.DEFAULT_NUMBERING_SCHEME]['table']:
                                    default_label = schemes[settings.DEFAULT_NUMBERING_SCHEME]['table'][dump_oliveira] + bulge_prime
                                    try:
                                        def_gpcrdb = ResidueGenericNumber.objects.get(label=default_label,
                                            scheme=schemes[settings.DEFAULT_NUMBERING_SCHEME]['obj'])
                                    except ResidueGenericNumber.DoesNotExist as e:
                                        def_gpcrdb = ResidueGenericNumber()
                                        def_gpcrdb.label = default_label
                                        def_gpcrdb.scheme = schemes[settings.DEFAULT_NUMBERING_SCHEME]['obj']
                                        def_gpcrdb.protein_segment = r.protein_segment
                                        def_gpcrdb.save()
                                        self.logger.info('Created generic number {:s} in numbering scheme {:s}'
                                            .format(default_label,
                                            schemes[settings.DEFAULT_NUMBERING_SCHEME]['obj'].short_name))
                                
                                # if default number was found/added successfully, process the alternative numbers
                                if def_gpcrdb:
                                    # add default generic number to residue record
                                    r.generic_number = def_gpcrdb

                                    # dict of sequence-based numbers, for use in structure-based numbers (5.46x461)
                                    seq_based_labels = {}

                                    # sequence-based schemes first (the sequence-based numbers are needed for the
                                    # structure based schemes)
                                    for scheme_name, scheme in schemes.items():
                                        if scheme['type'] == 'sequence':
                                            # is this number in the scheme defined for this protein?
                                            if scheme_name == schemes[protein.residue_numbering_scheme.slug]['seq_based']:
                                                seq_based_label = dump_seq_based
                                            # if not convert the number to the correct scheme
                                            else:
                                                for d, c in schemes[schemes[protein.residue_numbering_scheme.slug]['seq_based']]['table'].items():
                                                    if c == dump_seq_based:
                                                        seq_based_label = scheme['table'][d]
                                                        break

                                            # fetch/insert the number
                                            try:
                                                seq_based = ResidueGenericNumber.objects.get(label=seq_based_label,
                                                    scheme=scheme['obj'])
                                            except ResidueGenericNumber.DoesNotExist as e:
                                                seq_based = ResidueGenericNumber()
                                                seq_based.label = seq_based_label
                                                seq_based.scheme = scheme['obj']
                                                seq_based.protein_segment = r.protein_segment
                                                seq_based.save()
                                            r.alternative_generic_number.add(seq_based)

                                            # add added number to the dict for later use
                                            seq_based_labels[scheme_name] = seq_based_label
                                            
                                    # structure-based numbers
                                    for scheme_name, scheme in schemes.items():
                                        if scheme['type'] == 'structure':
                                            # is this number in the scheme defined for this protein?
                                            if scheme_name == protein.residue_numbering_scheme.slug:
                                                struct_based_label = dump_gpcrdb + bulge_prime
                                            # if not convert the number to the correct scheme
                                            else:
                                                for d, c in schemes[protein.residue_numbering_scheme.slug]['table'].items():
                                                    if c == dump_gpcrdb:
                                                        struct_based_label = scheme['table'][d] + bulge_prime
                                                        break

                                            # add the sequence-based label (5x461 -> 5.46x461)
                                            split_struct_based_label = struct_based_label.split('x')
                                            struct_based_label = seq_based_labels[scheme['seq_based']] + 'x' + split_struct_based_label[1]

                                            # fetch/insert the number
                                            try:
                                                struct_based = ResidueGenericNumber.objects.get(
                                                    label=struct_based_label, scheme=scheme['obj'])
                                            except ResidueGenericNumber.DoesNotExist as e:
                                                struct_based = ResidueGenericNumber()
                                                struct_based.label = struct_based_label
                                                struct_based.scheme = scheme['obj']
                                                struct_based.protein_segment = r.protein_segment
                                                struct_based.save()
                                            
                                            # add to residue as a display number or alternative number?
                                            if scheme_name == protein.residue_numbering_scheme.slug:
                                                r.display_generic_number = struct_based
                                            else:
                                                r.alternative_generic_number.add(struct_based)
                        try:
                            r.save()
                            self.logger.info('Added generic numbers for residue {:n}{!s} for protein {!s}'.format(i, aa,
                                protein.entry_name))
                        except Exception as msg:
                            print(msg)
                            self.logger.error(
                                'Failed to create generic numbers for residue {:n}{!s} for protein {!s}'.format(i, aa,
                                protein.entry_name))
                self.logger.info('Completed creating residues for protein {}'.format(protein.entry_name))
            self.logger.info('COMPLETED CREATING RESIDUES')

    def parse_residue_data_file(self, file_name):
        print('Parsing residue data from {}'.format(file_name))
        residue_data = {}
        residue_data_fh = open(file_name, 'r')

        for line in residue_data_fh:
            id,num,res_name,oli,gpcrdb,bw,bw2,bs,prot_name,sec_str_name = [x.strip().strip('"') for x in line.split(',')] #double strip due to some weird bug...
            #the data will be in dict of lists
            if prot_name not in residue_data.keys():
                residue_data[prot_name] = []
            residue_data[prot_name].append([num, res_name, oli, gpcrdb, bw, bs, sec_str_name])

        print('done')
        return residue_data