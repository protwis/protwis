from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import connection

from protein.models import Protein
from protein.models import ProteinSegment
from residue.models import Residue
from residue.models import ResidueGenericNumber
from residue.models import ResidueNumberingScheme

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
                    'gpcrdb_display': {'type': False},
                    'gpcrdba': {
                        'type': 'structure',
                        'seq_based': 'bw',
                    },
                    'gpcrdba': {
                        'type': 'structure',
                        'seq_based': 'woot',
                    },
                    'gpcrdba': {
                        'type': 'structure',
                        'seq_based': 'pin',
                    },
                    'gpcrdba': {
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
                    if scheme['type']:
                        with open(os.sep.join([self.generic_numbers_source_dir, 'mapping_' + scheme_name + '.txt']), "r", encoding='UTF-8') as scheme_table_file:
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
                            self.logger.info('Created residue {:n}{!s} for protein {!s}'.format(i, aa, protein.name))
                        except Exception as msg:
                            print(msg)
                            self.logger.error('Failed to create residue {:n}{!s} for protein {!s}'.format(i, aa, protein.name))
                  
                        for res_record in residue_data[protein.entry_name]:
                            if int(res_record[0]) == r.sequence_number and res_record[1] == r.three_letter():
                                dump_gpcrdb = res_record[3]
                                dump_seq_based = res_record[4]
                                dump_segment = res_record[6]
                                r.protein_segment = ProteinSegment.objects.get(slug=dump_segment)

                                # default gpcrdb number
                                def_gpcrdb = False
                                for d, c in schemes[protein.residue_numbering_scheme.slug]['table'].items():
                                    if c == dump_gpcrdb:
                                        try:
                                            def_gpcrdb = ResidueGenericNumber.objects.get(label=d,
                                                scheme=schemes[settings.DEFAULT_NUMBERING_SCHEME]['obj'])
                                        except ResidueGenericNumber.DoesNotExist as e:
                                            def_gpcrdb = ResidueGenericNumber()
                                            def_gpcrdb.label = d
                                            def_gpcrdb.scheme = schemes[settings.DEFAULT_NUMBERING_SCHEME]['obj']
                                            def_gpcrdb.protein_segment = r.protein_segment
                                            def_gpcrdb.save()
                                            self.logger.info('Created generic number {:s} in numbering scheme {:s}'.format(d, schemes[settings.DEFAULT_NUMBERING_SCHEME]['obj'].short_name))
                                        r.generic_number = def_gpcrdb
                                        break
                                if def_gpcrdb:
                                    split_gpcrdb = def_gpcrdb.label.split('x')
                                    def_gpcrdb_res_pos = split_gpcrdb[1]

                                    for scheme_name, scheme in schemes.items():
                                        # class specific sequence based number
                                        if scheme['type'] == 'sequence':
                                            # is this number in the scheme defined for this protein?
                                            if scheme_name == schemes[protein.residue_numbering_scheme.slug]['seq_based']:
                                                label = dump_seq_based
                                            # if not convert the number to the correct scheme
                                            else:
                                                for d, c in schemes[schemes[protein.residue_numbering_scheme.slug]['seq_based']]['table'].items():
                                                    if c == dump_seq_based:
                                                        label = scheme['table'][d]
                                                        break
                                            try:
                                                seq_based = ResidueGenericNumber.objects.get(label=label,
                                                    scheme=scheme['obj'])
                                            except ResidueGenericNumber.DoesNotExist as e:
                                                seq_based = ResidueGenericNumber()
                                                seq_based.label = label
                                                seq_based.scheme = scheme['obj']
                                                seq_based.protein_segment = r.protein_segment
                                                seq_based.save()
                                            r.alternative_generic_number.add(seq_based)
                                            
                                            # display number
                                            if scheme_name == schemes[protein.residue_numbering_scheme.slug]['seq_based']:
                                                display_label = seq_based.label + 'x' + def_gpcrdb_res_pos
                                                try:
                                                    display = ResidueGenericNumber.objects.get(label=display_label,
                                                        scheme=schemes['gpcrdb_display']['obj'])
                                                except ResidueGenericNumber.DoesNotExist as e:
                                                    display = ResidueGenericNumber()
                                                    display.label = display_label
                                                    display.scheme = schemes['gpcrdb_display']['obj']
                                                    display.protein_segment = r.protein_segment
                                                    display.save()
                                                r.display_generic_number = display
                                        # class specific gpcrdb number
                                        elif scheme['type'] == 'structure':
                                            if scheme_name == protein.residue_numbering_scheme.slug:
                                                label = dump_gpcrdb
                                            else:
                                                for d, c in schemes[protein.residue_numbering_scheme.slug]['table'].items():
                                                    if c == dump_gpcrdb:
                                                        label = scheme['table'][d]
                                                        break
                                            try:
                                                gpcrdb = ResidueGenericNumber.objects.get(label=label,
                                                    scheme=scheme['obj'])
                                            except ResidueGenericNumber.DoesNotExist as e:
                                                gpcrdb = ResidueGenericNumber()
                                                gpcrdb.label = label
                                                gpcrdb.scheme = scheme['obj']
                                                gpcrdb.protein_segment = r.protein_segment
                                                gpcrdb.save()
                                            r.alternative_generic_number.add(gpcrdb)
                        try:
                            r.save()
                            self.logger.info('Added generic numbers for residue {:n}{!s}for protein {!s}'.format(i, aa, protein.name))
                        except Exception as msg:
                            print(msg)
                            self.logger.error('Failed to create residue {:n}{!s}for protein {!s}'.format(i, aa, protein.name))
                self.logger.info('COMPLETED CREATING RESIDUES FROM FILE {}'.format(os.sep.join([self.dump_source_dir, arg])))
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