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


class Command(BaseCommand):
    

    logger = logging.getLogger(__name__)

    #avoiding pathing shenanigans
    generic_numbers_source_dir = os.sep.join([settings.DATA_DIR, 'residue_data'])
    help = 'Creates residues from protein records from the filenames specifeid as arguments. The files are looked up in the directory {}'.format(generic_numbers_source_dir)
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
            'residue_generic_number',
            'residue_set',
            'residue',
                
        ]

        for table in tables_to_truncate:
            cursor.execute("TRUNCATE TABLE {!s} CASCADE".format(table))

    def create_residues(self, args):
        self.logger.info('CREATING RESIDUES')  
        for arg in args:
            residue_data = {}
            if os.path.exists(os.sep.join([self.generic_numbers_source_dir, arg])):
                residue_data = self.parse_residue_data_file(os.sep.join([self.generic_numbers_source_dir, arg]))
                self.logger.info('USING DATA FROM OLD GPCRB')
                if len(ResidueNumberingScheme.objects.all()) == 0:
                    self.add_numbering_schemes()
                #just a shortcut to prevent gazilion of subqueries
                oliveira_id = ResidueNumberingScheme.objects.get(slug='oliveira')
                bw_id = ResidueNumberingScheme.objects.get(slug='bw')
                gpcrdb_id = ResidueNumberingScheme.objects.get(slug='gpcrdb')
                baldwin_id = ResidueNumberingScheme.objects.get(slug='baldwin')
            else:
                print("Can't find {!s}".format(os.sep.join([self.generic_numbers_source_dir, arg])))
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
                            self.logger.info('Created residue {:n}{!s}for protein {!s}'.format(i, aa, protein.name))
                        except Exception as msg:
                            print(msg)
                            self.logger.error('Failed to create residue {:n}{!s}for protein {!s}'.format(i, aa, protein.name))
                  
                        for res_record in residue_data[protein.entry_name]:
                            if int(res_record[0]) == r.sequence_number and res_record[1] == r.three_letter():
                                r.protein_segment = ProteinSegment.objects.get(slug=res_record[6])

                                try:
                                    oliveira = ResidueGenericNumber.objects.get(label=res_record[2], scheme=oliveira_id)
                                except ResidueGenericNumber.DoesNotExist as e:
                                    oliveira = ResidueGenericNumber(label=res_record[2], scheme=oliveira_id)
                                    oliveira.save()
                                r.generic_number.add(oliveira)
                                try:
                                    bw = ResidueGenericNumber.objects.get(label=res_record[3], scheme=bw_id)
                                except ResidueGenericNumber.DoesNotExist as e:
                                    bw = ResidueGenericNumber(label=res_record[3], scheme=bw_id)
                                    bw.save()
                                r.generic_number.add(bw)
                                try:
                                    gpcrdb = ResidueGenericNumber.objects.get(label=res_record[4], scheme=gpcrdb_id)
                                except ResidueGenericNumber.DoesNotExist as e:
                                    gpcrdb = ResidueGenericNumber(label=res_record[4], scheme=gpcrdb_id)
                                    gpcrdb.save()
                                r.generic_number.add(gpcrdb)
                                try:
                                    baldwin = ResidueGenericNumber.objects.get(label=res_record[5], scheme=baldwin_id)
                                except ResidueGenericNumber.DoesNotExist as e:
                                    baldwin = ResidueGenericNumber(label=res_record[5], scheme=baldwin_id)
                                    baldwin.save()
                                r.generic_number.add(baldwin)

                        try:
                            r.save()
                            self.logger.info('Added generic numbers for residue {:n}{!s}for protein {!s}'.format(i, aa, protein.name))
                        except Exception as msg:
                            print(msg)
                            self.logger.error('Failed to create residue {:n}{!s}for protein {!s}'.format(i, aa, protein.name))
                self.logger.info('COMPLETED CREATING RESIDUES FROM FILE {}'.format(os.sep.join([self.generic_numbers_source_dir, arg])))
            self.logger.info('COMPLETED CREATING RESIDUES')

    def parse_residue_data_file(self, file_name):
        print('Parsing residue data from {}'.format(file_name))
        residue_data = {}
        residue_data_fh = open(file_name, 'r')

        for line in residue_data_fh:
            id,num,res_name,family,oli,gpcrdb,bw,bs,prot_name,sec_str_name = [x.strip('"') for x in line.split(',')]
            #the data will be in dict of lists
            if prot_name not in residue_data.keys():
                residue_data[prot_name] = []
            residue_data[prot_name].append([num, res_name, oli, gpcrdb, bw, bs, sec_str_name])

        print('done')
        return residue_data

    def add_numbering_schemes(self):
        #FIXME temporary workaround, will be (?) in a separate file
        rns = ResidueNumberingScheme.objects.create(slug="oliveira", name="Oliveira")
        rns = ResidueNumberingScheme.objects.create(slug="bw", name="Ballesteros-Weinstein")
        rns = ResidueNumberingScheme.objects.create(slug="gpcrdb", name="GPCRdb generic")
        rns = ResidueNumberingScheme.objects.create(slug="baldwin", name="Baldwin-Schwartz")
