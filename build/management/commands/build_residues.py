from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import connection

from protein.models import Protein
from residue.models import Residue
from residue.models import ResidueGenericNumber
from residue.models import ResidueNumberingScheme

from optparse import make_option
from os import path
import logging


class Command(BaseCommand):
    

    logger = logging.getLogger(__name__)

    generic_numbers_source_file = settings.BASE_DIR + '\data\generic_residue_numbers_dump.csv'
    help = 'Creates residues from protein records, if the file {} is found, the generic numbers are also added (if present)'.format(generic_numbers_source_file)
    #option_list = BaseCommand.option_list + (
    #    make_option('--update-generic-numbers',
    #                action='store_true',
    #                dest='generic',
    #                default=False,
    #                help='Update the residue records with generic numbers extracted from the old gpcrdb'),
    #    )
    

    def handle(self, *args, **options):

        # delete any existing residue data
        try:
            self.truncate_residue_tables()
        except Exception as msg:
            print(msg)
            self.logger.error(msg)

        # create residue records for all proteins
        try:
            self.create_residues()
        except Exception as msg:
            print(msg)
            self.logger.error(msg)

    def truncate_residue_tables(self):
        cursor = connection.cursor()
        
        tables_to_truncate = [
            #Following the changes in the models - SM
            'residue',
            'residue_set',
            'generic_number'
        ]

        for table in tables_to_truncate:
            cursor.execute("TRUNCATE TABLE {!s} CASCADE".format(table))

    def create_residues(self):
        self.logger.info('CREATING RESIDUES')
        
        residue_data = {}
        if path.exists(self.generic_numbers_source_file):
            residue_data = self.parse_residue_data_file()
            self.logger.info('USING DATA FROM OLD GPCRB')
            if len(ResidueNumberingScheme.objects.all()) == 0:
                self.add_numbering_schemes()
                #just a shortcut to prevent gazilion of subqueries
            oliveira_id = ResidueNumberingScheme.objects.get(slug='oliveira')
            bw_id = ResidueNumberingScheme.objects.get(slug='bw')
            gpcrdb_id = ResidueNumberingScheme.objects.get(slug='gpcrdb')
            baldwin_id = ResidueNumberingScheme.objects.get(slug='baldwin')
        else:
            print("Can't find {!s}".format(self.generic_numbers_source_file))
        proteins = Protein.objects.all()

        for protein in proteins:
            for i, aa in enumerate(protein.sequence):
                r = Residue()
                r.protein = protein
                r.sequence_number = i+1
                r.amino_acid = aa
                try:
                    r.save()
                    self.logger.info('Created residue {:n}{!s}for protein {!s}'.format(i, aa, protein.name))
                except Exception as msg:
                    print(msg)
                    self.logger.error('Failed to create residue {:n}{!s}for protein {!s}'.format(i, aa, protein.name))

                generic_numbers = []
                
                if protein.entry_name in residue_data.keys():                    
                    for res_record in residue_data[protein.entry_name]:
                        if int(res_record[0]) == r.sequence_number and res_record[1] == r.three_letter():
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

        self.logger.info('COMPLETED CREATING RESIDUES')

    def parse_residue_data_file(self):
        print('Parsing residue data')
        residue_data = {}
        residue_data_fh = open(self.generic_numbers_source_file, 'r')

        for line in residue_data_fh:
            id,num,oli,gpcrdb,bw,bs,res_name,prot_name,sec_str_id = [x.strip('"') for x in line.split(',')]
            #the data will be in dict of lists
            if prot_name not in residue_data.keys():
                residue_data[prot_name] = []
            residue_data[prot_name].append([num, res_name, oli, gpcrdb, bw, bs])

        print('done')
        return residue_data

    def add_numbering_schemes(self):
        #FIXME temporary workaround, will be (?) in a separate file
        rns = ResidueNumberingScheme(slug="oliveira", name="Oliveira")
        rns.save()
        rns = ResidueNumberingScheme(slug="bw", name="Ballesteros-Weinstein")
        rns.save()
        rns = ResidueNumberingScheme(slug="gpcrdb", name="GPCRdb generic")
        rns.save()
        rns = ResidueNumberingScheme(slug="baldwin", name="Baldwin-Schwartz")
        rns.save()
