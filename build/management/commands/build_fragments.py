from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import connection

from fragment.models import Fragment, Interaction
from residue.models import Residue
from protein.models import Protein
from structure.models import Structure

from optparse import make_option
import logging, os

class Command(BaseCommand):

    logger = logging.getLogger(__name__)
    
    fragments_dir = os.sep.join([settings.DATA_DIR, 'fragment_data'])

    def handle(self):
        # delete any existing fragment data
        try:
            self.truncate_protein_tables()
        except Exception as msg:
            # print(msg)
            self.logger.error(msg)
        fragments = self.read_fragments_dir()
        for fragment_fname in fragments:
            self.create_fragment(fragment_fname)

        self.logger.info("DONE CREATING FRAGMENTS")
        
    def truncate_fragment_tables(self):
        cursor = connection.cursor()
        
        tables_to_truncate = [
            #Following the changes in the models - SM
            'fragment',                
        ]

        for table in tables_to_truncate:
            cursor.execute("TRUNCATE TABLE {!s} CASCADE".format(table))

    def read_fragments_dir(self):
        return [f for f in os.listdir(self.fragments_dir) if os.path.isfile(os.path.join(self.fragments_dir, f))]

    def create_fragment(self, fragment_file_name):
        tokens = fragment_file_name.strip().replace('.pdb', '').split('_')
        #Not the most efficient way, but gives the vie on what is going on
        generic_num = float("%s.%s" %(tokens[0], tokens[1]))
        res_name = tokens[2]
        protein_entry_name = tokens[3]
        pdb_code = tokens[4]

        #Fetching needed objects
        try:
            p = Protein.objects.get(entry_name=protein_entry_name)
        except Protein.DoesNotExist as e:
            print("Can't find the protein: {} {}".format(protein_entry_name, e))
            self.logger.error("Failed to create fragment {} due to PROTEIN lookup error.".format(fragment_file_name.replace('.pdb', '')))
            return
        try:
            r = Residue.objects.get(protein=p, amino_acid=res_name, generic_number__label=str(generic_num))
        except Residue.DoesNotExist as e:
            print("Can't find the residue: {!s} {}: {}".format(generic_num, res_name, e))
            self.logger.error("Failed to create fragment {} due to RESIDUE lookup error.".format(fragment_file_name.replace('.pdb', '')))
            return
        try:
            s = Structure.objects.get(pdb_code__index=pdb_code)
        except Structure.DoesNotExist as e:
            print("Can't find the structure {}: {}".format(pdb_code, e))
            self.logger.error("Failed to create fragment {} due to STRUCTURE lookup error.".format(fragment_file_name.replace('.pdb', '')))
            return
        if len(tokens) > 5:
            if len(tokens) == 7:
                feature = '_'.join(tokens[5], tokens[6])
            elif len(tokens) == 6:
                feature = token[5]
            try:
                i = Interaction.objects.get(slug=feature)
            except Intraction.DoesNotExist:
                i - Interaction(slug=feature, description='')
                                
            Fragment.objects.create(residue=r, protein=p, structure=s, interaction=i)
            self.logger.info("Created interacting moiety: {}".format(fragment_file_name))
        else:
            Fragment.objects.create(residue=r, protein=p, structure=s)
            self.logger.info("Created rotamer: {}".format(fragment_file_name))

