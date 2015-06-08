from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import connection

from interaction.models import *
from residue.models import Residue
from protein.models import Protein
from structure.models import Structure,Rotamer,Fragment

from Bio.PDB import PDBParser
from optparse import make_option
import logging, os

class Command(BaseCommand):

    logger = logging.getLogger(__name__)
    
    #fragments_dir = os.sep.join([settings.DATA_DIR, 'fragment_data'])
    fragments_dir = 'C:\\Users\clz967\documents\protwis_vagrant\shared\protwis\data\fragments'

    #Interactions definitions to be uploaded to the db
    interactions = {
        'HBD': 'Hydrogen bond donor',
        'HBD-ICT': 'Hurr durr',
        'HBD_H2O': 'Hydrogen bond, water proxy',
        'HBA': 'Hydrogen bond acceptor',
        'HBA-IAN': 'Hurr durr',
        'AFE': 'Aromatic face to edge',
        'AFF': 'Aromatic face to face',
        'IAN': 'Hurr durr',
        'ICT': 'Hurr durr',
        }

    def handle(self):
        # delete any existing fragment data
        fragments = self.read_fragments_dir()
        for fragment_fname in fragments:
            self.create_fragment(fragment_fname)

        self.logger.info("DONE CREATING FRAGMENTS")
        

    def read_fragments_dir(self):
        return [f for f in os.listdir(self.fragments_dir) if os.path.isfile(os.path.join(self.fragments_dir, f))]


    def create_fragment(self, fragment_file_name):

        tokens = fragment_file_name.strip().replace('.pdb', '').split('_')
        #Not the most efficient way, but gives the overview on what is going on
        generic_num = float("%s.%s" %(tokens[0], tokens[1]))
        res_name = tokens[2]
        protein_entry_name = tokens[3]
        pdb_code = tokens[4]

        if len(tokens) > 5:
            if len(tokens) == 7:
                feature = '_'.join(tokens[5], tokens[6])
            elif len(tokens) == 6:
                feature = token[5]

        #Checking the if the crystal is in the database
        try:
            s = Structure.objects.get(pdb_code__index=pdb_code)
        except Structure.DoesNotExist:
            self.logger.warning('Cannot find the structure {} in the database. Skipping the fragment {}'.foramat(pdb_code, fragment_file_name.strip().replace('.pdb', '')))
            return

        #ResidueFragmentInteractionType
        try:
            i = ResidueFragmentInteractionType.objects.get_or_create(slug=feature, name=self.interactions[feature])
        except ResidueFragmentInteractionType.DoesNotExist:
            self.logger.info("Failed to find or create feature {}...".format(feature))
        #Rotamer and Fragment
        try:
            fragment_struct = PDBParser(PERMISSIVE=True).get_structure('frag', fragment_file_name)[0]
            fragment_pdb_data = ''
            for residue in fragment_struct.get_residues():
                hetfield, resseq, icode=residue.get_id()
                if hetfield == '': #Amino acid
                    try:
                        r = Residue.objects.get()
                        rot = Rotamer(residue=r.id, Structure=s.id, pdbdata=extract_pdb_data(residue))
                        rot.save()
                    except Exception as msg:
                        self.logger.error('Failed to add rotamer {}:{}{}\n'.format(pdb_code, residue.resnum, resseq, msg))
                        return
                else:
                    fragment_pdb_data += extract_pdb_data(residue)
            try:
                f = Fragment(residue=r.id, ligand='', structure=s.id, pdbdata=fragment_pdb_data)
                f.save()
            except Exception as msg:
                self.logger.error('Failed to add fragment {}'.format(fragment_file_name))

