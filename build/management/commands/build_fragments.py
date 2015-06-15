from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import connection

from interaction.models import *
from ligand.models import LigandRole
from residue.models import Residue
from protein.models import Protein
from structure.models import Structure,Rotamer,Fragment,PdbData
from structure.functions import extract_pdb_data

from Bio.PDB import PDBParser
import Bio.PDB.Polypeptide as polypeptide

import logging, os, sys

class Command(BaseCommand):

    logger = logging.getLogger(__name__)
    
    fragments_dir = os.sep.join([settings.DATA_DIR, 'fragment_data'])
    #Interactions definitions to be uploaded to the db
    interactions = {
        'HBD': 'Hydrogen bond donor',
        'HBD-ICT': 'Hurr durr',
        'HBD_H2O': 'Hydrogen bond, water proxy',
        'HBA': 'Hydrogen bond acceptor',
        'HBA-IAN': 'Hurr durr',
        'AFE': 'Aromatic face to edge',
        'AEF': 'Aromatic edge to face',
        'AFF': 'Aromatic face to face',
        'IAN': 'Hurr durr',
        'ICT': 'Hurr durr',
        'KPI': 'Hurr durr',
        }

    def handle(self, *args, **options):
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
                feature = '_'.join([tokens[5], tokens[6]])
            elif len(tokens) == 6:
                feature = tokens[5]

        #Checking the if the crystal is in the database
        try:
            s = Structure.objects.get(pdb_code__index=pdb_code)
        except Structure.DoesNotExist:
            self.logger.warning('Cannot find the structure {} in the database. Skipping the fragment {}'.format(pdb_code, fragment_file_name.strip().replace('.pdb', '')))
            return

        #ResidueFragmentInteractionType
        try:
            i, created = ResidueFragmentInteractionType.objects.get_or_create(slug=feature, name=self.interactions[feature])
        except Exception:
            self.logger.info("Failed to find or create feature {}...".format(feature))
        #Rotamer and Fragment
        try:
            fragment_struct = PDBParser(PERMISSIVE=True).get_structure('frag', os.sep.join([self.fragments_dir, fragment_file_name]))[0]
            fragment_pdb_data = ''
            r = None
            for residue in fragment_struct.get_residues():
                hetfield, resseq, icode=residue.get_id()
                if hetfield == ' ': #Amino acid
                    try:
                        r = Residue.objects.get(sequence_number=int(resseq), amino_acid=polypeptide.three_to_one(residue.resname),protein_conformation=s.protein_conformation)
                        d, created = PdbData.objects.get_or_create(pdb=extract_pdb_data(residue))
                        rot, created = Rotamer.objects.get_or_create(residue=r, structure=s, pdbdata=d)
                        #rot.save()
                    except Exception as msg:
                        self.logger.error('Failed to add rotamer {}:{}{}\n'.format(pdb_code, resseq, msg))
                        return
                else:
                    fragment_pdb_data += extract_pdb_data(residue)
            try:
                fd, created = PdbData.objects.get_or_create(pdb=fragment_pdb_data)
                #Taking the first ligand from the list, since existing fragments do not contain the ligand info
                f, created = Fragment.objects.get_or_create(residue=r, ligand=s.ligands.all()[0], structure=s, pdbdata=fd)
                #f.save()
            except Exception as msg:
                self.logger.error('Failed to add fragment {}\n{}'.format(fragment_file_name, msg))
        except Exception as msg:
            self.logger.error('Failed to add fragment {} to the db\n{}'.format(fragment_file_name, msg))
        #StructureLigandInteraction
        try:
            lr, created = LigandRole.objects.get_or_create(name='unknown',slug='unknown')
            sli, created = StructureLigandInteraction.objects.get_or_create(structure=s, ligand=s.ligands.all()[0], ligand_role=lr)
        except Exception as msg:
            self.logger.error("Failed to add fragment {} to the db\n{}".format(fragment_file_name, msg))
        try:
            rfi, created = ResidueFragmentInteraction.objects.get_or_create(structure_ligand_pair=sli, rotamer=rot, fragment=f, interaction_type=i)
            self.logger.info("Successfully added interacting fragment {}".format(fragment_file_name))
        except Exception as msg:
            self.logger.error("Failed to add fragment {} to the db\n{}".format(fragment_file_name, msg))
