from django.core.management.base import BaseCommand, CommandError
from django.core.management import call_command
from django.conf import settings
from django.db import connection

from Bio import SeqIO

import logging

class Command(BaseCommand):

    help = "Function evaluating sequence compatibility with wildtype protein."

    logger = logging.getLogger(__name__)

    def add_arguments(self, parser):
        parser.add_argument('pdb_file')

    def handle(self, *args, **options):

        print(options['pdb_file'])
        seq_atom = SeqIO.parse(options['pdb_file'], 'pdb-atom')
        print('\n'.join([str(x.seq) for x in seq_atom]))
        print("SEQRES:")
        seq_seqres = SeqIO.parse(options['pdb_file'], 'pdb-seqres')
        print('\n'.join([str(x.seq) for x in seq_seqres]))

    #def parse_structure(self, pdb_struct):
    #    """
    #    extracting sequence and preparing dictionary of residues
    #    bio.pdb reads pdb in the following cascade: model->chain->residue->atom
    #    """
    #    for chain in pdb_struct:
    #        self.residues[chain.id] = {}
    #        self.pdb_seq[chain.id] = Seq('')
            
    #        for res in chain:
    #        #in bio.pdb the residue's id is a tuple of (hetatm flag, residue number, insertion code)
    #            if res.resname == "HID":
    #                resname = polypeptide.three_to_one('HIS')
    #            else:
    #                if res.resname not in self.residue_list:
    #                    continue
    #                self.residues[chain.id][res.id[1]] = MappedResidue(res.id[1], polypeptide.three_to_one(res.resname))
    
    #        self.pdb_seq[chain.id] = ''.join([self.residues[chain.id][x].name for x in sorted(self.residues[chain.id].keys())])
            
    #        for pos, res in enumerate(sorted(self.residues[chain.id].keys()), start=1):
                #self.residues[chain.id][res].pos_in_aln = pos