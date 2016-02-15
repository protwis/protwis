from django.core.management.base import BaseCommand, CommandError
from django.core.management import call_command
from django.conf import settings
from django.db import connection
from structure.sequence_parser import SequenceParser

from Bio import SeqIO, pairwise2
from Bio.PDB import PDBParser, PPBuilder, parse_pdb_header

import logging

class Command(BaseCommand):

    help = "Function evaluating sequence compatibility with wildtype protein."

    logger = logging.getLogger(__name__)

    def add_arguments(self, parser):
        parser.add_argument('pdb_file')

    def handle(self, *args, **options):
        print("Working on file {}".format(options['pdb_file']))
        
        seq_seqres = list(SeqIO.parse(options['pdb_file'], 'pdb-seqres'))
        seq_pdb = SeqIO.parse(options['pdb_file'], 'pdb-atom')
        #print("\n".join([str(x.seq) for x in seq_pdb]))
        #print(seq_seqres)
        #print([x for x in seq_seqres])
        #print('\n'.join([str(x.seq) for x in seq_seqres]))

        print("PDB structure:")
        pdb_header = parse_pdb_header(options['pdb_file'])
        print(pdb_header['name'])        
        pdb_struct = PDBParser(PERMISSIVE=True, QUIET=True).get_structure('pdb', options['pdb_file'])
        for chain in pdb_struct[0]:
            print(chain.id)
            poly = PPBuilder()
            for pp in poly.build_peptides(chain):
                print(pp.get_sequence())
                pwa = pairwise2.align.localms(seq_seqres[0], pp.get_sequence(), 2, -4, -4, -.1, one_alignment_only=True)[0]
                print(pwa)
              

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