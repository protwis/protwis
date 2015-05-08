from subprocess import Popen, PIPE
from io import StringIO
from Bio.Blast import NCBIXML

from django.conf import settings
from common.selection import SimpleSelection
from protein.models import ProteinSegment

import os,sys,tempfile,logging

logger = logging.getLogger("structure_gpcr")

#==============================================================================
# I have put it into separate class for the sake of future uses
class BlastSearch(object):
    
    logger = logging.getLogger("structure_gpcr")

    def __init__ (self, blast_path = 'blastp', blastdb = os.sep.join([settings.STATICFILES_DIRS[0], 'blast', 'protwis_blastdb']), top_results = 1):
  
        self.blast_path = blast_path
        self.blastdb = blastdb
        #print(blastdb)
        #typicaly top scored result is enough, but for sequences with missing residues it is better to use more results to avoid getting sequence of e.g. different species
        self.top_results = top_results
      
    #takes Bio.Seq sequence as an input and returns a list of tuples with the alignments 
    def run (self, input_seq):
    
        output = []
        #Windows has problems with Popen and PIPE
        if sys.platform == 'win32':
            tmp = tempfile.NamedTemporaryFile()
            tmp.write(bytes(input_seq.seq+'\n', 'latin1'))
            tmp.seek(0)
            blast = Popen('%s -db %s -outfmt 5' %(self.blast_path, self.blastdb), universal_newlines=True, shell=True, stdin=tmp, stdout=PIPE, stderr=PIPE)
            (blast_out, blast_err) = blast.communicate()
        else:
        #Rest of the world:
            blast = Popen('%s -db %s -outfmt 5' %(self.blast_path, self.blastdb), universal_newlines=True, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
            (blast_out, blast_err) = blast.communicate(input=input_seq.seq)
        if len(blast_err) != 0:
            self.logger.debug(blast_err)

        result = NCBIXML.read(StringIO(blast_out))
        for aln in result.alignments[:self.top_results]:         
            self.logger.debug("Looping over alignments, current hit: {}".format(aln.hit_id))
            output.append((aln.hit_id, aln))
        return output

#==============================================================================

#stores information about alignments and b-w numbers
class MappedResidue(object):
  
    def __init__(self, res_num, res_name):
          
        self.number = res_num
        self.name = res_name
        self.pos_in_aln = 0
        self.mapping = {}
        self.bw = 0.
        self.gpcrdb = 0.       
  
    def add_bw_number(self, bw_number=''):
    
        self.bw = bw_number


    def add_gpcrdb_number(self, gpcrdb_number=''):

        #PDB format does not allow fractional part longer than 2 digits
        #so numbers x.xx1 are negative
        if len(gpcrdb_number) > 4:
          self.gpcrdb = '-' + gpcrdb_number[:4].replace('x', '.')
        else:
          self.gpcrdb = gpcrdb_number.replace('x', '.')

#==============================================================================

#turns selection into actual residues
class SelectionParser(object):

    def __init__(self, selection):
    
        self.generic_numbers = []
        self.helices = []
        
        self.db_segments = [x.slug for x in ProteinSegment.objects.all()]
        
        for segment in selection.segments:
            if 'TM' in segment.item.slug:
                self.helices.append(int(segment.item.slug[-1]))
            elif segment not in self.db_segments:
                self.residues.append(segment)
        logger.debug("Helices selected: {}; Residues: {}".format(self.helices, self.generic_numbers))
    
#==============================================================================


def check_gn(pdb_struct):
        
    gn_list = []
    for chain in pdb_struct:
        for residue in chain:
            try:
                if -8.1 < residue['CA'].get_bfactor() < 8.1:
                    gn_list.append(residue['CA'])
                    return True
            except:
                continue
    return False
