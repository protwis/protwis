from subprocess import Popen, PIPE
from StringIO import StringIO
from Bio.Blast import NCBIXML


#==============================================================================
# I have put it into separate class for the sake of future uses
class BlastSearch(object):
  
    def __init__ (self, blast_path = 'blastp', blastdb = '../protected/data/gpcrs_iuphar_blastdb', top_results = 1):
  
        self.blast_path = blast_path
        self.blastdb = blastdb
        #typicaly top scored result is enough, but for sequences with missing residues it is better to use more results to avoid getting sequence of e.g. different species
        self.top_results = top_results
      
    #takes Bio.Seq sequence as an input and returns a list of tuples with the alignments 
    def run (self, input_seq):
    
        output = []
        
        #invoking blast with default settings
        blast = Popen('%s -db %s -outfmt 5' %(self.blast_path, self.blastdb), universal_newlines=True, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        blast_out, blast_err = blast.communicate(input=input_seq.seq)
        #print blast_err
        result = NCBIXML.parse(StringIO(blast_out)).next()
        
        for aln in result.alignments[:self.top_results]:
            seq_id = aln.hit_id.split("|")
            #first condition is for working with "default" databases, where SwissProt ids come with some junk
            if 'sp' in seq_id:
                upid = seq_id[seq_id.index('sp')+1]
            else:
                #0 or 1, the index actualy depends on blast version used
                upid = seq_id[0]
                
            if upid is None:
                continue
            
            #print len(aln.hsps)

            #output.append(((SeqRecord(Seq(aln.hsps[0].sbjct), id=upid), aln.hsps[0].sbjct_start), (SeqRecord(Seq(aln.hsps[0].query), id='pdb'), aln.hsps[0].query_start)))
            output.append((upid, aln))
        return output
  
      
    def run_online (self, input_seq = ''):
        #TODO
        pass

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
          
      
    def add_mapping(self, uprot_code, uprot_num):
    
        self.mapping[uprot_code] = uprot_num
  
  
    def add_bw_number(self, bw_number=''):
    
        self.bw = bw_number


    def add_gpcrdb_number(self, gpcrdb_number=''):
        #PDB format does not allow fractional part longer than 2 digits
        #so numbers x.xx1 are negative
        if len(gpcrdb_number) > 4:
          self.gpcrdb = '-' + gpcrdb_number.replace('x', '.')
        else:
          self.gpcrdb = gpcrdb_number.replace('x', '.')
    
    
    def get_mapping(self, uprot_code):
    
        if self.mapping.has_key(uprot_code):
            return self.mapping[uprot_code]
    
        return 0