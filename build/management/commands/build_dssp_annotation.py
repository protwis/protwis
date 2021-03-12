from django.core.management.base import BaseCommand
from build.management.commands.base_build import Command as BaseBuild

from residue.models import ResidueDataType, ResidueDataPoint
from protein.models import *
from structure.models import *


from common.tools import fetch_from_cache, save_to_cache, fetch_from_web_api

import Bio.PDB as PDB

import logging
from urllib import request, parse
import json,time

class Command(BaseBuild):
    help = 'Add dssp annotations to structures.'

    logger = logging.getLogger(__name__)
    def add_arguments(self, parser):
        parser.add_argument('-p', '--proc',
            type=int,
            action='store',
            dest='proc',
            default=1,
            help='Number of processes to run')

    def handle(self, *args, **options):

        # All human proteins and xtaled
        self.structures = Structure.objects.all()

        self.dssp_type, created = ResidueDataType.objects.get_or_create(slug=slugify('dssp'), name='DSSP')


        self.prepare_input(options['proc'], self.structures)
        self.logger.info('Finished dssp annotations')

    def dssp_dict(self,handle, chain): 
      """Internal function used by mask_dssp_dict (PRIVATE). 
   
      Return a DSSP dictionary that maps (chainid, resid) to an amino acid, 
      secondary structure symbol, solvent accessibility value, and hydrogen bond 
      information (relative dssp indices and hydrogen bond energies) from an open 
      DSSP file object. 
   
      Parameters 
      ---------- 
      handle : file 
          the open DSSP output file handle 
   
      """ 
      dssp = {} 
      start = 0 
      for l in handle.split('\n'): 
          sl = l.split() 
          if len(sl) < 2: 
              continue 
          if sl[1] == "RESIDUE": 
              # Start parsing from here 
              start = 1 
              continue 
          if not start: 
              continue 
          if l[9] == " ": 
              # Skip -- missing residue 
              continue 
   
          dssp_index = int(l[:5]) 
          resseq = int(l[5:10]) 
          icode = l[10] 
          chainid = l[11] 
          if chainid!=chain: 
            continue
          aa = l[13] 
          ss = l[16] 
          if ss == " ": 
              ss = "-" 
          try: 
              NH_O_1_relidx = int(l[38:45]) 
              NH_O_1_energy = float(l[46:50]) 
              O_NH_1_relidx = int(l[50:56]) 
              O_NH_1_energy = float(l[57:61]) 
              NH_O_2_relidx = int(l[61:67]) 
              NH_O_2_energy = float(l[68:72]) 
              O_NH_2_relidx = int(l[72:78]) 
              O_NH_2_energy = float(l[79:83]) 
   
              acc = int(l[34:38]) 
              phi = float(l[103:109]) 
              psi = float(l[109:115]) 
          except ValueError as exc: 
              # DSSP output breaks its own format when there are >9999 
              # residues, since only 4 digits are allocated to the seq num 
              # field.  See 3kic chain T res 321, 1vsy chain T res 6077. 
              # Here, look for whitespace to figure out the number of extra 
              # digits, and shift parsing the rest of the line by that amount. 
              if l[34] != ' ': 
                  shift = l[34:].find(' ') 
   
                  NH_O_1_relidx = int(l[38 + shift:45 + shift]) 
                  NH_O_1_energy = float(l[46 + shift:50 + shift]) 
                  O_NH_1_relidx = int(l[50 + shift:56 + shift]) 
                  O_NH_1_energy = float(l[57 + shift:61 + shift]) 
                  NH_O_2_relidx = int(l[61 + shift:67 + shift]) 
                  NH_O_2_energy = float(l[68 + shift:72 + shift]) 
                  O_NH_2_relidx = int(l[72 + shift:78 + shift]) 
                  O_NH_2_energy = float(l[79 + shift:83 + shift]) 
   
                  acc = int((l[34 + shift:38 + shift])) 
                  phi = float(l[103 + shift:109 + shift]) 
                  psi = float(l[109 + shift:115 + shift]) 
              else: 
                  raise ValueError(exc) 
          # dssp[resseq] = (aa, ss, acc, phi, psi, dssp_index, 
          #         NH_O_1_relidx, NH_O_1_energy, O_NH_1_relidx, O_NH_1_energy, 
          #         NH_O_2_relidx, NH_O_2_energy, O_NH_2_relidx, O_NH_2_energy) 
          dssp[resseq] = ss
      return dssp 
   

    ## This will be saved into DB
    # H   Alpha helix
    # B   Beta bridge
    # E   Strand
    # G   Helix-3
    # I   Helix-5
    # T   Turn
    # S   Bend

    # @transaction.atomic
    def main_func(self, positions, iteration,count,lock):
        while count.value<len(self.structures):
            with lock:
                s = self.structures[count.value]
                count.value +=1 
                self.logger.info('Generating DSSP data for \'{}\'... ({} out of {})'.format(s, count.value, len(self.structures)))
            print(s)

            pdbcode = s.pdb_code.index.lower()
            chain = s.preferred_chain

            # Grab DSSP db index number
            url = 'http://mrs.cmbi.ru.nl/search?db=dssp&q=%s&count=3' % (pdbcode)
            r = request.urlopen(url)
            t = r.geturl()
            d_id = t.split('=')[2][:-3]

            # Grab DSSP file
            url = 'http://mrs.cmbi.ru.nl/download?db=dssp&nr=$index'
            cache_dir = ['dssp', 'id']
            dssp = fetch_from_web_api(url, d_id, cache_dir, raw=True)

            # Parse file
            dssp = self.dssp_dict(dssp,chain)

            rs = Residue.objects.filter(protein_conformation=s.protein_conformation).all()

            for r in rs:
                if r.sequence_number in dssp:
                    point, created = ResidueDataPoint.objects.get_or_create(data_type=self.dssp_type, residue=r, value_text=dssp[r.sequence_number])
