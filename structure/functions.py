from Bio.Blast import NCBIXML, NCBIWWW
from Bio.PDB import PDBParser
from Bio.PDB.PDBIO import Select
import Bio.PDB.Polypeptide as polypeptide
from Bio.PDB.AbstractPropertyMap import AbstractPropertyMap
from Bio.PDB.Polypeptide import CaPPBuilder, is_aa
try:
    from Bio.PDB.vectors import rotaxis
except:
    from Bio.PDB import rotaxis

from django.conf import settings
from common.selection import SimpleSelection
from common.alignment import Alignment
from protein.models import Protein, ProteinSegment, ProteinConformation, ProteinState
from residue.functions import dgn, ggn
from residue.models import Residue, ResidueGenericNumberEquivalent
from structure.models import Structure, Rotamer

from subprocess import Popen, PIPE
from io import StringIO
import os
import sys
import tempfile
import logging
import math
import urllib
from collections import OrderedDict
import Bio.PDB as PDB
import csv
# from openpyxl import Workbook
import numpy
import zipfile
import pprint
import json
import yaml


logger = logging.getLogger("protwis")

ATOM_FORMAT_STRING="%s%5i %-4s%c%3s %c%4i%c   %8.3f%8.3f%8.3f%s%6.2f      %4s%2s%2s\n"

#==============================================================================
# I have put it into separate class for the sake of future uses
class BlastSearch(object):


    def __init__ (self, blast_path='blastp',
        blastdb=os.sep.join([settings.STATICFILES_DIRS[0], 'blast', 'protwis_blastdb']), top_results=1):

        self.blast_path = blast_path
        self.blastdb = blastdb
        #typicaly top scored result is enough, but for sequences with missing
        #residues it is better to use more results to avoid getting sequence of
        #e.g.  different species
        self.top_results = top_results

    #takes Bio.Seq sequence as an input and returns a list of tuples with the
    #alignments
    def run (self, input_seq):

        output = []
        #Windows has problems with Popen and PIPE
        if sys.platform == 'win32':
            tmp = tempfile.NamedTemporaryFile()
            logger.debug("Running Blast with sequence: {}".format(input_seq))
            tmp.write(bytes(str(input_seq) + '\n', 'latin1'))
            tmp.seek(0)
            blast = Popen('%s -db %s -outfmt 5' % (self.blast_path, self.blastdb), universal_newlines=True, stdin=tmp,
                stdout=PIPE, stderr=PIPE)
            (blast_out, blast_err) = blast.communicate()
        else:
        #Rest of the world:
            blast = Popen('%s -db %s -outfmt 5' % (self.blast_path, self.blastdb), universal_newlines=True, shell=True,
                stdin=PIPE, stdout=PIPE, stderr=PIPE)
            (blast_out, blast_err) = blast.communicate(input=str(input_seq))

        if len(blast_err) != 0:
            logger.debug(blast_err)
        if blast_out!='\n':
            result = NCBIXML.read(StringIO(blast_out))
            for aln in result.alignments[:self.top_results]:
                logger.debug("Looping over alignments, current hit: {}".format(aln.hit_id))
                output.append((aln.hit_id, aln))
        return output
#==============================================================================

class BlastSearchOnline(object):

    def __init__(self, blast_program='blastp', db='swissprot', top_results=1):
        self.blast_program = blast_program
        self.db = db
        self.top_results = top_results
        pass

    def run(self, input_seq):
        output = []

        result = NCBIXML.read(NCBIWWW.qblast(self.blast_program, self.db, input_seq, auto_format='xml'))
        for aln in result.alignments[:self.top_results]:
            logger.debug("Looping over alignments, current hit: {}".format(aln.hit_id))
            output.append((aln.hit_id, aln))
        return output

    def get_uniprot_entry_name (self, up_id):

        #get the 'entry name' field for given uniprot id
        #'entry name' is a protein id in gpcrdb
        url = "http://www.uniprot.org/uniprot/?query=accession:%s&columns=entry name&format=tab" %up_id

        #used urllib, urllib2 throws an error here for some reason
        try:
            response = urllib.urlopen(url)
            page = response.readlines()
            return page[1].strip().lower()

        except urllib.HTTPError as error:
            print(error)
            return ''

#==============================================================================

#stores information about alignments and b-w numbers
class MappedResidue(object):

    def __init__ (self, res_num, res_name):

        self.number = res_num
        self.name = res_name
        self.pos_in_aln = 0
        self.mapping = {}
        self.bw = 0.
        self.gpcrdb = 0.
        self.gpcrdb_id = 0
        self.segment = ''
        self.display = ''
        self.residue_record = None

    def add_bw_number (self, bw_number=''):

        self.bw = bw_number

    def add_segment (self, segment=''):

        self.segment = segment

    def add_display_number (self, display = ''):

        self.display = display

    def add_gpcrdb_number_id (self, gpcrdb_number_id=''):

        self.gpcrdb_id = gpcrdb_number_id

    def add_residue_record (self, residue_record= None):

        self.residue_record = residue_record

    def add_gpcrdb_number (self, gpcrdb_number=''):

        #PDB format does not allow fractional part longer than 2 digits
        #so numbers x.xx1 are negative
        if len(gpcrdb_number.split('.')[1]) > 2:
          self.gpcrdb = '-' + gpcrdb_number[:4].replace('x', '.')
        else:
          self.gpcrdb = gpcrdb_number.replace('x', '.')

#==============================================================================

#turns selection into actual residues
class SelectionParser(object):

    def __init__ (self, selection):

        self.generic_numbers = []
        self.helices = []
        self.substructures = []

        for segment in selection.segments:
            logger.debug('Segments in selection: {}'.format(segment))
            if segment.type == 'helix':
                self.helices.append(int(segment.item.slug[-1]))
            elif segment.type == 'residue':
                self.generic_numbers.append(segment.item.label.replace('x','.'))
            else:
                self.substructures.append(segment.item.slug)
        logger.debug("Helices selected: {}; Residues: {}; Other substructures:{}".format(self.helices, self.generic_numbers, self.substructures))


#==============================================================================
class GenericNumbersSelector(Select):

    def __init__(self, generic_numbers=[], helices=[], parsed_selection=None):

        self.generic_numbers = generic_numbers
        self.helices = helices
        if parsed_selection:
            self.generic_numbers=parsed_selection.generic_numbers
            self.helices = parsed_selection.helices

    def accept_residue(self, residue):

        try:
            if "{:.2f}".format(residue['CA'].get_bfactor()) in self.generic_numbers:
                return 1
            if -8.1 < residue['CA'].get_bfactor() < 0 and "{:.3f}".format(-residue['CA'].get_bfactor() + 0.001) in self.generic_numbers:
                return 1
            if -8.1 < residue['CA'].get_bfactor() < 8.1 and int(math.floor(abs(residue['CA'].get_bfactor()))) in self.helices:
                return 1
        except:
            return 0
        return 0

#==============================================================================
class SubstructureSelector(Select):

    def __init__(self, segment_mapping, parsed_selection=None):

        self.residues = []

        for tm in parsed_selection.helices:
            self.residues.extend(segment_mapping['TM{}'.format(tm)])
        for substr in parsed_selection.substructures:
            self.residues.extend(segment_mapping[substr])

    def accept_residue(self, residue):

        try:
            if int(residue.id[1]) in self.residues:
                return 1
        except:
            return 0
        return 0

#==============================================================================
class CASelector(object):

    def __init__ (self, parsed_selection, ref_pdbio_struct, alt_structs):

        self.ref_atoms = []
        self.alt_atoms = {}

        self.selection = parsed_selection
        try:
            self.ref_atoms.extend(self.select_generic_numbers(ref_pdbio_struct))
            self.ref_atoms.extend(self.select_helices(ref_pdbio_struct))
        except Exception as msg:
            logger.warning("Can't select atoms from the reference structure!\n{!s}".format(msg))

        for alt_struct in alt_structs:
            try:
                self.alt_atoms[alt_struct.id] = []
                self.alt_atoms[alt_struct.id].extend(self.select_generic_numbers(alt_struct))
                self.alt_atoms[alt_struct.id].extend(self.select_helices(alt_struct))

            except Exception as msg:
                logger.warning("Can't select atoms from structure {!s}\n{!s}".format(alt_struct.id, msg))


    def select_generic_numbers (self, structure):

        if self.selection.generic_numbers == []:
            return []

        atom_list = []

        for chain in structure:
            for res in chain:
                try:
                    if 0 < res['CA'].get_bfactor() < 8.1 and "{:.2f}".format(res['CA'].get_bfactor()) in self.selection.generic_numbers:
                        atom_list.append(res['CA'])
                    if -8.1 < res['CA'].get_bfactor() < 0 and "{:.3f}".format(-res['CA'].get_bfactor() + 0.001) in self.selection.generic_numbers:
                        atom_list.append(res['CA'])
                except :
                    continue

        if atom_list == []:
            logger.warning("No atoms with given generic numbers {} for  {!s}".format(self.selection.generic_numbers, structure.id))
        return atom_list


    def select_helices (self, structure):

        if self.selection.helices == []:
            return []

        atom_list = []
        for chain in structure:
            for res in chain:
                try:
                    if -8.1 < res['CA'].get_bfactor() < 8.1 and int(math.floor(abs(res['CA'].get_bfactor()))) in self.selection.helices:
                        atom_list.append(res['CA'])
                except Exception as msg:
                    continue

        if atom_list == []:
            logger.warning("No atoms with given generic numbers for  {!s}".format(structure.id))

        return atom_list


    def get_consensus_atom_sets (self, alt_id):

        tmp_ref = []
        tmp_alt = []

        for ref_at in self.ref_atoms:
            for alt_at in self.alt_atoms[alt_id]:
                if ref_at.get_bfactor() == alt_at.get_bfactor():
                    tmp_ref.append(ref_at)
                    tmp_alt.append(alt_at)

        if len(tmp_ref) != len(tmp_alt):
            return ([], [])

        return (tmp_ref, tmp_alt)


    def get_consensus_gn_set (self):

        gn_list = []
        for alt_id in self.alt_atoms.keys():
            tmp_ref, tmp_alt = self.get_consensus_atom_sets(alt_id)

            for ref_ca in tmp_ref:
                for alt_ca in tmp_alt:
                    if ref_ca.get_bfactor() == alt_ca.get_bfactor():
                        if 0 < ref_ca.get_bfactor() < 8.1 and "{:.2f}".format(ref_ca.get_bfactor()) in self.selection.generic_numbers:
                            gn_list.append("{:.2f}".format(ref_ca.get_bfactor()))
                        if -8.1 < ref_ca.get_bfactor() < 0 and "{:.3f}".format(-ref_ca.get_bfactor() + 0.001) in self.selection.generic_numbers:
                            gn_list.append("{:.3f}".format(-ref_ca.get_bfactor() + 0.001))
                        if 0 < ref_ca.get_bfactor() < 8.1 and int(math.floor(abs(ref_ca.get_bfactor()))) in self.selection.helices:
                            gn_list.append("{:.2f}".format(ref_ca.get_bfactor()))
                        if -8.1 < ref_ca.get_bfactor() < 0 and int(math.floor(abs(ref_ca.get_bfactor()))) in self.selection.helices:
                            gn_list.append("{:.3f}".format(-ref_ca.get_bfactor() + 0.001))
        return gn_list


    def get_ref_atoms (self):

        return self.ref_atoms


    def get_alt_atoms (self, alt_id):

        try:
            return self.alt_atoms[alt_id]
        except:
            return []


    def get_alt_atoms_all (self):

        return self.alt_atoms


#==============================================================================
class BackboneSelector():
    """
    Selects backbone atoms from reference structure and rotamer for Superposition
    """

    similarity_dict = {
        "fragment_residue" : 0,
        "interaction_type" : 1,
        "target_residue" : 2
        }

    #similarity_rules = [[['H', 'F', 'Y', 'W'], ['AEF', 'AFF'], ['H', 'F', 'Y', 'W']],
    #    [['Y'], ['AFE'], ['F']],
    #    [['S', 'T'], ['HBA', 'HBD'], ['S', 'T']],]

    #interaction_type slugs changed along the way
    similarity_rules = [[['H', 'F', 'Y', 'W'], ['aro_ef_protein', 'aro_ff'], ['H', 'F', 'Y', 'W']],
        [['Y'], ['aro_fe_protein'], ['F']],
        [['S', 'T'], ['polar_acceptor_protein', 'polar_donor_protein'], ['S', 'T']],]

    def __init__(self, ref_pdbio_struct, fragment, use_similar=False):

        self.ref_atoms = []
        self.alt_atoms = []

        self.ref_atoms = self.select_ref_atoms(fragment, ref_pdbio_struct, use_similar)
        self.alt_atoms = self.select_alt_atoms(PDBParser(PERMISSIVE=True, QUIET=True).get_structure('ref', StringIO(str(fragment.rotamer.pdbdata)))[0])


    def select_ref_atoms(self, fragment, ref_pdbio_struct, use_similar=False):

        for chain in ref_pdbio_struct:
            for res in chain:
                try:
                    gn = self.get_generic_number(res)
                    if gn == fragment.rotamer.residue.display_generic_number.label:
                        # logger.info("Ref {}:{}\tFragment {}:{}".format(polypeptide.three_to_one(res.resname), self.get_generic_number(res), fragment.rotamer.residue.amino_acid, fragment.rotamer.residue.display_generic_number.label))
                        if polypeptide.three_to_one(res.resname) == fragment.rotamer.residue.amino_acid:
                            return [res['CA'], res['N'], res['O']]
                        else:
                            if use_similar:
                                for rule in self.similarity_rules:
                                    if polypeptide.three_to_one(res.resname) in rule[self.similarity_dict["target_residue"]] and fragment.rotamer.residue.amino_acid in rule[self.similarity_dict["target_residue"]] and fragment.interaction_type.slug in rule[self.similarity_dict["interaction_type"]]:
                                        return [res['CA'], res['N'], res['O']]
                        # else:
                        #     if fragment.interaction_type.slug not in ['acc', 'hyd']:
                        #         return [res['CA'], res['N'], res['O']]
                except Exception as msg:
                    continue
        return []


    def select_alt_atoms(self, rotamer_pdbio_struct):

        for chain in rotamer_pdbio_struct:
            for res in chain:
                try:
                    return [res['CA'], res['N'], res['O']]
                except:
                    continue
        return []


    def get_generic_number(self, res):

        if 'CA' not in res:
            return 0.0
        if 0 < res['CA'].get_bfactor() < 8.1:
            return "{:.2f}x{!s}".format(res['N'].get_bfactor(), self._get_fraction_string(res['CA'].get_bfactor()))
        if -8.1 < res['CA'].get_bfactor() < 0:
            return "{:.2f}x{!s}".format(res['N'].get_bfactor(),  self._get_fraction_string(res['CA'].get_bfactor() - 0.001))
        return 0.0

    #TODO: Is this function really neccessary?
    def _get_fraction_string(self, number):

        if number > 0:
            return "{:.2f}".format(number).split('.')[1]
        else:
            return "{:.3f}".format(number).split('.')[1]


    def get_ref_atoms(self):

        return self.ref_atoms


    def get_alt_atoms(self):

        return self.alt_atoms


#==============================================================================
def check_gn (pdb_struct):

    for chain in pdb_struct:
        for residue in chain:
            try:
                if -8.1 < residue['CA'].get_bfactor() < 8.1:
                    return True
            except:
                continue
    return False


#==============================================================================
def get_segment_template (protein, segments=['TM1', 'TM2', 'TM3', 'TM4','TM5','TM6', 'TM7'], state=None):

    a = Alignment()
    a.load_reference_protein(protein)
    #You are so gonna love it...
    if state:
        a.load_proteins([x.protein_conformation.protein.parent for x in list(Structure.objects.order_by('protein_conformation__protein__parent','resolution').exclude(protein_conformation__protein=protein.id, protein_conformation__state=state))])
    else:
        a.load_proteins([x.protein_conformation.protein.parent for x in list(Structure.objects.order_by('protein_conformation__protein__parent','resolution').exclude(protein_conformation__protein=protein.id))])
    a.load_segments(ProteinSegment.objects.filter(slug__in=segments))
    a.build_alignment()
    a.calculate_similarity()

    return a.proteins[1]


#==============================================================================
def fetch_template_structure (template_protein):

    return Structure.objects.get(protein_conformation__protein__parent=template_protein.entry_name)


#==============================================================================
def extract_pdb_data(residue):
    """Returns PDB string of a given residue"""
    pdb_string = ''
    hetfield, resseq, icode=residue.get_id()
    resname=residue.get_resname()
    segid=residue.get_segid()
    atom_number = 1
    for atom in residue:
        pdb_string += get_atom_line(atom, hetfield, segid, atom_number, resname, resseq, icode, residue.get_parent().get_id())
        atom_number += 1
    return pdb_string


#==============================================================================
# def convert_csv_to_xlsx(self, csv, separator):
#     wb = Workbook()
#     sheet = wb.active

#     CSV_SEPARATOR = separator

#     with open(csv) as f:
#         reader = csv.reader(f)
#         for r, row in enumerate(reader):
#             for c, col in enumerate(row):
#                 for idx, val in enumerate(col.split(CSV_SEPARATOR)):
#                     cell = sheet.cell(row=r+1, column=idx+1)
#                     cell.value = val

#     wb.save("my_file.xlsx")



#==============================================================================
def get_atom_line(atom, hetfield, segid, atom_number, resname, resseq, icode, chain_id, charge="  "):
    """Returns an ATOM PDB string."""
    if hetfield!=" ":
        record_type="HETATM"
    else:
        record_type="ATOM  "
    if atom.element:
        element = atom.element.strip().upper()
        element = element.rjust(2)
    else:
        element = "  "
    name=atom.get_fullname()
    altloc=atom.get_altloc()
    x, y, z=atom.get_coord()
    bfactor=atom.get_bfactor()
    occupancy=atom.get_occupancy()
    try:
        occupancy_str = "%6.2f" % occupancy
    except TypeError:
        if occupancy is None:
            occupancy_str = " " * 6
            import warnings
            from Bio import BiopythonWarning
            warnings.warn("Missing occupancy in atom %s written as blank" % repr(atom.get_full_id()), BiopythonWarning)
        else:
            raise TypeError("Invalid occupancy %r in atom %r" % (occupancy, atom.get_full_id()))
        pass
    args=(record_type, atom_number, name, altloc, resname, chain_id, resseq, icode, x, y, z, occupancy_str, bfactor, segid, element, charge)
    return ATOM_FORMAT_STRING % args

#==============================================================================
class HSExposureCB(AbstractPropertyMap):
    """
    Abstract class to calculate Half-Sphere Exposure (HSE). #GP: Biopython class repurposed
    The HSE can be calculated based on the CA-CB vector, or the pseudo CB-CA
    vector based on three consecutive CA atoms. This is done by two separate
    subclasses.
    """
    def __init__(self, model, radius, offset=0, hse_up_key='HSE_U', hse_down_key='HSE_D', angle_key=None, check_chain_breaks=False, 
                 check_knots=False, receptor=None, signprot=None,  restrict_to_chain=[], check_hetatoms=False):
        """
        @param model: model
        @type model: L{Model}

        @param radius: HSE radius
        @type radius: float

        @param offset: number of flanking residues that are ignored in the calculation
        of the number of neighbors
        @type offset: int

        @param hse_up_key: key used to store HSEup in the entity.xtra attribute
        @type hse_up_key: string

        @param hse_down_key: key used to store HSEdown in the entity.xtra attribute
        @type hse_down_key: string

        @param angle_key: key used to store the angle between CA-CB and CA-pCB in
        the entity.xtra attribute
        @type angle_key: string
        """
        assert(offset>=0)
        # For PyMOL visualization
        self.ca_cb_list=[]
        ppb=CaPPBuilder()
        ppl=ppb.build_peptides(model)
        hse_map={}
        hse_list=[]
        hse_keys=[]
        ### GP
        if model.get_id()!=0:
            model = model[0]
        residues_in_pdb,residues_with_proper_CA=[],[]
        if check_chain_breaks==True:
            # for m in model:
                for chain in model:
                    for res in chain:
                        # try:
                            if is_aa(res):
                                residues_in_pdb.append(res.get_id()[1])
                        # except:
                        #     if is_aa(chain):
                        #         residues_in_pdb.append(chain.get_id()[1])
                        #         print('chain', chain, res)
                        #         break
        het_resis, het_resis_close = [], []
        for chain in model:
            for res in chain:
                if res.get_id()[0]!=' ':
                    het_resis.append(res)
        self.clash_pairs = []
        self.chain_breaks = []
        
        if check_knots:
            possible_knots = PossibleKnots(receptor, signprot)
            knot_resis = possible_knots.get_resnums()
            self.remodel_resis = {}
        if len(restrict_to_chain)>0:
            restricted_ppl = []
            for p in ppl:
                if p[0].get_parent().get_id() in restrict_to_chain:
                    restricted_ppl.append(p)
            ppl = restricted_ppl
        ###########
        for pp1 in ppl:
            for i in range(0, len(pp1)):
                residues_with_proper_CA.append(pp1[i].get_id()[1])
                if i==0:
                    r1=None
                else:
                    r1=pp1[i-1]
                r2=pp1[i]
                if i==len(pp1)-1:
                    r3=None
                else:
                    r3=pp1[i+1]
                # This method is provided by the subclasses to calculate HSE
                result=self._get_cb(r1, r2, r3)
                if result is None:
                    # Missing atoms, or i==0, or i==len(pp1)-1
                    continue
                pcb, angle=result
                hse_u=0
                hse_d=0
                ca2=r2['CA'].get_vector()
                residue_up=[]   ### GP
                residue_down=[] ### GP
                for pp2 in ppl:
                    for j in range(0, len(pp2)):
                        try:
                            if r2.get_id()[1]-1!=r1.get_id()[1] or r2.get_id()[1]+1!=r3.get_id()[1]:
                                pass
                            else:
                                raise Exception
                        except:
                            if pp1 is pp2 and abs(i-j)<=offset:
                            # neighboring residues in the chain are ignored
                                continue
                        ro=pp2[j]
                        if not is_aa(ro) or not ro.has_id('CA'):
                            continue
                        cao=ro['CA'].get_vector()
                        d=(cao-ca2)
                        if d.norm()<radius:
                            if d.angle(pcb)<(math.pi/2):
                                hse_u+=1
                                ### GP
                                # Puts residues' names in a list that were found in the upper half sphere
                                residue_up.append(ro)

                                ### end of GP code
                            else:
                                hse_d+=1
                                ### GP
                                # Puts residues' names in a list that were found in the lower half sphere
                                residue_down.append(ro)
                                ### end of GP code
                res_id=r2.get_id()
                chain_id=r2.get_parent().get_id()
                # Fill the 3 data structures
                hse_map[(chain_id, res_id)]=(hse_u, hse_d, angle)
                hse_list.append((r2, (residue_up, residue_down, hse_u, hse_d, angle)))
                ### GP residue_up and residue_down added to hse_list
                hse_keys.append((chain_id, res_id))
                # Add to xtra
                r2.xtra[hse_up_key]=hse_u
                r2.xtra[hse_down_key]=hse_d
                if angle_key:
                    r2.xtra[angle_key]=angle

                ### GP checking for knots
                if check_knots:
                    for knot in knot_resis:
                        if knot[0][1]==pp1[i].get_id()[1] and knot[0][0]==pp1[i].get_parent().get_id():
                            # print(pp1[i].get_parent().get_id(),pp1[i]) #print reference
                            for r in residue_up:
                                if r.get_parent().get_id()==knot[1][0] and r.get_id()[1] in knot[1][1]:
                                    # print('close: ', r.get_parent().get_id(),r) #print res within radius
                                    resi_range = [knot[1][1][0], knot[1][1][-1]]
                                    if knot[1][0] not in self.remodel_resis:
                                        self.remodel_resis[knot[1][0]] = [resi_range]
                                    else:
                                        if resi_range not in self.remodel_resis[knot[1][0]]:
                                            self.remodel_resis[knot[1][0]].append(resi_range)

                ### GP checking for atom clashes
                include_prev, include_next = False, False
                try:
                    if pp1[i].get_id()[1]-1!=pp1[i-1].get_id()[1]:
                        include_prev = True
                except:
                    include_prev = False
                try:
                    if pp1[i].get_id()[1]+1!=pp1[i+1].get_id()[1]:
                        include_next = True
                except:
                    include_next = False
                for atom in pp1[i]:
                    ref_vector = atom.get_vector()
                    for other_res in residue_up:
                        try:
                            if other_res==pp1[i-1] and include_prev==False:
                                continue
                            elif len(pp1)>=i+1 and other_res==pp1[i+1] and include_next==False:
                                continue
                            else:
                                raise Exception
                        except:
                            for other_atom in other_res:
                                other_vector = other_atom.get_vector()
                                d = other_vector-ref_vector
                                if d.norm()<2:
                                    if len(str(pp1[i]['CA'].get_bfactor()).split('.')[1])==1:
                                        clash_res1 = float(str(pp1[i]['CA'].get_bfactor())+'0')
                                    else:
                                        clash_res1 = pp1[i]['CA'].get_bfactor()
                                    if len(str(other_res['CA'].get_bfactor()).split('.')[1])==1:
                                        clash_res2 = float(str(other_res['CA'].get_bfactor())+'0')
                                    else:
                                        clash_res2 = other_res['CA'].get_bfactor()
                                    self.clash_pairs.append([(clash_res1, pp1[i].get_id()[1]), (clash_res2, other_res.get_id()[1])])
                    if check_hetatoms:
                        for het_res in het_resis:
                            for het_atom in het_res:
                                het_atom_vector = het_atom.get_vector()
                                d = het_atom_vector-ref_vector
                                if d.norm()<6:
                                    het_resis_close.append(het_res)
        ### GP checking HETRESIS to remove if not interacting with AAs
        self.hetresis_to_remove = [i for i in het_resis if i not in het_resis_close]
        if check_chain_breaks:
            for r in residues_in_pdb:
                if r not in residues_with_proper_CA:
                    self.chain_breaks.append(r)


    def _get_cb(self, r1, r2, r3):
        """
        Method to calculate CB-CA vector.

        @param r1, r2, r3: three consecutive residues (only r2 is used)
        @type r1, r2, r3: L{Residue}
        """
        if r2.get_resname()=='GLY':
            return self._get_gly_cb_vector(r2), 0.0
        else:
            if r2.has_id('CB') and r2.has_id('CA'):
                vcb=r2['CB'].get_vector()
                vca=r2['CA'].get_vector()
                return (vcb-vca), 0.0
        return None

    def _get_gly_cb_vector(self, residue):
        """
        Return a pseudo CB vector for a Gly residue.
        The pseudoCB vector is centered at the origin.

        CB coord=N coord rotated over -120 degrees
        along the CA-C axis.
        """
        try:
            n_v=residue["N"].get_vector()
            c_v=residue["C"].get_vector()
            ca_v=residue["CA"].get_vector()
        except:
            return None
        # center at origin
        n_v=n_v-ca_v
        c_v=c_v-ca_v
        # rotation around c-ca over -120 deg
        rot=rotaxis(-math.pi*120.0/180.0, c_v)
        cb_at_origin_v=n_v.left_multiply(rot)
        # move back to ca position
        cb_v=cb_at_origin_v+ca_v
        # This is for PyMol visualization
        self.ca_cb_list.append((ca_v, cb_v))
        return cb_at_origin_v


class PossibleKnots():
    def __init__(self, receptor, signprot):
        self.receptor = Protein.objects.get(entry_name=receptor)
        self.signprot = Protein.objects.get(entry_name=signprot)
        self.possible_knots = {'ICL3-H4':[['R','A'],['G.H4.11','G.H4.14','G.H4.15','G.h4s6.01']],
                               'h1ha-hehf':[['A','A'],['H.HE.08', 'H.hdhe.05']]}
        self.output = []

    def get_resnums(self):
        if self.receptor and self.signprot:
            for knot_label, values in self.possible_knots.items():
                chain1, chain2 = values[0]
                region1 = list(Residue.objects.filter(protein_conformation__protein=self.receptor, protein_segment__slug=knot_label.split('-')[0]).values_list('sequence_number', flat=True))
                if len(region1)==0:
                    region1 = list(Residue.objects.filter(protein_conformation__protein=self.signprot, protein_segment__slug=knot_label.split('-')[0]).values_list('sequence_number', flat=True))
                if len(region1)==0:
                    raise AssertionError('Protein segment slug error for loop knot: No residues found for {} in {}'.format(knot_label.split('-')[0], self.receptor, self.signprot))
                if knot_label=='h1ha-hehf':
                    for i in range(0,3):
                        region1.append(region1[-1]+1)
                for r in values[1]:
                    region2 = Residue.objects.get(protein_conformation__protein=self.signprot, display_generic_number__label=r)
                    self.output.append([[chain2,region2.sequence_number],[chain1,region1]])
        return self.output


class PdbChainSelector():
    def __init__(self, pdb_code, protein):
        self.pdb_code = pdb_code
        try:
            protein.entry_name
        except:
            protein = Protein.objects.get(entry_name=protein)
        self.protein = protein
        self.chains = []
        self.dssp_dict = OrderedDict()
        self.dssp_info = OrderedDict()
        self.aux_residues = []

    def run_dssp(self):
        pdb = PDB.PDBList()
        pdb.retrieve_pdb_file(self.pdb_code, pdir='./', file_format="pdb")
        p = PDB.PDBParser()
        f = 'pdb{}.ent'.format(self.pdb_code.lower())
        wt_residues = [i for i in Residue.objects.filter(protein_conformation__protein=self.protein).exclude(protein_segment__slug__in=['N-term','C-term'])]
        gn_residues = [i.sequence_number for i in wt_residues if i.generic_number and i.protein_segment.slug not in ['ECL1','ECL2','ICL3','ECL3']]
        structure = p.get_structure(self.pdb_code, f)
        for chain in structure[0]:
            ch = chain.get_id()
            self.chains.append(ch)
            self.dssp_dict[ch] = OrderedDict()
            self.dssp_info[ch] = OrderedDict([('H',0),('B',0),('E',0),('G',0),('I',0),('T',0),('S',0),('-',0)])
        if len(self.dssp_dict)>1:
            dssp = PDB.DSSP(structure[0], f, dssp='/env/bin/dssp')
            for key in dssp.keys():
                if int(key[1][1]) in gn_residues:
                    self.dssp_dict[key[0]][key[1][1]] = dssp[key]
                    self.dssp_info[key[0]][dssp[key][2]] = self.dssp_info[key[0]][dssp[key][2]]+1
        os.remove(f)

    def get_seqnums_by_secondary_structure(self, chain, secondary_structure):
        ''' Returns list of sequence numbers that match the secondary structural property. \n
            H: alpha helix, B: isolated beta-bridge, E: strand, G: 3-10 helix, I: pi-helix, T: turn, S: bend, -: Other
        '''
        if secondary_structure not in ['H','B','E','G','I','T','S','-']:
            raise ValueError('Incorrect secondary_structure input. Input can be H, B, E, G, I, T, S, -')
        output = []
        for seqnum, val in self.dssp_dict[chain].items():
            if val[2]==secondary_structure:
                output.append(seqnum)
        return output

    def select_chain(self):
        num_helix_res = []
        seq_lengths = []
        for c, val in self.dssp_info.items():
            seq_length = 0
            num_helix_res.append(val['H']+val['G']+val['I'])
            for s, num in val.items():
                seq_length+=num
            seq_lengths.append(seq_length)

        max_res = num_helix_res[0]
        max_i = 0
        for i in range(1,len(num_helix_res)):
            if num_helix_res[i]>max_res:
                if num_helix_res[i]-max_res>=seq_lengths[max_i]-seq_lengths[i]:
                    max_res = num_helix_res[i]
                    max_i = i
            elif num_helix_res[i]<max_res:
                if max_res-num_helix_res[i]<seq_lengths[i]-seq_lengths[max_i]:
                    max_res = num_helix_res[i]
                    max_i = i
            elif num_helix_res[i]==max_res:
                if seq_lengths[max_i]<seq_lengths[i]:
                    max_res = num_helix_res[i]
                    max_i = i
        return self.chains[max_i]


class PdbStateIdentifier():
    def __init__(self, structure, tm2_gn='2x41', tm6_gn='6x38', tm3_gn='3x44', tm7_gn='7x52', inactive_cutoff=2, intermediate_cutoff=7.15):
        self.structure_type = None

        try:
            if structure.protein_conformation.protein.parent==None:
                raise Exception
            self.structure = structure
            self.structure_type = 'structure'
            family = structure.protein_conformation.protein.family
        except:
            try:
                structure.protein_conformation.protein
                self.structure = structure
                self.structure_type = 'refined'
                family = structure.protein_conformation.protein.family
            except:
                try:
                    structure.protein
                    self.structure = structure
                    self.structure_type = 'hommod'
                    family = structure.protein.family
                except:
                    structure.receptor_protein
                    self.structure = structure
                    self.structure_type = 'complex'
                    family = structure.receptor_protein.family
        if tm2_gn=='2x41' and tm6_gn=='6x38' and tm3_gn=='3x44' and tm7_gn=='7x52' and inactive_cutoff==2 and intermediate_cutoff==7.15:
            if family.slug.startswith('002') or family.slug.startswith('003'):
                tm6_gn, tm7_gn = '6x33', '7x51'
                inactive_cutoff, intermediate_cutoff = 2.5, 5.75
            elif family.slug.startswith('004'):
                inactive_cutoff, intermediate_cutoff = 5, 7.15
        self.tm2_gn, self.tm6_gn, self.tm3_gn, self.tm7_gn = tm2_gn, tm6_gn, tm3_gn, tm7_gn
        self.inactive_cutoff = inactive_cutoff
        self.intermediate_cutoff = intermediate_cutoff
        self.state = None
        self.activation_value = None
        self.line = False

    def run(self):
        if self.structure_type=='structure':
            self.parent_prot_conf = ProteinConformation.objects.get(protein=self.structure.protein_conformation.protein.parent)
            ssno = StructureSeqNumOverwrite(self.structure)
            ssno.seq_num_overwrite('pdb')
        elif self.structure_type=='refined':
            self.parent_prot_conf = ProteinConformation.objects.get(protein=self.structure.protein_conformation.protein)
        elif self.structure_type=='hommod':
            self.parent_prot_conf = ProteinConformation.objects.get(protein=self.structure.protein)
        elif self.structure_type=='complex':
            self.parent_prot_conf = ProteinConformation.objects.get(protein=self.structure.receptor_protein)
        # class A and T
        if self.parent_prot_conf.protein.family.slug.startswith('001') or self.parent_prot_conf.protein.family.slug.startswith('007'):
            tm6 = self.get_residue_distance(self.tm2_gn, self.tm6_gn)
            tm7 = self.get_residue_distance(self.tm3_gn, self.tm7_gn)
            print(tm6, tm7, tm6-tm7)
            if tm6!=False and tm7!=False:
                self.activation_value = tm6-tm7
                if self.activation_value<self.inactive_cutoff:
                    self.state = ProteinState.objects.get(slug='inactive')
                elif self.inactive_cutoff<=self.activation_value<=self.intermediate_cutoff:
                    self.state = ProteinState.objects.get(slug='intermediate')
                elif self.activation_value>self.intermediate_cutoff:
                    self.state = ProteinState.objects.get(slug='active')
        # class B
        elif self.parent_prot_conf.protein.family.slug.startswith('002') or self.parent_prot_conf.protein.family.slug.startswith('003'):
            tm2_gn_b = ResidueGenericNumberEquivalent.objects.get(default_generic_number__label=self.tm2_gn, scheme__short_name='GPCRdb(B)').label
            tm6_gn_b = ResidueGenericNumberEquivalent.objects.get(default_generic_number__label=self.tm6_gn, scheme__short_name='GPCRdb(B)').label
            tm3_gn_b = ResidueGenericNumberEquivalent.objects.get(default_generic_number__label=self.tm3_gn, scheme__short_name='GPCRdb(B)').label
            tm7_gn_b = ResidueGenericNumberEquivalent.objects.get(default_generic_number__label=self.tm7_gn, scheme__short_name='GPCRdb(B)').label

            tm6 = self.get_residue_distance(tm2_gn_b, tm6_gn_b)
            tm7 = self.get_residue_distance(tm3_gn_b, tm7_gn_b)
            if tm6!=False and tm7!=False:
                self.activation_value = tm6-tm7
                if self.activation_value<self.inactive_cutoff:
                    self.state = ProteinState.objects.get(slug='inactive')
                elif self.inactive_cutoff<=self.activation_value<=self.intermediate_cutoff:
                    self.state = ProteinState.objects.get(slug='intermediate')
                elif self.activation_value>self.intermediate_cutoff:
                    self.state = ProteinState.objects.get(slug='active')
        # class C
        elif self.parent_prot_conf.protein.family.slug.startswith('004'):
            tm2_gn_c = ResidueGenericNumberEquivalent.objects.get(default_generic_number__label=self.tm2_gn, scheme__short_name='GPCRdb(C)').label
            tm6_gn_c = ResidueGenericNumberEquivalent.objects.get(default_generic_number__label=self.tm6_gn, scheme__short_name='GPCRdb(C)').label
            tm3_gn_c = ResidueGenericNumberEquivalent.objects.get(default_generic_number__label=self.tm3_gn, scheme__short_name='GPCRdb(C)').label
            tm7_gn_c = ResidueGenericNumberEquivalent.objects.get(default_generic_number__label=self.tm7_gn, scheme__short_name='GPCRdb(C)').label

            tm6 = self.get_residue_distance(tm2_gn_c, tm6_gn_c)
            tm7 = self.get_residue_distance(tm3_gn_c, tm7_gn_c)
            if tm6!=False and tm7!=False:
                self.activation_value = tm6-tm7
                if self.activation_value<self.inactive_cutoff:
                    self.state = ProteinState.objects.get(slug='inactive')
                elif self.inactive_cutoff<=self.activation_value<=self.intermediate_cutoff:
                    self.state = ProteinState.objects.get(slug='intermediate')
                elif self.activation_value>self.intermediate_cutoff:
                    self.state = ProteinState.objects.get(slug='active')
        # class D
        #########
        # class F
        elif self.parent_prot_conf.protein.family.slug.startswith('006'):
            tm2_gn_f = ResidueGenericNumberEquivalent.objects.get(default_generic_number__label=self.tm2_gn, scheme__short_name='GPCRdb(F)').label
            tm6_gn_f = ResidueGenericNumberEquivalent.objects.get(default_generic_number__label=self.tm6_gn, scheme__short_name='GPCRdb(F)').label
            tm3_gn_f = ResidueGenericNumberEquivalent.objects.get(default_generic_number__label=self.tm3_gn, scheme__short_name='GPCRdb(F)').label
            tm7_gn_f = ResidueGenericNumberEquivalent.objects.get(default_generic_number__label=self.tm7_gn, scheme__short_name='GPCRdb(F)').label

            tm6 = self.get_residue_distance(tm2_gn_f, tm6_gn_f)
            tm7 = self.get_residue_distance(tm3_gn_f, tm7_gn_f)
            if tm6!=False and tm7!=False:
                self.activation_value = tm6-tm7
                if self.activation_value<0:
                    self.state = ProteinState.objects.get(slug='inactive')
                elif 0<=self.activation_value<=2:
                    self.state = ProteinState.objects.get(slug='intermediate')
                elif self.activation_value>2:
                    self.state = ProteinState.objects.get(slug='active')
        else:
            print('{} is not class A,B,C,F'.format(self.structure))
        if self.structure_type=='structure':
            ssno.seq_num_overwrite('pdb')

    def get_residue_distance(self, residue1, residue2):
        try:
            res1 = Residue.objects.get(protein_conformation__protein=self.structure.protein_conformation.protein.parent, display_generic_number__label=dgn(residue1, self.parent_prot_conf))
            res2 = Residue.objects.get(protein_conformation__protein=self.structure.protein_conformation.protein.parent, display_generic_number__label=dgn(residue2, self.parent_prot_conf))
            print(res1, res1.id, res2, res2.id)
            try:
                rota1 = Rotamer.objects.filter(structure=self.structure, residue__sequence_number=res1.sequence_number)
                if len(rota1)==0:
                    raise Exception
            except:
                rota1 = Rotamer.objects.filter(structure=self.structure, residue__display_generic_number__label=dgn(residue1, self.structure.protein_conformation))
            rota1 = right_rotamer_select(rota1, self.structure.preferred_chain[0])
            try:
                rota2 = Rotamer.objects.filter(structure=self.structure, residue__sequence_number=res2.sequence_number)
                if len(rota2)==0:
                    raise Exception
            except:
                rota2 = Rotamer.objects.filter(structure=self.structure, residue__display_generic_number__label=dgn(residue2, self.structure.protein_conformation))
            rota2 = right_rotamer_select(rota2, self.structure.preferred_chain[0])
            rotas = [rota1, rota2]
            io1 = StringIO(rotas[0].pdbdata.pdb)
            rota_struct1 = PDB.PDBParser(QUIET=True).get_structure('structure', io1)[0]
            io2 = StringIO(rotas[1].pdbdata.pdb)
            rota_struct2 = PDB.PDBParser(QUIET=True).get_structure('structure', io2)[0]

            for chain1, chain2 in zip(rota_struct1, rota_struct2):
                for r1, r2 in zip(chain1, chain2):
                    # print(self.structure, r1.get_id()[1], r2.get_id()[1], self.calculate_CA_distance(r1, r2), self.structure.state.name)
                    line = '{},{},{},{},{}\n'.format(self.structure, self.structure.state.name, round(self.calculate_CA_distance(r1, r2), 2), r1.get_id()[1], r2.get_id()[1])
                    self.line = line
                    return self.calculate_CA_distance(r1, r2)
        except:
            try:
                res1 = Residue.objects.get(protein_conformation=self.parent_prot_conf, display_generic_number__label=dgn(residue1, self.parent_prot_conf))
                res2 = Residue.objects.get(protein_conformation=self.parent_prot_conf, display_generic_number__label=dgn(residue2, self.parent_prot_conf))
                if self.structure_type=='refined':
                    pdb_data = self.structure.pdb_data.pdb
                elif self.structure_type=='hommod':
                    pdb_data = self.structure.pdb_data.pdb
                io = StringIO(pdb_data)
                struct = PDB.PDBParser(QUIET=True).get_structure('structure', io)[0]
                for chain in struct:
                    r1 = chain[res1.sequence_number]
                    r2 = chain[res2.sequence_number]
                    print(self.structure, r1.get_id()[1], r2.get_id()[1], self.calculate_CA_distance(r1, r2), self.structure.state.name)
                    line = '{},{},{},{},{}\n'.format(self.structure, self.structure.state.name, round(self.calculate_CA_distance(r1, r2), 2), r1.get_id()[1], r2.get_id()[1])
                    self.line = line
                    return self.calculate_CA_distance(r1, r2)

            except:
                print('Error: {} no matching rotamers ({}, {})'.format(self.structure.pdb_code.index, residue1, residue2))
                return False   

    def calculate_CA_distance(self, residue1, residue2):
        diff_vector = residue1['CA'].get_coord()-residue2['CA'].get_coord()
        return numpy.sqrt(numpy.sum(diff_vector * diff_vector))


class StructureSeqNumOverwrite():
    def __init__(self, structure):
        self.structure = structure
        path = os.sep.join([settings.DATA_DIR, 'structure_data','wt_pdb_lookup', '{}.json'.format(self.structure.pdb_code.index)])
        if os.path.isfile(path):
            with open(path, 'r') as lookup_file:
                self.lookup = json.load(lookup_file)
            self.wt_pdb_table, self.pdb_wt_table = OrderedDict(), OrderedDict()
            for i in self.lookup:
                self.wt_pdb_table[i['WT_POS']] = i['PDB_POS']
                self.pdb_wt_table[i['PDB_POS']] = i['WT_POS']
        else:
            self.lookup = OrderedDict()
            self.wt_pdb_table = OrderedDict()
            self.pdb_wt_table = OrderedDict()
            
    def seq_num_overwrite(self, overwrite_target):
        ''' Overwrites Residue object sequence numbers in GPCRDB
            @param overwrite_target: 'pdb' if converting pdb to wt, 'wt' if the other way around 
        '''
        resis = Residue.objects.filter(protein_conformation=self.structure.protein_conformation)
        if overwrite_target=='pdb':
            target_dict = self.pdb_wt_table
        elif overwrite_target=='wt':
            target_dict = self.wt_pdb_table
        for r in resis:
            if r.sequence_number in target_dict:
                r.sequence_number = int(target_dict[r.sequence_number])
                r.save()


class StructureBuildCheck():
    def __init__(self):
        with open(os.sep.join([settings.DATA_DIR, 'structure_data', 'annotation', 'xtal_segends.yaml']), 'r') as f:
            self.annotation_data = yaml.load(f, Loader=yaml.FullLoader)

    def check_rotamers(self, pdb):
        structure = Structure.objects.get(pdb_code__index=pdb)
        structure_residues = Residue.objects.filter(protein_conformation=structure.protein_conformation)
        wt_residues = Residue.objects.filter(protein_conformation__protein=structure.protein_conformation.protein.parent)
        key = structure.protein_conformation.protein.parent.entry_name+'_'+pdb
        try:
            annotation = self.annotation_data[key]
        except:
            raise Exception('Warning: {} not annotated'.format(pdb))
        errors = []
        for i in range(1,9):
            if annotation[str(i)+'e']!='-':
                anno_len = int(annotation[str(i)+'e']-annotation[str(i)+'b']+1)
                if i==8:
                    struct_len = len(structure_residues.filter(protein_segment__slug='H8'))
                    wt_len = len(wt_residues.filter(protein_segment__slug='H8'))
                else:
                    struct_len = len(structure_residues.filter(protein_segment__slug='TM'+str(i)))
                    wt_len = len(wt_residues.filter(protein_segment__slug='TM'+str(i)))
                if anno_len!=struct_len:
                    if wt_len!=struct_len:
                        errors.append(['H'+str(i), anno_len, struct_len])
        return errors


def update_template_source(template_source, keys, struct, segment, just_rot=False):
    ''' Update the template_source dictionary with structure info for backbone and rotamers.
    '''
    for k in keys:
        if just_rot==True:
            try:
                template_source[segment][k][1] = struct
            except:
                pass
        else:
            try:
                template_source[segment][k][0] = struct
                template_source[segment][k][1] = struct
            except:
                pass
    return template_source

def compare_and_update_template_source(template_source, segment, signprot_pdb_array, i, cgn, template_source_key, segs_for_alt_complex_struct, alt_complex_struct, main_structure):
    if cgn in signprot_pdb_array[segment]:
        if segment in segs_for_alt_complex_struct and signprot_pdb_array[segment][cgn]!='x':
            update_template_source(template_source, [str(template_source_key)], alt_complex_struct, segment)
        elif signprot_pdb_array[segment][cgn]!='x':
            update_template_source(template_source, [str(template_source_key)], main_structure, segment)
        else:
            update_template_source(template_source, [str(template_source_key)], None, segment)
    else:
        update_template_source(template_source, [str(template_source_key)], None, segment)
    return template_source


def right_rotamer_select(rotamer, chain=None):
    ''' Filter out compound rotamers.
    '''
    if len(rotamer)>1:
        for i in rotamer:
            if i.pdbdata.pdb.startswith('COMPND')==False:
                if chain!=None and chain==i.pdbdata.pdb[21]:
                    rotamer = i
                    break
    else:
        rotamer=rotamer[0]
    return rotamer



