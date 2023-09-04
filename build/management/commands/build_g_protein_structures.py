from django.core.management.base import BaseCommand, CommandError
from build.management.commands.base_build import Command as BaseBuild
from django.conf import settings
from django.db import connection
from django.db import IntegrityError

from protein.models import (Protein, ProteinConformation, ProteinState, ProteinFamily, ProteinAlias,
        ProteinSequenceType, Species, Gene, ProteinSource, ProteinSegment)
from residue.models import (ResidueNumberingScheme, ResidueGenericNumber, Residue, ResidueGenericNumberEquivalent)
from signprot.models import SignprotComplex, SignprotStructure, SignprotStructureExtraProteins
from common.models import WebResource, WebLink, Publication
from structure.models import StructureType, StructureStabilizingAgent, PdbData, Rotamer
from structure.functions import get_pdb_ids, create_structure_rotamer, fetch_signprot_data, build_signprot_struct
from common.tools import test_model_updates
from protein.management.commands.blastp import CustomBlast

import re
from Bio import pairwise2, SeqIO
from collections import OrderedDict
import logging
import shlex, subprocess
from io import StringIO
from Bio.PDB import PDBParser, PPBuilder, PDBIO, Polypeptide
from Bio import pairwise2
import pprint
import json
import yaml
import urllib
import django.apps
import traceback
import sys, os
import datetime


AA = {"ALA":"A", "ARG":"R", "ASN":"N", "ASP":"D",
     "CYS":"C", "GLN":"Q", "GLU":"E", "GLY":"G",
     "HIS":"H", "ILE":"I", "LEU":"L", "LYS":"K",
     "MET":"M", "PHE":"F", "PRO":"P", "SER":"S",
     "THR":"T", "TRP":"W", "TYR":"Y", "VAL":"V",
     "YCM":"C", "CSD":"C", "TYS":"Y", "SEP":"S"} #non-standard AAs


class Command(BaseBuild):

    #Setting the variables for the test tracking of the model upadates
    tracker = {}
    all_models = django.apps.apps.get_models()[6:]
    test_model_updates(all_models, tracker, initialize=True)
    local_uniprot_dir = os.sep.join([settings.DATA_DIR, "g_protein_data", "uniprot"])
    local_uniprot_beta_dir = os.sep.join([settings.DATA_DIR, "g_protein_data", "uniprot_beta"])
    local_uniprot_gamma_dir = os.sep.join([settings.DATA_DIR, "g_protein_data", "uniprot_gamma"])
    with open(os.sep.join([settings.DATA_DIR, "g_protein_data", "g_protein_display_names.yaml"]), "r") as y:
        display_name_lookup = yaml.load(y, Loader=yaml.FullLoader)

    def add_arguments(self, parser):
        parser.add_argument('-p', '--proc', type=int, action='store', dest='proc', default=1, help='Number of processes to run')
        parser.add_argument("--purge_complex", default=False, action="store_true", help="Purge G protein complex structures from database")
        parser.add_argument("--purge_non_complex", default=False, action="store_true", help="Purge G protein non-complex structures from database")
        parser.add_argument("--only_signprot_structures", default=False, action="store_true", help="Only build SignprotStructure objects")
        parser.add_argument("-s", default=False, type=str, action="store", nargs="+", help="PDB codes to build")
        parser.add_argument("--debug", default=False, action="store_true", help="Debug mode")

    def handle(self, *args, **options):
        startTime = datetime.datetime.now()
        self.options = options
        if self.options["purge_complex"]:
            Residue.objects.filter(protein_conformation__protein__entry_name__endswith="_a", protein_conformation__protein__family__parent__parent__name="Alpha").delete()
            ProteinConformation.objects.filter(protein__entry_name__endswith="_a", protein__family__parent__parent__name="Alpha").delete()
            Protein.objects.filter(entry_name__endswith="_a", family__parent__parent__name="Alpha").delete()
            self.tracker = {}
            test_model_updates(self.all_models, self.tracker, initialize=True)
        if self.options["purge_non_complex"]:
            SignprotStructureExtraProteins.objects.all().delete()
            SignprotStructure.objects.all().delete()
            self.tracker = {}
            test_model_updates(self.all_models, self.tracker, initialize=True)

        if not self.options["only_signprot_structures"]:
            if self.options["s"]:
                self.scs = SignprotComplex.objects.filter(structure__pdb_code__index__in=[i.upper() for i in self.options["s"]])
            else:
                self.scs = SignprotComplex.objects.filter(protein__family__slug__startswith='100_001')
            self.prepare_input(self.options['proc'], self.scs, 1)
            test_model_updates(self.all_models, self.tracker, check=True)
        if not self.options["s"]:
            ### Build SignprotStructure objects from non-complex signprots
            g_prot_alphas = Protein.objects.filter(family__slug__startswith="100_001", accession__isnull=False)#.filter(entry_name="gnai1_human")
            complex_structures = SignprotComplex.objects.filter(protein__family__slug__startswith="100_001").values_list("structure__pdb_code__index", flat=True)
            for a in g_prot_alphas:
                pdb_list = get_pdb_ids(a.accession)
                for pdb in pdb_list:
                    if pdb not in complex_structures:
                        try:
                            data = fetch_signprot_data(pdb, a, os.listdir(self.local_uniprot_beta_dir), os.listdir(self.local_uniprot_gamma_dir))
                            if data:
                                ss = build_signprot_struct(a, pdb, data)
                                self.build_gprot_extra_proteins(a, ss, data)
                        except Exception as msg:
                            self.logger.error("SignprotStructure of {} {} failed\n{}: {}".format(a.entry_name, pdb, type(msg), msg))
            test_model_updates(self.all_models, self.tracker, check=True)
        if self.options["debug"]:
            print(datetime.datetime.now() - startTime)


    def main_func(self, positions, iterations, count, lock):
        # setting up processes
        if not positions[1]:
            scs = self.scs[positions[0]:]
        else:
            scs = self.scs[positions[0]:positions[1]]
        # while count.value<len(positions):
        #     with lock:
        #         sc = positions[count.value]
        #         count.value +=1
        for sc in scs:
            # Building protein and protconf objects for g protein structure in complex
            self.logger.info("Protein, ProteinConformation and Residue build for alpha subunit of {} is building".format(sc))
            try:
                # Alpha subunit
                try:
                    alpha_protein = Protein.objects.get(entry_name=sc.structure.pdb_code.index.lower()+"_a")
                except:
                    alpha_protein = Protein()
                    alpha_protein.entry_name = sc.structure.pdb_code.index.lower()+"_a"
                    alpha_protein.accession = None
                    alpha_protein.name = sc.structure.pdb_code.index.lower()+"_a"
                    alpha_protein.sequence = sc.protein.sequence
                    alpha_protein.family = sc.protein.family
                    alpha_protein.parent = sc.protein
                    alpha_protein.residue_numbering_scheme = sc.protein.residue_numbering_scheme
                    alpha_protein.sequence_type = ProteinSequenceType.objects.get(slug="mod")
                    alpha_protein.source = ProteinSource.objects.get(name="OTHER")
                    alpha_protein.species = sc.protein.species
                    alpha_protein.save()

                try:
                    alpha_protconf = ProteinConformation.objects.get(protein__entry_name=sc.structure.pdb_code.index.lower()+"_a")
                except:
                    alpha_protconf = ProteinConformation()
                    alpha_protconf.protein = alpha_protein
                    alpha_protconf.state = ProteinState.objects.get(slug="active")
                    alpha_protconf.save()

                pdbp = PDBParser(PERMISSIVE=True, QUIET=True)
                s = pdbp.get_structure("struct", StringIO(sc.structure.pdb_data.pdb))
                chain = s[0][sc.alpha]
                nums = []
                structure_seq = ''
                for res in chain:
                    if "CA" in res and res.id[0]==" ":
                        if sc.structure.pdb_code.index=='7RYC':
                            if res.get_id()[1]<1005:
                                continue
                        nums.append(res.get_id()[1])
                        structure_seq+=Polypeptide.three_to_one(res.get_resname())

                if self.options['debug']:
                    print('Annotated protein:')
                    print(sc.protein)
                    print('Structure seq:')
                    print(structure_seq)

                resis = Residue.objects.filter(protein_conformation__protein=sc.protein)
                num_i = 0
                temp_seq2 = ""
                pdb_num_dict = OrderedDict()
                # Create first alignment based on sequence numbers
                try:
                    for n in nums:
                        if sc.structure.pdb_code.index=="6OIJ" and n<30:
                            nr = n+6
                        elif sc.structure.pdb_code.index in ['7MBY', '7F9Y', '7F9Z'] and n>58:
                            nr = n-35
                        elif sc.structure.pdb_code.index in ['7EIB', '7F2O']:
                            nr = n-2
                        elif sc.structure.pdb_code.index in ['7P00'] and n>52:
                            nr = n-18
                        elif sc.structure.pdb_code.index=='7RYC':
                            nr = n-994
                        elif sc.structure.pdb_code.index in ['7W53','7W55','7W56','7W57','7WKD']:
                            nr = n-2
                        elif sc.structure.pdb_code.index=='7X9Y' and n>58:
                            nr = n-408
                        elif sc.structure.pdb_code.index in ['7WXU','7WY5'] and n>63:
                            nr = n-35
                        elif sc.structure.pdb_code.index=='7WY0':
                            if n<67:
                                nr = n+8
                            elif n>66:
                                nr = n-1
                        elif sc.structure.pdb_code.index=='7XW9' and n>52:
                            nr = n-50
                        elif sc.structure.pdb_code.index=='8H8J':
                            nr = n-3
                        else:
                            nr = n
                        pdb_num_dict[n] = [chain[n], resis.get(sequence_number=nr)]
                except Residue.DoesNotExist:
                    nr = resis[0].sequence_number
                    for n in nums:
                        nr = n-(n-nr)
                        pdb_num_dict[n] = [chain[n], resis.get(sequence_number=nr)]

                # Find mismatches
                mismatches = []
                for n, res in pdb_num_dict.items():
                    if AA[res[0].get_resname()]!=res[1].amino_acid:
                        mismatches.append(res)

                pdb_lines = sc.structure.pdb_data.pdb.split("\n")
                seqadv = []
                for l in pdb_lines:
                    if l.startswith("SEQADV"):
                        seqadv.append(l)
                mutations, shifted_mutations = OrderedDict(), OrderedDict()
                # Search for annotated engineered mutations in pdb SEQADV
                for s in seqadv:
                    line_search = re.search("SEQADV\s{1}[A-Z\s\d]{4}\s{1}([A-Z]{3})\s{1}([A-Z]{1})\s+(\d+)[\s\S\d]{5}([\s\S\d]{12})([A-Z]{3})\s+(\d+)(\s\S+)",s)
                    if line_search!=None:
                        if line_search.group(2)==sc.alpha:
                            if line_search.group(4).strip()==sc.protein.accession:
                                if line_search.group(3)==line_search.group(6):
                                    mutations[int(line_search.group(3))] = [line_search.group(1), line_search.group(5)]
                                else:
                                    shifted_mutations[int(line_search.group(3))] = [line_search.group(1), line_search.group(5), int(line_search.group(6))]
                            else:
                                # Exception for 6G79
                                if line_search.group(3)!=line_search.group(6) and "CONFLICT" in line_search.group(7):
                                    mutations[int(line_search.group(3))] = [line_search.group(1), line_search.group(5)]
                                # Exception for 5G53
                                if line_search.group(4).strip()!=sc.protein.accession:
                                    mutations[int(line_search.group(3))] = [line_search.group(1), line_search.group(5)]
                remaining_mismatches = []

                # Check and clear mismatches that are registered in pdb SEQADV as engineered mutation
                for m in mismatches:
                    num = m[0].get_id()[1]
                    if num in mutations:
                        if m[0].get_resname()!=mutations[num][0] and m[1].amino_acid!=AA[mutations[num][1]]:
                            remaining_mismatches.append(m)
                    elif num in shifted_mutations:
                        remaining_mismatches.append(m)
                    else:
                        remaining_mismatches.append(m)

                if self.options["debug"]:
                    print(sc)
                    print(mutations)
                    print(shifted_mutations)
                    print(mismatches)
                    print("======")
                    print(remaining_mismatches)
                    pprint.pprint(pdb_num_dict)

                no_seqnum_shift = ['6OY9', '6OYA', '6LPB', '6WHA', '7D77', '6XOX', '7L1U', '7L1V']

                # Check if HN is mutated to GNAI1 for the scFv16 stabilizer
                if sc.protein.entry_name!='gnai1_human' and len(remaining_mismatches)>0:
                    target_HN = resis.filter(protein_segment__slug='G.HN')
                    gnai1_HN = Residue.objects.filter(protein_conformation__protein__entry_name='gnai1_human', protein_segment__slug='G.HN')
                    pdb_HN_seq = ''
                    for num, val in pdb_num_dict.items():
                        if num<=target_HN.reverse()[0].sequence_number:
                            pdb_HN_seq+=Polypeptide.three_to_one(val[0].get_resname())
                    if self.options['debug']:
                        print('Checking if HN is gnai1_human')
                        print(pdb_HN_seq)
                        print(''.join(gnai1_HN.values_list('amino_acid', flat=True)))
                    gnai1_HN_seq = ''.join(gnai1_HN.values_list('amino_acid', flat=True))
                    if len(pdb_HN_seq)>0:
                        pw2 = pairwise2.align.localms(gnai1_HN_seq, pdb_HN_seq, 3, -4, -3, -1)
                        ref_seq, temp_seq = str(pw2[0][0]), str(pw2[0][1])
                        length, match = 0,0
                        for r, t in zip(ref_seq, temp_seq):
                            if self.options['debug']:
                                print(r,t)
                            if t!='-' and r!='-':
                                if r==t:
                                    match+=1
                                length+=1
                        identity = match/length*100
                        if self.options['debug']:
                            print(identity)
                        if identity>85 and length/len(temp_seq)>.5:
                            if sc.structure.pdb_code.index not in ['7DFL','7S8L','7S8P','7S8N','7MBY','7AUE','7XW9']:
                                no_seqnum_shift.append(sc.structure.pdb_code.index)
                            if self.options['debug']:
                                print('INFO: HN has {}% with gnai1_human HN, skipping seqnum shift correction'.format(round(identity)))
                        elif sc.structure.pdb_code.index in ['7KH0']:
                            no_seqnum_shift.append(sc.structure.pdb_code.index)

                # Mismatches remained possibly to seqnumber shift, making pairwise alignment to try and fix alignment
                if len(remaining_mismatches)>0 and sc.structure.pdb_code.index not in no_seqnum_shift:
                    ppb = PPBuilder()
                    seq = ""
                    for pp in ppb.build_peptides(chain, aa_only=False):
                        seq += str(pp.get_sequence())
                    seq = structure_seq
                    if sc.structure.pdb_code.index in ['7JVQ','7L1U','7L1V','7D68','7EZK']:
                        pw2 = pairwise2.align.localms(sc.protein.sequence, seq, 3, -4, -3, -1)
                    else:
                        pw2 = pairwise2.align.localms(sc.protein.sequence, seq, 2, -1, -.5, -.1)
                    ref_seq, temp_seq = str(pw2[0][0]), str(pw2[0][1])
                    if self.options['debug']:
                        print('default pw')
                        for i,j in zip(ref_seq, temp_seq):
                            print(i,j)
                        print('==========')

                    alignment_fragments = [i for i in str(pw2[0][1]).split('-') if i!='' and len(i)<10]

                    ### Try constricted alignment with higher gap penalty for H5 fragment structures
                    if len(alignment_fragments)>0 and len(seq)<50:
                        pw2 = pairwise2.align.localms(sc.protein.sequence, seq, 3, -4, -3, -1)
                        ref_seq, temp_seq = str(pw2[0][0]), str(pw2[0][1])
                        if self.options['debug']:
                            print('strict pw')
                            for i,j in zip(ref_seq, temp_seq):
                                print(i,j)
                            print('==========')
                        alignment_fragments = [i for i in str(pw2[0][1]).split('-') if i!='' and len(i)<10]

                    ### Chimera mapping to wt ###
                    if len(alignment_fragments)>0:
                        # blast chimeras to find best chimera match
                        cb = CustomBlast(os.sep.join([settings.STATICFILES_DIRS[0], 'blast', 'g_protein_chimeras']))
                        matched_chimera_key = cb.run(temp_seq)

                        # parsing gapped chimera fasta
                        chimeras = SeqIO.to_dict(SeqIO.parse(open(os.sep.join([settings.DATA_DIR, 'g_protein_data', 'g_protein_chimeras_gapped.fasta'])), "fasta"))
                        if self.options['debug']:
                            print('matched chimera', matched_chimera_key, chimeras[matched_chimera_key].seq)

                        # pairwise alignment to best match
                        chimera_pw2 = pairwise2.align.localms(chimeras[matched_chimera_key].seq, seq, 3, -4, -3, -.5)
                        ref_seq_chim, temp_seq_chim = str(chimera_pw2[0][0]), str(chimera_pw2[0][1])
                        if self.options['debug']:
                            print('chimera pw')
                            for i,j in zip(ref_seq_chim, temp_seq_chim):
                                print(i,j)
                            print('==========')
                        chimera_wt_key = matched_chimera_key.split('|')[0]

                        # reassign wt best match as ref
                        ref_seq = chimeras[chimera_wt_key].seq
                        temp_seq = temp_seq_chim
                    ##############################

                    wt_pdb_dict = OrderedDict()
                    pdb_wt_dict = OrderedDict()
                    j, k = 0, 0
                    if self.options["debug"]:
                        print('wt pw')
                        print(len(ref_seq), len(temp_seq))
                    for i, ref, temp in zip(range(0,len(ref_seq)), ref_seq, temp_seq):
                        if self.options["debug"]:
                            print(i, ref, temp) # alignment check
                        if ref!="-" and temp!="-":
                            wt_pdb_dict[resis[j]] = pdb_num_dict[nums[k]]
                            pdb_wt_dict[pdb_num_dict[nums[k]][0]] = resis[j]
                            j+=1
                            k+=1
                        elif ref=="-" and temp=="-":
                            pass
                        elif ref=="-":
                            wt_pdb_dict[i] = pdb_num_dict[nums[k]]
                            pdb_wt_dict[pdb_num_dict[nums[k]][0]] = i
                            k+=1
                        elif temp=="-":
                            wt_pdb_dict[resis[j]] = i
                            pdb_wt_dict[i] = resis[j]
                            j+=1

                    # Custom fix for 7JJO isoform difference
                    if sc.structure.pdb_code.index in ['7JJO', '7JOZ', '7AUE', '7EZK']:
                        pdb_num_dict = OrderedDict()
                        for wt_res, st_res in wt_pdb_dict.items():
                            if type(st_res)==type([]):
                                pdb_num_dict[wt_res.sequence_number] = [st_res[0], wt_res]
                    else:
                        for i, r in enumerate(remaining_mismatches):
                            # Adjust for shifted residue when residue is a match
                            if r[0].get_id()[1]-remaining_mismatches[i-1][0].get_id()[1]>1:
                                try:
                                    pdb_num_dict[r[0].get_id()[1]-1][1] = pdb_wt_dict[chain[r[0].get_id()[1]-1]]
                                except:
                                    print('Warning: Resnum {} not in structure {}'.format(r[0].get_id()[1]-1, sc.structure.pdb_code.index))
                            # Adjust for shifted residue when residue is mutated and it's logged in SEQADV
                            if r[0].get_id()[1] in shifted_mutations:
                                pdb_num_dict[r[0].get_id()[1]][1] = resis.get(sequence_number=shifted_mutations[r[0].get_id()[1]][2])
                            # Adjust for shift
                            else:
                                pdb_num_dict[r[0].get_id()[1]][1] = pdb_wt_dict[r[0]]
                        if sc.structure.pdb_code.index=='7JVQ':
                            pdb_num_dict[198][1] = Residue.objects.get(protein_conformation__protein=sc.protein, sequence_number=346)
                            pdb_num_dict[235][1] = Residue.objects.get(protein_conformation__protein=sc.protein, sequence_number=383)
                        elif sc.structure.pdb_code.index=='6PB0':
                            pdb_num_dict[205][1] = Residue.objects.get(protein_conformation__protein=sc.protein, sequence_number=205)
                        elif sc.structure.pdb_code.index=='7RYC':
                            pdb_num_dict[1198][1] = Residue.objects.get(protein_conformation__protein=sc.protein, sequence_number=314)
                            pdb_num_dict[1110][1] = Residue.objects.get(protein_conformation__protein=sc.protein, sequence_number=230)
                        elif sc.structure.pdb_code.index=='7P00':
                            pdb_num_dict[272][1] = Residue.objects.get(protein_conformation__protein=sc.protein, sequence_number=271)
                        elif sc.structure.pdb_code.index in ['7EIB','7F2O']:
                            pdb_num_dict[256][1] = Residue.objects.get(protein_conformation__protein=sc.protein, sequence_number=271)
                        elif sc.structure.pdb_code.index in ['7F9Y','7F9Z','7MBY']:
                            pdb_num_dict[289][1] = Residue.objects.get(protein_conformation__protein=sc.protein, sequence_number=271)

                ### Custom alignment fix for 6WHA, 7MBY mini-Gq/Gi/Gs chimera
                elif sc.structure.pdb_code.index in ['6WHA']:
                    if sc.structure.pdb_code.index=='6WHA':
                        ref_seq  = "MTLESIMACCLSEEAKEARRINDEIERQLRRDKRDARRELKLLLLGTGESGKSTFIKQMRIIHGSGYSDEDKRGFTKLVYQNIFTAMQAMIRAMDTLKIPYKYEHNKAHAQLVREVDVEKVSAFENPYVDAIKSLWNDPGIQECYDRRREYQLSDSTKYYLNDLDRVADPAYLPTQQDVLRVRVPTTGIIEYPFDLQSVIFRMVDVGGQRSERRKWIHCFENVTSIMFLVALSEYDQVLVESDNENRMEESKALFRTIITYPWFQNSSVILFLNKKDLLEEKIM--YSHLVDYFPEYDGP----QRDAQAAREFILKMFVDL---NPDSDKIIYSHFTCATDTENIRFVFAAVKDTILQLNLKEYNLV"
                        temp_seq = "----------VSAEDKAAAERSKMIDKNLREDGEKARRTLRLLLLGADNSGKSTIVK----------------------------------------------------------------------------------------------------------------------------------GIFETKFQVDKVNFHMFDVG-----RRKWIQCFNDVTAIIFVVDSSDYNR----------LQEALNDFKSIWNNRWLRTISVILFLNKQDLLAEKVLAGKSKIEDYFPEFARYTTPDPRVTRAKY-FIRKEFVDISTASGDGRHICYPHFTC-VDTENARRIFNDCKDIILQMNLREYNLV"



                    pdb_num_dict = OrderedDict()
                    temp_resis = [res for res in chain]
                    temp_i = 0
                    mapped_cgns = []
                    for i, aa in enumerate(temp_seq):
                        if aa!="-":
                            ref_split_on_gaps = ref_seq[:i+1].split("-")
                            ref_seqnum = i-(len(ref_split_on_gaps)-1)+1
                            res = resis.get(sequence_number=ref_seqnum)
                            if res.display_generic_number.label in mapped_cgns:
                                next_presumed_cgn = self.get_next_presumed_cgn(res)
                                if next_presumed_cgn:
                                    res = next_presumed_cgn
                                    while res and res.display_generic_number.label in mapped_cgns:
                                        res = self.get_next_presumed_cgn(res)
                                else:
                                    print("Warning: {} CGN does not exist. Incorrect mapping of {} in {}".format(next_presumed_cgn, chain[nums[temp_i]], sc.structure))
                            if res:
                                mapped_cgns.append(res.display_generic_number.label)
                            pdb_num_dict[nums[temp_i]] = [chain[nums[temp_i]], res]
                            temp_i+=1

                bulked_rotamers = []
                for key, val in pdb_num_dict.items():
                    # print(key, val) # sanity check
                    if not isinstance(val[1], int):
                        res_obj = Residue()
                        res_obj.sequence_number = val[0].get_id()[1]
                        res_obj.amino_acid = AA[val[0].get_resname()]
                        res_obj.display_generic_number = val[1].display_generic_number
                        res_obj.generic_number = val[1].generic_number
                        res_obj.protein_conformation = alpha_protconf
                        res_obj.protein_segment = val[1].protein_segment
                        res_obj.save()
                        rot = create_structure_rotamer(val[0], res_obj, sc.structure)
                        bulked_rotamers.append(rot)
                    else:
                        self.logger.info("Skipped {} as no annotation was present, while building for alpha subunit of {}".format(val[1], sc))
                if self.options["debug"]:
                    pprint.pprint(pdb_num_dict)
                Rotamer.objects.bulk_create(bulked_rotamers)
                self.logger.info("Protein, ProteinConformation and Residue build for alpha subunit of {} is finished".format(sc))
            except Exception as msg:
                if self.options["debug"]:
                    print("Error: ", sc, msg)
                self.logger.info("Protein, ProteinConformation and Residue build for alpha subunit of {} has failed".format(sc))

    @staticmethod
    def get_next_presumed_cgn(res):
        try:
            next_num = str(int(res.display_generic_number.label[-2:])+1)
            if len(next_num)==1:
                next_num = "0"+next_num
            next_cgn = res.display_generic_number.label[:-2]+next_num
            presumed_cgn = ResidueGenericNumber.objects.get(label=next_cgn)
            res = Residue.objects.filter(display_generic_number=presumed_cgn)[0]
            return res
        except ResidueGenericNumber.DoesNotExist:
            return False


    def build_gprot_extra_proteins(self, alpha_prot, signprot_structure, data):
        # Extra proteins
        # Alpha - ### A bit redundant, consider changing this in the future
        if data["alpha"]:
            alpha_sep = SignprotStructureExtraProteins()
            alpha_sep.wt_protein = alpha_prot
            alpha_sep.structure = signprot_structure
            alpha_sep.protein_conformation = ProteinConformation.objects.get(protein=alpha_prot)
            alpha_sep.display_name = self.display_name_lookup[alpha_prot.family.name]
            alpha_sep.note = None
            alpha_sep.chain = data["alpha_chain"]
            alpha_sep.category = "G alpha"
            cov = round(data["alpha_coverage"]/len(alpha_prot.sequence)*100)
            if cov>100:
                self.logger.warning("SignprotStructureExtraProtein Alpha subunit sequence coverage of {} is {}% which is longer than 100% in structure {}".format(alpha_sep, cov, signprot_structure))
                cov = 100
            alpha_sep.wt_coverage = cov
            alpha_sep.save()
        # Beta
        if data["beta"]:
            beta_prot = Protein.objects.get(accession=data["beta"])
            beta_sep = SignprotStructureExtraProteins()
            beta_sep.wt_protein = beta_prot
            beta_sep.structure = signprot_structure
            beta_sep.protein_conformation = ProteinConformation.objects.get(protein=beta_prot)
            beta_sep.display_name = self.display_name_lookup[beta_prot.name]
            beta_sep.note = None
            beta_sep.chain = data["beta_chain"]
            beta_sep.category = "G beta"
            beta_sep.wt_coverage = None
            beta_sep.save()
        # Gamma
        if data["gamma"]:
            gamma_prot = Protein.objects.get(accession=data["gamma"])
            gamma_sep = SignprotStructureExtraProteins()
            gamma_sep.wt_protein = gamma_prot
            gamma_sep.structure = signprot_structure
            gamma_sep.protein_conformation = ProteinConformation.objects.get(protein=gamma_prot)
            gamma_sep.display_name = self.display_name_lookup[gamma_prot.name]
            gamma_sep.note = None
            gamma_sep.chain = data["gamma_chain"]
            gamma_sep.category = "G gamma"
            gamma_sep.wt_coverage = None
            gamma_sep.save()
        self.logger.info("Created SignprotStructure: {}".format(signprot_structure.pdb_code))
