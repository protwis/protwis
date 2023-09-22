from build.management.commands.base_build import Command as BaseBuild
from django.conf import settings

from signprot.models import SignprotComplex
from structure.models import StructureExtraProteins, Structure
from structure.functions import ParseStructureCSV
from protein.models import Protein, ProteinConformation
from residue.models import Residue
from common.tools import test_model_updates
from common.definitions import G_PROTEIN_DISPLAY_NAME as g_prot_dict, ARRESTIN_DISPLAY_NAME as arrestin_dict
from interaction.models import StructureLigandInteraction

from Bio.PDB import PDBParser
from io import StringIO
import os
import math
import yaml
import django.apps

class Command(BaseBuild):
    help = "Build StructureExtraProteins data"

    #Setting the variables for the test tracking of the model upadates
    tracker = {}
    all_models = django.apps.apps.get_models()[6:]
    test_model_updates(all_models, tracker, initialize=True)

    def add_arguments(self, parser):
        super(Command, self).add_arguments(parser=parser)
        parser.add_argument('--purge', help='Purge existing entries in GPCRdb', default=False, action='store_true')

    def handle(self, *args, **options):
        if options['purge']:
            StructureExtraProteins.objects.all().delete()
            self.tracker = {}
            test_model_updates(self.all_models, self.tracker, initialize=True)

        self.psc = ParseStructureCSV()
        self.psc.parse_g_proteins()
        self.psc.parse_arrestins()

        self.build_g_protein_heterotrimer()
        test_model_updates(self.all_models, self.tracker, check=True)
        self.build_arrestin_extra_protein()
        test_model_updates(self.all_models, self.tracker, check=True)

    def build_g_protein_heterotrimer(self):
        sc = SignprotComplex.objects.filter(protein__family__slug__startswith='100')
        for s in sc:
            # Alpha
            sep = StructureExtraProteins()
            sep.wt_protein = s.protein
            sep.structure = s.structure
            sep.protein_conformation = ProteinConformation.objects.get(protein__entry_name=s.structure.pdb_code.index.lower()+'_a')
            sep.display_name = g_prot_dict[s.protein.family.name]
            if 'note'!='' in self.psc.structures[s.structure.pdb_code.index]['g_protein']['note']:
                sep.note = self.psc.structures[s.structure.pdb_code.index]['g_protein']['note']
            else:
                sep.note = None
            sep.chain = s.alpha
            sep.category = 'G alpha'
            wt_resis_len = len(Residue.objects.filter(protein_conformation__protein=s.protein))
            struct_resis_len = len(Residue.objects.filter(protein_conformation=sep.protein_conformation))
            if s.structure.pdb_code.index=='6E67':
                struct_resis_len = 175
            sep.wt_coverage = round(struct_resis_len/wt_resis_len*100)
            sep.save()
            # Beta
            if s.beta_protein:
                beta_sep = StructureExtraProteins()
                beta_sep.wt_protein = s.beta_protein
                beta_sep.structure = s.structure
                beta_sep.protein_conformation = ProteinConformation.objects.get(protein=s.beta_protein)
                beta_sep.display_name = g_prot_dict[s.beta_protein.name]
                beta_sep.note = None
                beta_sep.chain = s.beta_chain
                beta_sep.category = 'G beta'
                beta_sep.wt_coverage = None
                beta_sep.save()
            # Gamma
            if s.gamma_protein:
                gamma_sep = StructureExtraProteins()
                gamma_sep.wt_protein = s.gamma_protein
                gamma_sep.structure = s.structure
                gamma_sep.protein_conformation = ProteinConformation.objects.get(protein=s.gamma_protein)
                gamma_sep.display_name = g_prot_dict[s.gamma_protein.name]
                gamma_sep.note = None
                gamma_sep.chain = s.gamma_chain
                gamma_sep.category = 'G gamma'
                gamma_sep.wt_coverage = None
                gamma_sep.save()

    def build_arrestin_extra_protein(self):
        for pdb, vals in self.psc.structures.items():
            if 'arrestin' in vals:
                sep = StructureExtraProteins()
                sep.wt_protein = Protein.objects.get(entry_name=vals['arrestin']['protein'])
                sep.structure = Structure.objects.get(pdb_code__index=pdb)
                sep.protein_conformation = None
                sep.note = vals['arrestin']['note']
                sep.chain = vals['arrestin']['chain']
                sep.category = 'Arrestin'
                sep.display_name = arrestin_dict[vals['arrestin']['protein']]
                struct_resis_len = len([i for i in PDBParser(PERMISSIVE=True, QUIET=True).get_structure('struct', StringIO(sep.structure.pdb_data.pdb))[0][sep.chain]])
                if pdb in ['4ZWJ', '5DGY']:
                    struct_resis_len = 347
                elif pdb=='5W0P':
                    struct_resis_len = 364
                wt_resis_len = len(Residue.objects.filter(protein_conformation__protein=sep.wt_protein))
                sep.wt_coverage = round(struct_resis_len/wt_resis_len*100)
                sep.save()

    # Deprecated
    def build_from_notes(self):
        for struct, vals in self.notes.items():
            if 'category' in vals:
                sep = StructureExtraProteins()
                if vals['category']=='G alpha':
                    try:
                        wt_protein = Protein.objects.get(entry_name=vals['prot'].lower()+'_bovin')
                    except Protein.DoesNotExist:
                        try:
                            wt_protein = Protein.objects.get(entry_name=vals['prot'].lower()+'_human')
                        except Protein.DoesNotExist:
                            wt_protein = vals['prot']
                    sep.display_name = g_prot_dict[vals['prot']]
                elif vals['category']=='Arrestin':
                    wt_protein = Protein.objects.get(entry_name=vals['prot'].lower())
                    sep.display_name = arrestin_dict[vals['prot']]
                elif vals['category']=='Antibody':
                    wt_protein = None
                    sep.display_name = 'Antibody'

                sep.wt_protein = wt_protein
                sep.structure = Structure.objects.get(pdb_code__index=struct)
                sep.protein_conformation = None
                sep.note = vals['note']
                sep.chain = vals['chain']
                sep.category = vals['category']
                if vals['category']=='Antibody':
                    sep.wt_coverage = None
                else:
                    try:
                        wt_resis = Residue.objects.filter(protein_conformation__protein=wt_protein)
                        sep.wt_coverage = round(vals['length']/len(wt_resis)*100)
                    except Protein.DoesNotExist:
                        sep.wt_coverage = None

                sep.save()
