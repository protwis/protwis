from django.core.management.base import BaseCommand
from protein.models import Protein
from ligand.models import AssayExperiment
from interaction.models import StructureLigandInteraction

from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem

# Disable the RDkit verbosity
RDLogger.DisableLog('rdApp.*')

class Command(BaseCommand):
    def handle(self, *args, **options):
        SIM_CUTOFF = 0.5
        PACT_CUTOFF = 7
        states = ["active", "inactive"]
        ligand_modalities = {
            "active": ["agonist", "agonist-partial"],
            "inactive": ["antagonist", "inverse-agonist"]
        }

        result_dictionary = {}
        # Process existing structures per ligand modality
        for state in states:
            result_dictionary[state] = {}
            ligand_structures = StructureLigandInteraction.objects.filter(ligand__ligand_type__slug='small-molecule', ligand_role__slug__in=ligand_modalities[state]).exclude(ligand__smiles = None).exclude(ligand__pdbe = None)\
                                .exclude(structure__structure_type__slug__startswith='af-').values_list('ligand_id', 'structure__protein_conformation__protein__parent__entry_name', 'structure__protein_conformation__protein__parent__family__slug', 'ligand__pdbe', 'structure__pdb_code__index')

            # Loop over all structures/ligands
            for pair in ligand_structures:
                # Add to structure-receptor slug to list
                if pair[2] not in result_dictionary[state]:
                    result_dictionary[state][pair[2]] = {}

                if pair[3] not in result_dictionary[state][pair[2]]:
                    result_dictionary[state][pair[2]][pair[3]] = [[pair[4]], "structure_match"]
                else:
                    result_dictionary[state][pair[2]][pair[3]][0].append(pair[4])

            # Loop over all structures/ligands again and match ligand to other receptors based on bioactivity data
            for pair in ligand_structures:
                data = list(AssayExperiment.objects.filter(ligand_id=pair[0], p_activity_value__gte=PACT_CUTOFF).values_list('protein__family__slug', flat=True).distinct())
                for slug in data:
                    # Skip same receptor as in structure
                    if slug == pair[2]:
                        continue

                    # Add to structure-receptor slug to list
                    if slug not in result_dictionary[state]:
                        result_dictionary[state][slug] = {}
                    elif slug in result_dictionary[state] and len(result_dictionary[state][slug]) > 0:
                        continue # skip receptors with actual structure matches

                    if pair[3] not in result_dictionary[state][slug]:
                        result_dictionary[state][slug][pair[3]] = [[pair[4]], "ligand_match"]
                    else:
                        result_dictionary[state][slug][pair[3]][0].append(pair[4])

            # Grab all unique smiles from the structure ligands
            structure_smiles = StructureLigandInteraction.objects.filter(ligand__ligand_type__slug='small-molecule', ligand_role__slug__in=ligand_modalities[state]).exclude(ligand__smiles = None).exclude(ligand__pdbe = None)\
                                .exclude(structure__structure_type__slug__startswith='af-').values_list('ligand__pdbe', 'structure__pdb_code__index', 'ligand__smiles').distinct()

            # Create Morgan Fingerprints from the SMILES
            structure_ligands = {}
            for ligand in structure_smiles:
                if ligand[0] not in structure_ligands:
                    structure_ligands[ligand[0]] = {}
                    #inactives[ligand[0]]["smiles"] = ligand[2]
                    mol = Chem.MolFromSmiles(ligand[2])
                    structure_ligands[ligand[0]]["mol"] = mol
                    structure_ligands[ligand[0]]["smiles"] = ligand[2]
                    structure_ligands[ligand[0]]["fingerprint"] = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
                    structure_ligands[ligand[0]]["structures"] = [ligand[1]]
                else:
                    structure_ligands[ligand[0]]["structures"].append(ligand[1])

            # For each receptor try to find ligands that are highly similar to ligands with a structure
            all_gpcrs = Protein.objects.filter(family__slug__startswith="00", sequence_type__slug='wt', species__id=1).prefetch_related("family")
            for receptor in all_gpcrs:
                # Do similarity search on bioactivity data ONLY if we don't have structure or ligand matches
                if receptor.family.slug not in result_dictionary[state] or len(result_dictionary[state][receptor.family.slug]) == 0:
                    if receptor.entry_name.startswith("ta2") or receptor.entry_name.startswith("t2r"):
                        continue
                    # Search ligands for a species-specific receptor
                    #bioactives = list(AssayExperiment.objects.filter(protein=receptor, p_activity_value__gte=PACT_CUTOFF).values_list('ligand__smiles', flat=True).distinct())
                    # Search ligands for receptor accross all species
                    bioactives = list(AssayExperiment.objects.filter(protein__family__slug = receptor.family.slug, p_activity_value__gte=PACT_CUTOFF).values_list('ligand__smiles', flat=True).distinct())
                    if len(bioactives) > 0:
                        # ENCODING the molecules
                        fingerprints = []
                        for ligand in bioactives:
                            try:
                                mol = Chem.MolFromSmiles(ligand)
                                fingerprints.append([AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048), ligand])
                            except TypeError:
                                continue

                        # Perform similarity search (NOTE - not optimized for batch processing)
                        if len(fingerprints) > 0:
                            current_pair = []
                            max_sim = 0
                            for struct_lig in structure_ligands:
                                for fp in fingerprints:
                                    sim = DataStructs.TanimotoSimilarity(structure_ligands[struct_lig]["fingerprint"], fp[0])
                                    if sim > max_sim:
                                        max_sim = sim
                                        current_pair = [struct_lig, fp]

                            if max_sim > SIM_CUTOFF:
                                struct_lig = current_pair[0]
                                fp = current_pair[1]
                                if receptor.family.slug not in result_dictionary[state]:
                                    result_dictionary[state][receptor.family.slug] = {}
                                result_dictionary[state][receptor.family.slug][struct_lig] = [structure_ligands[struct_lig]["structures"], fp[1], max_sim, "similarity_match"]
                                #DEBUG
                                #print(receptor, receptor.family.slug, result_dictionary[state][slug][struct_lig])
        return result_dictionary
        # CHECK RESULT
        #print(json.dumps(result_dictionary))

        # FEW STATS
        # matches = {"structure_match": 0, "ligand_match": 0, "similarity_match" : 0}
        # for state in result_dictionary:
        #     local_matches = {"structure_match": 0, "ligand_match": 0, "similarity_match" : 0}
        #     for slug in result_dictionary[state]:
        #         first_ligand = list(result_dictionary[state][slug].keys())[0]
        #         matches[result_dictionary[state][slug][first_ligand][-1]] += 1
        #         local_matches[result_dictionary[state][slug][first_ligand][-1]] += 1
        #     print(state, len(result_dictionary[state].keys()))
        #     print(state, local_matches)
        # print(matches)
