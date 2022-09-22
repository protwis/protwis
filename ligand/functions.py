#mccabe complexity: ["error", 31]
import math

from django.utils.text import slugify
from django.db import IntegrityError
from django.db.models import Q

#from chembl_webresource_client import new_client
from common.models import WebResource, Publication
from ligand.models import Ligand, LigandType, BiasedData, Endogenous_GTP, BalancedLigands
from protein.models import Protein

# def get_or_make_ligand(ligand_id, type_id, name = None, pep_or_prot = None):
#     if type_id=='PubChem CID' or type_id=='SMILES':
#         if type_id=='PubChem CID':
#             pubchem_lookup_value = 'cid'
#         elif type_id=='SMILES':
#             pubchem_lookup_value = 'smiles'
#
#         try:
#             web_resource = WebResource.objects.get(slug='pubchem')
#         except:
#             # abort if pdb resource is not found
#             raise Exception('PubChem resource not found, aborting!')
#         if name:
#             ligand_name = name
#         else:
#             ligand_name = False
#
#         try:
#             # if this name is canonical and it has a ligand record already
#             if (ligand_name==False):
#
#                 l = None
#                 ls = Ligand.objects.filter(canonical=True,
#                    ids__web_resource=web_resource,
#                    ids__index=ligand_id)
#
#                 for ligand in ls:
#                     l = ligand
#                     #print (l)
#                     break
#                 if l == None:
#                     l = Ligand.objects.get(canonical=True,
#                     ids__web_resource=web_resource,
#                     ids__index=ligand_id)
#
#             else:
#                l = Ligand.objects.get(name=ligand_name, canonical=True,
#                    ids__web_resource=web_resource,
#                    ids__index=ligand_id)
#
#             #l = Ligand.objects.get(name=ligand_name, canonical=True,
#             #    ids__web_resource=web_resource,
#             #    ids__index=ligand_id)
#             #
#         except Ligand.DoesNotExist:
#             try:
#                 # if exists under different name
#                 l_canonical = Ligand.objects.get(ids__web_resource=web_resource,
#                     ids__index=ligand_id, canonical=True)
#                 #print (created)
#                 try:
#                     l, created = Ligand.objects.get_or_create(properities = l_canonical.properities,
#                         name = ligand_name, canonical = False)
#                 except IntegrityError:
#                     l = Ligand.objects.get(properities = l_canonical.properities,
#                         name = ligand_name, canonical = False)
#             except Ligand.DoesNotExist:
#                 # fetch ligand from pubchem
#                 default_ligand_type = 'Small molecule'
#                 lt, created = LigandType.objects.get_or_create(slug=slugify(default_ligand_type),
#                     defaults={'name': default_ligand_type})
#                 l = Ligand()
#                 #print (ligand_name)
#                 l = l.load_from_pubchem(pubchem_lookup_value, ligand_id, lt, ligand_name)
#                 #print (l)
#                 if l == None and type_id=='SMILES': #insert manually if smiles and unfound in pubchem
#                     try:
#                         l = Ligand.objects.get(name=ligand_name, canonical=True,
#                                                 smiles=ligand_id)
#                     except Ligand.DoesNotExist:
#                         try:
#                             l = Ligand.objects.get(name__startswith=ligand_name, canonical=True,properities__smiles=ligand_id) #if no properities exist
#                         except Ligand.DoesNotExist:
#                             try:
#                                 l = Ligand.objects.get(name=ligand_name, canonical=True,properities__smiles=None) #if no properities exist
#                                 l.properities.smiles = ligand_id
#                                 l.properities.save()
#                                 l.save()
#                             except Ligand.DoesNotExist:
#                                 ## now insert a new ligand, but first make sure name is unique
#                                 if Ligand.objects.filter(name=ligand_name).exists():
#                                     ls = Ligand.objects.filter(name__startswith=ligand_name, canonical=True).order_by("pk")
#                                     last = ""
#                                     for l_temp in ls:
#                                         if last == "":
#                                             last = l_temp.name
#                                         try:
#                                             last = int(l_temp.name.split("_")[-1])
#                                         except ValueError:
#                                             continue
#
#                                     # TODO: better solution
#                                     if last == ligand_name or last == "": #no addition yet
#                                         ligand_name = ligand_name +"_1"
#                                     else:
#                                         ligand_name = ligand_name +"_"+str(int(last)+1)
#                                 l = Ligand()
#                                 l.name = ligand_name
#                                 # lp = LigandProperities()
#                                 # lp.smiles = ligand_id
#                                 # lp.ligand_type = lt
#                                 # lp.save()
#                                 # l.properities = lp
#                                 l.canonical = True #maybe false, but that would break stuff.
#                                 l.ambigious_alias = False
#                                 try:
#                                     l.save()
#                                 except IntegrityError:
#                                     l = Ligand.objects.get(name=ligand_name, canonical=True)
#
#     elif name:
#
#         # if this name is canonical and it has a ligand record already
#         if Ligand.objects.filter(name=name, canonical=True).exists():
#             l = Ligand.objects.get(name=name, canonical=True)
#
#         # if this matches an alias that only has "one" parent canonical name - eg distinct
#         elif Ligand.objects.filter(name=name, canonical=False,
#             ambigious_alias=False).exists():
#             l = Ligand.objects.get(name=name, canonical=False, ambigious_alias=False)
#
#         # if this matches an alias that only has several canonical parents, must investigate, start
#         # with empty.
#         elif Ligand.objects.filter(name=name, canonical=False,
#             ambigious_alias=True).exists():
#             # lp = LigandProperities()
#             # lp.save()
#             l = Ligand()
#             # l.properities = lp
#             l.name = name
#             l.canonical = False
#             l.ambigious_alias = True
#             l.save()
#             l.load_by_name(name)
#
#         # if neither a canonical or alias exists, create the records. Remember to check for
#         # canonical / alias status.
#         else:
#             # lp = LigandProperities()
#             # lp.save()
#             l = Ligand()
#             # l.properities = lp
#             l.name = str(name)
#             l.canonical = True
#             l.ambigious_alias = False
#             try:
#                 l.save()
#                 l.load_by_name(str(name))
#             except IntegrityError:
#                 l = Ligand.objects.get(name=str(name), canonical=True)
#             #if provided, update the ligand_type field of properities
#             #with correct labeling as peptide or protein
#             if pep_or_prot:
#                 l.properities.ligand_type = LigandType.objects.get(name = pep_or_prot)
#     else:
#         l = None
#
#     return l

#def fetch_chembl_refs(lig_chembl_id, target_accesion):

#    target_id = new_client.target.filter(accession=target_accesion)

#    assay = new_client.assay.filter(target_id=target_id, compound=lig_chembl_id)

#    refs = [x['document_chembl_id'] for x in assay]
#    #https://www.ebi.ac.uk/chembl/doc/inspect/CHEMBL2766014

#    return refs

################################################################################
#################### On The Fly Calculation BLOCK ##############################
################################################################################

def assess_reference(data_dict, user=False):
    gtp_info = {}
    reference = {}
    tested = {}
    skip = False
    reference_ligand = None
    checks = ['Principal', 'Secondary', 'None', None]
    #if reference ligand is not provided by user
    if user == False:
        receptor_id = data_dict[list(data_dict.keys())[0]]['receptor_id']
        lig_ids = set([data_dict[assay]['ligand_id'] for assay in data_dict])
        for status in checks:
            endo_objs = list(Endogenous_GTP.objects.filter(Q(ligand_species_id__common_name="Human") | Q(ligand_species_id__isnull=True),
                                                           ligand__in = lig_ids,
                                                           receptor = receptor_id,
                                                           endogenous_status=status).values_list("ligand", flat=True).distinct())

            if len(endo_objs) > 0:
                if status not in ['None', None]:
                    #if either principal or secondary check numerosity
                    data = list(Endogenous_GTP.objects.filter(Q(ligand_species_id__common_name="Human") | Q(ligand_species_id__isnull=True),
                                                         receptor=receptor_id,
                                                         ligand__in = endo_objs,
                                                         endogenous_status=status,
                                                        ).values_list("ligand_id", "pec50", "pKi").distinct())
                    if len(data) == 1:
                        #if single principal ligand break the cycle and define the reference ligand
                        reference_ligand = data[0][0]
                        break
                    else:
                        #if multiple principals, assess best based on pEC50, than pKi when no pEC50 is provided
                        pec50 = -100
                        pki = -100
                        for entry in data:
                            try:
                                if float(entry[1].split(' | ')[-1]) > pec50:
                                    reference_ligand = entry[0]
                            except ValueError:
                                try:
                                    if float(entry[2].split(' | ')[-1]) > pki:
                                        reference_ligand = entry[0]
                                except ValueError:
                                    continue
                        if status == "Principal" and reference_ligand:
                            break
                else:
                #check for highest ranking of the endogenous ligands
                    for endo_id in endo_objs:
                        data = Endogenous_GTP.objects.filter(Q(ligand_species_id__common_name="Human") | Q(ligand_species_id__isnull=True),
                                                            receptor=receptor_id,
                                                            ligand=endo_id,
                                                            ).values_list("potency_ranking").distinct()

                        try:
                            gtp_info[endo_id] = sum(data, ())[0]
                        except IndexError:
                            pass

                    if len(gtp_info) != 0: #check for actual endogenous ligands, otherwise skip pub
                        try:
                            reference_ligand = min(gtp_info, key=gtp_info.get)
                        except TypeError:
                            print("You have encountered the only exception, hooray!")
                            #grabbing the first ligand among the available as reference
                            reference_ligand = list(gtp_info.keys())[0]

        if not reference_ligand:
            skip = True
            # print("Skipping publication, no reference ligand tested")
            return reference, tested, skip

    #if reference ligand is provided by the user
    else:
        reference_ligand = int(user)


    for assay in data_dict.keys():
        if data_dict[assay]['ligand_id'] == reference_ligand:
            reference[assay] = data_dict[assay]
            reference_name = data_dict[assay]['ligand_name']
        else:
            tested[assay] =  data_dict[assay]

    for assay in tested.keys():
        tested[assay]['Reference_ligand'] = reference_name

    return reference, tested, skip

#check if that publication has balanced ligands, otherwise skip and go to next one
def retrieve_balanced_references(receptor, pub, subtype):
    skip = False
    hits = BalancedLigands.objects.filter(receptor_id=receptor,
                                          publication_id=pub,
                                          subtype_balanced=subtype).values_list(
                                                                    "first_pathway",
                                                                    "second_pathway",
                                                                    "ligand_id",
                                                                    "ligand_id__name")
    if len(hits) == 0:
        skip = True
    balanced_refs = {}
    for drug in hits:
        cmp = drug[0] + ' - ' + drug[1]
        balanced_refs[cmp] = [drug[2], drug[3]]

    return balanced_refs, skip

def assess_pathway_pairs(ranking, data, balanced_refs, subtype):
        #assess the pathway pairs in the tested ligands
        pairs = {}
        level = 'primary_effector_family'
        if subtype:
            level = 'primary_effector_subtype'
        #make all the comparison to strongest ranked pathway
        for drug in ranking:
            for i in range(1, len(ranking[drug])):
                #we don't want empty comparisons (manly in case of subtype)
                if (data[ranking[drug][0]][level] != None) and (data[ranking[drug][i]][level] != None):
                    comparison = data[ranking[drug][0]][level] + ' - ' + data[ranking[drug][i]][level]
                    if comparison not in pairs.keys():
                        pairs[comparison] = []
                    pairs[comparison].append(drug)
        matched_keys = set(list(balanced_refs.keys())).intersection(set(list(pairs.keys())))
        #now i need to create the reference and tested dicts based on the entries
        #and also the comparisons that are tied to the balanced reference and the pathway pair
        reference = {}
        tested = {}
        for key in matched_keys:
            balanced_id = balanced_refs[key][0]
            balanced_name = balanced_refs[key][1]
            if balanced_id in pairs[key]:
                pairs[key].remove(balanced_id)
            for entry in ranking[balanced_id]:
                reference[entry] = data[entry]
            for drug in pairs[key]:
                for entry in ranking[drug]:
                    tested[entry] = data[entry]
                    tested[entry]['Reference_ligand'] = balanced_name

        return reference, tested

def assess_comparisons(reference, tested):
    # Do we need to add these fields for comparison's sake?
    # 'receptor_isoform', #this can be null
    # 'active_receptor_complex', #this can be null
    common_traits = ['tissue',
                     'species',
                     'primary_effector_family',
                     'primary_effector_subtype',
                     'experiment',
                     'molecule_1',
                     'molecule_2',
                     'measured_process',
                     'assay_type']

    comparisons = {}
    for assay in reference:
        comparisons[assay] = []
        for test in tested:
            shared = list({x:reference[assay][x] for x in reference[assay] if x in tested[test] and reference[assay][x] == tested[test][x]}.keys())
            if set(common_traits).issubset(shared):
                comparisons[assay].append(test)
    return comparisons

def calculate_first_delta(comparisons, reference, tested, subtype=False):
    #Delta = log(Emax/EC50) of tested - log(Emax/EC50) of reference
    skip = False
    if subtype == False:
        find_best_subtype(comparisons, reference, tested)
    for ref in list(comparisons.keys()):
        r_emax, r_ec50, r_logtauka = reference[ref]['Emax'], reference[ref]['EC50'], reference[ref]['Tau_KA']
        try:
            r_logemaxec50 = math.log((r_emax/r_ec50),10)
        except (ValueError, TypeError, ZeroDivisionError):
            #excepting situations where one value is None or 0.0 [ValueError, TypeError, ZeroDivisionError]
            r_logemaxec50 = None
        #Check for lacking of ALL values, then break cycle and skip publication
        if (r_logemaxec50 == None) and (r_logtauka == None):
            skip = True
            return comparisons, tested, skip

        for test in comparisons[ref]:
            t_emax, t_ec50, t_logtauka = tested[test]['Emax'], tested[test]['EC50'], tested[test]['Tau_KA']
            try:
                t_logemaxec50 = math.log((t_emax/t_ec50),10)
            except (ValueError, TypeError, ZeroDivisionError):
                #excepting situations where one value is None or 0.0 [ValueError, TypeError, ZeroDivisionError]
                t_logemaxec50 = None
            #Here assess if both values of tested ligand are negative
            #and then check for qualitative values and report it
            if (t_logtauka == None) and (t_logemaxec50 == None):
                delta_logtauka = None
                delta_logemaxec50 = None
                continue
            #Calculate delta
            try:
                delta_logtauka = round((t_logtauka - r_logtauka), 3)
            except TypeError:
                delta_logtauka = None
            tested[test]['Delta_log(Tau/KA)'] = delta_logtauka
            try:
                delta_logemaxec50 = round((t_logemaxec50 - r_logemaxec50), 3)
            except TypeError:
                delta_logemaxec50 = None
            tested[test]['Delta_log(Emax/EC50)'] = delta_logemaxec50
            #Adding the section ragarding reference values
            tested[test]['Reference_Emax'] = r_emax
            tested[test]['Reference_EC50'] = r_ec50
            tested[test]['Reference_Tau/KA'] = r_logtauka
    return comparisons, tested, skip

def find_best_subtype(comparisons, reference, tested):
    families = {}
    to_be_deleted = []
    for ref in list(comparisons.keys()):
        #assessing the values
        r_emax, r_ec50 = reference[ref]['Emax'], reference[ref]['EC50']
        try:
            r_logemaxec50 = math.log((r_emax/r_ec50),10)
        except (ValueError, TypeError, ZeroDivisionError):
            #setting a largely negative value insted of None for comparisons sake
            #excepting situations where one value is None or 0.0 [ValueError, TypeError, ZeroDivisionError]
            r_logemaxec50 = -100

        if reference[ref]['primary_effector_family'] not in families.keys():
            #adding new family and values
            families[reference[ref]['primary_effector_family']] = [r_logemaxec50, ref]
        else:
            #updating with most relevant subtype
            if r_logemaxec50 > families[reference[ref]['primary_effector_family']][0]:
                #Adding the obsolete tests by comparing old key to updated key
                to_be_deleted += list(set(comparisons[families[reference[ref]['primary_effector_family']][1]]) - set(comparisons[ref]))
                #Delete obsolete key from comparisons register
                del comparisons[families[reference[ref]['primary_effector_family']][1]]
                #update reference key
                families[reference[ref]['primary_effector_family']] = [r_logemaxec50, ref]
            else:
                #Adding the obsolete tests by comparing old key to updated key
                to_be_deleted += list(set(comparisons[ref]) - set(comparisons[families[reference[ref]['primary_effector_family']][1]]))
                #Delete obsolete key from comparisons register
                del comparisons[ref]

    #get unique obsolete test keys and remove them
    to_be_deleted = list(set(to_be_deleted))
    for test in to_be_deleted:
        del tested[test]

def calculate_second_delta(comparisons, tested, rank_method, subtype=False, pathway=False):
    ranking = assess_pathway_preferences(comparisons, tested, rank_method, subtype, pathway)
    #STEPS
    for drug in ranking.keys():
        #perform analysis only when we have multiple pathways to actually compare
        if len(ranking[drug]) > 1:
        #Set reference pathway (first in list)
            path1 = ranking[drug][0]
            tested[path1]['Pathway Rank'] = 'P1'
            reference_emax = tested[path1]['Emax']
            reference_ec50 = tested[path1]['EC50']
            path_count = 1
            for test in ranking[drug][1:]:

                if subtype:
                    tested[test]['P1'] = tested[path1]['primary_effector_subtype']
                    #Match pathway levels and primary family for sake of subtype investigation
                    conditions = (tested[test]['pathway_level'] == tested[path1]['pathway_level']) and (tested[test]['primary_effector_family'] == tested[path1]['primary_effector_family'])
                else:
                    tested[test]['P1'] = tested[path1]['primary_effector_family']
                    #Match pathway levels + skip matching if Arrestin involved
                    conditions = (tested[test]['pathway_level'] == tested[path1]['pathway_level']) or ('Arrestin' in [tested[path1]['primary_effector_family'], tested[test]['primary_effector_family']])

                if conditions:
                    path_count +=1
                    tested[test]['Pathway Rank'] = 'P'+str(path_count)
                    if subtype:
                        tested[test]['P1'] = tested[path1]['primary_effector_subtype']
                    if pathway:
                        try:
                            delta_logtauka = round(tested[path1]['Tau_KA'] - tested[test]['Tau_KA'], 3)
                        except TypeError:
                            if tested[test]['qualitative_activity'] == 'No activity':
                                delta_logtauka = 'Full Bias'
                            else:
                                delta_logtauka = None
                        try:
                            tested[path1]['log(Emax/EC50)'] = round(math.log((tested[path1]['Emax']/tested[path1]['EC50']),10), 3)
                        except (TypeError, ValueError):
                            tested[path1]['log(Emax/EC50)'] = None
                        try:
                            tested[test]['log(Emax/EC50)'] = round(math.log((tested[test]['Emax']/tested[test]['EC50']),10), 3)
                        except (TypeError, ValueError):
                            tested[test]['log(Emax/EC50)'] = None
                        try:
                            delta_logemaxec50 = round(tested[path1]['log(Emax/EC50)'] - tested[test]['log(Emax/EC50)'], 3)
                        except TypeError:
                            if tested[test]['qualitative_activity'] == 'No activity':
                                delta_logemaxec50 = 'Full Bias'
                            elif tested[path1]['log(Emax/EC50)'] == None:
                                delta_logemaxec50 = None
                            elif tested[test]['log(Emax/EC50)'] == None:
                                delta_logemaxec50 = None

                        tested[test]['Delta_log(Tau/KA)'] = delta_logtauka
                        tested[test]['Delta_log(Emax/EC50)'] = delta_logemaxec50
                        tested[test]['Reference_Emax'] = reference_emax
                        tested[path1]['Reference_Emax'] = reference_emax
                        tested[test]['Reference_EC50'] = reference_ec50
                        tested[path1]['Reference_EC50'] = reference_ec50
                    else:
                        try:
                            deltadelta_logtauka = round(tested[path1]['Delta_log(Tau/KA)'] - tested[test]['Delta_log(Tau/KA)'], 3)
                        except (TypeError, KeyError):
                            if tested[test]['qualitative_activity'] == 'No activity':
                                deltadelta_logtauka = 'Full Bias'
                            elif 'Delta_log(Tau/KA)' not in tested[path1] or tested[path1]['Delta_log(Tau/KA)'] == None:
                                deltadelta_logtauka = None
                            elif 'Delta_log(Tau/KA)' not in tested[test] or tested[test]['Delta_log(Tau/KA)'] == None:
                                deltadelta_logtauka = None
                        try:
                            deltadelta_logemaxec50 = round(tested[path1]['Delta_log(Emax/EC50)'] - tested[test]['Delta_log(Emax/EC50)'], 3)
                            bias_factor = round(10**deltadelta_logemaxec50, 3)
                        except (TypeError, KeyError):
                            bias_factor = None
                            if tested[test]['qualitative_activity'] == 'No activity':
                                deltadelta_logemaxec50 = 'Full Bias'
                            elif 'Delta_log(Emax/EC50)' not in tested[path1] or tested[path1]['Delta_log(Emax/EC50)'] == None:
                                deltadelta_logemaxec50 = None
                            elif 'Delta_log(Emax/EC50)' not in tested[test] or tested[test]['Delta_log(Emax/EC50)'] == None:
                                deltadelta_logemaxec50 = None
                        tested[test]['DeltaDelta_log(Tau/KA)'] = deltadelta_logtauka
                        tested[test]['DeltaDelta_log(Emax/EC50)'] = deltadelta_logemaxec50
                        tested[test]['Bias factor'] = bias_factor
                else:
                    del tested[test]
    return tested

def assess_pathway_preferences(comparisons, tested, rank_method, subtype=False, pathway=False):
    pathway_preference = {}
    if pathway:
        comparisons, tested = find_best_pathway_family(comparisons, tested)
    #calculate values for ranking (or replace with qualitative activity when missing)
    for assay in comparisons.keys():
        for test in comparisons[assay]:
            if tested[test]['ligand_id'] not in pathway_preference.keys():
                pathway_preference[tested[test]['ligand_id']] = {}
            if pathway:
                path_label = tested[test]['primary_effector_family']
                Tau_KA = tested[test]['Tau_KA']
                try:
                    Emax_EC50 = round(math.log((tested[test]['Emax']/tested[test]['EC50']),10), 3)
                except (TypeError, ValueError):
                    Emax_EC50 = None
            elif subtype:
                path_label = str(tested[test]['primary_effector_family']) + ' - ' + str(tested[test]['primary_effector_subtype'])
                try:
                    Tau_KA = tested[test]['Delta_log(Tau/KA)']
                except KeyError:
                    Tau_KA = tested[test]['qualitative_activity']
                try:
                    Emax_EC50 = tested[test]['Delta_log(Emax/EC50)']
                except KeyError:
                    Emax_EC50 = tested[test]['qualitative_activity']
            else:
                path_label = tested[test]['primary_effector_family']
                try:
                    Tau_KA = tested[test]['Delta_log(Tau/KA)']
                except KeyError:
                    Tau_KA = tested[test]['qualitative_activity']
                try:
                    Emax_EC50 = tested[test]['Delta_log(Emax/EC50)']
                except KeyError:
                    Emax_EC50 = tested[test]['qualitative_activity']
            #will overwrite if multiple experiments,
            #then when parsing down will fetch ONLY the first match
            #To quick fix, don't allow overwrite:
            if path_label not in pathway_preference[tested[test]['ligand_id']].keys():
                #saving also the record ID so it will be easy to match the values
                pathway_preference[tested[test]['ligand_id']][path_label] = [Tau_KA, Emax_EC50, test]
            else:
                #Overwrite ONLY if both values are higher
                #because we want to keep the best outcome
                #from the different experiments performed
                try:
                    if (Tau_KA > pathway_preference[tested[test]['ligand_id']][path_label][0]) and (Emax_EC50 > pathway_preference[tested[test]['ligand_id']][path_label][1]):
                        pathway_preference[tested[test]['ligand_id']][path_label] = [Tau_KA, Emax_EC50, test]
                except TypeError:
                #Nedd to assess when comparing str to float
                    if isinstance(Tau_KA, float) and isinstance(Emax_EC50, float):
                        pathway_preference[tested[test]['ligand_id']][path_label] = [Tau_KA, Emax_EC50, test]
                    #only option is to overwrite when both new values are floats and we have the error

    #ranking accordingly to ∆Tau/KA or ∆Emax/EC50 (depending on how many missing values are)
    for val in pathway_preference.keys():
        temp = []
        if rank_method == 'tau':
            tmp_dict = dict(sorted(pathway_preference[val].items(), key=lambda item: -1000 if (item[1][0] == 'No activity' or item[1][0] == None or item[1][0] == 'Inverse agonism/antagonism') else item[1][0], reverse=True))
        elif rank_method == 'emax':
            tmp_dict = dict(sorted(pathway_preference[val].items(), key=lambda item: -1000 if (item[1][1] == 'No activity' or item[1][1] == None or item[1][1] == 'Inverse agonism/antagonism') else item[1][1], reverse=True))
        else: #for the build procedure
            none_tau = len([pathway_preference[val][key][0] for key in pathway_preference[val] if pathway_preference[val][key][0] is None])
            none_emax = len([pathway_preference[val][key][1] for key in pathway_preference[val] if pathway_preference[val][key][1] is None])
            if none_tau <= none_emax:
                tmp_dict = dict(sorted(pathway_preference[val].items(), key=lambda item: -1000 if (item[1][0] == 'No activity' or item[1][0] == None or item[1][0] == 'Inverse agonism/antagonism') else item[1][0], reverse=True))
            else:
                tmp_dict = dict(sorted(pathway_preference[val].items(), key=lambda item: -1000 if (item[1][1] == 'No activity' or item[1][1] == None or item[1][1] == 'Inverse agonism/antagonism') else item[1][1], reverse=True))
        #provide ranked keys
        for key in tmp_dict:
            temp.append(tmp_dict[key][2])
        pathway_preference[val] = temp

    return pathway_preference

def find_best_pathway_family(comparisons, tested):
    for ligand in list(comparisons.keys()):
        families = {}
        to_be_deleted = []
        for data in comparisons[ligand]:
        #assessing the values
            emax, ec50 = tested[data]['Emax'], tested[data]['EC50']
            try:
                logemaxec50 = math.log((emax/ec50),10)
            except (ValueError, TypeError, ZeroDivisionError):
                #setting a largely negative value insted of None for comparisons sake
                #excepting situations where one value is None or 0.0 [ValueError, TypeError, ZeroDivisionError]
                logemaxec50 = -100
            if tested[data]['primary_effector_family'] not in families.keys():
                #adding new family and values
                families[tested[data]['primary_effector_family']] = [logemaxec50, data]
            else:
                #updating with most relevant value
                if logemaxec50 > families[tested[data]['primary_effector_family']][0]:
                    #Adding the obsolete old data to the lists of tests to be deleted
                    to_be_deleted.append(families[tested[data]['primary_effector_family']][1])
                    #The data is a better duplicate and thus need to overwrite previous version
                    families[tested[data]['primary_effector_family']] = [logemaxec50, data]
                else:
                    #Adding the obsolete tests to the lists of tests to be deleted
                    to_be_deleted.append(data)
        #Delete obsolete key from comparisons register
        comparisons[ligand] = [item for item in comparisons[ligand] if item not in to_be_deleted]
        #get unique obsolete test keys and remove them
        for test in to_be_deleted:
            del tested[test]
    return comparisons, tested

def define_ligand_pathways(master):
    ligands = {}
    for key in master:
        if master[key]['ligand_id'] not in ligands.keys():
            ligands[master[key]['ligand_id']] = []
        ligands[master[key]['ligand_id']].append(key)
    return ligands

def OnTheFly(receptor_id, rank_method='Default', subtype=False, pathway=False, user=False, balanced=False):
    receptor_name = list(Protein.objects.filter(id=receptor_id).values_list("name", flat=True))[0]

    #fetching data given the receptor id
    if user == False:
        test_data = BiasedData.objects.filter(receptor=receptor_id) #37 for AGTR1
        pub_ids = list(BiasedData.objects.filter(receptor=receptor_id).values_list("publication", flat=True).distinct())
        lig_ids = list(BiasedData.objects.filter(receptor=receptor_id).values_list("ligand_id", flat=True).distinct())
    else:
        pub_ids = list(BiasedData.objects.filter(receptor=receptor_id, ligand=user).values_list("publication", flat=True).distinct())
        test_data = BiasedData.objects.filter(receptor=receptor_id, publication__in=pub_ids)
        lig_ids = list(BiasedData.objects.filter(receptor=receptor_id, publication__in=pub_ids).values_list("ligand_id", flat=True).distinct())

    # Performance: first collect all publication and ligand data
    pub_objs = Publication.objects.filter(id__in=pub_ids).values_list("id", "web_link_id__index", "year", "journal_id__name", "authors")
    pub_objs_dict = {pub_obj[0]:pub_obj[1:] for pub_obj in pub_objs}

    lig_objs = Ligand.objects.filter(id__in=lig_ids).values_list("id", "name")
    lig_objs_dict = {lig_obj[0]:lig_obj[1:] for lig_obj in lig_objs}

    publications = {}
    for entry in test_data:
        if entry.publication_id not in publications.keys():
            publications[entry.publication_id] = {}
        if entry.id not in publications[entry.publication_id].keys():
            if entry.publication_id in pub_objs_dict:
                pub_data = pub_objs_dict[entry.publication_id]
            else:
                pub_data = Publication.objects.filter(id=entry.publication_id).values_list("web_link_id__index",
                                                                                        "year",
                                                                                        "journal_id__name",
                                                                                        "authors")

            if entry.ligand_id in lig_objs_dict:
                ligand_data = lig_objs_dict[entry.ligand_id]
            else:
                ligand_data = Ligand.objects.filter(id=entry.ligand_id).values_list("name")

            publications[entry.publication_id][entry.id] = {'experiment': entry.experiment,
                                                            'endogenous_status': entry.endogenous_status,
                                                            'publication': entry.publication_id,
                                                            'doi': pub_data[0],
                                                            'pub_year': pub_data[1],
                                                            'journal': pub_data[2],
                                                            'authors': pub_data[3],
                                                            'receptor': receptor_name,
                                                            'receptor_id': entry.receptor_id,
                                                            'ligand_id': entry.ligand_id,
                                                            'ligand_name': ligand_data[0],
                                                            'receptor_isoform': entry.receptor_isoform,
                                                            'active_receptor_complex': entry.active_receptor_complex,
                                                            'cell_line': entry.cell_line,
                                                            'tissue': entry.tissue,
                                                            'species': entry.species,
                                                            'primary_effector_family': entry.primary_effector_family,
                                                            'primary_effector_subtype': entry.primary_effector_subtype,
                                                            'molecule_1': entry.molecule_1,
                                                            'molecule_2': entry.molecule_2,
                                                            'measured_process': entry.measured_process,
                                                            'pathway_level': entry.pathway_level,
                                                            'assay_type': entry.assay_type,
                                                            'EC50': entry.EC50,
                                                            'EC50_sign': entry.EC50_sign,
                                                            'qualitative_activity': entry.qualitative_activity,
                                                            'Emax': entry.Emax,
                                                            'Emax_sign': entry.Emax_sign,
                                                            'Tau_KA': entry.Tau_KA,
                                                            'delta_Tau_KA': entry.delta_Tau_KA,
                                                            'time_resolved': entry.time_resolved}


    for pub in list(publications.keys()):
        #Calculation branch 1 for Biased ligands (Regular, Balanced, Subtype)
        if pathway == False:
            #If we use Balanced referece, we need a set of preliminary steps
            if balanced == True:
                balanced_refs, skip = retrieve_balanced_references(receptor_id, pub, subtype)
                if skip == True:
                    del publications[pub]
                    continue
                #define ligand pathway and assess the pathway ranking for each ligand
                ligands = define_ligand_pathways(publications[pub])
                ranking = assess_pathway_preferences(ligands, publications[pub], rank_method, subtype)
                reference, tested = assess_pathway_pairs(ranking, publications[pub], balanced_refs, subtype)
            else:
                #Assess reference and tested ligands
                reference, tested, skip = assess_reference(publications[pub], user)
                if skip == True:
                    del publications[pub]
                    continue

            #Assess available comparisons
            comparisons = assess_comparisons(reference, tested)

            #calculate the first delta, remove excess data if needed
            comparisons, tested, skip = calculate_first_delta(comparisons, reference, tested, subtype)
            #all tested have calculated Δlog(Emax/EC50) and Δlog(Tau/KA)
            #this check might be added somewhere else
            if skip == True:
                del publications[pub]
                continue

            #calculate the second delta (across pathways)
            tested = calculate_second_delta(comparisons, tested, rank_method, subtype, pathway) #need subtype (defined by button)

            #save the updated data in the original key
            reference.update(tested)
            publications[pub] = reference

        #Calculation branch 2 for Pathway preference
        else:
            ligands = define_ligand_pathways(publications[pub])
            publications[pub] = calculate_second_delta(ligands, publications[pub], rank_method, subtype, pathway)

    # Cleaning time, removing empty publication records
    publications = {pub:value_id for pub,value_id in publications.items() if value_id != {}}
    return publications

def AddPathwayData(master, data, rank, pathway=False):
    master[rank+' - Pathway'] = data['primary_effector_family']
    try:
        master[rank+' - EC50'] = -math.log(data['EC50'],10)
    except (TypeError, ValueError):
        master[rank+' - EC50'] = data['EC50']
    master[rank+' - Emax'] = data['Emax']
    master[rank+' - Measured molecule 1'] = data['molecule_1']
    master[rank+' - Measured molecule 2'] = data['molecule_1']
    master[rank+' - Biological process'] = data['measured_process']
    master[rank+' - Cell line'] = data['cell_line']
    #This is should also check for the Pathway Preferred columns
    if set(['DeltaDelta_log(Emax/EC50)','DeltaDelta_log(Tau/KA)']).issubset(set(data.keys())):
        master['P1-'+rank+' - ΔΔLog(Emax/EC50)'] = data['DeltaDelta_log(Emax/EC50)']
        master['P1-'+rank+' - ΔΔLog(Tau/KA)'] = data['DeltaDelta_log(Tau/KA)']
        try:
            master[rank+' - Bias factor'] = round(data['Bias factor'], 1)
        except:
            master[rank+' - Bias factor'] = 'Full Bias'
    if pathway:
        master[rank+' - log(Tau/KA)'] = data['Tau_KA']
        try:
            master[rank+' - log(Emax/EC50)'] = data['log(Emax/EC50)']
        except KeyError:
            master[rank+' - log(Emax/EC50)'] = None
        if set(['Delta_log(Emax/EC50)','Delta_log(Tau/KA)']).issubset(set(data.keys())):
            master['P1-'+rank+' - ΔLog(Emax/EC50)'] = data['Delta_log(Emax/EC50)']
            master['P1-'+rank+' - ΔLog(Tau/KA)'] = data['Delta_log(Tau/KA)']
    else:
        master[rank+' - Subtype'] = data['primary_effector_subtype']
        master[rank+' - Δlog(Tau/KA)'] = data['delta_Tau_KA']
        try:
            master[rank+' - Δlog(Emax/EC50)'] = data['Delta_log(Emax/EC50)']
        except KeyError: #Delta_log(Emax/EC50) was not calculated
            master[rank+' - Δlog(Emax/EC50)'] = None
