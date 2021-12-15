##### Function for all the stuff starting point

# inputs required:
# - Receptor
# -	Endogenous ligand
# -	User defined ligand
# -	Effector family calculation
# -	Pathway-Preferring calculation
# -	Subtype calculation

from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import connection
from django.db import IntegrityError
from django.utils.text import slugify
from django.http import HttpResponse, JsonResponse
from decimal import Decimal
from build.management.commands.base_build import Command as BaseBuild
from common.tools import fetch_from_cache, save_to_cache, fetch_from_web_api
from protein.models import Protein, ProteinCouplings
from ligand.models import  Ligand, LigandProperities, LigandType, Endogenous_GTP, BiasedData, LigandVendorLink

from ligand.functions import get_or_make_ligand
from common.models import WebLink, WebResource, Publication
import logging
import math
import pandas as pd
import os
import traceback
import time
import requests
import timeit

def assess_reference(data_dict, user=False):
    gtp_info = {}
    reference = {}
    tested = {}
    #if reference ligand is not provided by user
    if user == False:
        #check for highest ranking
        for assay in data_dict.keys():
            #grab the endogenous data for the ligand if present
            data = Endogenous_GTP.objects.filter(receptor=data_dict[assay]['receptor_id'],
                                                ligand=data_dict[assay]['ligand_id'],
                                                # ligand_specie__latin_name=data_dict[assay]['specie'],
                                                ).values_list("potency_ranking").distinct()
            try:
                gtp_info[data_dict[assay]['ligand_id']] = sum(data, ())[0]
            except IndexError:
                pass
        reference_ligand = min(gtp_info, key=gtp_info.get)
    else:
        reference_ligand = user
    for assay in data_dict.keys():
        if data_dict[assay]['ligand_id'] == reference_ligand:
            reference[assay] =  data_dict[assay]
        else:
            tested[assay] =  data_dict[assay]
    return reference, tested

def assess_comparisons(reference, tested):
    common_traits = ['tissue',
                     'specie',
                     'primary_effector_family',
                     # 'receptor_isoform', #this can be null
                     # 'active_receptor_complex', #this can be null
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
            r_logemaxec50 = math.log(r_emax/r_ec50)
        except TypeError:
            r_logemaxec50 = None
        #Check for lacking of ALL values, then break cycle and skip publication
        if (r_logemaxec50 == None) and (r_logtauka == None):
            skip = True
            return comparisons, tested, skip

        for test in comparisons[ref]:
            t_emax, t_ec50, t_logtauka = tested[test]['Emax'], tested[test]['EC50'], tested[test]['Tau_KA']
            try:
                t_logemaxec50 = math.log(t_emax/t_ec50)
            except TypeError:
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
            tested[test]['Reference_ligand'] = reference[ref]['ligand_name']
            tested[test]['Reference_Emax'] = r_emax
            tested[test]['Reference_EC50'] = r_ec50
            tested[test]['Reference_Tau/KA'] = r_logtauka
    return comparisons, tested, skip

def find_best_subtype(comparisons, reference, tested):
    families = {}
    for ref in list(comparisons.keys()):
        #assessing the values
        r_emax, r_ec50, r_logtauka = reference[ref]['Emax'], reference[ref]['EC50'], reference[ref]['Tau_KA']
        try:
            r_logemaxec50 = math.log(r_emax/r_ec50)
        except TypeError:
            r_logemaxec50 = None

        if reference[ref]['primary_effector_family'] not in families.keys():
            #adding new family and values
            families[reference[ref]['primary_effector_family']] = [r_logemaxec50, ref]
        else:
            #updating with most relevant subtype
            if r_logemaxec50 > families[reference[ref]['primary_effector_family']][0]:
                families[reference[ref]['primary_effector_family']] = [r_logemaxec50, ref]
                for test in comparisons[families[reference[ref]['primary_effector_family']][1]]:
                    del tested[test]
                del comparisons[families[reference[ref]['primary_effector_family']][1]]
            else:
                for test in comparisons[ref]:
                    del tested[test]
                del comparisons[ref]

def calculate_second_delta(comparisons, tested):
    ranking = assess_pathway_preferences(comparisons, tested)
    #STEPS
    for drug in ranking.keys():
        #perform analysis only when we have multiple pathways to actually compare
        if len(ranking[drug]) > 1:
        #Set reference pathway (first in list)
            path1 = ranking[drug][0]
            tested[path1]['Pathway Rank'] = 'P1'
            path_count = 1
            for test in ranking[drug][1:]:
                path_count +=1
                tested[test]['Pathway Rank'] = 'P'+str(path_count)
                #Match pathway levels + skip matching if Arrestin involved
                if (tested[test]['pathway_level'] == tested[path1]['pathway_level']) or ('Arrestin' in [tested[path1]['primary_effector_family'], tested[test]['primary_effector_family']]):
                    tested[test]['P1'] = tested[path1]['primary_effector_family']
                    try:
                        deltadelta_logemaxec50 = round(tested[path1]['Delta_log(Emax/EC50)'] - tested[test]['Delta_log(Emax/EC50)'], 3)
                    except (TypeError, KeyError):
                        if tested[test]['qualitative_activity'] == 'No activity':
                            deltadelta_logemaxec50 = 'Full Bias'
                        elif tested[path1]['Delta_log(Emax/EC50)'] == None:
                            deltadelta_logemaxec50 = None
                    tested[test]['DeltaDelta_log(Emax/EC50)'] = deltadelta_logemaxec50
                    try:
                        deltadelta_logtauka = round(tested[path1]['Delta_log(Tau/KA)'] - tested[test]['Delta_log(Tau/KA)'], 3)
                    except (TypeError, KeyError):
                        if tested[test]['qualitative_activity'] == 'No activity':
                            deltadelta_logtauka = 'Full Bias'
                        elif tested[path1]['Delta_log(Tau/KA)'] == None:
                            deltadelta_logtauka = None
                    tested[test]['DeltaDelta_log(Tau/KA)'] = deltadelta_logtauka
    return tested, ranking

def assess_pathway_preferences(comparisons, tested):
    pathway_preference = {}
    #calculate values for ranking (or replace with qualitative activity when missing)
    for assay in comparisons.keys():
        for test in comparisons[assay]:
            if tested[test]['ligand_id'] not in pathway_preference.keys():
                pathway_preference[tested[test]['ligand_id']] = {}
            try:
                pathway_preference[tested[test]['ligand_id']][tested[test]['primary_effector_family']] = [tested[test]['Delta_log(Tau/KA)'], tested[test]['Delta_log(Emax/EC50)']]
            except KeyError:
                try:
                    pathway_preference[tested[test]['ligand_id']][tested[test]['primary_effector_family']] = [tested[test]['Delta_log(Tau/KA)'], tested[test]['qualitative_activity']]
                except KeyError:
                    try:
                        pathway_preference[tested[test]['ligand_id']][tested[test]['primary_effector_family']] = [tested[test]['qualitative_activity'], tested[test]['Delta_log(Emax/EC50)']]
                    except KeyError:
                        pathway_preference[tested[test]['ligand_id']][tested[test]['primary_effector_family']] = [tested[test]['qualitative_activity'], tested[test]['qualitative_activity']]
    #ranking accordingly to ∆Tau/KA or ∆Emax/EC50 (depending on how many missing values are)
    for id in pathway_preference.keys():
        temp = []
        none_tau = len([pathway_preference[id][key][0] for key in pathway_preference[id] if pathway_preference[id][key][0] is None])
        none_emax = len([pathway_preference[id][key][1] for key in pathway_preference[id] if pathway_preference[id][key][1] is None])
        if none_tau <= none_emax:
            pathway_preference[id] = list(dict(sorted(pathway_preference[id].items(), key=lambda item: -1000 if (item[1][0] == 'No activity' or item[1][1] == None) else item[1][0], reverse=True)).keys())
        else:
            pathway_preference[id] = list(dict(sorted(pathway_preference[id].items(), key=lambda item: -1000 if (item[1][1] == 'No activity' or item[1][1] == None) else item[1][1], reverse=True)).keys())
    #provide ranked keys
        for path in pathway_preference[id]:
            temp.append([key for key in tested.keys() if (tested[key]['ligand_id'] == id and tested[key]['primary_effector_family'] == path)][0])
        pathway_preference[id] = temp
    return pathway_preference

#list of all the data from receptor 114 (CCR1)
test_data = BiasedData.objects.filter(receptor=114)
#list of all the data from receptor 37 (AGTR1) (subtype test)
test_data = BiasedData.objects.filter(receptor=37)

#gather all the data into a single dictionary
#split entries by publication for further analysis and into ids per publication
#add a filter option for the user defined ligand as reference
publications = {}
for entry in test_data:
    if entry.publication_id not in publications.keys():
        publications[entry.publication_id] = {}
    if entry.id not in publications[entry.publication_id].keys():
        pub_data = Publication.objects.filter(id=entry.publication_id).values_list("web_link_id__index",
                                                                                   "year",
                                                                                   "journal_id__name",
                                                                                   "authors")
        ligand_data = Ligand.objects.filter(id=entry.ligand_id).values_list("name",
                                                                            "properities_id")
        publications[entry.publication_id][entry.id] = {'experiment': entry.experiment,
                                                        'endogenous_status': entry.endogenous_status,
                                                        'publication': entry.publication_id,
                                                        'doi': pub_data[0][0],
                                                        'pub_year': pub_data[0][1],
                                                        'journal': pub_data[0][2],
                                                        'authors': pub_data[0][3],
                                                        'receptor': entry.receptor.name,
                                                        'receptor_id': entry.receptor_id,
                                                        'ligand_id': entry.ligand_id,
                                                        'ligand_name': ligand_data[0][0],
                                                        'ligand__properities_id': ligand_data[0][1], #for browser vendors
                                                        'receptor_isoform': entry.receptor_isoform,
                                                        'active_receptor_complex': entry.active_receptor_complex,
                                                        'cell_line': entry.cell_line,
                                                        'tissue': entry.tissue,
                                                        'specie': entry.specie,
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
                                                        'delta_Tau_KA': entry.delta_Tau_KA}

for pub in publications.keys():
    for assay in publications[pub].keys():
        if publications[pub][assay]['ligand_id'] in [278403, 284573]:
            publications[pub][assay]['ligand_id'] = 273200

#find the endogenous ligand for each publication based on endogenous_status (make a call)
#TESTING pub 1135
for pub in publications:
    print("Performing analysis on publication id: " + str(pub))
    #Assess reference and tested ligands
    reference, tested = assess_reference(publications[pub])
    #9 reference dicts, 27 tested dicts
    #Assess available comparisons
    comparisons = assess_comparisons(reference, tested)
    #OUTPUT comparisons:
    # Reference id --> tested ids
    # {90178: [90151, 90160, 90169],
    #  90179: [90152, 90161, 90170], #different subtype
    #  90180: [90153, 90162, 90171],
    #  90181: [90154, 90163, 90172], #different subtype
    #  90182: [90155, 90164, 90173],
    #  90183: [90156, 90165, 90174],
    #  90184: [90157, 90166, 90175], #different subtype
    #  90185: [90158, 90167, 90176], #different subtype
    #  90186: [90159, 90168, 90177]}
    #calculate the first delta, remove excess data if needed
    print("Calculating the first Delta")
    comparisons, tested, skip = calculate_first_delta(comparisons, reference, tested)
    #OUTPUT comparisons (shaved):
    # {90178: [90151, 90160, 90169],
    #  90180: [90153, 90162, 90171],
    #  90182: [90155, 90164, 90173],
    #  90183: [90156, 90165, 90174],
    #  90186: [90159, 90168, 90177]}
    #9 reference dicts, 15 tested dicts
    #all tested have calculated Δlog(Emax/EC50) and Δlog(Tau/KA)
    #this check might be added somewhere else
    if skip == True:
        continue
    #calculate the second delta (across pathways)
    print("Test that can be assessed for second delta: " + str(len(tested)))
    tested, ranking = calculate_second_delta(comparisons, tested) #need more?
    reference.update(tested)
    publications[pub] = reference)
    #OUTPUT Ranking:
    #ligand id --> ranking by pathway (Δlog(Emax/EC50) and Δlog(Tau/KA))
    # {1824: [90153, 90155, 90151, 90156, 90159],
    #  1815: [90162, 90165, 90168, 90160, 90164],
    #  1806: [90173, 90169, 90171, 90174, 90177]}


###########################
# Quindi per ogni ligando comparato devo avere
# la lista di dati qui sotto descritta
#
#
#RANK ORDER (Scatter Plot)
'Gi/o',        #Highest biased pathway for ligand
'Rotigotine',  #ligand name                                                     #CHECK
'10.1124/jpet.114.220293', #publication DOI                                     #CHECK
2014,          #publication year                                                #CHECK
'Journal of Pharmacology and Experimental Therapeutics', #publication journal   #CHECK
'Brust TF, Hayes MP, Roman DL, Burris KD, Watts VJ',     #publication authors   #CHECK
1635,          # ligand ID                                                      #CHECK
'dopamine',    # reference ligand name                                          #CHECK
 None,         # qualitative activity --> this should exclude the exclude list  #CHECK
 1.8,          # ΔLog(Emax/EC50)                                                #CHECK
 None,         # ΔΔLog(Emax/EC50)                                               #CHECK
 0,            # ranking                                                        #what is this?
 2.1,          # ΔLog(Tau/KA)                                                   #CHECK
 None,         # ΔΔLog(Tau/KA)                                                  #CHECK
 0.0,          # EC50                                                           #CHECK
 95.0,         # Emax                                                           #CHECK
 1.2e-09,      # EC50 of reference                                              #CHECK
 100.0,        # Emax of reference                                              #CHECK
 'Gi/o',       # Pathway of Reference                                           #CHECK
 'Dopamine',   # Name reference ligand                                          #CHECK
 None,         # subtype                                                        #CHECK
 None          # compared subtype                                               #CHECK
###########################
exclude_list = ["Agonism", "Partial agonism", "Medium activity"] #full agonism should be removed
publications = list(AnalyzedAssay.objects.filter(
                family__isnull=False,
                experiment__receptor=receptor,
                assay_description=self.assay,
                experiment__source=self.source).exclude(
                qualitative_activity__in=exclude_list,
                ).values_list(
                "family", #pathway                                              0
                "experiment__ligand__name", #name                               1
                "experiment__publication__web_link__index", # DOI               2
                "experiment__publication__year", #year                          3
                "experiment__publication__journal__name", #journal              4
                "experiment__publication__authors",  #authors                   5
                "experiment__ligand",    #ligand_id for hash                    6
                "experiment__endogenous_ligand__name", #endogenous              7
                "qualitative_activity",  #activity values                       8
                "delta_emax_ec50", #ΔLog(Emax/EC50)                           9
                "log_bias_factor",  #ΔΔLog(Emax/EC50)                           10
                "order_no",         #ranking                                    11
                "relative_transduction_coef",    #ΔLog(TAU/Ka)                               12
                "delta_relative_transduction_coef",         #ΔΔLog(TAU/Ka)                              13
                "quantitive_activity",   #EC 50                                 14
                "quantitive_efficacy",   #Emax                                  15
                "reference_assay_initial_id__quantitive_activity", #            16 EC 50 compared drug
                "reference_assay_initial_id__quantitive_efficacy", #            17 Emax compared drug
                "reference_assay_initial_id__family", # Pathway compared drug   18
                "reference_assay_initial_id__biased_experiment__ligand__name", #19 Name compared drug
                "signalling_protein",   #subtype                                20
                "reference_assay_initial_id__signalling_protein", #             21 Compared signalling protein
                ).distinct()) #check
################################################################################
#Post Parsing data for Browser section
#Setting the function for adding based on pathway rank
def add_pathway_data(master, data, rank):
    master[rank+' - Pathway'] = data['primary_effector_family']
    try:
        master[rank+' - Delta_log(Emax/EC50)'] = data['Delta_log(Emax/EC50)']
    except KeyError: #Delta_log(Emax/EC50) was not calculated
        master[rank+' - Delta_log(Emax/EC50)'] = None
    master[rank+' - EC50'] = data['EC50']
    master[rank+' - Emax'] = data['Emax']
    master[rank+' - Delta_log(Tau/KA)'] = data['delta_Tau_KA']
    master[rank+' - Measured molecule 1'] = data['molecule_1']
    master[rank+' - Measured molecule 2'] = data['molecule_1']
    master[rank+' - Biological process'] = data['measured_process']
    if set(['DeltaDelta_log(Emax/EC50)','DeltaDelta_log(Tau/KA)']).issubset(set(data.keys())):
        master['P1-'+rank+' - ΔΔLog(Emax/EC50)'] = data['DeltaDelta_log(Emax/EC50)']
        master['P1-'+rank+' - ΔΔLog(Tau/KA)'] = data['DeltaDelta_log(Tau/KA)']
#Setting the final dict
table = pd.DataFrame()
#receptor_id
receptor_info = Protein.objects.filter(id=37).values_list("family__parent__parent__parent__name",
                                                          "family__parent__name",
                                                          "entry_name",
                                                          "name")
#Parsing through the pubs
#in a one big open dict
#over the data we go
#swearing all the way
for pub in publications:
    ligands = {}
    slice = pd.DataFrame()
    for key in publications[pub]:
        if 'Pathway Rank' in publications[pub][key].keys():
            if publications[pub][key]['ligand_id'] not in ligands.keys():
                labs = []
                ligands[publications[pub][key]['ligand_id']] = {}
                ligands[publications[pub][key]['ligand_id']]['Name'] = publications[pub][key]['ligand_name']
                ligands[publications[pub][key]['ligand_id']]['Species'] = publications[pub][key]['specie']
                ligands[publications[pub][key]['ligand_id']]['Cell line'] = publications[pub][key]['cell_line']
                ligands[publications[pub][key]['ligand_id']]['Authors'] = publications[pub][key]['authors']
                ligands[publications[pub][key]['ligand_id']]['DOI/Reference'] = publications[pub][key]['doi']
                ligands[publications[pub][key]['ligand_id']]['#Vendors'] = LigandVendorLink.objects.filter(lp_id=publications[pub][key]['ligand__properities_id']).values_list("vendor_id").distinct().count()
                ligands[publications[pub][key]['ligand_id']]['#Articles'] = BiasedData.objects.filter(ligand_id=publications[pub][key]['ligand_id']).values_list("publication_id").distinct().count()
                for authors in BiasedData.objects.filter(ligand_id=publications[pub][key]['ligand_id']).values_list("publication_id__authors").distinct():
                    if authors[0].split(',')[-1] not in labs:
                        labs.append(authors[0].split(',')[-1])
                ligands[publications[pub][key]['ligand_id']]['#Labs'] = len(labs)
            add_pathway_data(ligands[publications[pub][key]['ligand_id']], publications[pub][key], publications[pub][key]['Pathway Rank'])
        else:
            Reference_ligand = publications[pub][key]['ligand_name']
    for drug in ligands:
        slice = slice.append(ligands[drug], ignore_index=True)
    slice['Class'] = receptor_info[0][0]
    slice['Receptor family'] = receptor_info[0][1].strip('receptors')
    slice['Uniprot'] = receptor_info[0][2].split('_')[0].upper()
    slice['IUPHAR'] = receptor_info[0][3]
    slice['Reference ligand'] = Reference_ligand

    table = table.append(slice, ignore_index=True)
