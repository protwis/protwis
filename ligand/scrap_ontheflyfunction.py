##### Function for all the stuff starting point

# inputs required:
# - Receptor
# -	Endogenous ligand
# -	User defined ligand
# -	Effector family calculation
# -	Pathway-Preferring calculation
# -	Subtype calculation

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
            for test in ranking[drug][1:]:
                #Match pathway levels + skip matching if Arrestin involved
                if (tested[test]['pathway_level'] == tested[path1]['pathway_level']) or ('Arrestin' in [tested[path1]['primary_effector_family'], tested[test]['primary_effector_family']]):
                    tested[path1]['preferred_pathway'] = True
                    try:
                        deltadelta_logemaxec50 = round(tested[path1]['Delta_log(Emax/EC50)'] - tested[test]['Delta_log(Emax/EC50)'], 3)
                    except (TypeError, KeyError):
                        if tested[test]['qualitative_activity'] == 'No activity':
                            deltadelta_logemaxec50 = 'High Bias'
                        elif tested[path1]['Delta_log(Emax/EC50)'] == None:
                            deltadelta_logemaxec50 = None
                    tested[test]['DeltaDelta_log(Emax/EC50)'] = deltadelta_logemaxec50
                    try:
                        deltadelta_logtauka = round(tested[path1]['Delta_log(Tau/KA)'] - tested[test]['Delta_log(Tau/KA)'], 3)
                    except (TypeError, KeyError):
                        if tested[test]['qualitative_activity'] == 'No activity':
                            deltadelta_logtauka = 'High Bias'
                        elif tested[path1]['Delta_log(Tau/KA)'] == None:
                            deltadelta_logtauka = None
                    tested[test]['DeltaDelta_log(Tau/KA)'] = deltadelta_logtauka
    return tested

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
        publications[entry.publication_id][entry.id] = {'experiment': entry.experiment,
                                                        'endogenous_status': entry.endogenous_status,
                                                        'publication': entry.publication_id,
                                                        'receptor': entry.receptor.name,
                                                        'receptor_id': entry.receptor_id,
                                                        'ligand_id': entry.ligand_id,
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
for pub in publications:
    print("Performing analysis on publication id: " + str(pub))
    #Assess reference and tested ligands
    reference, tested = assess_reference(publications[pub])
    #Assess available comparisons
    comparisons = assess_comparisons(reference, tested)
    #calculate the first delta, remove excess data if needed
    print("Calculating the first Delta")
    comparisons, tested, skip = calculate_first_delta(comparisons, reference, tested)
    #this check might be added somewhere else
    if skip == True:
        continue
    #calculate the second delta (across pathways)
    print("Test that can be assessed for second delta: " + str(len(tested)))
    tested = calculate_second_delta(comparisons, tested) #need more?
