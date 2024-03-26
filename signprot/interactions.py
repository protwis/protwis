import json
import re
import string
import random

from collections import Counter
from decimal import Decimal

from residue.models import ResidueGenericNumberEquivalent
from signprot.models import SignprotComplex
from protein.models import Protein, ProteinCouplings
from common.definitions import AMINO_ACID_GROUPS

from django.core.exceptions import ObjectDoesNotExist


def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))

def get_entry_names(request):
    """Extract a list of entry names from the post request"""
    prot_confs = request.POST.getlist("pos[]")
    complex_objs = SignprotComplex.objects.prefetch_related('structure__protein_conformation__protein').filter(structure__protein_conformation__in=prot_confs)
    entry_names = [complex_obj.structure.protein_conformation.protein.entry_name for complex_obj in complex_objs]

    return entry_names


def get_ignore_info(request):
    """Extract a dict of residues to ignore from the post request"""
    ignore_dict = request.POST.get("ignore")
    return json.loads(ignore_dict)

def get_class_slug(common_class):
    # access the list of the most common element and get the value
    value = common_class[0][0]
    # extract the class character, e.g. Class B1 Receptor -> B
    value = value.split(' ')[1][0]
    # return the lowercase character
    return value.lower()


def get_protein_segments(request):
    """From a list of given generic numbers (3x50), query appropriate generic residue
    number objects"""
    segments = []
    segment_raw = request.POST.getlist("seg[]")
    selected_receptor_classes = request.POST.getlist("selectedreceptorclasses[]")

    for s in segment_raw:
        try:
            gen_object = ResidueGenericNumberEquivalent.objects.filter(
                label=s, scheme__slug__in=['gpcrdba']
            ).get()
            segments.append(gen_object)
        except ObjectDoesNotExist as e:
            print("For {} a {} ".format(s, e))
            continue
    return segments


def get_generic_numbers(signature_data):
    """Parse the generic numbers in the signature data"""
    generic_numbers = []
    for _, segments in signature_data["common_generic_numbers"].items():
        for elem, num in segments.items():
            gnl = []
            for x, dn in num.items():
                gnl.append(x)
            generic_numbers.append(gnl)

    return generic_numbers


def get_signature_features(signature_data, generic_numbers, feats):
    """Extract the signature features and prepare for visualization"""
    signature_features = []
    x = 0

    tmp = list()
    for segment, s in signature_data["a_pos"].consensus.items():
        for p, r in s.items():
            tmp.append({
                "aa": r[0],
                "aa_cons": r[2]
                })

    for i, feature in enumerate(signature_data["a_pos"].feature_stats):
        for j, segment in enumerate(feature):
            for k, freq in enumerate(segment):
                # freq0: score
                # freq1: level of conservation
                # freq2: a - b explanation
                try:
                    if int(freq[0]) > 0:
                        dkey = int(x)
                        dfeature = str(feats[i][0])
                        dfeature_code = str(feats[i][1])
                        dlength = str(feats[i][2])
                        dgn = str(generic_numbers[j][k])
                        dfreq = int(freq[0])
                        dcons = int(freq[1])

                        sort_code = dfeature_code + "_" + dlength
                        if sort_code in AMINO_ACID_GROUPS:
                            sort_score = len(AMINO_ACID_GROUPS[sort_code])
                        elif dfeature_code == 'Y':
                            print('Y_?')
                            sort_code = dfeature_code + "_" + '?'
                            sort_score = len(AMINO_ACID_GROUPS[sort_code])
                        else:
                            sort_score = 99

                        signature_features.append(
                            {
                                "key": dkey,
                                "feature": dfeature,
                                "feature_code": dfeature_code,
                                "length": dlength,
                                "gn": dgn,
                                "freq": dfreq,
                                "cons": dcons,
                                "sort_score": sort_score,
                                # 'expl': str(freq[2]),
                                "aa": str(tmp[k]["aa"]),
                                "aa_cons": int(tmp[k]["aa_cons"]),
                                "sort_code": str(sort_code),
                            }
                        )
                    x += 1
                except Exception as e:
                    print(e)
                    continue
    return signature_features


def group_signature_features(signature_features):
    """Further prepare signature feature dict for visualization"""
    grouped_features = {}
    for feature in signature_features:
        if feature["gn"] not in grouped_features:
            grouped_features[feature["gn"]] = []
        grouped_features[feature["gn"]].append(feature)

    for key in grouped_features:
        curr_group = grouped_features[key]
        grouped_features[key] = sorted(
            curr_group, key=lambda feature: feature["freq"], reverse=True
        )
    return grouped_features


def get_signature_consensus(signature_data, generic_numbers):
    """Extract the signature consensus and prepare for visualization"""
    sigcons = []
    x = 0
    for segment, cons in signature_data["feats_cons_pos"].items():
        for pos in cons:
            # pos0: Code
            # pos1: Name
            # pos2: Score
            # pos3: level of conservation

            # res0: Property Abbreviation
            # res1: Feature Score
            # res2: Conservation Level
            try:
                sigcons.append(
                    {
                        "key": int(x),
                        "gn": str(generic_numbers[x]),
                        "code": str(pos[0]),
                        "feature": str(pos[1]),
                        "score": int(pos[2]),
                        "cons": int(pos[3]),
                        "length": str(pos[4]),
                    }
                )
                x += 1
            except Exception as e:
                print(e)
    return sigcons


def prepare_signature_match(signature_match, effector):
    repl_str = id_generator(6)
    sign_true_1 = '<div class="{}">'.format(repl_str)
    sign_true_2 = '{}</div>'
    sign_false = '<div></div>'
    class_coupling = 'coupling '
    coupling_data = prepare_coupling_data_container()
    coupling_data = process_coupling_data(coupling_data, effector)

    coupling_data_dict = {}
    for entry in coupling_data:
        coupling_data_dict[entry['rec_obj'].entry_name] = entry
    out = {}
    for elem in signature_match["scores"].items():
        entry = elem[0].protein.entry_name
        out[entry] = {
            "entry": elem[0].protein.entry_short(),
            "prot": elem[0].protein.name,
            "score": elem[1][0],
            "nscore": round(elem[1][1], 0),
            "class": elem[0].protein.get_protein_class().strip().split(' ')[1],
            "family": elem[0].protein.get_protein_family(),
            "subfamily": elem[0].protein.get_protein_subfamily(),
        }
        coupling_entry = coupling_data_dict.get(entry)
        for prot in coupling_entry['coupling'].keys():
            out[entry][prot] = {}
            if coupling_entry['coupling'][prot] == 'primary':
                html_val = sign_true_1.replace(repl_str, class_coupling+'prim') + sign_true_2.format("1'")
                text_val = "1'"
            elif coupling_entry['coupling'][prot] == 'secondary':
                html_val = sign_true_1.replace(repl_str, class_coupling+'seco') + sign_true_2.format("2'")
                text_val = "2'"
            elif coupling_entry['coupling'][prot] == 'No Data':
                html_val = sign_true_1.replace(repl_str, class_coupling+'nodata') + sign_true_2.format('No Data')
                text_val = 'No Data'
            elif coupling_entry['coupling'][prot] == '':
                html_val = sign_false
                text_val = ''
            else:
                html_val = sign_true_1.replace(repl_str, class_coupling+'value') + sign_true_2.format(str(coupling_entry['coupling'][prot]))
                text_val = str(coupling_entry['coupling'][prot])
            out[entry][prot]['html'] = html_val
            out[entry][prot]['text'] = text_val
    return out


def prepare_coupling_data_container():
    class_names = {}
    data = {}
    prot_names = []
    # Calling the db for all the GPCR proteins
    proteins = (
        Protein.objects.filter(
            sequence_type__slug="wt",
            family__slug__startswith="00",
            species__common_name="Human")
        .prefetch_related("family",
                          "family__parent__parent",
                          "family__parent__parent__parent")
        )
    # Calling the db for all the coupling data available
    couplings = ProteinCouplings.objects.prefetch_related(
        "protein", "g_protein_subunit", "g_protein"
    )
    # Setting up the data dictionary with the human GPCR proteins
    for protein in proteins:
        if protein.entry_short() not in prot_names:
            prot_names.append(protein.entry_short())
            protein_class = protein.family.slug.split("_")[0]
            if protein_class not in class_names:
                class_names[protein_class] = re.sub(
                    r"\([^)]*\)", "", protein.family.parent.parent.parent.name
                )
            protein_class_name = class_names[protein_class].strip()

            data[protein.entry_short()] = {
                "class": protein_class_name,
                "pretty": protein.short()[:15],
                "rec_class": protein.family.parent.parent.name,
                "rec_obj": protein,
                "rec_pdb": protein.entry_short(),
                "sources": {},
            }
    # Parsing through the coupling data and filling up the data dictionary
    for record in couplings:
        protein_name = record.protein.entry_short()
        # Skip entries without any annotation
        if protein_name not in prot_names:
            continue

        source = record.source
        transduction = record.transduction
        emax_val = record.logemaxec50
        g_protein_family = record.g_protein.name
        g_protein_family = g_protein_family.replace(" family", "")
        # Adding the information about subunit
        if record.g_protein_subunit:
            g_protein_subunit = record.g_protein_subunit.entry_name
            g_protein_subunit = g_protein_subunit.replace("_human", "")
        # Adding all the sources of coupling data
        if source not in data[protein_name]["sources"]:
            data[protein_name]["sources"][source] = {}
        # Adding the information about g protein coupling family
        if g_protein_family not in data[protein_name]["sources"][source]:
            data[protein_name]["sources"][source][g_protein_family] = {}
        # Adding transduction data (based on source information)
        if transduction:
            data[protein_name]["sources"][source][g_protein_family] = transduction
        else:
            if "subunits" not in data[protein_name]["sources"][source][g_protein_family]:
                data[protein_name]["sources"][source][g_protein_family] = {"subunits": {}, "best": -2.00}
            data[protein_name]["sources"][source][g_protein_family]["subunits"][g_protein_subunit] = round(Decimal(emax_val), 2)
            if round(Decimal(emax_val), 2) == -0.00:
                data[protein_name]["sources"][source][g_protein_family]["subunits"][g_protein_subunit] = 0.00
            # get the higher number into 'best'
            if emax_val > data[protein_name]["sources"][source][g_protein_family]["best"]:
                data[protein_name]["sources"][source][g_protein_family]["best"] = round(Decimal(emax_val), 2)
    # Purge the record without any coupling data
    # purged_data = {key: value for key, value in data.items() if len(data[key]['sources'])>0 }
    # Highlight data with no sources somehow (like add a No Data field further down)
    # return purged_data #or data if we want to keep blank data
    return data

def process_coupling_data(data, effector):
    results = []
    # For the g proteins use the consensus options from the coupling browser
    # (either GtP or at least two sources)
    # For the arrestins use the Bouvier only data)
    for entry in data.keys():
        record = data[entry]
        temp_dict = {}
        temp_dict["rec_class"] = record["rec_class"]
        temp_dict["rec_obj"] = record["rec_obj"]
        temp_dict["pretty"] = record["pretty"]
        temp_dict["rec_pdb"] = record["rec_pdb"]
        temp_dict['coupling'] = extract_coupling(record, effector)
        results.append(temp_dict)
    return results

def extract_coupling(entry, effector):
# For the g proteins use the consensus options from the coupling browser
# (either GtP or at least two sources)
# For the arrestins use the Bouvier only data)
    sources = list(entry['sources'].keys())
    couplings = {}
    # not needed and to be better worded
    if effector == 'G alpha':
        capital = 'G'
        refined_coupling = {'Gs': "",
                            'Gi/o': "",
                            'Gq/11': "",
                            'G12/13': "",
                            'Gs_emax': "",
                            'Gi/o_emax': "",
                            'Gq/11_emax': "",
                            'G12/13_emax': ""}
        for source in sources:
            # The actual primary and secondary, based on GtoP
            if source == 'GuideToPharma':
                for subunit in entry['sources'][source].keys():
                    refined_coupling[subunit] = entry['sources'][source][subunit]
            else:
                # or if GtoP is not available, based on two sources
                couplings[source] = {}
                transducers = list(entry['sources'][source].keys())
                # Keeping only G Protein values given the effector selection
                transducers = [item for item in transducers if item.startswith(capital)]
                for transducer in transducers:
                    couplings[source][transducer] = float(entry['sources'][source][transducer]['best'])

    elif effector == 'A':
        capital = 'B'
        refined_coupling = {'arrb1': "",
                            'arrb2': "",
                            'arrb1_emax': "",
                            'arrb2_emax': ""}
        for source in sources:
            if source == 'Bouvier':
                transducers = list(entry['sources'][source].keys())
                # Keeping only Arrestin values given the effector selection
                transducers = [item for item in transducers if item.startswith(capital)]
                # Then parsing if we have at least one Arrestin (which should be labeled as 'Beta')
                if len(transducers) == 1:
                    couplings[source] = {}
                    transducer = transducers[0]
                    for subunit in entry['sources'][source][transducer]['subunits']:
                        couplings[source][subunit] = float(entry['sources'][source][transducer]['subunits'][subunit])
        # Here we need at least one source that provides info (Bouvier Data)
        if len(couplings.keys()) == 1:
            source = list(couplings.keys())[0]
            filtered_dict = {arrestin:value for arrestin,value in couplings[source].items() if value != 0.0}
            if len(filtered_dict) > 0:
                sorted_list = list({arrestin: value for arrestin, value in sorted(filtered_dict.items(), key=lambda item: item[1], reverse=True)})
                refined_coupling[sorted_list[0]] = filtered_dict[sorted_list[0]]
                try:
                    refined_coupling[sorted_list[1]] = filtered_dict[sorted_list[1]]
                except IndexError:
                    pass
    # If more sources are provided, we define the mean values as in the browser
    if len(couplings.keys()) > 1:
        means = {}
        for source in couplings.keys():
            for transducer in couplings[source]:
              if transducer not in means:
                  means[transducer] = [couplings[source][transducer]]
              else:
                  means[transducer].append(couplings[source][transducer])
                  # Filtering out the 0.0 values
                  means[transducer] = list(filter((0.0).__ne__, means[transducer]))
        for transducer in means.keys():
            if len(means[transducer]) > 1:
                refined_coupling[transducer+'_emax'] = round(sum(means[transducer])/len(means[transducer]),2)
    # If all the data is empty because we do not have data
    # Add the label "No Data" that it will be displayed in the table
    refined_values = list(set(list(refined_coupling.values())))
    if len(refined_values) == 1 and refined_values[0] == '':
        for item in refined_coupling.keys():
            refined_coupling[item] = "No Data"

    return refined_coupling
