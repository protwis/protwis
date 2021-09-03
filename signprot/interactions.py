import json
import re
import time
from itertools import chain
import string
import random

from collections import Counter

from residue.models import ResidueGenericNumberEquivalent
from signprot.models import SignprotComplex
from protein.models import Protein, ProteinCouplings
from common.definitions import *

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
    most_common_class = Counter(selected_receptor_classes).most_common(1)
    slug_ending = get_class_slug(most_common_class)

    for s in segment_raw:
        try:
            gen_object = ResidueGenericNumberEquivalent.objects.filter(
                # label=s, scheme__slug__in=['gpcrdb' + slug_ending]
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


def prepare_signature_match(signature_match):
    repl_str = id_generator(6)
    sign_true_1 = '<div class="{}">'.format(repl_str)
    sign_true_2 = '{}</div>'
    sign_false = '<div></div>'
    gprots = ['Gs','Gi/o','Gq/11','G12/13']
    class_coupling = 'coupling '

    coupling_data = prepare_coupling_data_container()
    coupling_data = fill_coupling_data_container(coupling_data)
    coupling_data = process_coupling_data(coupling_data)

    coupling_data_dict = {}
    for e in coupling_data:
        coupling_data_dict[e['rec_obj'].entry_name] = e

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

    for elem in signature_match["scores"].items():
        entry = elem[0].protein.entry_name
        coupling_entry = coupling_data_dict.get(entry)

        sources = ["GuideToPharma", "Aska", "Merged"]

        for source in sources:
            out[entry][source] = {}
            for gprot in gprots:
                out[entry][source][gprot] = {}
                if coupling_entry:
                    ce = coupling_entry
                    cl = ce['coupling'][source].get(gprot, '')
                    if ce[source][gprot]:
                        if cl[:4] == 'prim':
                            html_val = sign_true_1.replace(repl_str, class_coupling+cl[:4]) + sign_true_2.format(cl)
                            text_val = cl
                        elif cl[:4] == 'seco':
                            html_val = sign_true_1.replace(repl_str, class_coupling+cl[:4]) + sign_true_2.format(cl)
                            text_val = cl
                        elif cl[:2] == 'no':
                            html_val = sign_true_1.replace(repl_str, class_coupling+cl[:2]) + sign_true_2.format(cl)
                            text_val = cl
                        else:
                            html_val = sign_false
                            text_val = ''
                        out[entry][source][gprot]['html'] = html_val
                        out[entry][source][gprot]['text'] = text_val
                    else:
                        out[entry][source][gprot]['html'] = sign_false
                        out[entry][source][gprot]['text'] = ''
                else:
                    out[entry][source][gprot]['html'] = sign_false
                    out[entry][source][gprot]['text'] = ''

    return out


def prepare_coupling_data_container():
    class_names = {}
    data = {}

    complex_objs = SignprotComplex.objects.prefetch_related('structure__protein_conformation__protein').values_list('structure__protein_conformation__protein__parent_id', flat=True)
    proteins = (
        Protein.objects.filter(
            sequence_type__slug="wt",
            family__slug__startswith="00",
            species__common_name="Human")
        # .exclude(id__in=complex_objs)
        .prefetch_related("family")
    )

    for p in proteins:
        p_class = p.family.slug.split("_")[0]
        if p_class not in class_names:
            class_names[p_class] = re.sub(
                r"\([^)]*\)", "", p.family.parent.parent.parent.name
            )
        p_class_name = class_names[p_class].strip()

        data[p.entry_short()] = {
            "class": p_class_name,
            "pretty": p.short()[:15],
            "GuideToPharma": {},
            "Aska": {},
            "rec_class": p.get_protein_class(),
            "rec_obj": p,
            "rec_pdb": p.entry_short(),
        }

    return data


def fill_coupling_data_container(data):
    threshold_primary = -0.1
    threshold_secondary = -1

    distinct_sources = ["GuideToPharma", "Aska"]

    couplings = ProteinCouplings.objects.filter(source__in=distinct_sources).prefetch_related(
        "protein", "g_protein_subunit", "g_protein"
    )

    for c in couplings:
        p = c.protein.entry_short()
        # Skip entries without any annotation
        if p not in data:
            continue

        s = c.source
        t = c.transduction
        m = c.logmaxec50
        gf = c.g_protein.name
        gf = gf.replace(" family", "")

        if c.g_protein_subunit:
            g = c.g_protein_subunit.entry_name
            g = g.replace("_human", "")

        if s not in data[p]:
            data[p][s] = {}

        if gf not in data[p][s]:
            data[p][s][gf] = {}

        # If transduction in GuideToPharma data
        if t:
            data[p][s][gf] = t
        else:
            if "subunits" not in data[p][s][gf]:
                data[p][s][gf] = {"subunits": {}, "best": -2.00}
            data[p][s][gf]["subunits"][g] = round(Decimal(m), 2)
            if round(Decimal(m), 2) == -0.00:
                data[p][s][gf]["subunits"][g] = 0.00
            # get the lowest number into 'best'
            if m > data[p][s][gf]["best"]:
                data[p][s][gf]["best"] = round(Decimal(m), 2)


    return data


def process_coupling_data(data):
    res = []

    for entry in data.keys():
        i = data[entry]

        e = {}

        c_gtop = extract_coupling_bool(i, "GuideToPharma")
        p_gtop = extract_coupling_primary(i["GuideToPharma"])

        c_aska = extract_coupling_bool(i, "Aska")
        p_aska = extract_coupling_primary(c_aska[1])

        c_merg = extract_coupling_bool(i, "Merged")
        p_merg = extract_coupling_primary(c_merg[1])

        e['coupling'] = {}
        e["GuideToPharma"] = {}
        e["Aska"] = {}
        e["Merged"] = {}

        e["rec_class"] = i["rec_class"]
        e["rec_obj"] = i["rec_obj"]
        e["key"] = entry

        e["coupling"]["GuideToPharma"] = i["GuideToPharma"]
        e["coupling"]["Aska"] = c_aska[1]
        e["coupling"]["Merged"] = c_merg[1]

        for x in ["Gs", "Gi/o", "Gq/11", "G12/13"]:
            e["GuideToPharma"][x] = c_gtop[x]
            e["Aska"][x] = c_aska[0][x]
            e["Merged"][x] = c_merg[0][x]

        e["GuideToPharma"]["gprot"] = p_gtop
        e["Aska"]["gprot"] = p_aska
        e["Merged"]["gprot"] = p_merg

        res.append(e)

    return res


def extract_coupling_bool(gp, source):
    distinct_g_families = ['Gs','Gi/o', 'Gq/11', 'G12/13', ]
    threshold_primary = -0.1
    threshold_secondary = -1

    if source == 'GuideToPharma':
        gp = gp[source]
        c = {"Gi/o": False, "Gs": False, "Gq/11": False, "G12/13": False}
        for key in c:
            if key in gp:
                c[key] = True
        return c

    elif source == 'Aska':
        gp = gp[source]
        c = {"Gi/o": False, "Gs": False, "Gq/11": False, "G12/13": False}
        c_levels = {}

        for gf in distinct_g_families:
            if gf in gp:
                if gp[gf]['best']>threshold_primary:
                    c[gf] = True
                    c_levels[gf] = "primary"
                elif gp[gf]['best']>threshold_secondary:
                    c[gf] = True
                    c_levels[gf] = "secondary"
                else:
                    c[gf] = True
                    c_levels[gf] = "no coupling"
        return (c, c_levels)

    elif source == 'Merged':
        c = {"Gi/o": False, "Gs": False, "Gq/11": False, "G12/13": False}
        c_levels = {}
        v = gp

        for gf in distinct_g_families:
            values = []
            if 'GuideToPharma' in v and gf in v['GuideToPharma']:
                values.append(v['GuideToPharma'][gf])
            if 'Aska' in v and gf in v['Aska']:
                best = v['Aska'][gf]['best']
                if best > threshold_primary:
                    values.append('primary')
                elif best > threshold_secondary:
                    values.append('secondary')
                else:
                    values.append("no coupling")
            if 'primary' in values:
                c[gf] = True
                c_levels[gf] = "primary"
            elif 'secondary' in values:
                c[gf] = True
                c_levels[gf] = "secondary"
            elif 'no coupling' in values:
                c[gf] = True
                c_levels[gf] = "no coupling"

        return (c, c_levels)

def extract_coupling_primary(gp):
    p = []
    for key in gp:
        if gp[key] == "primary":
            p.append(key)
    return p
