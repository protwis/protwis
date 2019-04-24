import json
import re
import time
from itertools import chain

from residue.models import ResidueGenericNumberEquivalent
from common.definitions import *

from django.core.exceptions import ObjectDoesNotExist


def get_entry_names(request):
    """Extract a list of entry names from the post request"""
    return request.POST.getlist("pos[]")


def get_ignore_info(request):
    """Extract a dict of residues to ignore from the post request"""
    ignore_dict = request.POST.get("ignore")
    return json.loads(ignore_dict)


def get_protein_segments(request):
    """From a list of given generic numbers (3x50), query appropriate generic residue
    number objects"""
    segments = []
    segment_raw = request.POST.getlist("seg[]")
    for s in segment_raw:
        try:
            gen_object = ResidueGenericNumberEquivalent.objects.filter(
                label=s, scheme__slug__in=["gpcrdba"]
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
                if dn != "":
                    rexp = r"(?<=<b>)\d{1,}|\.?\d{2,}[\-?\d{2,}]*|x\d{2,}"
                    gn = re.findall(rexp, dn)
                else:
                    gn = "".join([str(trans[elem]), ".", str(x)])
                gnl.append("".join(gn))
            generic_numbers.append(gnl)
    return generic_numbers


def get_signature_features(signature_data, generic_numbers, feats):
    """Extract the signature features and prepare for visualization"""
    signature_features = []
    x = 0
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
