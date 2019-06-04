import re
import django
from itertools import combinations, chain
import pickle
from decimal import *
from pathlib import Path

from IPython.display import display
import pandas as pd

from protein.models import Protein, ProteinSegment, ProteinFamily, ProteinGProteinPair


def get_receptor_classes():
    return ProteinFamily.objects.filter(parent=1).exclude(slug=100)


def get_gprot_classes():
    return ProteinFamily.objects.filter(parent=533)


def get_receptor_segments():
    return ProteinSegment.objects.filter(proteinfamily="GPCR")


def get_gprot_segments():
    return ProteinSegment.objects.filter(proteinfamily="Alpha")


def extract_per_attr(data, attr, name):
    """Return elements of data for which attr is equal to name"""
    return [elem for elem in data if name in elem[attr]]


def calc_consensus_from_signature(signature_dict):
    from signprot.interactions import get_generic_numbers, get_signature_consensus
    
    signature = signature_dict["signature"]

    sig_data = signature.prepare_display_data()
    gn = get_generic_numbers(sig_data)
    gn_flat = list(chain.from_iterable(gn))

    signature_dict["consensus"] = get_signature_consensus(sig_data, gn_flat)
    return signature_dict


def aggregate_consensus_data(entry, origin=None):
    data = []
    consensus = entry["consensus"]
    rec_class = entry["rec_class"]
    gprot = entry["gprot"]
    while len(consensus) > 0:
        a = consensus.pop()
        a["gprot"] = gprot
        a["rec_class"] = rec_class
        a["origin"] = origin
        data.append(a)
    return data


def summarize_df(df):
    print("Dataframe description:")
    display(df.describe())
    print("\n")

    print("Dataframe size:")
    print(df.shape)
    print("\n")

    display(df.head())


def show_group_top_n(df, group='feature', n=5):
    count_col = "{}_count".format(group)
    print(count_col)
    d1 = df.groupby(group)[group].agg(
        {count_col: len}).sort_values(
        count_col, ascending=False).head(n).reset_index()

    display(d1)


def compare_sets(
    df1, df2, method=set.intersection, drop_list=["origin", "key", "score", "cons"]
):
    v_ds1 = df1.drop(drop_list, 1)
    v_ds2 = df2.drop(drop_list, 1)
    colnames = v_ds1.columns

    summarize_df(v_ds1)
    summarize_df(v_ds2)

    ds1 = set([tuple(line) for line in v_ds1.values])
    ds2 = set([tuple(line) for line in v_ds2.values])

    comp = pd.DataFrame(list(method(ds2, ds1)))
    try:
        comp.columns = colnames
        return comp
    except ValueError as e:
        print("Value Error\n{}:\nNo entries overlap between the two sets.".format(e))


def prepare_coupling_data_container():
    class_names = {}
    data = {}

    proteins = (
        Protein.objects.filter(
            sequence_type__slug="wt",
            family__slug__startswith="00",
            species__common_name="Human",
        )
        .all()
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


def fill_coupling_data_container(data, sources=["GuideToPharma", "Aska"]):
    distinct_g_families = []
    distinct_g_subunit_families = {}
    distinct_sources = sources

    couplings = ProteinGProteinPair.objects.all().prefetch_related(
        "protein", "g_protein_subunit", "g_protein"
    )

    for c in couplings:
        p = c.protein.entry_short()
        s = c.source
        t = c.transduction
        m = c.log_rai_mean
        gf = c.g_protein.name
        # print(gf)
        gf = gf.replace(" family", "")

        if s not in distinct_sources:
            continue

        if gf not in distinct_g_families:
            distinct_g_families.append(gf)
            distinct_g_subunit_families[gf] = []

        if c.g_protein_subunit:
            g = c.g_protein_subunit.entry_name
            g = g.replace("_human", "")
            # print("g",g)
            if g not in distinct_g_subunit_families[gf]:
                distinct_g_subunit_families[gf].append(g)
                distinct_g_subunit_families[gf] = sorted(
                    distinct_g_subunit_families[gf]
                )

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

        c = extract_coupling_bool(i["GuideToPharma"])
        p = extract_coupling_primary(i["GuideToPharma"])

        e["rec_class"] = i["rec_class"]
        e["rec_obj"] = i["rec_obj"]
        e["key"] = entry
        e["coupling"] = i["GuideToPharma"]
        e["Gi/Go"] = c["Gi/Go"]
        e["Gs"] = c["Gs"]
        e["Gq/G11"] = c["Gq/G11"]
        e["G12/G13"] = c["G12/G13"]
        e["gprot"] = p

        res.append(e)

    return res


def extract_coupling_bool(gp):
    c = {"Gi/Go": False, "Gs": False, "Gq/G11": False, "G12/G13": False}
    for key in c:
        if key in gp:
            c[key] = True
    return c


def extract_coupling_primary(gp):
    p = []
    i = 0
    for key in gp:
        if gp[key] == "primary":
            p.append(key)
    return p


# import cPickle as pickle
def pickle_signature(data, path, filename):


    with open(path+filename, 'wb+') as out_file:

        p = pickle.Pickler(out_file)
        p.fast = True
        p.dump(data) # d is your dictionary
        p.clear_memo()

        # pickle.dump(data, out_file)


def load_pickle_signature(path, result_file, which, with_type):
    if with_type == 0:
        wt = 'file_with'
    elif with_type == 1:
        wt = 'file_wo'
    else:
        print('Pick either 0 for "file_with" or 1 for "file_wo" column')
        return

    file_path = Path(path + result_file.iloc[which][wt])
    with file_path.open('rb') as f:
        obj = pickle.load(f)

    return obj


def consensus_to_dataframe(entry):
    data = []
    consensus = entry['consensus']
    rec_class = entry['rec_class']
    gprot = entry['gprot']
    while len(consensus) > 0:
        a = consensus.pop()
        a['gprot'] = gprot
        a['rec_class'] = rec_class
        data.append(a)

    ds1 = pd.DataFrame(data)
    return ds1


