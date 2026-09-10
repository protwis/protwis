"""Find candidate per-TM anchor generic-number positions for the snake plot.

Standalone analysis script (not a management command) for finding candidate
per-TM "anchor" generic-number positions used by DrawSnakePlot
(common/diagrams_gpcr.py) as a static, per-render lookup instead of a live
zero-crossing calculation.

Run from inside the app container:
    docker compose exec app python common/diagrams_gpcr_constants_calc.py

An anchor is only useful for aligning TMs to the membrane midplane if it is
BOTH:
  (a) present in virtually every GPCR's sequence (coverage), and
  (b) actually close to midplane_distance == 0 (otherwise averaging raw
      anchor row-positions across TMs doesn't approximate the true membrane
      midplane at all - each TM's anchor would just be some other, unrelated,
      well-conserved-but-off-center position).

For each TM: compute cross-class mean midplane_distance (from the
highest-resolution structure per protein) for every generic-number position,
restrict to those within NEAR_ZERO_THRESHOLD, and take up to TOP_N_ANCHORS of
those ranked by closeness to zero (best alignment first) as an ordered
fallback list - DrawSnakePlot tries them in order and uses whichever is
present in a given receptor's own sequence.

Prints a ready-to-paste TM_ANCHOR_CANDIDATES block for common/diagrams_gpcr.py.
"""
import os
import sys
import statistics
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'protwis.settings')
import django
django.setup()

from residue.models import Residue
from protein.models import Protein
from structure.models import Structure
from angles.models import ResidueAngle

TMS = ['TM1', 'TM2', 'TM3', 'TM4', 'TM5', 'TM6', 'TM7']
CLASS_SLUG_PATH = 'structure__protein_conformation__protein__parent__family__parent__parent__parent__slug'
NEAR_ZERO_THRESHOLD = 3.0  # max |mean midplane_distance| (cross-class pooled) to be anchor-eligible
TOP_N_ANCHORS = 3  # how many eligible candidates (closest-to-zero first) to keep as fallbacks per TM


def stage_a_coverage():
    """Per TM: {generic_number_label: set(protein_id)} covering every WT protein with a residue there."""
    total_wt = Protein.objects.filter(sequence_type__slug='wt').count()

    covering = {tm: defaultdict(set) for tm in TMS}
    rows = Residue.objects.filter(
        protein_conformation__protein__sequence_type__slug='wt',
        protein_segment__slug__in=TMS,
        generic_number__isnull=False,
    ).values_list('protein_segment__slug', 'generic_number__label', 'protein_conformation__protein_id')

    for tm, gn_label, protein_id in rows:
        covering[tm][gn_label].add(protein_id)

    return covering, total_wt


def build_highest_resolution_pool():
    all_qs = Structure.objects.filter(refined=False, structure_type__origin='experiment')
    best_by_protein = {}
    for s in all_qs.select_related('protein_conformation__protein__parent').order_by('resolution'):
        pid = s.protein_conformation.protein.parent_id
        if pid not in best_by_protein:
            best_by_protein[pid] = s.pk
    return Structure.objects.filter(pk__in=list(best_by_protein.values()))


def compute_distance_means(pool_qs):
    """Per TM: {generic_number_label: mean_midplane_distance}, cross-class pooled, for every position with data."""
    rows = ResidueAngle.objects.filter(
        structure__in=pool_qs,
        residue__generic_number__isnull=False,
        residue__protein_segment__slug__in=TMS,
    ).values_list('residue__protein_segment__slug', 'residue__generic_number__label', 'midplane_distance')

    values_by_tm_label = {tm: defaultdict(list) for tm in TMS}
    for tm, label, value in rows:
        if value is not None:
            values_by_tm_label[tm][label].append(value)

    return {
        tm: {label: statistics.mean(values) for label, values in labels.items()}
        for tm, labels in values_by_tm_label.items()
    }


def main():
    print("Stage A: computing per-TM generic-number coverage across all WT GPCR proteins...")
    covering, total_wt = stage_a_coverage()

    print("Building highest_resolution structure pool (one best-resolution structure per protein)...")
    pool_qs = build_highest_resolution_pool()
    print("  highest_resolution: {} structures".format(pool_qs.count()))

    print("Computing cross-class mean midplane_distance for every candidate generic number...")
    distance_means = compute_distance_means(pool_qs)

    anchors_per_tm = {}
    for tm in TMS:
        all_candidates = covering[tm]
        means_for_tm = distance_means.get(tm, {})

        print("\n{} candidates (|mean| <= {} required to be anchor-eligible):".format(tm, NEAR_ZERO_THRESHOLD))
        ranked = sorted(
            all_candidates.keys(),
            key=lambda l: abs(means_for_tm[l]) if l in means_for_tm else float('inf'))
        for label in ranked[:15]:
            mean = means_for_tm.get(label)
            coverage_pct = 100 * len(all_candidates[label]) / total_wt if total_wt else 0
            eligible = mean is not None and abs(mean) <= NEAR_ZERO_THRESHOLD
            mean_str = "{:+.2f}".format(mean) if mean is not None else "no data"
            print("    [{}] {}: mean={}  coverage={:.1f}%".format(
                "OK" if eligible else "  ", label, mean_str, coverage_pct))

        eligible_labels = [
            label for label in all_candidates
            if label in means_for_tm and abs(means_for_tm[label]) <= NEAR_ZERO_THRESHOLD
        ]
        if not eligible_labels:
            print("  WARNING: no candidates for {} met |mean| <= {}; falling back to unfiltered "
                  "closest-to-zero selection (fix likely needs a wider threshold or more structure data)".format(
                      tm, NEAR_ZERO_THRESHOLD))
            eligible_labels = list(all_candidates.keys())

        eligible_labels.sort(key=lambda l: abs(means_for_tm.get(l, float('inf'))))
        chosen_labels = eligible_labels[:TOP_N_ANCHORS]
        anchors_per_tm[tm] = (chosen_labels, all_candidates, means_for_tm)

    print("\n" + "=" * 60)
    print("ANCHOR SELECTION REPORT (highest_resolution pool)")
    print("=" * 60)

    for tm in TMS:
        chosen_labels, all_candidates, means_for_tm = anchors_per_tm[tm]
        print("\n{} anchor set (top {} closest-to-zero, near-zero-eligible candidates):".format(tm, TOP_N_ANCHORS))
        if not chosen_labels:
            print("  no eligible generic numbers found for this TM")
            continue
        cumulative_covered = set()
        for i, label in enumerate(chosen_labels, start=1):
            cumulative_covered |= all_candidates[label]
            pct = 100 * len(cumulative_covered) / total_wt if total_wt else 0
            mean = means_for_tm.get(label)
            print("  {}. {}: mean midplane_distance={:+.2f}  ({}/{}, {:.1f}% cumulative coverage)".format(
                i, label, mean, len(cumulative_covered), total_wt, pct))

    print("\n" + "=" * 60)
    print("---- copy into common/diagrams_gpcr.py ----")
    print("=" * 60)
    print("TM_ANCHOR_CANDIDATES = {")
    for tm in TMS:
        chosen_labels, _, means_for_tm = anchors_per_tm[tm]
        means_comment = ', '.join("{}: {:+.2f}".format(l, means_for_tm[l]) for l in chosen_labels)
        tm_num = tm.replace('TM', '')
        print("    {}: {},  # mean midplane_distance: {}".format(tm_num, chosen_labels, means_comment))
    print("}")
    print("---- end ----")


if __name__ == '__main__':
    main()
