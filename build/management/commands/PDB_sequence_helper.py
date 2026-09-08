"""
PDB_sequence_helper.py

This module contains helper functions for extracting, aligning, and analyzing
GPCR structure sequences (PDB vs. wild-type) in the GPCRdb/Protwis project.
It includes:

1. **Sequence Extraction & Distance Calculation**
   - `generate_seq_and_distances_from_pdb_text()`: Parse raw PDB text (from the DB),
     build the amino acid sequence for a specified chain, and compute CA–CA distances.
   - `distances_stats()`: Compute mean, standard deviation, and identify outlier CA–CA
     distances, which often reveal potential misalignment or unusual structural
     features.

2. **Alignment & Alignment Corrections**
   - `pre_align_modifications()`: Applies custom sequence slicing/trimming for
     certain known tricky PDB codes. We shouldn't need this for the PairwiseAligner.
   - `decide_penalty()`: Chooses gap/mismatch scores for local alignment, depending
     on the PDB code's known quirks. We got these values from previous code.
   - `run_pairwisealigner()`: Runs a Biopython PairwiseAligner for local alignment
     of the wild-type and extracted PDB sequences.
   - `format_local_alignment_pairwise2_blocks()`: Reconstructs pairwise2‐like
     alignment strings from the PairwiseAligner blocks.

3. **Misalignment Detection & Fixing**
   - `detect_alignment_mistakes_and_reposition()`: Identifies suspicious “misplaced”
     chunks near gaps and moves them to a more appropriate place in the alignment.
     This is a one-pass fix that can address repeated motifs and duplicated segments
     that appear out of order.

Typical Workflow
----------------
1) **Sequence + Distance Extraction**
   - Use `generate_seq_and_distances_from_pdb_text(pdb_text, preferred_chain)`
     to get `pdb_seq` plus a list of CA–CA distances and author residue numbers.

2) **Outlier Detection**
   - Pass the distances to `distances_stats()` to identify positions that deviate
     significantly from expected values (the "outlier_indexes").

3) **Alignment**
   - Run `run_pairwisealigner(pdb_code, wt_seq, pdb_seq)` to obtain two alignment
     strings (`ref_seq`, `temp_seq`) and a mapping (`pdb_map`) from the raw PDB
     indices to alignment positions.

4) **Misalignment Fix**
   - Call `detect_alignment_mistakes_and_reposition(...)` with the alignment plus
     the outlier indexes. This will detect a suspicious chunk near a gap and move
     it to the left edge of the gap if appropriate.

"""

import json
import os
from io import StringIO
from django.conf import settings
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

# from Bio.PDB.Polypeptide import PPBuilder
from Bio.Align import PairwiseAligner
import numpy as np

# Pre-computed alignments (pdb_code -> [ref_seq, temp_seq]) for structures where
# the automatic PairwiseAligner can't be trusted to reproduce the correct result.
# Same file/format used by add_annotation.py's Command.custom_mappings.
with open(
    os.sep.join(
        [settings.DATA_DIR, "structure_data", "annotation", "custom_mappings.json"]
    )
) as _cm:
    CUSTOM_MAPPINGS = json.load(_cm)


# preferred_chain = structure.preferred_chain
# if len(preferred_chain.split(','))>1: #if A,B
#     preferred_chain = preferred_chain.split(',')[0]


def generate_seq_and_distances_from_pdb_text(
    pdb_text, preferred_chain="A", residues_to_remove=None, aa_map=None
):
    """
    Extracts a protein sequence and computes C-alpha (CA) distances from a raw PDB string.
    Optionally removes specified residues before processing.

    This function:
      1. Parses the provided PDB text with Biopython's PDBParser.
      2. Selects the specified chain (`preferred_chain`).
      3. If `residues_to_remove` is provided, detaches those residues from the chain.
      4. Iterates over residues, selecting which to include per `aa_map` (see below).
      5. Builds the amino-acid sequence (single-letter codes).
      6. Computes successive CA–CA distances, or inserts `None` if CA is missing.

    Parameters
    ----------
    pdb_text : str
        Full PDB content as a single string (e.g., from `structure.pdb_data.pdb`).
    preferred_chain : str, optional
        The chain identifier to parse, by default 'A'.
    residues_to_remove : list of int, optional
        A list of PDB residue sequence numbers to remove from the chain before
        generating the sequence and distances.
    aa_map : dict, optional
        Maps PDB three-letter residue names to one-letter codes (e.g. the
        `AA` dict in build_structures.py, which also covers modified
        residues like YCM/CSD/TYS/SEP/TPO). When provided, a residue is
        included iff its name is a key in this map -- regardless of its
        hetero-atom flag -- and the caller MUST use the identical map to
        build any other residue index (e.g. build_structures.py's
        `pdbseq`) it intends to look up against this function's output by
        position, since the two indices only stay aligned if both make the
        same inclusion decision for every residue. When omitted, falls back
        to the original behavior: skip all hetero-atom-flagged residues,
        translate the rest via `Bio.SeqUtils.seq1` (falling back to "X" for
        an unrecognized standard-flagged residue).

    Returns
    -------
    tuple
        (pdb_seq, distances, resnums) where:
          - pdb_seq (str): the one-letter amino acid sequence extracted from the chain
          - distances (list of float or None): CA–CA distances for each successive pair,
            with None inserted if no distance can be computed at that position.
          - resnums (list of int): the author (PDB) residue sequence number for each
            residue in `pdb_seq`, in chain order. A non-consecutive jump between two
            entries indicates a real gap (e.g. a construct-level deletion) even when
            the corresponding CA–CA distance looks like a normal bonded distance --
            see `find_chain_breaks`.

    Raises
    ------
    ValueError
        If the specified chain is not found in the parsed structure.
    """
    # 1. Parse the PDB data using Biopython
    parser = PDBParser(PERMISSIVE=True, QUIET=True)
    biopython_structure = parser.get_structure("ref", StringIO(pdb_text))[0]

    # 2. Access the preferred chain
    if preferred_chain not in biopython_structure:
        raise ValueError(f"Preferred chain '{preferred_chain}' not found.")

    chain = biopython_structure[preferred_chain]

    # 3. *** NEW: Remove specified residues ***
    if residues_to_remove:
        # PDB residue ID is a tuple: (hetero_flag, sequence_number, insertion_code)
        # We only care about the sequence_number, assuming standard residues.
        ids_to_remove = [res_id for res_id in residues_to_remove]

        # We must collect IDs first then detach, to avoid modifying the list while iterating
        residue_ids_in_chain = [res.get_id() for res in chain]

        for res_id_tuple in residue_ids_in_chain:
            # res_id_tuple is like (' ', 25, ' ')
            if res_id_tuple[1] in ids_to_remove:
                chain.detach_child(res_id_tuple)

    # 4. Iterate over residues in the modified chain
    sequence = ""
    distances = []
    resnums = []
    prev_ca = None

    for residue in chain:
        resname = residue.get_resname()

        if aa_map is not None:
            aa = aa_map.get(resname)
            if aa is None:
                continue
        else:
            # Skip heteroatoms or waters
            if residue.id[0] != " ":
                continue
            try:
                aa = seq1(resname)
            except Exception:
                aa = "X"  # Unknown amino acid

        sequence += aa
        resnums.append(residue.id[1])

        # First residue has no preceding one, so distance is None
        if prev_ca is None:
            distances.append(None)
        # For subsequent residues, compute distance or use placeholder
        elif "CA" in residue:
            current_ca = residue["CA"].get_vector()
            distances.append((current_ca - prev_ca).norm())
        else:
            # Insert placeholder for missing distance
            distances.append(None)

        # Update prev_ca for the next iteration
        if "CA" in residue:
            prev_ca = residue["CA"].get_vector()
        else:
            prev_ca = None

    return sequence, distances, resnums


def distances_stats(distances, threshold=3, debug=False):
    """
    Computes basic statistics on CA–CA distances and identifies outliers.

    Parameters
    ----------
    distances : list of float or None
        The list of CA–CA distances from a protein structure sequence.
    threshold : int or float, optional
        The multiplier of standard deviation used to define outliers
        (default is 3, i.e., values beyond mean ± 3*std).
    debug : bool, optional
        If True, prints a summary of mean, std, and outlier positions for
        reference, by default False.

    Returns
    -------
    outlier_indexes : list of int
        A list of zero-based positions where distances are flagged as outliers.

    Notes
    -----
    - Ignores `None` entries in the distance list when computing stats.
    - Prints a summary of mean, std, and outlier positions for reference.
    """

    # Filter out None values
    filtered_distances = [d for d in distances if d is not None]

    # Compute mean and standard deviation using the full dataset
    mean_all = np.mean(filtered_distances)
    std_dev_all = np.std(filtered_distances, ddof=1)

    # Set threshold (you can adjust k, default is 3)
    k = threshold
    lower_bound = mean_all - k * std_dev_all
    upper_bound = mean_all + k * std_dev_all

    # Select inliers (values within the acceptable range)
    filtered_data = [d for d in filtered_distances if lower_bound <= d <= upper_bound]

    # Identify outliers
    outliers = [d for d in filtered_distances if d < lower_bound or d > upper_bound]

    # Print results
    if debug:
        print(f"Mean of all distances: {mean_all:.4f}")
        print(f"Standard deviation of all distances: {std_dev_all:.4f}")
        print(f"Lower bound for normal values: {lower_bound:.4f}")
        print(f"Upper bound for normal values: {upper_bound:.4f}")
        print(f"Filtered mean (after removing true outliers): {np.mean(filtered_data):.4f}")
        print(f"Filtered variance: {np.var(filtered_data, ddof=1):.4f}")
        print(f"Filtered standard deviation: {np.std(filtered_data, ddof=1):.4f}")
        print(f"True outliers: {outliers}")
    outlier_indexes = [i for i, d in enumerate(distances) if d in outliers]
    for idx, d in enumerate(distances):
        if d is None:
            continue
        # if d < lower_bound or d > upper_bound:
        #     print(f"Outlier at position {idx}")
    return outlier_indexes


def find_chain_breaks(distances, threshold=4.5, resnums=None):
    """
    Flags real backbone chain breaks using an absolute CA-CA distance
    threshold, rather than `distances_stats`'s whole-structure relative
    z-score, optionally combined with author residue numbering.

    A real peptide-bonded CA-CA distance is tightly clustered around 3.8 A
    regardless of local backbone conformation, so any missing residue pushes
    the effective distance far beyond that. A fixed threshold catches this
    reliably even when the rest of the structure's distance distribution has
    enough natural spread to hide a real break from a `mean +/- 3*std` test
    (as happens with `distances_stats`).

    Distance alone still misses one real case: a construct-level deletion
    (e.g. a flexible loop deliberately excised for crystallization/cryo-EM)
    can leave the two flanking residues directly, genuinely bonded to each
    other at a completely normal ~3.8 A distance, even though a real gap in
    the wild-type sequence separates them. When `resnums` is given, any
    non-consecutive jump in author residue numbering between two
    consecutively-modeled residues is flagged as a break too, independent of
    distance -- this is the only reliable signal for that case.

    Parameters
    ----------
    distances : list of float or None
        CA-CA distances, as returned by
        `generate_seq_and_distances_from_pdb_text` (distances[i] is the
        distance between residue i-1 and i).
    threshold : float, optional
        Distance in Angstrom beyond which a gap is considered a real chain
        break, by default 4.5.
    resnums : list of int, optional
        Author (PDB) residue sequence numbers, as returned by
        `generate_seq_and_distances_from_pdb_text`, parallel to `distances`.
        When provided, a break is also flagged wherever `resnums[i] -
        resnums[i-1] != 1`.

    Returns
    -------
    list of int
        Indexes into `distances`/the sequence where a chain break occurs
        (index i means the break falls between residue i-1 and i).
    """
    breaks = {i for i, d in enumerate(distances) if d is not None and d > threshold}
    if resnums is not None:
        for i in range(1, len(resnums)):
            if resnums[i] - resnums[i - 1] != 1:
                breaks.add(i)
    return sorted(breaks)


# Author residue number after which a real construct-level deletion exists but
# was renumbered continuously by the depositor (no numbering jump) across a
# normal ~3.8 A CA-CA bond -- the one case `find_chain_breaks` cannot detect on
# its own (see its docstring). Confirmed by aligning each PDB code's sequence
# against its WT reference: the initial alignment already places the flanking
# residues correctly on either side of a real gap, but `consolidate_structural_islands`
# then merges them into one "island" for lack of a detected break, destroying
# that correct alignment. Forcing the break here lets it recognize the two
# sides as separate, already-correctly-aligned islands and leave them alone.
MANUAL_CHAIN_BREAK_AFTER_RESNUM = {
    "9JQZ": 214,
    "9JQY": 214,
}


def manual_chain_break_indexes(resnums, pdb_code):
    """
    Extra `find_chain_breaks`-style break indexes for a real construct
    deletion that a depositor renumbered straight through, based on
    `MANUAL_CHAIN_BREAK_AFTER_RESNUM`.

    Parameters
    ----------
    resnums : list of int
        Author (PDB) residue sequence numbers, parallel to `pdb_seq`.
    pdb_code : str
        The PDB code, used to look up a manual break point.

    Returns
    -------
    list of int
        Indexes into `resnums`/`pdb_seq` where a break should be forced
        (index i means the break falls between residue i-1 and i), suitable
        for merging into `find_chain_breaks`'s output.
    """
    after_resnum = MANUAL_CHAIN_BREAK_AFTER_RESNUM.get(pdb_code)
    if after_resnum is None:
        return []
    return [i for i in range(1, len(resnums)) if resnums[i - 1] == after_resnum]


def pre_align_modifications(pdb_code, seq):
    # this function is used to modify the extracted sequence before alignment
    # only for old pairwise2 alignment, not for pairwise aligner

    if pdb_code == "6U1N":
        return seq[:265]
    elif pdb_code in ["1GZM", "3C9L"]:
        return seq[:-3]
    elif pdb_code == "8WU1":
        return seq[:-13]
    else:
        return seq


def decide_penalty(pdb_code):
    """
    Chooses gap and mismatch penalties for the local alignment, based on prior
    experience with certain structures. Some PDB codes require non-default
    values to avoid pathological alignments.

    Parameters
    ----------
    pdb_code : str
        The PDB code, used to look up special penalty values.

    Returns
    -------
    tuple of four floats
        These are the values used for PairwiseAligner scoring.
    """
    if pdb_code in [
        "6NBI",
        "6NBF",
        "6NBH",
        "6U1N",
        "6M1H",
        "6PWC",
        "7JVR",
        "7SHF",
        "7EJ0",
        "7EJ8",
        "7EJA",
        "7EJK",
        "7VVJ",
        "7TS0",
        "7W6P",
        "7W7E",
        "8IRS",
        "8FLQ",
        "8FLR",
        "8FLS",
        "8FLU",
        "8FU6",
        "8IRU",
        "7Y35",
        "7Y36",
        "8TB7",
        "8SZI",
        "8TZQ",
        "8U02",
        "8W8Q",
        "8W8R",
        "8W8S",
        "9N29",
        "9PEE",
        "9N09",
        "9UAP",
        "9UCP",
        "8ZWF",
        "9ISI",
        "9PLO",
        "9PLN",
        "9PQD",
        "9LFA",
        "9UAZ"
    ]:
        return 3, -4, -3, -1
    elif pdb_code in ["5VEW","5VEX","6LN2","6KJV","6KK7","6KK1","8HN8","8HOC"]:
        # 5VEW's engineered disulfide (I317C/G361C, see structures.tsv) puts two
        # point mutations one residue apart. With cheap gaps (open<=-5), the
        # aligner "wobbles" -- opens a 1-residue gap on each side right next to
        # each other instead of accepting two plain substitutions -- and since
        # only one side's gap count feeds the WT-position bookkeeping further
        # down in create_rotamers(), that balanced-looking wobble permanently
        # shifts every residue after it by one WT position for the rest of the
        # chain. A steeper gap-open cost makes the straight, no-gap alignment
        # win instead. 5VEW previously shared the cheap-gap group above, which
        # was the actual cause here, not a fix for it.
        return 3, -4, -8, -2
    elif pdb_code in ["6KUX", "6KUY", "6KUW"]:
        return 3, -4, -4, -1.5
    elif pdb_code in ["7YMJ"]:
        return 3, -5, -4, -4
    else:
        return 3, -4, -5, -2


def run_pairwisealigner(pdb_code, wt_seq, pdb_seq):
    """
    Runs a local alignment (PairwiseAligner) between a wild-type sequence
    and a PDB sequence, applying code-specific gap penalties.

    Parameters
    ----------
    pdb_code : str
        The PDB code for which we might choose specialized penalty settings.
    wt_seq : str
        The reference (wild-type) sequence.
    pdb_seq : str
        The extracted PDB sequence (ungapped).

    Returns
    -------
    tuple
        (ref_seq, temp_seq, pdb_map)

        - ref_seq (str): Gapped alignment of the wild-type sequence.
        - temp_seq (str): Gapped alignment of the PDB sequence.
        - pdb_map (dict): raw PDB index -> alignment index, used for
          mapping outlier positions back into 'temp_seq'.

    Notes
    -----
    If `pdb_code` has an entry in CUSTOM_MAPPINGS (loaded from
    gpcr/structure_data/annotation/custom_mappings.json), that pre-computed
    (ref_seq, temp_seq) pair is used directly and the PairwiseAligner is not
    run at all -- mirrors the override in add_annotation.py's handle().
    """
    if pdb_code in CUSTOM_MAPPINGS:
        ref_seq, temp_seq = CUSTOM_MAPPINGS[pdb_code]
        pdb_map = {}
        pdb_idx = 0
        for aln_idx, ch in enumerate(temp_seq):
            if ch != "-":
                pdb_map[pdb_idx] = aln_idx
                pdb_idx += 1
        return ref_seq, temp_seq, pdb_map, True

    a1, a2, a3, a4 = decide_penalty(pdb_code)
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = a1
    aligner.mismatch_score = a2
    aligner.open_gap_score = a3
    aligner.extend_gap_score = a4
    alignments = aligner.align(wt_seq, pdb_seq)
    best = alignments[0]
    parent_blocks, pdb_blocks = best.aligned
    ref_seq, temp_seq, pdb_map = format_local_alignment_pairwise2_blocks(
        wt_seq, pdb_seq, best
    )

    return ref_seq, temp_seq, pdb_map, False


def format_local_alignment_pairwise2_blocks(
    wt_seq, pdb_seq, best_alignment, return_pdb_map=True
):
    """
    Reconstructs a pairwise2-like alignment from Biopython PairwiseAligner output.

    The PairwiseAligner object returns 'aligned' blocks that specify index ranges
    in each sequence. We expand them into full strings with gaps for unaligned
    regions, reproducing a classic "pairwise2" style representation.

    Parameters
    ----------
    wt_seq : str
        Wild-type (target) sequence, ungapped.
    pdb_seq : str
        The PDB sequence being aligned, ungapped.
    best_alignment : Bio.Align.Alignment
        The best (highest score) alignment object from PairwiseAligner.align().
    return_pdb_map : bool, optional
        If True, build and return a dictionary mapping each residue index in
        `pdb_seq` to its alignment index in the returned 'temp_seq' string.

    Returns
    -------
    ref_str : str
        The final gapped alignment of the wild-type sequence.
    tmp_str : str
        The final gapped alignment of the PDB sequence (parallel to ref_str).
    pdb_map : dict, optional
        Maps raw PDB sequence indices (0-based) -> alignment indices. Only
        returned if `return_pdb_map` is True.

    Notes
    -----
    - If no aligned blocks exist, we simply pad the shorter sequence with dashes.
    - We preserve the original ordering of alignment blocks, inserting gaps
      where needed. This is reminiscent of what `pairwise2.align.localxx` might
      produce.
    """
    # The .aligned property has two arrays of shape (N,2):
    # best_alignment.aligned = (
    #    [(t_start, t_end), (t_start2, t_end2), ...],
    #    [(q_start, q_end), (q_start2, q_end2), ...]
    # )
    target_blocks, query_blocks = best_alignment.aligned

    # Convert them to Python lists of tuples
    target_blocks = [tuple(x) for x in target_blocks]
    query_blocks = [tuple(x) for x in query_blocks]

    # If there are no blocks, just align the entire sequences with dashes
    if len(target_blocks) == 0:
        length = max(len(wt_seq), len(pdb_seq))
        ref_str = wt_seq.ljust(length, "-")
        tmp_str = pdb_seq.ljust(length, "-")
        if return_pdb_map:
            # Build a trivial map: each residue in pdb_seq is at the same index
            pdb_map = {}
            for i, c in enumerate(pdb_seq):
                pdb_map[i] = i  # or min(i, length-1)
            return ref_str, tmp_str, pdb_map
        else:
            return ref_str, tmp_str

    # We'll build these as lists of strings and join at the end
    out_ref = []
    out_tmp = []

    # Optionally keep a map from raw PDB index -> alignment index
    pdb_map = {}
    pdb_curr = 0  # pointer into the original PDB sequence

    # We'll track how many chars we've output so far, so we can compute alignment indexes
    aln_pos = 0

    # Helper to append a chunk to out_tmp and update pdb_map
    def append_pdb_chunk(pdb_chunk):
        """
        pdb_chunk: string chunk that goes into out_tmp
        For each real residue (non-dash) in this chunk, map raw PDB index -> alignment index.
        """
        nonlocal pdb_curr, aln_pos

        out_tmp.append(pdb_chunk)
        # For each char in this chunk:
        for i, ch in enumerate(pdb_chunk):
            if ch != "-":
                pdb_map[pdb_curr] = aln_pos + i
                pdb_curr += 1
        # Advance the alignment position
        aln_pos += len(pdb_chunk)

    # 1) Sort blocks by the target (WT) start index
    blocks = sorted(zip(target_blocks, query_blocks), key=lambda x: x[0][0])

    # 2) Leading expansions
    (t0_start, t0_end), (q0_start, q0_end) = blocks[0]

    # We have some unaligned region at the front
    leadA = wt_seq[:t0_start]
    leadB = pdb_seq[:q0_start]

    if len(leadA) > len(leadB):
        out_ref.append(leadA)
        out_tmp_chunk = "-" * (len(leadA) - len(leadB)) + leadB
        append_pdb_chunk(out_tmp_chunk)
    elif len(leadB) > len(leadA):
        out_ref.append("-" * (len(leadB) - len(leadA)) + leadA)
        append_pdb_chunk(leadB)
    else:
        # same length
        out_ref.append(leadA)
        append_pdb_chunk(leadB)

    prev_t_end = t0_start
    prev_q_end = q0_start

    # 3) Iterate over each aligned block
    for i, ((t_start, t_end), (q_start, q_end)) in enumerate(blocks):
        # Gap between the previous block's end and this block's start
        gapA = wt_seq[prev_t_end:t_start]
        gapB = pdb_seq[prev_q_end:q_start]

        # Append gap region
        if len(gapA) > len(gapB):
            out_ref.append(gapA)
            append_pdb_chunk("-" * (len(gapA) - len(gapB)) + gapB)
        elif len(gapB) > len(gapA):
            out_ref.append("-" * (len(gapB) - len(gapA)) + gapA)
            append_pdb_chunk(gapB)
        else:
            # same length
            out_ref.append(gapA)
            append_pdb_chunk(gapB)

        # Now the aligned block
        blockA = wt_seq[t_start:t_end]
        blockB = pdb_seq[q_start:q_end]

        # If there's net gap difference, pad the shorter one
        lenA = len(blockA)
        lenB = len(blockB)
        if lenA > lenB:
            blockB += "-" * (lenA - lenB)
        elif lenB > lenA:
            blockA += "-" * (lenB - lenA)

        # Append them
        out_ref.append(blockA)
        append_pdb_chunk(blockB)

        # Update
        prev_t_end = t_end
        prev_q_end = q_end

    # 4) Trailing expansions
    last_t_end = blocks[-1][0][1]  # the target end of the last block
    last_q_end = blocks[-1][1][1]  # the query end of the last block

    tailA = wt_seq[last_t_end:]
    tailB = pdb_seq[last_q_end:]

    if len(tailA) > len(tailB):
        out_ref.append(tailA)
        append_pdb_chunk(tailB + "-" * (len(tailA) - len(tailB)))
    elif len(tailB) > len(tailA):
        out_ref.append("-" * (len(tailB) - len(tailA)) + tailA)
        append_pdb_chunk(tailB)
    else:
        # same length
        out_ref.append(tailA)
        append_pdb_chunk(tailB)

    # 5) Join them
    ref_str = "".join(out_ref)
    tmp_str = "".join(out_tmp)

    if return_pdb_map:
        return ref_str, tmp_str, pdb_map
    else:
        return ref_str, tmp_str


def detect_alignment_mistakes_and_reposition(
    pdb_code,
    wt_seq,
    pdb_seq,
    ref_seq,
    temp_seq,
    pdb_map,
    distances,
    outlier_indexes,
    aanumber=3,
    debug=False,
):
    """
    For each outlier distance, detect suspicious chunk near a large gap and "move" it
    to the start of that gap in the PDB alignment.

    Steps:
      1) Identify a chunk of residues in 'temp_seq' that appears right before a run of
         dashes ("the gap") and matches the segment in 'ref_seq' at the gap position.
      2) Remove the chunk from its old location (turning those chars into '-').
      3) Insert that chunk at the left edge of the gap, presumably where it
         belongs, to mirror the wild-type sequence.

    Parameters
    ----------
    pdb_code : str
        The PDB code (used for logging/printing).
    wt_seq : str
        The wild-type (ungapped) sequence.
    pdb_seq : str
        The PDB (ungapped) sequence.
    ref_seq : str
        The gapped alignment string for the wild-type sequence.
    temp_seq : str
        The gapped alignment string for the PDB sequence.
    pdb_map : dict
        A map from raw PDB residue index (0-based in pdb_seq) to alignment
        index in 'temp_seq'.
    distances : list of float or None
        The CA–CA distances for the PDB sequence, used to identify outliers.
    outlier_indexes : list of int
        The indexes in `pdb_seq` flagged as outliers.
    aanumber : int, optional
        How many residues to look backward from the outlier alignment position
        for a suspicious chunk, by default 3.
    debug : bool, optional
        If True, prints detailed progress and diagnostic information about
        each outlier and any repositioning performed, by default False.

    Returns
    -------
    str
        The updated 'temp_seq' alignment with suspicious chunks moved to a
        more likely position in the gap. This is a single-pass fix.

    NOTE:
    -----
    - Does not recalculate the alignment or 'pdb_map' after making changes.
      It's purely a quick fix for repeated motifs or duplications.
    - If the suspicious chunk is found at multiple positions, each might be
      moved in turn during this single pass. If deeper realignment is needed,
      you may want to recalculate after applying these fixes.
    """
    label_width = 50
    temp_seq_list = list(temp_seq)  # so we can mutate temp_seq characters

    if debug:
        print(f"\n=== Detecting alignment mistakes for {pdb_code} ===")
        print(f"Number of outliers: {len(outlier_indexes)}\n")

    for outlier_pdb_idx in outlier_indexes:
        if outlier_pdb_idx not in pdb_map:
            # print(
            #     f"  - Outlier PDB index {outlier_pdb_idx} is not in pdb_map; skipping."
            # )
            continue

        # Find where this residue is in the aligned PDB sequence
        temp_idx = pdb_map[outlier_pdb_idx]
        snippet_temp_seq_after_gap = "".join(temp_seq_list[temp_idx : temp_idx + 10])
        snippet_pdb_seq_after_gap = pdb_seq[outlier_pdb_idx : outlier_pdb_idx + 10]
        if debug:
            print(f"Outlier at raw PDB index {outlier_pdb_idx} (aligned pos {temp_idx})")
            print(f"Seq after the gap (temp_seq): '{snippet_temp_seq_after_gap}'")
            print(f"Seq after the gap (pdb_seq):  '{snippet_pdb_seq_after_gap}'")

        found_something = False

        # Look backward up to `aanumber` residues for suspicious chunk
        for block_size in range(1, aanumber + 1):
            block_start = temp_idx - block_size
            block_end = temp_idx
            if block_start < 0:
                break  # No room to look

            preceding_char_index = block_start - 1
            if preceding_char_index < 0:
                break

            # If the preceding char in temp_seq is a dash, that suggests a gap
            if temp_seq_list[preceding_char_index] == "-":
                # The suspicious chunk is what's in temp_seq right before the gap
                suspicious_chunk = "".join(temp_seq_list[block_start:block_end])
                if suspicious_chunk.startswith("-"):
                    # Already in a gap, skip
                    continue

                found_something = True
                if debug:
                    print(
                        f"  Possible mispositioned block (length={len(suspicious_chunk)}): '{suspicious_chunk}'"
                    )

                # Count how many dashes in that gap going further left
                gap_length = 0
                gap_pos = preceding_char_index
                while gap_pos >= 0 and temp_seq_list[gap_pos] == "-":
                    gap_length += 1
                    gap_pos -= 1

                # Show up to 10 chars before the gap for context
                context_start = max(0, gap_pos - 9)
                context_end = gap_pos + 1
                alignment_context = "".join(temp_seq_list[context_start:context_end])

                ref_context = ref_seq[context_start:context_end]
                if debug:
                    print(
                        f"Gap found at alignment pos {gap_pos + 1} (length={gap_length} dashes)"
                    )
                    print(
                        f"{'Alignment context before gap (temp_seq):':<{label_width}} {alignment_context}"
                    )
                    print(f"{'Corresponding WT context:':<{label_width}} {ref_context}")

                # Instead of placing chunk AFTER the entire gap,
                # we place it at the START of the gap:
                ref_gap_start = gap_pos + 1
                ref_gap_end = ref_gap_start + len(suspicious_chunk)

                # Confirm it doesn't go out of bounds
                if ref_gap_end <= len(ref_seq):
                    ref_gap_chunk = ref_seq[ref_gap_start:ref_gap_end]
                    if suspicious_chunk == ref_gap_chunk:
                        if debug:
                            print(
                                f"  >>> The suspicious chunk '{suspicious_chunk}' "
                                f"matches the ref_seq chunk '{ref_gap_chunk}' right after the gap!"
                            )
                            print(
                                "      This strongly suggests a misalignment. Moving chunk to the gap.\n"
                            )

                        # 1) Remove from old place
                        for i in range(block_start, block_end):
                            temp_seq_list[i] = "-"

                        # 2) Place it at the start of the gap
                        if ref_gap_end <= len(temp_seq_list):
                            for i, ch in enumerate(suspicious_chunk):
                                temp_seq_list[ref_gap_start + i] = ch
                        # else:
                        #     print(
                        #         "    [Warning] Not enough space to move chunk in temp_seq.\n"
                        #     )

                # Print alignment after fix attempt
                updated_temp_seq = "".join(temp_seq_list)
                if debug:
                    print("  Full alignment after fix attempt:")
                    print("  ref_seq: ", ref_seq)
                    print("  temp_seq:", updated_temp_seq)

        if not found_something:
            if debug:
                print("  No suspicious gap found nearby.\n")

    fixed_temp_seq = "".join(temp_seq_list)
    return fixed_temp_seq


def consolidate_structural_islands(pdb_seq, temp_seq, break_indexes):
    """
    Detects and fixes PDB residues that are scattered across a gap in the
    alignment even though they are physically bonded (no chain break between
    them) and therefore must occupy adjacent alignment columns.

    Unlike `detect_alignment_mistakes_and_reposition`, which looks for
    sequence-level clues (a chunk matching the ref sequence near a gap),
    this uses chain-break data directly (see `find_chain_breaks`, which --
    when given `resnums` -- flags a break on either an abnormal CA-CA
    distance or a non-consecutive author residue number, so a genuine
    construct-level deletion is never mistaken for physical contiguity just
    because its flanking residues happen to sit at a normal bond distance):
    any run of `pdb_seq` residues not separated by a break is a "structural
    island" that is physically contiguous in 3D, so it must map to
    contiguous alignment columns. A local aligner can still coincidentally
    match a few of an island's residues to unrelated WT positions purely by
    letter similarity (e.g. inside a disordered loop where only a handful of
    residues are resolved) -- this splits one physical island across
    several isolated, gap-separated alignment columns, which is a
    structural impossibility rather than a legitimate deletion.

    An island can straddle a mix of already-correct and scattered residues
    (e.g. a short stray cluster directly bonded to the start of an otherwise
    long, already well-aligned helix). Only the single longest
    already-contiguous run within the island is trusted as a fixed anchor;
    every other run -- regardless of its own length -- is treated as
    displaced and relocated as one block, snapped immediately next to its
    neighboring anchor. This is always safe: because an island is defined
    by the absence of a real chain break between its members (see
    `find_chain_breaks`), any split among its members can only be an
    alignment artifact, never a legitimate deletion, so islands are always
    fully consolidated into one contiguous run.

    Parameters
    ----------
    pdb_seq : str
        The raw (ungapped) PDB sequence, in chain order.
    temp_seq : str
        The current gapped PDB alignment string (e.g. the output of
        `detect_alignment_mistakes_and_reposition`).
    break_indexes : list of int
        Indexes into `pdb_seq` flagged as chain breaks (e.g. by
        `find_chain_breaks`; index i means the break falls between residue
        i-1 and i).

    Returns
    -------
    str
        The consolidated alignment string for the PDB sequence.
    """
    n = len(pdb_seq)
    if n == 0:
        return temp_seq

    breaks = set(break_indexes)
    temp_list = list(temp_seq)

    # 1) Partition pdb_seq into islands: runs with no real chain break between members.
    islands = []
    current = [0]
    for i in range(1, n):
        if i in breaks:
            islands.append(current)
            current = [i]
        else:
            current.append(i)
    islands.append(current)

    # 2) Recompute raw-pdb-index -> alignment-column map fresh from temp_seq. Any
    # prior repair pass may have moved characters without updating a pdb_map, so we
    # can't trust one passed in -- this is always correct given the current string.
    fresh_map = {}
    pdb_idx = 0
    for aln_idx, ch in enumerate(temp_list):
        if ch != "-":
            fresh_map[pdb_idx] = aln_idx
            pdb_idx += 1

    for island in islands:
        if len(island) < 2 or not all(idx in fresh_map for idx in island):
            continue

        positions = [fresh_map[idx] for idx in island]

        # Find every maximal already-contiguous run within the island.
        runs = []  # (start index within island/positions, length)
        run_start = 0
        for k in range(1, len(positions)):
            if positions[k] != positions[k - 1] + 1:
                runs.append((run_start, k - run_start))
                run_start = k
        runs.append((run_start, len(positions) - run_start))

        # Trust only the single longest already-contiguous run as the anchor
        # (ties go to the earliest-occurring run in chain order); every other
        # run in the island, regardless of its own length, is relocated next
        # to it below.
        anchors = [max(runs, key=lambda r: r[1])]

        anchor_set = set()
        for start, length in anchors:
            anchor_set.update(range(start, start + length))

        # Walk the island, relocating each maximal run of non-anchor ("loose")
        # members as one block, snapped next to its neighboring anchor.
        i = 0
        while i < len(island):
            if i in anchor_set:
                i += 1
                continue
            j = i
            while j < len(island) and j not in anchor_set:
                j += 1

            loose_island_idxs = island[i:j]
            loose_positions = positions[i:j]
            chars = [pdb_seq[idx] for idx in loose_island_idxs]

            prev_anchor_end_pos = positions[i - 1] if i > 0 else None
            next_anchor_start_pos = positions[j] if j < len(island) else None
            if prev_anchor_end_pos is not None:
                new_start = prev_anchor_end_pos + 1
            elif next_anchor_start_pos is not None:
                new_start = next_anchor_start_pos - len(chars)
            else:
                new_start = loose_positions[0]
            new_end = new_start + len(chars) - 1

            # Bail out (leave this loose run untouched) rather than clobber
            # unrelated data if there isn't clean room to relocate it.
            if new_start < 0 or new_end >= len(temp_list):
                i = j
                continue
            old_positions = set(loose_positions)
            destination_clear = all(
                p in old_positions or temp_list[p] == "-"
                for p in range(new_start, new_end + 1)
            )
            if not destination_clear:
                i = j
                continue

            for p in loose_positions:
                temp_list[p] = "-"
            for k, ch in enumerate(chars):
                temp_list[new_start + k] = ch

            i = j

    return "".join(temp_list)
