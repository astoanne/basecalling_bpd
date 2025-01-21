import numpy as np
import pysam
from Bio import pairwise2
import streamlit as st
def get_local_sequence(seq, query_to_signal, signal_index, radius=5):
    """
    Given a DNA sequence, a query_to_signal array (base-to-signal mapping),
    a signal index, and a radius (± number of bases to include around the found base),
    return the local substring of seq surrounding the found base.
    """
    base_index = np.searchsorted(query_to_signal, signal_index, side='right') - 1
    
    # Guard for boundaries
    if base_index < 0:
        base_index = 0
    elif base_index >= len(seq):
        base_index = len(seq) - 1

    start = max(0, base_index - radius)
    end = min(len(seq), base_index + radius + 1)
    local_sequence = seq[start:end]
    return local_sequence

def align_local_sequences(seqA, seqB):
    """
    Perform a local alignment between seqA and seqB using custom scoring.
    Returns (aligned_seqA, aligned_seqB, score, start, end) or None if no alignment found.
    """
    # Example scoring: match=2, mismatch=-2, gap opening=-5, gap extension=-1
    alignments = pairwise2.align.localms(seqA, seqB, 2, -2, -5, -1)
    #alignments = pairwise2.align.localms(seqA, seqB, 1, -1, -1, -1)
    if not alignments:
        return None

    best = alignments[0]
    aligned_seqA, aligned_seqB, score, start, end = best
    return (aligned_seqA, aligned_seqB, score, start, end)

def get_center_index_in_aligned(aligned_seq, original_center=5):
    """
    Given an aligned sequence (with possible gaps), find the position
    in the aligned sequence corresponding to the 'original_center'-th real base.
    """
    real_bases_count = 0
    for i, base in enumerate(aligned_seq):
        if base != '-':
            if real_bases_count == original_center:
                return i
            real_bases_count += 1
    return None

def analyze_center_5(aligned_seq1, aligned_seq2, original_center=5):
    """
    Extract ±2 bases around the 'original_center' in both aligned sequences,
    strip gaps, and compare. If they differ, return (old_5, new_5).
    Otherwise, return None.
    """
    c1 = get_center_index_in_aligned(aligned_seq1, original_center)
    c2 = get_center_index_in_aligned(aligned_seq2, original_center)
    if c1 is None or c2 is None:
        return None

    # Define a 5-base window around the center (±2)
    window_1 = aligned_seq1[max(0, c1 - 2): c1 + 3]
    window_2 = aligned_seq2[max(0, c2 - 2): c2 + 3]

    # Strip gaps
    stripped_1 = window_1.replace('-', '')
    stripped_2 = window_2.replace('-', '')

    if stripped_1 == stripped_2:
        return None  # no edit needed
    else:
        return (stripped_1, stripped_2)

def process_signal_index(seq1, seq2, query_to_signal1, query_to_signal2, signal_index, radius=5):
    """
    Extract local sequences from seq1 and seq2 using the same signal_index,
    perform local alignment, and return alignment data plus the base_index
    in seq1 that corresponds to this signal_index.
    """
    # 1) Extract local windows
    local_seq_1 = get_local_sequence(seq1, query_to_signal1, signal_index, radius)
    local_seq_2 = get_local_sequence(seq2, query_to_signal2, signal_index, radius)

    # 2) Do the local alignment
    alignment = align_local_sequences(local_seq_1, local_seq_2)
    if alignment is None:
        return None  # No alignment

    aligned_seq1, aligned_seq2, score, local_start, local_end = alignment

    # 3) (Optional) analyze the center 5 bases difference
    local_diff = analyze_center_5(aligned_seq1, aligned_seq2, original_center=radius)
    
    # 4) Map the alignment window back to seq1
    #    np.searchsorted(...) finds which base in seq1 aligns to signal_index
    base_index = np.searchsorted(query_to_signal1, signal_index, side='right') - 1

    return {
        "aligned_seq1": aligned_seq1,
        "aligned_seq2": aligned_seq2,
        "score": score,
        "local_start": local_start,  # alignment start in local_seq_1 coords
        "local_end": local_end,      # alignment end (exclusive) in local_seq_1 coords
        "base_index": base_index,
        "local_diff": local_diff
    }


def batch_edit_bam(
    input_path,
    edited_path,
    read_id,
    seq1,
    seq2,
    query_to_signal1,
    query_to_signal2,
    signal_indices,
    radius=5
):
    """
    1. For each signal_index in signal_indices, get alignment data from seq1 vs seq2.
       Use the alignment's start/end to plan an edit.
    2. Sort all edits in descending order by the 'start_index'.
    3. Apply all edits on seq1. Update the read's query_sequence with the new string in the BAM.
    """

    MIN_ALIGNMENT_SCORE = 0
    planned_edits = []

    for s_idx in signal_indices:
        alignment_data = process_signal_index(
            seq1=seq1,
            seq2=seq2,
            query_to_signal1=query_to_signal1,
            query_to_signal2=query_to_signal2,
            signal_index=s_idx,
            radius=radius
        )

        # If no alignment returned, skip
        if alignment_data is None:
            continue

        aligned_seq1 = alignment_data["aligned_seq1"]
        aligned_seq2 = alignment_data["aligned_seq2"]
        score        = alignment_data["score"]
        local_start  = alignment_data["local_start"]
        local_end    = alignment_data["local_end"]  # typically 'end' is exclusive
        base_index   = alignment_data["base_index"]
        local_diff   = alignment_data["local_diff"]

        # Skip if alignment is too weak
        if score < MIN_ALIGNMENT_SCORE:
            continue

        # If you specifically only want to edit if there's a difference in the 5-base center:
        if local_diff is None:
            continue

        # Convert the alignment's local_start/local_end to coordinates in the full seq1.
        # The local_seq_1 is extracted from: (base_index - radius) .. (base_index + radius)
        global_offset = max(0, base_index - radius)
        # local_start/local_end refer to slices in local_seq_1
        global_start = global_offset + local_start
        # 'local_end' is exclusive, so the last included base is local_end - 1
        global_end = global_offset + local_end - 1
        
        # If the alignment was from local_seq_1[2:7], that means
        # in seq1, the region is [global_start : global_end].

        # The old substring in seq1 that is being replaced:
        old_substring = seq1[global_start : global_end + 1]

        # The new substring is the aligned region from seq2, minus any gaps:
        new_substring = aligned_seq2[local_start : local_end].replace('-', '')

        # For reference, local_diff might still be "the center 5" difference,
        # but here we’re applying the entire alignment region as an edit.
        # If you truly only want to replace the center 5 bases, you can adapt the logic.

        # Basic boundary checks
        if global_start < 0 or global_end >= len(seq1) or global_end < global_start:
            continue

        planned_edits.append(
            (global_start, global_end, old_substring, new_substring, score)
        )

    # --- Now apply edits in descending order of start_index ---
    planned_edits.sort(key=lambda x: x[0], reverse=True)

    seq1_list = list(seq1)  # Convert seq1 to mutable list
    # for (start_index, end_index, old_sub, new_sub, alignment_score) in planned_edits:
    #     # Optional: check the expected old_substring
    #     # actual_old = ''.join(seq1_list[start_index:end_index+1])
    #     # if actual_old != old_sub:
    #     #     print(f"Warning: mismatch. Expected {old_sub}, found {actual_old}")

    #     # Do the replacement (handles expansions or contractions)
    #     seq1_list[start_index:end_index+1] = list(new_sub)

    #     # print(
    #     #     f"Proposed edit (score={alignment_score:.2f}): "
    #     #     f"{old_sub} -> {new_sub} at base_range={start_index}-{end_index}"
    #     # )
    #  # Initialize a list to store all edits
    edits_content = []

    # Loop through planned edits and build the content
    for (start_index, end_index, old_sub, new_sub, alignment_score) in planned_edits:
        seq1_list[start_index:end_index+1] = list(new_sub)
        edits_content.append(
            f"**Proposed edit** (score={alignment_score:.2f}):\n"
            f"`{old_sub}` -> `{new_sub}`\n"
            f"\n**Base Range**: `{start_index}-{end_index}`"
        )

    # Display all edits inside one expander
    with st.expander("Proposed Edits"):
        for edit in edits_content:
            st.markdown(edit)   

    # Create the final edited sequence
    new_seq1 = ''.join(seq1_list)
    

    # Write updated read to BAM
    with pysam.AlignmentFile(input_path, "rb", check_sq=False) as in_bam, \
         pysam.AlignmentFile(edited_path, "wb", header=in_bam.header) as out_bam:

        for read in in_bam:
            if read.query_name == read_id:
                read.query_sequence = new_seq1
            out_bam.write(read)

    return new_seq1
