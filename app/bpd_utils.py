from matplotlib import pyplot as plt
import numpy as np
import streamlit as st
from Bio.Align import PairwiseAligner

def filter_signal_indices(signal_indices, minimum_spacing=50):
    """
    Given a sorted list of signal indices, remove any indices
    that fall within 'minimum_spacing' of the previously accepted index.

    Parameters:
    -----------
    signal_indices: list of int
        A list of signal indices (need not be sorted, but will be sorted here).
    minimum_spacing: int
        The minimum number of bases that must separate consecutive signal indices.

    Returns:
    --------
    filtered: list of int
        A filtered list of signal indices such that no two indices
        are within 'minimum_spacing' bases of each other.
    """
    # First, ensure the signal_indices are sorted
    sorted_indices = sorted(signal_indices)

    filtered = []
    last_accepted = None

    for idx in sorted_indices:
        # If we haven't accepted any index yet, accept the first one
        if last_accepted is None:
            filtered.append(idx)
            last_accepted = idx
        else:
            # Check distance from the last accepted index
            if idx - last_accepted >= minimum_spacing:
                filtered.append(idx)
                last_accepted = idx

    return filtered

def signal_index_to_seq_index(query_to_signal, signal_index):
    """
    Given a mapping from bases to signal indices (query_to_signal)
    and a particular signal_index, return the corresponding
    base index in the sequence.

    :param query_to_signal: A sorted array-like that indicates, for each base,
                            the earliest signal index covering that base.
                            (Typically something like: [0, 5, 10, 15, ...])
    :param signal_index:    An integer representing a position within the signal.
    :return:                The 0-based index into seq that covers signal_index.
    """

    # Use searchsorted to find where 'signal_index' would be inserted
    # so that query_to_signal is kept sorted. Then subtract 1 to get
    # the actual base index that covers this signal index.
    base_index = np.searchsorted(query_to_signal, signal_index, side='right') - 1

    return base_index

def match_breakpoints(algorithm, method, tolerance=5):
    matched_algorithm = []
    matched_method = []
    false_positives = []
    false_negatives = []

    # Create copies to avoid modifying the original arrays
    algorithm_copy = list(algorithm)
    method_copy = list(method)

    for point in algorithm:
        # Find the closest point in the method within the tolerance
        close_points = [m for m in method_copy if abs(m - point) <= tolerance]
        if close_points:
            closest_point = min(close_points, key=lambda x: abs(x - point))
            matched_algorithm.append(point)
            matched_method.append(closest_point)
            method_copy.remove(closest_point)
        else:
            false_negatives.append(point)

    # Remaining points in method_copy are false positives
    false_positives = method_copy

    return matched_algorithm, matched_method, false_positives, false_negatives

def plot_breakpoints_with_labels(matched_bonito, matched_method, false_positives, false_negatives, method_name):
    # Create the plot
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Define y-axis positions with more spacing
    y_positions = {
        'Matched (Original)': 1.5,
        'Matched (Method)': 2.0,
        'False Positives': 2.5,
        'False Negatives': 1.0
    }
    
    # Scatter plots for matched breakpoints and false positives/negatives
    ax.scatter(matched_bonito, [y_positions['Matched (Original)']]*len(matched_bonito), color='green', label=f'Matched {method_name} (Original)')
    ax.scatter(matched_method, [y_positions['Matched (Method)']]*len(matched_method), color='orange', label=f'Matched {method_name} (Method)')
    ax.scatter(false_positives, [y_positions['False Positives']]*len(false_positives), color='red', label='False Positives')
    ax.scatter(false_negatives, [y_positions['False Negatives']]*len(false_negatives), color='purple', label='False Negatives')

    # Draw lines between matched breakpoints
    for a, m in zip(matched_bonito, matched_method):
        ax.plot([a, m], [y_positions['Matched (Original)'], y_positions['Matched (Method)']], color='black')

    # Add text labels for false positives and false negatives
    for fp in false_positives:
        ax.text(fp, y_positions['False Positives'], str(fp), fontsize=9, ha='right')
    for fn in false_negatives:
        ax.text(fn, y_positions['False Negatives'], str(fn), fontsize=9, ha='right')

    # Customize the plot
    ax.set_yticks(list(y_positions.values()))
    ax.set_yticklabels(list(y_positions.keys()))
    ax.set_xlabel('signal index')
    # ax.set_ylabel('Type')
    ax.legend()
    ax.set_title(f'Breakpoints Comparison - {method_name}')

    # Display the plot in Streamlit
    st.pyplot(fig)

def build_alignment_strings_with_index_maps(seq, ref_seq, alignment):
    """
    Given:
      seq       : the 'target' sequence used in aligner.align(seq, ref_seq)
      ref_seq   : the 'query' sequence used in aligner.align(seq, ref_seq)
      alignment : a single Alignment object from the new PairwiseAligner (e.g. alignments[0])

    Returns:
      (aligned_seq_str, aligned_ref_str, seq_index_map, ref_index_map)

      where aligned_seq_str and aligned_ref_str are the same length,
      containing letters (A/C/G/T/etc.) or '-' for gaps.

      seq_index_map[i] == j means aligned_seq_str[i] came from seq[j]
        (or None if aligned_seq_str[i] == '-')
      ref_index_map[i] == k means aligned_ref_str[i] came from ref_seq[k]
        (or None if aligned_ref_str[i] == '-')
    """
    # alignment.aligned is a list of coordinate-block pairs for (target, query)
    # Example: [((0, 4), (0, 4)), ((5, 8), (4, 7)), ...]
    # Each block says that seq[t_start:t_end] aligns with ref_seq[r_start:r_end].
    # Between blocks, there may be unaligned (gap) regions.
    blocks = alignment.aligned

    # blocks[0] -> target blocks (shape: (num_blocks, 2))
    # blocks[1] -> query blocks  (shape: (num_blocks, 2))
    target_blocks, query_blocks = blocks

    # print(alignment.aligned)
    aligned_seq_fragments = []
    aligned_ref_fragments = []
    seq_index_map = []
    ref_index_map = []

    # We'll track the previous block end to detect gaps
    prev_t_end = 0
    prev_r_end = 0

    for (t_start, t_end), (r_start, r_end) in zip(target_blocks, query_blocks):
        
        # 1) If there's a gap between prev_t_end and t_start in 'seq'
        gap_in_target = t_start - prev_t_end
        # 2) If there's a gap between prev_r_end and r_start in 'ref_seq'
        gap_in_query  = r_start - prev_r_end

        # Handle gap in 'seq' (extra bases in 'ref_seq' or just unaligned)
        if gap_in_target > 0:
            # The portion in seq that isn't in the aligned block
            unaligned_seq = seq[prev_t_end : t_start]
            aligned_seq_fragments.append(unaligned_seq)
            # We must put '-' in the ref for these unaligned seq bases? Actually the reverse:
            # If the target advanced but the query didn't, that means there's a gap in the REF.
            # So let's see if gap_in_target > 0 but gap_in_query == 0 => we have an insertion in 'seq' or unaligned region?
            # However, in local alignment, it's possible these bases simply aren't aligned at all.
            # We'll treat them as "gap in ref".
            aligned_ref_fragments.append('-' * gap_in_target)
            # Build index maps
            for j in range(gap_in_target):
                seq_index_map.append(prev_t_end + j)
                ref_index_map.append(None)

        # Handle gap in 'ref_seq'
        if gap_in_query > 0:
            unaligned_ref = ref_seq[prev_r_end : r_start]
            aligned_ref_fragments.append(unaligned_ref)
            aligned_seq_fragments.append('-' * gap_in_query)
            # Build index maps
            for j in range(gap_in_query):
                seq_index_map.append(None)
                ref_index_map.append(prev_r_end + j)

        # Now the actual aligned block
        block_seq = seq[t_start:t_end]       # substring in seq
        block_ref = ref_seq[r_start:r_end]   # substring in ref_seq

        # Usually these blocks are the same length, but let's handle the general case:
        len_s = len(block_seq)
        len_r = len(block_ref)

        if len_s == len_r:
            aligned_seq_fragments.append(block_seq)
            aligned_ref_fragments.append(block_ref)
            for j in range(len_s):
                seq_index_map.append(t_start + j)
                ref_index_map.append(r_start + j)
        elif len_s > len_r:
            # More bases in seq => pad ref
            aligned_seq_fragments.append(block_seq)
            aligned_ref_fragments.append(block_ref + '-'*(len_s - len_r))
            for j in range(len_s):
                seq_index_map.append(t_start + j)
                if j < len_r:
                    ref_index_map.append(r_start + j)
                else:
                    ref_index_map.append(None)
        else:
            # More bases in ref => pad seq
            aligned_seq_fragments.append(block_seq + '-'*(len_r - len_s))
            aligned_ref_fragments.append(block_ref)
            for j in range(len_r):
                if j < len_s:
                    seq_index_map.append(t_start + j)
                else:
                    seq_index_map.append(None)
                ref_index_map.append(r_start + j)

        prev_t_end = t_end
        prev_r_end = r_end

    # If there's leftover in seq beyond the last block
    if prev_t_end < len(seq):
        leftover_seq = seq[prev_t_end:]
        aligned_seq_fragments.append(leftover_seq)
        aligned_ref_fragments.append('-' * len(leftover_seq))
        for j in range(len(leftover_seq)):
            seq_index_map.append(prev_t_end + j)
            ref_index_map.append(None)

    # If there's leftover in ref_seq beyond the last block
    if prev_r_end < len(ref_seq):
        leftover_ref = ref_seq[prev_r_end:]
        aligned_ref_fragments.append(leftover_ref)
        aligned_seq_fragments.append('-' * len(leftover_ref))
        for j in range(len(leftover_ref)):
            seq_index_map.append(None)
            ref_index_map.append(prev_r_end + j)

    # Combine all fragments
    aligned_seq_str = ''.join(aligned_seq_fragments)
    aligned_ref_str = ''.join(aligned_ref_fragments)

    return aligned_seq_str, aligned_ref_str, seq_index_map, ref_index_map


def mismatch_regions_with_pairwisealigner(
    seq,
    ref_seq,
    query_to_signal,
    mode="local",
    treat_unaligned_ends_as_mismatch=False
):
    """
    1) Align seq (target) and ref_seq (query) using Bio.Align.PairwiseAligner.
    2) Build fully gapped alignment strings + index maps for seq and ref_seq.
    3) Identify mismatches in 'seq' (including places where ref has '-').
    4) Return a list of (signal_start, signal_end) intervals for mismatches.

    :param seq:  The primary sequence for which query_to_signal is defined (length N).
    :param ref_seq: The other sequence. Possibly different length from seq.
    :param query_to_signal: array of length N+1, so base i in seq => [q2s[i], q2s[i+1]).
    :param mode: "global" or "local" alignment.
    :param treat_unaligned_ends_as_mismatch:
        If True (and using local alignment), any portion of seq not aligned
        to ref_seq is considered mismatch (e.g. leading/trailing unaligned part).
    :return: list of (signal_start, signal_end) intervals for each contiguous mismatch region
    """
    if len(query_to_signal) != len(seq) + 1:
        raise ValueError(
            "query_to_signal must be length len(seq)+1. "
            f"Got len(seq)={len(seq)}, len(query_to_signal)={len(query_to_signal)}."
        )

    aligner = PairwiseAligner()
    aligner.mode = mode
    # Example scoring (tweak as needed):
    aligner.match_score = 2
    aligner.mismatch_score = -2
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -1

    # Perform alignment
    # Note: 'seq' is the 'target', 'ref_seq' is the 'query' in PairwiseAligner terminology
    all_alignments = aligner.align(seq, ref_seq)
    best_aln = all_alignments[0]
    # print("score:"+ str(best_aln.score))

    # Build fully gapped strings
    aligned_seq_str, aligned_ref_str, seq_index_map, ref_index_map = \
        build_alignment_strings_with_index_maps(seq, ref_seq, best_aln)

    # Now detect mismatches in aligned_seq_str
    # We'll track contiguous mismatch blocks in terms of original seq indices
    mismatch_regions = []
    mismatch_start = None

    for i in range(len(aligned_seq_str)):
        q_char = aligned_seq_str[i]   # from seq
        r_char = aligned_ref_str[i]   # from ref_seq
        q_idx  = seq_index_map[i]     # None if gap in seq
        # ref_index_map[i] is not crucial for marking mismatch in seq

        if q_char == '-':
            # Gap in seq => no actual base in seq to mismatch
            # If we were in a mismatch region, close it
            if mismatch_start is not None:
                mismatch_regions.append((mismatch_start, q_idx - 1))
                mismatch_start = None
            continue

        if r_char == '-':
            # Gap in ref => definitely a mismatch in seq at q_idx
            if mismatch_start is None:
                mismatch_start = q_idx
            # Continue, do not close yet, might be contiguous
        else:
            # Real base vs base
            if q_char != r_char:
                # mismatch
                if mismatch_start is None:
                    mismatch_start = q_idx
            else:
                # match
                if mismatch_start is not None:
                    mismatch_regions.append((mismatch_start, q_idx - 1))
                    mismatch_start = None

    # If we ended in a mismatch region, close it
    if mismatch_start is not None:
        mismatch_regions.append((mismatch_start, len(seq) - 1))

    # Optionally, treat unaligned ends (outside best scoring region) as mismatch
    if mode == "local" and treat_unaligned_ends_as_mismatch:
        # best_aln.aligned might look like [((t0_start, t0_end), (r0_start, r0_end)), ...]
        # We'll see how much of seq is unaligned at the start or end
        blocks = best_aln.aligned[0]  # first block
        t_start_aln = blocks[0][0]    # alignment start in seq
        blocks_last = best_aln.aligned[-1]
        t_end_aln   = blocks_last[0][1]  # alignment end in seq (exclusive)

        # front unaligned region
        if t_start_aln > 0:
            mismatch_regions.append((0, t_start_aln - 1))
        # back unaligned region
        if t_end_aln < len(seq):
            mismatch_regions.append((t_end_aln, len(seq) - 1))

        # Sort & merge any overlap
        mismatch_regions.sort()
        merged = []
        for region in mismatch_regions:
            if not merged or region[0] > merged[-1][1] + 1:
                merged.append(region)
            else:
                # overlap or adjacent => extend
                merged[-1] = (merged[-1][0], max(merged[-1][1], region[1]))
        mismatch_regions = merged

    # Convert mismatch indices [start, end] to [signal_start, signal_end)
    mismatch_signal_regions = []
    for (start_idx, end_idx) in mismatch_regions:
        signal_start = query_to_signal[start_idx]
        signal_end   = query_to_signal[end_idx + 1]
        mismatch_signal_regions.append((signal_start, signal_end))

    return mismatch_signal_regions


def filter_indices_by_mismatch_regions(to_edit_signal_indices, mismatched_regions):
    # Flatten mismatched_regions for easier comparison
    mismatched_intervals = np.array(mismatched_regions)
    
    # Define a helper function to check if an index is in any interval
    def is_in_mismatched_regions(index):
        # Use binary search to efficiently determine if index is within any region
        low = 0
        high = len(mismatched_intervals) - 1
        
        while low <= high:
            mid = (low + high) // 2
            start, end = mismatched_intervals[mid]
            if start <= index <= end:
                return True
            elif index < start:
                high = mid - 1
            else:
                low = mid + 1
        return False
    
    # Filter indices based on the mismatched regions
    return [index for index in to_edit_signal_indices if is_in_mismatched_regions(index)]