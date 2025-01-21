import math
import re
import streamlit as st
def calculate_accuracy_error_qscore_from_md_and_cigar(cigar: str, md_tag: str):
    """
    Parse the CIGAR string and MD tag, then compute:
      - Number of matches, mismatches, insertions, deletions
      - Alignment accuracy (matches / (matches + mismatches + insertions + deletions))
      - Error rate (1 - accuracy)
      - Phred-like Q-score = -10 * log10(error_rate)

    Returns a dictionary with these statistics.
    """

    # --------------------------------------------------------------------------
    # 1. Parse the CIGAR string to count:
    #    - total 'M' (which includes both real matches + mismatches),
    #    - total 'I' (insertions),
    #    - total 'D' (deletions).
    # --------------------------------------------------------------------------
    cigar_pattern = re.compile(r"(\d+)([MIDNSHP=XB])")
    total_cigar_m = 0  # sum of 'M' (and '='/'X' if present)
    total_cigar_i = 0  # sum of 'I'
    total_cigar_d = 0  # sum of 'D'

    for length_str, op in cigar_pattern.findall(cigar):
        length = int(length_str)
        if op in ('M', '=', 'X'):
            total_cigar_m += length
        elif op == 'I':
            total_cigar_i += length
        elif op == 'D':
            total_cigar_d += length
        # We ignore S, H, P, N, B, etc. for counting error stats

    # --------------------------------------------------------------------------
    # 2. Parse the MD tag to get:
    #    - total_md_matches: how many reference-matching bases
    #    - total_md_deletions: how many deleted bases (^AAA => 3 deletions)
    #    - total_md_mismatches: how many single-base mismatches
    #
    #   The MD tag is typically of the form:  MD:Z:<numbers and symbols>
    #   Where:
    #       - Runs of digits (e.g. '18') => that many matches
    #       - A single letter (e.g. 'A') => 1 mismatch
    #       - A caret plus letters (e.g. '^AC') => a deletion of those letters
    # --------------------------------------------------------------------------
    if md_tag.startswith("MD:Z:"):
        md_tag = md_tag[5:]  # remove 'MD:Z:' prefix if present

    # Find all tokens which can be:
    #   - a run of digits (\d+)
    #   - a caret plus letters (\^[A-Za-z]+), i.e. a deletion
    #   - a single mismatch letter ([A-Za-z])
    md_tokens = re.findall(r'(\d+|\^[A-Za-z]+|[A-Za-z])', md_tag)

    total_md_matches = 0
    total_md_deletions = 0
    total_md_mismatches = 0

    for token in md_tokens:
        if token.isdigit():
            # A run of matching bases
            total_md_matches += int(token)
        elif token.startswith('^'):
            # This indicates a deletion in the reference
            # e.g. ^AC => 2 deleted bases
            deleted_bases = token[1:]  # remove '^'
            total_md_deletions += len(deleted_bases)
        else:
            # A single mismatch (A, T, G, or C, etc.)
            total_md_mismatches += 1

    # --------------------------------------------------------------------------
    # 3. Combine MD and CIGAR info to get final stats:
    #
    #   - "matches" (true matches) comes from the total_md_matches
    #   - "mismatches" comes from total_md_mismatches
    #   - "insertions" comes from total_cigar_i
    #   - "deletions" comes from total_md_deletions
    #
    #  The sum of M in the CIGAR = matches + mismatches.
    #  So we expect total_cigar_m ≈ total_md_matches + total_md_mismatches.
    #
    #  Our denominator for "accuracy" will be:
    #     matches + mismatches + insertions + deletions
    #  i.e. we treat any mismatch or indel as an error.
    # --------------------------------------------------------------------------
    matches = total_md_matches
    mismatches = total_md_mismatches
    insertions = total_cigar_i
    deletions = total_md_deletions

    total_bases_aligned = matches + mismatches + insertions + deletions
    if total_bases_aligned == 0:
        # edge case: no aligned bases
        return {
            "accuracy": 0.0,
            "error_rate": 1.0,
            "qscore": 0.0,
            "matches": 0,
            "mismatches": 0,
            "insertions": 0,
            "deletions": 0,
        }

    accuracy = matches / total_bases_aligned
    error_rate = 1.0 - accuracy

    # Compute a simple Phred-like Q-score
    if error_rate > 0:
        qscore = -10.0 * math.log10(error_rate)
    else:
        # If error_rate == 0, we assign an arbitrarily large Q-score
        # (some tools cap Q-scores at 60 or similar).
        qscore = 60.0

    return {
        "accuracy": accuracy,
        "error_rate": error_rate,
        "qscore": qscore,
        "matches": matches,
        "mismatches": mismatches,
        "insertions": insertions,
        "deletions": deletions,
    }
def pretty_print_acc(acc_dict, acc_dict_new=None):
    if acc_dict_new is None:
        # Old behavior: Single dictionary display
        st.table([
            {"Metric": "Accuracy", "Value": acc_dict["accuracy"]},
            {"Metric": "Error Rate", "Value": acc_dict["error_rate"]},
            {"Metric": "Q-Score", "Value": acc_dict["qscore"]},
            {"Metric": "Matches", "Value": str(int(acc_dict["matches"]))},  # Convert to string
            {"Metric": "Mismatches", "Value": str(int(acc_dict["mismatches"]))},  # Convert to string
            {"Metric": "Insertions", "Value": str(int(acc_dict["insertions"]))},  # Convert to string
            {"Metric": "Deletions", "Value": str(int(acc_dict["deletions"]))},  # Convert to string
        ])
    else:
        # New behavior: Comparison display
        st.table([
            {"Metric": "Accuracy", "Original": acc_dict["accuracy"], "New": acc_dict_new["accuracy"]},
            {"Metric": "Error Rate", "Original": acc_dict["error_rate"], "New": acc_dict_new["error_rate"]},
            {"Metric": "Q-Score", "Original": acc_dict["qscore"], "New": acc_dict_new["qscore"]},
            {"Metric": "Matches", "Original": str(int(acc_dict["matches"])), "New": str(int(acc_dict_new["matches"]))},
            {"Metric": "Mismatches", "Original": str(int(acc_dict["mismatches"])), "New": str(int(acc_dict_new["mismatches"]))},
            {"Metric": "Insertions", "Original": str(int(acc_dict["insertions"])), "New": str(int(acc_dict_new["insertions"]))},
            {"Metric": "Deletions", "Original": str(int(acc_dict["deletions"])), "New": str(int(acc_dict_new["deletions"]))},
        ])



def compute_accuracy(sam_file):
    total_bases = 0
    total_matches = 0
    total_mismatches = 0
    total_insertions = 0
    total_deletions = 0

    with open(sam_file, 'r') as f:
        for line in f:
            if line.startswith('@'):  # Skip header lines
                continue

            fields = line.strip().split('\t')
            cigar = fields[5]  # CIGAR string
            optional_fields = fields[11:]  # Optional fields

            # Skip malformed or invalid CIGAR strings
            if cigar == '*' or not cigar:
                continue

            try:
                # Parse the CIGAR string
                matches = sum(int(x) for x in re.findall(r'(\d+)M', cigar))
                insertions = sum(int(x) for x in re.findall(r'(\d+)I', cigar))
                deletions = sum(int(x) for x in re.findall(r'(\d+)D', cigar))

                # Get mismatches from NM:i tag
                mismatches = 0
                for field in optional_fields:
                    if field.startswith('NM:i:'):
                        mismatches = int(field.split(':')[2]) - insertions - deletions
                        break

                # Update totals
                aligned_bases = matches + insertions + deletions + mismatches
                total_bases += aligned_bases
                total_matches += matches
                total_insertions += insertions
                total_deletions += deletions
                total_mismatches += mismatches

            except ValueError:
                print(f"Skipping malformed CIGAR string: {cigar}")
                continue

    # Compute metrics
    accuracy = (total_matches / total_bases) if total_bases > 0 else 0
    error_rate = 1 - accuracy  # Error rate is the complement of accuracy
    qscore = -10 * math.log10(error_rate) if accuracy < 1 else float('inf')  # Q-score

    return {
        "accuracy": accuracy,  # Full precision
        "error_rate": error_rate,  # Full precision
        "qscore": qscore,  # Full precision
        "matches": total_matches,
        "mismatches": total_mismatches,
        "insertions": total_insertions,
        "deletions": total_deletions,
    }
