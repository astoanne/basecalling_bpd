def diff_cigar(original_cigar, edited_cigar):
    """
    Computes the difference between two CIGAR strings.

    :param original_cigar: The original CIGAR string.
    :param edited_cigar: The edited CIGAR string.
    :return: A list of differences. Each difference is a tuple of the form:
             (original, edited), where `original` and `edited` are CIGAR
             elements from the original and edited strings, respectively.
             If one string is longer, unmatched elements will be included.
    """

    def _parse_cigar(cigar):
        """
        Parses a CIGAR string into a list of tuples (length, operation).
        For example: "10M5I2D" -> [(10, 'M'), (5, 'I'), (2, 'D')].

        :param cigar: A CIGAR string.
        :return: A list of tuples representing the parsed CIGAR.
        """
        return [(int(length), op) for length, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar)]

    # Parse both CIGAR strings
    parsed_original = _parse_cigar(original_cigar)
    parsed_edited = _parse_cigar(edited_cigar)

    # Compare the two lists element by element
    diff_result = []
    max_len = max(len(parsed_original), len(parsed_edited))

    for i in range(max_len):
        original_element = parsed_original[i] if i < len(parsed_original) else None
        edited_element = parsed_edited[i] if i < len(parsed_edited) else None

        if original_element != edited_element:
            diff_result.append((original_element, edited_element))

    return diff_result
def extract_cigar_strings(sam_file):
    """
    Extracts CIGAR strings from a SAM file.

    Args:
        sam_file (str): Path to the SAM file.

    Returns:
        list: A list of CIGAR strings extracted from the SAM file.
    """
    cigar_strings = []

    # Open the SAM file and extract CIGAR strings
    with open(sam_file, 'r') as f:
        for line in f:
            if line.startswith('@'):  # Skip header lines
                continue
            fields = line.strip().split('\t')
            if len(fields) > 5:  # Ensure the CIGAR field exists
                cigar = fields[5]  # CIGAR string is the 6th column (0-indexed)
                cigar_strings.append(cigar)

    return cigar_strings