import os
import sys
import pysam
import re
if not os.path.exists("minimap2"):
    os.system("git clone https://github.com/lh3/minimap2")
    os.system(f"{sys.executable} -m pip install -U setuptools wheel")  # Ensure build tools are installed
    os.system("cd minimap2 && python setup.py install")

# Now import the mappy module
try:
    import mappy as mp
    print("Successfully imported mappy!")
except ImportError:
    print("Failed to import mappy. Check if minimap2 was installed properly.")
def extract_md_tag(io_read):
    """
    Extract the MD tag value (excluding the 'MD:Z:' prefix) from the alignment data.

    Args:
        io_read: An object or dictionary containing the alignment data,
                          where `full_align` holds the SAM alignment string.

    Returns:
        str: The MD tag value if found (excluding 'MD:Z:'), otherwise a message indicating it is not found.
    """
    alignment_data = io_read.full_align  # Assuming 'full_align' contains the SAM alignment string
    tags = alignment_data.get('tags', [])
    md_tag = next((tag for tag in tags if tag.startswith("MD:Z:")), None)

    if md_tag:
        # Remove the 'MD:Z:' prefix before returning
        return md_tag[5:]
    else:
        return "MD tag not found."

def extract_read_info(bam_path, read_id):
    """
    Extracts the sequence for a specific read in the BAM file.

    Parameters:
        bam_path (str): Path to the BAM file.
        read_id (str): Query name of the read to find.

    Returns:
        str: Sequence of the read, or None if not found.
    """
    with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as bam_fh:
        for read in bam_fh:
            if read.query_name == read_id:
                return read.query_sequence
    return None

def get_cigar_and_md_tag(bam_path, read_id):
    """
    Extract the CIGAR string and MD tag for a specific read from a BAM file.

    Args:
        bam_path (str): Path to the BAM file.
        read_id (str): Query name (read ID) for which to extract the CIGAR and MD tag.

    Returns:
        tuple: A tuple containing the 'cigar' string and 'md_tag'.
               If the read or MD tag is not found, the respective values will be None.
    """
    cigar = None
    md_tag = None

    with pysam.AlignmentFile(bam_path, "rb") as bam_fh:
        for read in bam_fh:
            if read.query_name == read_id:
                # Extract the CIGAR string
                cigar = read.cigarstring

                # Extract the MD tag
                tags = dict(read.tags)
                md_tag = tags.get("MD", None)

                # Exit the loop once the read is found
                break

    return cigar, md_tag
def edit_bam(input_path, edited_path, read_id, start_index, end_index, replacement_sequence):
    with pysam.AlignmentFile(input_path, "rb", check_sq=False) as input_bam, \
         pysam.AlignmentFile(edited_path, "wb", header=input_bam.header) as output_bam:

        for read in input_bam:
            if read.query_name == read_id:
                ####sanity check###
                #absurd_replacement = 'K' * len(replacement_sequence)
                #################
                query_sequence = list(read.query_sequence)
                original_length = len(query_sequence)
                
                #print(f"Original sequence length: {original_length}")
                if 0 <= start_index <= end_index < original_length:
                    # Replace the sequence
                    original_seq = ''.join(query_sequence[start_index:end_index + 1])
                    #print(f"Original bases at {start_index}-{end_index}: {original_seq}")

                    ####sanity check###
                    #query_sequence[start_index:end_index + 1] = list(absurd_replacement)
                    query_sequence[start_index:end_index + 1] = list(replacement_sequence)
                    
                    new_seq = ''.join(query_sequence[start_index:start_index + len(replacement_sequence)])
                    #print(f"Modified bases at {start_index}-{start_index + len(replacement_sequence) - 1}: {new_seq}")
                    
                    # Update the read's query sequence
                    read.query_sequence = ''.join(query_sequence)
            
            output_bam.write(read)

            

def generate_fastq(sequence, output_file, read_id, quality_char="I"):
    """
    Generate a FASTQ file with the given sequence.

    Args:
        sequence (str): The nucleotide sequence (e.g., "ACGTTT").
        output_file (str): Path to the output .fq file.
        read_id (str): Identifier for the sequence 
        quality_char (str): Quality score character (default is "I", Phred score 40).
    """
    # Ensure the sequence is uppercase
    sequence = sequence.upper()
    
    # Generate quality scores (same length as the sequence)
    quality_scores = quality_char * len(sequence)

    # Write to the FASTQ file
    with open(output_file, "w") as fq:
        fq.write(f"@{read_id}\n")       # Sequence ID
        fq.write(f"{sequence}\n")      # Sequence
        fq.write("+\n")                # Separator
        fq.write(f"{quality_scores}\n")  # Quality scores

def generate_sam_file(reference_filepath, query_seq, sam_output_path, read_id="query1", quality_char="I"):
    # Initialize the aligner
    aligner = mp.Aligner(str(reference_filepath), preset="map-hifi")  # Preset for high-quality reads
    if not aligner:
        raise Exception("Failed to initialize the aligner.")

    # Open the SAM file for writing
    with open(sam_output_path, 'w') as sam_file:
        # Write the SAM header
        sam_file.write(f"@HD\tVN:1.6\tSO:unsorted\n")
        for ref_name in aligner.seq_names:  # seq_names is a list of reference names
            ref_seq = aligner.seq(ref_name)  # Retrieve the sequence for this reference
            ref_len = len(ref_seq)  # Get the length of the reference sequence
            sam_file.write(f"@SQ\tSN:{ref_name}\tLN:{ref_len}\n")

        # Perform alignment
        for aln in aligner.map(query_seq):
            if aln.is_primary:  # Check if it's the primary alignment
                # Generate the SAM format fields
                flag = 0  # Adjust flags based on read orientation, etc.
                rname = aln.ctg
                pos = aln.r_st + 1  # 1-based position
                mapq = aln.mapq
                cigar = aln.cigar_str
                rnext = "*"  # No reference for paired-end reads
                pnext = 0
                tlen = 0
                seq = query_seq
                qual = quality_char * len(query_seq)

                # Compute NM tag (number of mismatches)
                nm = aln.NM

                # Write the alignment line to the SAM file
                sam_file.write(f"{read_id}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}\tNM:i:{nm}\n")

    print(f"SAM file generated at: {sam_output_path}")