
import numpy as np
import base64
from pathlib import Path
def diff_str(str1, str2):
    """
    Compare two strings and return all differences as a list of tuples.
    Each tuple contains:
      - Index of the difference
      - Character in str1 (None if it doesn't exist)
      - Character in str2 (None if it doesn't exist)
    """
    differences = []
    max_len = max(len(str1), len(str2))

    for i in range(max_len):
        char1 = str1[i] if i < len(str1) else None
        char2 = str2[i] if i < len(str2) else None
        if char1 != char2:
            differences.append((i, char1, char2))

    return differences



def find_floor_index(arr, num):
    # Find the largest element in the array that is smaller than or equal to num
    floor_val = arr[arr <= num].max()
    # Get the index of that value
    return np.where(arr == floor_val)[0][0]

def clean_list(input_list):
    """
    Convert all NumPy integers (e.g., np.int64) in a list to Python int.
    Other data types in the list remain unchanged.

    Parameters
    ----------
    input_list : list
        A list that may contain NumPy integer types and/or other data.

    Returns
    -------
    list
        A new list where all NumPy integers have been converted to Python ints.
    """
    cleaned = []
    for item in input_list:
        # Check if item is a NumPy integer
        if isinstance(item, np.integer):
            cleaned.append(int(item))  # convert to Python int
        else:
            cleaned.append(item)
    return cleaned



def create_download_link(file_path: str, link_text: str = None) -> str:
    """
    Create a clickable download link using an HTML anchor tag with base64-encoded data.
    :param file_path: Local path to the file you want to enable for download.
    :param link_text: The text displayed in the link (defaults to the file name if not provided).
    :return: A string containing an HTML anchor tag.
    """
    file_path = Path(file_path)
    link_text = link_text or file_path.name

    # Read and encode the file
    with open(file_path, 'rb') as f:
        file_data = f.read()
    b64_data = base64.b64encode(file_data).decode()

    # You can customize the MIME type depending on your file.
    # For text-based files, "data:text/plain" might be more appropriate, or "data:application/octet-stream" for generic.
    return f'<a href="data:application/octet-stream;base64,{b64_data}" download="{file_path.name}">{link_text}</a>'
