from collections import Counter

import numpy as np

# Module imports
from primalscheme3.core.seq_functions import extend_ambiguous_base


def trucnate_msa(msa: np.ndarray, mapping_index: int = 0) -> np.ndarray:
    """
    Removes all the bases outside the reference genome
    """
    return msa[:, msa[mapping_index] != ""]


def create_mapping(
    msa: np.ndarray, mapping_index: int = 0
) -> tuple[np.ndarray, np.ndarray]:
    """
    This returns a tuple of two items:
        mapping_array: list[int | None]
        truncated_msa: np.ndarray

    mapping_array: Each position in the list corrasponds to the same index in the MSA, The value in the list is the position in the reference genome
    truncated_msa: The MSA with all the bases outside the reference genome removed
    """
    # As NP is modifited in place, returning is not nessisary but is done for clarity
    # Truncate the msa
    truncated_msa = trucnate_msa(msa, mapping_index)
    # Create the empty mapping array
    mapping_list = [None] * truncated_msa.shape[1]
    mapping_array = np.array(mapping_list)
    # Select the reference genome
    reference_genome = truncated_msa[mapping_index]
    # Iterate over the msa genome
    current_ref_index = 0
    for col_index in range(truncated_msa.shape[1]):
        # If the base is not a gap, assign the mapping
        if reference_genome[col_index] not in {"", "-"}:
            mapping_array[col_index] = current_ref_index
            # increase refence index
            current_ref_index += 1
    return (mapping_array, truncated_msa)


def generate_consensus(msa: np.ndarray) -> str:
    """
    Generates a consensus sequence from an msa
    """
    consensus = []
    # For each column in the msa
    for col in range(msa.shape[1]):
        # Create the counter
        col_counter = Counter()
        # For each row in the msa
        for row in range(msa.shape[0]):
            # Update the counter with the de-ambiguous bases
            col_counter.update(extend_ambiguous_base(msa[row, col]))

        # Remove invalid bases if other bases are available
        col_counter.pop("N", None)

        if len(col_counter) == 0:
            consensus.append("N")
        else:
            consensus.append(col_counter.most_common(1)[0][0])
    return "".join(consensus)


def generate_reference(msa: np.ndarray) -> str:
    """
    Generates a reference string from the first row of a multiple sequence alignment (MSA) array.

    Args:
    - msa (np.ndarray): A numpy array representing a multiple sequence alignment.

    Returns:
    - str: A string representing the reference sequence, obtained by joining the characters in the first row of the MSA and removing any gaps represented by hyphens ("-").
    """

    return "".join(msa[0]).replace("-", "")
