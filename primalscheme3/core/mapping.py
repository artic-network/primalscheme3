from collections import Counter

import numpy as np

# Module imports
from primalscheme3.core.seq_functions import extend_ambiguous_base


def create_mapping(
    msa: np.ndarray, mapping_index: int = 0
) -> tuple[np.ndarray, np.ndarray]:
    """
    This returns a tuple of two items:
        mapping_array: list[int | None]
        truncated_msa: np.ndarray
    mapping_array: Each position in the list corresponds to the same index in the MSA, The value in the list is the position in the reference genome
    """
    # As NP is modified in place, returning is not necessary but is done for clarity
    # Create the empty mapping array
    mapping_list = [None] * msa.shape[1]
    mapping_array = np.array(mapping_list)
    # Select the reference genome
    reference_genome = msa[mapping_index]
    # Iterate over the msa genome
    current_ref_index = 0
    for col_index in range(msa.shape[1]):
        # If the base is not a gap, assign the mapping
        if reference_genome[col_index] not in {"", "-"}:
            mapping_array[col_index] = current_ref_index
            # increase reference index
            current_ref_index += 1
    return (mapping_array, msa)


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


def ref_index_to_msa(mapping_array: np.ndarray) -> dict[int, int]:
    """
    Convert a reference index to an MSA index
    """
    ref_dict = {x: i for i, x in enumerate(list(mapping_array)) if x is not None}
    ref_dict[max(ref_dict.keys()) + 1] = (
        max(ref_dict.values()) + 1
    )  # This ensures that an fprimer with non-inclusive end will not cause key error.

    return ref_dict


def check_for_end_on_gap(ref_index_to_msa: dict[int, int], ref_index) -> bool:
    """
    Check if a slice of a mapping array ends on a gap
    Returns True if the slice ends on a gap, False otherwise

    # Example

           5' AGAGTGTGGGGGTAGTGTTACG          > MPXV_142_LEFT_1 170931:170953
    TTTTTTTTATAGAGTGTGGGGGTAGTGTTACG-------GAT >MT903345
    TTTTTTTTATAGAGTGT-GGGGTAGTGTTACGGATATCTGAT >KJ642613.1
                                           ^ 173598: MSA
                                           ^ 170953: ref

    In this example, the slice of the primer ends on a gap. So slicing the array with
    - array[:, ref_to_msa[170931]:ref_to_msa[170953]] will return "GGGGTAGTGTTACG-------" as the non exclusive end captures the gap
    - fix_end_on_gap() will return in indexes to slice the array without the gap ie "TGTGGGGGTAGTGTTACG"

    """
    exclusive_msa_end = ref_index_to_msa[ref_index]
    inclusive_msa_end = ref_index_to_msa[ref_index - 1]
    return exclusive_msa_end - inclusive_msa_end != 1


def fix_end_on_gap(ref_index_to_msa: dict[int, int], ref_index) -> int:
    """
    Returns the MSA index of the non-inclusive end of a slice with the gap removed
    """
    return ref_index_to_msa[ref_index - 1] + 1
