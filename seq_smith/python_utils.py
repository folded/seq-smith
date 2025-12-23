import numpy as np

from ._seq_smith import Alignment, AlignmentFragment, FragmentType


def make_score_matrix(
    alphabet: str,
    match_score: int,
    mismatch_score: int,
    *,
    ambiguous: str = "",
    ambiguous_match_score: int | None = None,
) -> np.ndarray:
    """
    Creates a scoring matrix for a given alphabet.

    Args:
        alphabet (str): A string containing all characters in the alphabet.
        match_score (int): The score for a match.
        mismatch_score (int): The score for a mismatch.
        ambiguous (str): Characters that represent an ambiguous match.
        ambiguous_match_score (int): The score for an ambiguous match.

    Returns:
        numpy.ndarray: A 2D numpy array representing the scoring matrix.
    """
    alpha_len = len(alphabet)
    score_matrix = np.full((alpha_len, alpha_len), mismatch_score, dtype=np.int32)
    np.fill_diagonal(score_matrix, match_score)
    if ambiguous:
        if ambiguous_match_score is None:
            raise ValueError(
                "ambiguous_match_score and ambiguous must both be provided when ambiguous characters are in use.",
            )
        if not set(alphabet).issuperset(set(ambiguous)):
            raise ValueError(
                "all ambiguous characters must be included in the alphabet.",
            )
        ambiguous_chars = set(ambiguous)
        for i, char in enumerate(alphabet):
            if char in ambiguous_chars:
                score_matrix[i, :] = ambiguous_match_score
                score_matrix[:, i] = ambiguous_match_score
    return score_matrix


def encode(seq: str, alphabet: str) -> bytes:
    """
    Encodes a sequence into a byte array using the provided alphabet.

    Args:
        seq (str): The sequence to encode.
        alphabet (str): A string containing all characters in the alphabet.

    Returns:
        bytes: The encoded sequence as a byte array.
    """
    char_to_index = {char: idx for idx, char in enumerate(alphabet)}
    return bytes(char_to_index[char] for char in seq)


def decode(encoded_seq: bytes, alphabet: str) -> str:
    """
    Decodes a byte-encoded sequence back to a string using the provided alphabet.

    Args:
        encoded_seq (bytes): The byte-encoded sequence.
        alphabet (str): A string containing all characters in the alphabet.

    Returns:
        str: The decoded sequence as a string.
    """
    return "".join(alphabet[b] for b in encoded_seq)


def format_alignment_ascii(
    seqa_bytes: bytes,
    seqb_bytes: bytes,
    fragments: list[AlignmentFragment],
    alphabet: str,
) -> tuple[str, str]:
    """
    Formats an alignment into a human-readable ASCII string.

    Args:
        seqa_bytes (bytes): The first sequence as a byte array.
        seqb_bytes (bytes): The second sequence as a byte array.
        fragments (list[AlignmentFragment]): A list of alignment fragments.
        alphabet (str): The alphabet used for encoding/decoding.

    Returns:
        tuple[str, str]: A tuple containing the aligned sequences as ASCII strings.
    """
    seqa = decode(seqa_bytes, alphabet)
    seqb = decode(seqb_bytes, alphabet)
    aligned_seqa_list = []
    aligned_seqb_list = []

    for frag in fragments:
        match frag.fragment_type:
            case FragmentType.Match:
                aligned_seqa_list.append(seqa[frag.sa_start : frag.sa_start + frag.len])
                aligned_seqb_list.append(seqb[frag.sb_start : frag.sb_start + frag.len])
            case FragmentType.AGap:
                aligned_seqa_list.append("-" * frag.len)
                aligned_seqb_list.append(seqb[frag.sb_start : frag.sb_start + frag.len])
            case FragmentType.BGap:
                aligned_seqa_list.append(seqa[frag.sa_start : frag.sa_start + frag.len])
                aligned_seqb_list.append("-" * frag.len)

    return "".join(aligned_seqa_list), "".join(aligned_seqb_list)


def generate_cigar(alignment: Alignment) -> str:
    """
    Generates a CIGAR string from an Alignment object.

    Args:
        alignment (Alignment): An Alignment object representing the alignment
          of seqa (the reference) and seqb (the query).

    Returns:
        str: The CIGAR string.
    """
    cigar_parts = []
    for frag in alignment.fragments:
        op = ""
        match frag.fragment_type:
            case FragmentType.Match:
                op = "M"
            case FragmentType.AGap:
                op = "I"  # Insertion in the query relative to the reference.
            case FragmentType.BGap:
                op = "D"  # Deletion in the query relative to the reference.
        cigar_parts.append(f"{frag.len}{op}")
    return "".join(cigar_parts)
