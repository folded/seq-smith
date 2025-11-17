import numpy as np

from ._seq_smith import AlignFrag, FragType


def make_score_matrix(alphabet: str, match_score: int, mismatch_score: int) -> np.ndarray:
    """
    Creates a scoring matrix for a given alphabet.

    Args:
        alphabet (str): A string containing all characters in the alphabet.
        match_score (int): The score for a match.
        mismatch_score (int): The score for a mismatch.

    Returns:
        numpy.ndarray: A 2D numpy array representing the scoring matrix.
    """
    alpha_len = len(alphabet)
    score_matrix = np.full((alpha_len, alpha_len), mismatch_score, dtype=np.int32)
    np.fill_diagonal(score_matrix, match_score)
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
    align_frags: list[AlignFrag],
    alphabet: str,
) -> tuple[str, str]:
    """
    Formats an alignment into a human-readable ASCII string.

    Args:
        seqa_bytes (bytes): The first sequence as a byte array.
        seqb_bytes (bytes): The second sequence as a byte array.
        align_frags (list[AlignFrag]): A list of alignment fragments.
        alphabet (str): The alphabet used for encoding/decoding.

    Returns:
        tuple[str, str]: A tuple containing the aligned sequences as ASCII strings.
    """
    seqa = decode(seqa_bytes, alphabet)
    seqb = decode(seqb_bytes, alphabet)
    aligned_seqa_list = []
    aligned_seqb_list = []

    for frag in align_frags:
        match frag.frag_type:
            case FragType.Match:
                aligned_seqa_list.append(seqa[frag.sa_start : frag.sa_start + frag.len])
                aligned_seqb_list.append(seqb[frag.sb_start : frag.sb_start + frag.len])
            case FragType.AGap:
                aligned_seqa_list.append("-" * frag.len)
                aligned_seqb_list.append(seqb[frag.sb_start : frag.sb_start + frag.len])
            case FragType.BGap:
                aligned_seqa_list.append(seqa[frag.sa_start : frag.sa_start + frag.len])
                aligned_seqb_list.append("-" * frag.len)

    return "".join(aligned_seqa_list), "".join(aligned_seqb_list)
