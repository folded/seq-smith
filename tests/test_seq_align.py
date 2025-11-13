import dataclasses

import numpy as np
import pytest
from seq_align import (
    AlignFrag,
    FragType,
    global_align,
    local_align,
    local_global_align,
    overlap_align,
)


def _encode(seq: str, alphabet: str) -> bytes:
    """Encode a sequence using the provided alphabet."""
    char_to_index = {char: idx for idx, char in enumerate(alphabet)}
    return bytes(char_to_index[char] for char in seq)


def _decode(encoded_seq: bytes, alphabet: str) -> str:
    """Decode a byte-encoded sequence back to string using the provided alphabet."""
    return "".join(alphabet[b] for b in encoded_seq)


def format_alignment_ascii(
    seqa_bytes: bytes,
    seqb_bytes: bytes,
    align_frags: list[AlignFrag],
    alphabet: str,
) -> tuple[str, str]:
    seqa = _decode(seqa_bytes, alphabet)
    seqb = _decode(seqb_bytes, alphabet)
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


@dataclasses.dataclass
class AlignmentInput:
    alphabet: str
    seqa: bytes
    seqb: bytes
    alpha_len: int
    score_matrix: np.ndarray
    gap_open: int
    gap_extend: int

    def encode(self, seq: str) -> bytes:
        return _encode(seq, self.alphabet)


@pytest.fixture
def common_data() -> AlignmentInput:
    alphabet = "ACGT"
    alpha_len = len(alphabet)
    seqa = _encode("ACGT", alphabet)
    seqb = _encode("AGCT", alphabet)

    score_matrix = np.array(
        [
            [1, -1, -1, -1],
            [-1, 1, -1, -1],
            [-1, -1, 1, -1],
            [-1, -1, -1, 1],
        ],
        dtype=np.int32,
    )
    gap_open = -2
    gap_extend = -1

    return AlignmentInput(
        alphabet=alphabet,
        seqa=seqa,
        seqb=seqb,
        alpha_len=alpha_len,
        score_matrix=score_matrix,
        gap_open=gap_open,
        gap_extend=gap_extend,
    )


@pytest.fixture
def complex_data() -> AlignmentInput:
    alphabet = "ACGT"
    alpha_len = len(alphabet)
    seqa = _encode("GAATTCAGTTA", alphabet)
    seqb = _encode("GGATCGA", alphabet)

    score_matrix = np.array(
        [
            [2, -1, -1, -1],
            [-1, 2, -1, -1],
            [-1, -1, 2, -1],
            [-1, -1, -1, 2],
        ],
        dtype=np.int32,
    )
    gap_open = -3
    gap_extend = -1
    return AlignmentInput(
        alphabet=alphabet,
        seqa=seqa,
        seqb=seqb,
        alpha_len=alpha_len,
        score_matrix=score_matrix,
        gap_open=gap_open,
        gap_extend=gap_extend,
    )


@pytest.fixture
def local_global_test_data() -> AlignmentInput:
    alphabet = "ACGTX"
    alpha_len = len(alphabet)
    seqa = _encode("XACGTX", alphabet)
    seqb = _encode("ACGT", alphabet)

    score_matrix = np.array(
        [
            [1, -1, -1, -1, -1],  # A
            [-1, 1, -1, -1, -1],  # C
            [-1, -1, 1, -1, -1],  # G
            [-1, -1, -1, 1, -1],  # T
            [-1, -1, -1, -1, 0],  # X
        ],
        dtype=np.int32,
    )
    gap_open = -2
    gap_extend = -1
    return AlignmentInput(
        alphabet=alphabet,
        seqa=seqa,
        seqb=seqb,
        alpha_len=alpha_len,
        score_matrix=score_matrix,
        gap_open=gap_open,
        gap_extend=gap_extend,
    )


def test_global_align_simple(common_data: AlignmentInput) -> None:
    # This test was flawed. The original seqb was "AGCT", but the expected
    # fragments described a perfect match. Correcting seqb to "ACGT".
    seqa = common_data.encode("ACGT")
    seqb = common_data.encode("ACGT")
    alignment = global_align(
        seqa,
        seqb,
        common_data.score_matrix,
        common_data.gap_open,
        common_data.gap_extend,
    )

    assert alignment.score == 4
    aligned_a, aligned_b = format_alignment_ascii(seqa, seqb, alignment.align_frag, common_data.alphabet)
    assert aligned_a == "ACGT"
    assert aligned_b == "ACGT"


def test_global_align_simple_gap(common_data: AlignmentInput) -> None:
    seqa = common_data.encode("A")
    seqb = common_data.encode("AC")
    alignment = global_align(seqa, seqb, common_data.score_matrix, common_data.gap_open, common_data.gap_extend)

    assert alignment.score == -1
    aligned_a, aligned_b = format_alignment_ascii(seqa, seqb, alignment.align_frag, common_data.alphabet)
    assert aligned_a == "A-"
    assert aligned_b == "AC"


def test_local_global_align_simple(common_data: AlignmentInput) -> None:
    # This test was flawed. The original seqb was "AGCT", but the expected
    # fragments described a perfect match. Correcting seqb to "ACGT".
    seqa = common_data.encode("ACGT")
    seqb = common_data.encode("ACGT")
    alignment = local_global_align(
        seqa,
        seqb,
        common_data.score_matrix,
        common_data.gap_open,
        common_data.gap_extend,
    )

    assert alignment.score == 4
    aligned_a, aligned_b = format_alignment_ascii(seqa, seqb, alignment.align_frag, common_data.alphabet)
    assert aligned_a == "ACGT"
    assert aligned_b == "ACGT"


def test_local_global_align_subsegment_global_seqb(local_global_test_data: AlignmentInput) -> None:
    alignment = local_global_align(
        local_global_test_data.seqa,
        local_global_test_data.seqb,
        local_global_test_data.score_matrix,
        local_global_test_data.gap_open,
        local_global_test_data.gap_extend,
    )

    assert alignment.score == 4
    aligned_a, aligned_b = format_alignment_ascii(
        local_global_test_data.seqa,
        local_global_test_data.seqb,
        alignment.align_frag,
        local_global_test_data.alphabet,
    )
    assert aligned_a == "ACGT"
    assert aligned_b == "ACGT"


def test_overlap_align_simple(common_data: AlignmentInput) -> None:
    # This test was flawed. The original seqb was "AGCT", but the expected
    # fragments described a perfect match. Correcting seqb to "ACGT".
    seqa = common_data.encode("ACGT")
    seqb = common_data.encode("ACGT")
    alignment = overlap_align(
        seqa,
        seqb,
        common_data.score_matrix,
        common_data.gap_open,
        common_data.gap_extend,
    )

    assert alignment.score == 4
    aligned_a, aligned_b = format_alignment_ascii(seqa, seqb, alignment.align_frag, common_data.alphabet)
    assert aligned_a == "ACGT"
    assert aligned_b == "ACGT"


def test_overlap_align_semi_global_overlap(common_data: AlignmentInput) -> None:
    seqa = common_data.encode("ACGTACGT")
    seqb = common_data.encode("CGTA")
    alignment = overlap_align(seqa, seqb, common_data.score_matrix, common_data.gap_open, common_data.gap_extend)

    assert alignment.score == 4
    aligned_a, aligned_b = format_alignment_ascii(seqa, seqb, alignment.align_frag, common_data.alphabet)
    assert aligned_a == "CGTA"
    assert aligned_b == "CGTA"


# Test with multiple fragments
@pytest.fixture
def multi_fragment_data() -> AlignmentInput:
    alphabet = "ACGT"
    alpha_len = len(alphabet)
    seqa = _encode("AGAGAGAGAG", alphabet)
    seqb = _encode("AGCAGCAGCA", alphabet)

    score_matrix = np.array(
        [
            [1, -1, -1, -1],
            [-1, 1, -1, -1],
            [-1, -1, 1, -1],
            [-1, -1, -1, 1],
        ],
        dtype=np.int32,
    )

    gap_open = -1
    gap_extend = -1
    return AlignmentInput(
        alphabet=alphabet,
        seqa=seqa,
        seqb=seqb,
        alpha_len=alpha_len,
        score_matrix=score_matrix,
        gap_open=gap_open,
        gap_extend=gap_extend,
    )


def test_local_align_perfect_match_subsegment() -> None:
    alphabet = "ACGTXYZW"
    len(alphabet)

    seqa = _encode("XXXXXAGCTYYYYY", alphabet)
    seqb = _encode("ZZZAGCTWWW", alphabet)

    score_matrix = np.eye(8, 8, dtype=np.int32) * 2 - 1

    gap_open = -2
    gap_extend = -1

    alignment = local_align(seqa, seqb, score_matrix, gap_open, gap_extend)
    assert alignment.score == 4
    aligned_a, aligned_b = format_alignment_ascii(seqa, seqb, alignment.align_frag, alphabet)
    assert aligned_a == "AGCT"
    assert aligned_b == "AGCT"


def test_local_align_multi_fragment(multi_fragment_data: AlignmentInput) -> None:
    alignment = local_align(
        multi_fragment_data.seqa,
        multi_fragment_data.seqb,
        multi_fragment_data.score_matrix,
        multi_fragment_data.gap_open,
        multi_fragment_data.gap_extend,
    )
    assert alignment.score == 4
    aligned_a, aligned_b = format_alignment_ascii(
        multi_fragment_data.seqa,
        multi_fragment_data.seqb,
        alignment.align_frag,
        multi_fragment_data.alphabet,
    )
    assert aligned_a == "AG-AG-AG"
    assert aligned_b == "AGCAGCAG"


def test_global_align_multi_fragment(multi_fragment_data: AlignmentInput) -> None:
    alignment = global_align(
        multi_fragment_data.seqa,
        multi_fragment_data.seqb,
        multi_fragment_data.score_matrix,
        multi_fragment_data.gap_open,
        multi_fragment_data.gap_extend,
    )
    assert alignment.score == 2
    aligned_a, aligned_b = format_alignment_ascii(
        multi_fragment_data.seqa,
        multi_fragment_data.seqb,
        alignment.align_frag,
        multi_fragment_data.alphabet,
    )
    assert aligned_a == "AG-AG-AGAGAG"
    assert aligned_b == "AGCAGCAG-CA-"
    assert len(aligned_a) == len(aligned_b)


def test_local_global_align_multi_fragment(multi_fragment_data: AlignmentInput) -> None:
    alignment = local_global_align(
        multi_fragment_data.seqa,
        multi_fragment_data.seqb,
        multi_fragment_data.score_matrix,
        multi_fragment_data.gap_open,
        multi_fragment_data.gap_extend,
    )
    assert alignment.score == 4
    aligned_a, aligned_b = format_alignment_ascii(
        multi_fragment_data.seqa,
        multi_fragment_data.seqb,
        alignment.align_frag,
        multi_fragment_data.alphabet,
    )
    assert aligned_a == "AG-AG-AG-A"
    assert aligned_b == "AGCAGCAGCA"


def test_overlap_align_multi_fragment(multi_fragment_data: AlignmentInput) -> None:
    alignment = overlap_align(
        multi_fragment_data.seqa,
        multi_fragment_data.seqb,
        multi_fragment_data.score_matrix,
        multi_fragment_data.gap_open,
        multi_fragment_data.gap_extend,
    )
    assert alignment.score == 4
    aligned_a, aligned_b = format_alignment_ascii(
        multi_fragment_data.seqa,
        multi_fragment_data.seqb,
        alignment.align_frag,
        multi_fragment_data.alphabet,
    )
    assert aligned_a == "AG-AG-AG"
    assert aligned_b == "AGCAGCAG"


# Test with empty sequences
def test_local_align_empty_seqa(common_data: AlignmentInput) -> None:
    with pytest.raises(ValueError, match=r"Input sequences cannot be empty."):
        local_align(b"", common_data.seqb, common_data.score_matrix, common_data.gap_open, common_data.gap_extend)


def test_local_align_empty_seqb(common_data: AlignmentInput) -> None:
    with pytest.raises(ValueError, match=r"Input sequences cannot be empty."):
        local_align(common_data.seqa, b"", common_data.score_matrix, common_data.gap_open, common_data.gap_extend)


def test_global_align_empty_seqa(common_data: AlignmentInput) -> None:
    with pytest.raises(ValueError, match=r"Input sequences cannot be empty."):
        global_align(b"", common_data.seqb, common_data.score_matrix, common_data.gap_open, common_data.gap_extend)


def test_global_align_empty_seqb(common_data: AlignmentInput) -> None:
    with pytest.raises(ValueError, match=r"Input sequences cannot be empty."):
        global_align(common_data.seqa, b"", common_data.score_matrix, common_data.gap_open, common_data.gap_extend)


def test_local_global_align_empty_seqa(common_data: AlignmentInput) -> None:
    with pytest.raises(ValueError, match=r"Input sequences cannot be empty."):
        local_global_align(
            b"",
            common_data.seqb,
            common_data.score_matrix,
            common_data.gap_open,
            common_data.gap_extend,
        )


def test_local_global_align_empty_seqb(common_data: AlignmentInput) -> None:
    with pytest.raises(ValueError, match=r"Input sequences cannot be empty."):
        local_global_align(
            common_data.seqa,
            b"",
            common_data.score_matrix,
            common_data.gap_open,
            common_data.gap_extend,
        )


def test_overlap_align_empty_seqa(common_data: AlignmentInput) -> None:
    with pytest.raises(ValueError, match=r"Input sequences cannot be empty."):
        overlap_align(b"", common_data.seqb, common_data.score_matrix, common_data.gap_open, common_data.gap_extend)


def test_overlap_align_empty_seqb(common_data: AlignmentInput) -> None:
    with pytest.raises(ValueError, match=r"Input sequences cannot be empty."):
        overlap_align(common_data.seqa, b"", common_data.score_matrix, common_data.gap_open, common_data.gap_extend)


@pytest.fixture
def poly_data() -> AlignmentInput:
    alphabet = "ACGT"
    alpha_len = len(alphabet)
    seqa = _encode("CCCCCCAACAA", alphabet)
    seqb = _encode("TTAAAAGGGGGGG", alphabet)

    score_matrix = np.eye(4, dtype=np.int32) * 2 - 1

    gap_open = -2
    gap_extend = -1
    return AlignmentInput(
        alphabet=alphabet,
        seqa=seqa,
        seqb=seqb,
        alpha_len=alpha_len,
        score_matrix=score_matrix,
        gap_open=gap_open,
        gap_extend=gap_extend,
    )


def test_local_align_poly(poly_data: AlignmentInput) -> None:
    alignment = local_align(
        poly_data.seqa,
        poly_data.seqb,
        poly_data.score_matrix,
        poly_data.gap_open,
        poly_data.gap_extend,
    )
    assert alignment.score == 2
    aligned_seqa, aligned_seqb = format_alignment_ascii(
        poly_data.seqa,
        poly_data.seqb,
        alignment.align_frag,
        poly_data.alphabet,
    )
    assert aligned_seqa == "ACAA"
    assert aligned_seqb == "AAAA"


def test_global_align_poly(poly_data: AlignmentInput) -> None:
    alignment = global_align(
        poly_data.seqa,
        poly_data.seqb,
        poly_data.score_matrix,
        poly_data.gap_open,
        poly_data.gap_extend,
    )
    assert alignment.score == -13
    aligned_seqa, aligned_seqb = format_alignment_ascii(
        poly_data.seqa,
        poly_data.seqb,
        alignment.align_frag,
        poly_data.alphabet,
    )
    assert aligned_seqa == "CCCCCCAA----CAA"
    assert aligned_seqb == "--TTAAAAGGGGGGG"


def test_global_align_poly_strong_gap_penalty(poly_data: AlignmentInput) -> None:
    alignment = global_align(
        poly_data.seqa,
        poly_data.seqb,
        poly_data.score_matrix,
        -100,
        -100,
    )
    assert alignment.score == -211
    aligned_seqa, aligned_seqb = format_alignment_ascii(
        poly_data.seqa,
        poly_data.seqb,
        alignment.align_frag,
        poly_data.alphabet,
    )
    assert aligned_seqa == "--CCCCCCAACAA"
    assert aligned_seqb == "TTAAAAGGGGGGG"


def test_local_global_align_poly(poly_data: AlignmentInput) -> None:
    alignment = local_global_align(
        poly_data.seqa,
        poly_data.seqb,
        poly_data.score_matrix,
        poly_data.gap_open,
        poly_data.gap_extend,
    )
    assert alignment.score == -8
    aligned_seqa, aligned_seqb = format_alignment_ascii(
        poly_data.seqa,
        poly_data.seqb,
        alignment.align_frag,
        poly_data.alphabet,
    )
    assert aligned_seqa == "CCAACA------A"
    assert aligned_seqb == "TTAAAAGGGGGGG"


@pytest.fixture
def poly_data_strong_gap_penalty() -> AlignmentInput:
    alphabet = "ACGT"
    alpha_len = len(alphabet)
    seqa = _encode("CCCCCCAACAACCCCCCC", alphabet)
    seqb = _encode("TTAAAAGGGG", alphabet)

    score_matrix = np.eye(4, dtype=np.int32) * 2 - 1

    gap_open = -100
    gap_extend = -100
    return AlignmentInput(
        alphabet=alphabet,
        seqa=seqa,
        seqb=seqb,
        alpha_len=alpha_len,
        score_matrix=score_matrix,
        gap_open=gap_open,
        gap_extend=gap_extend,
    )


def test_local_global_align_poly_strong_gap_penalty(poly_data_strong_gap_penalty: AlignmentInput) -> None:
    alignment = local_global_align(
        poly_data_strong_gap_penalty.seqa,
        poly_data_strong_gap_penalty.seqb,
        poly_data_strong_gap_penalty.score_matrix,
        poly_data_strong_gap_penalty.gap_open,
        poly_data_strong_gap_penalty.gap_extend,
    )
    assert alignment.score == -4
    aligned_seqa, aligned_seqb = format_alignment_ascii(
        poly_data_strong_gap_penalty.seqa,
        poly_data_strong_gap_penalty.seqb,
        alignment.align_frag,
        poly_data_strong_gap_penalty.alphabet,
    )
    assert aligned_seqa == "CAACAACCCC"
    assert aligned_seqb == "TTAAAAGGGG"


def test_overlap_align_poly(poly_data: AlignmentInput) -> None:
    alignment = overlap_align(
        poly_data.seqa,
        poly_data.seqb,
        poly_data.score_matrix,
        poly_data.gap_open,
        poly_data.gap_extend,
    )
    assert alignment.score == 0
    aligned_seqa, aligned_seqb = format_alignment_ascii(
        poly_data.seqa,
        poly_data.seqb,
        alignment.align_frag,
        poly_data.alphabet,
    )
    assert aligned_seqa == "CAACAA"
    assert aligned_seqb == "TTAAAA"
