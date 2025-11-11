import dataclasses

import numpy as np
import pytest
from seq_align import (
    AlignFrag,
    FragType,
    global_align,
    glocal_align,
    local_align,
    overlap_align,
)


def _encode(seq: str, alphabet: str) -> bytes:
    """Encode a sequence using the provided alphabet."""
    char_to_index = {char: idx for idx, char in enumerate(alphabet)}
    return bytes(char_to_index[char] for char in seq)


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
def glocal_test_data() -> AlignmentInput:
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
    alignment = global_align(
        common_data.seqa,
        common_data.seqb,
        common_data.score_matrix,
        common_data.gap_open,
        common_data.gap_extend,
    )

    assert alignment.score == 0
    assert alignment.frag_count == 1
    expected_frags = [
        AlignFrag(frag_type=FragType.Match, sa_start=0, sb_start=0, len=4),
    ]
    assert alignment.align_frag == expected_frags


def test_global_align_simple_gap(common_data: AlignmentInput) -> None:
    seqa = common_data.encode("A")
    seqb = common_data.encode("AC")
    alignment = global_align(seqa, seqb, common_data.score_matrix, common_data.gap_open, common_data.gap_extend)

    assert alignment.score == -1
    assert alignment.frag_count == 2
    expected_frags = [
        AlignFrag(frag_type=FragType.Match, sa_start=0, sb_start=0, len=1),
        AlignFrag(frag_type=FragType.AGap, sa_start=1, sb_start=1, len=1),
    ]
    assert alignment.align_frag == expected_frags


def test_glocal_align_simple(common_data: AlignmentInput) -> None:
    alignment = glocal_align(
        common_data.seqa,
        common_data.seqb,
        common_data.score_matrix,
        common_data.gap_open,
        common_data.gap_extend,
    )

    assert alignment.score == 0
    assert alignment.frag_count == 1
    expected_frags = [
        AlignFrag(frag_type=FragType.Match, sa_start=0, sb_start=0, len=4),
    ]
    assert alignment.align_frag == expected_frags


def test_glocal_align_subsegment_global_seqb(glocal_test_data: AlignmentInput) -> None:
    alignment = glocal_align(
        glocal_test_data.seqa,
        glocal_test_data.seqb,
        glocal_test_data.score_matrix,
        glocal_test_data.gap_open,
        glocal_test_data.gap_extend,
    )

    assert alignment.score == 4
    assert alignment.frag_count == 1

    expected_frags = [
        AlignFrag(frag_type=FragType.Match, sa_start=1, sb_start=0, len=4),
    ]
    assert alignment.align_frag == expected_frags


def test_overlap_align_simple(common_data: AlignmentInput) -> None:
    alignment = overlap_align(
        common_data.seqa,
        common_data.seqb,
        common_data.score_matrix,
        common_data.gap_open,
        common_data.gap_extend,
    )

    assert alignment.score == 0
    assert alignment.frag_count == 1

    expected_frags = [
        AlignFrag(frag_type=FragType.Match, sa_start=0, sb_start=0, len=4),
    ]
    assert alignment.align_frag == expected_frags


def test_overlap_align_semi_global_overlap(common_data: AlignmentInput) -> None:
    seqa = common_data.encode("ACGTACGT")
    seqb = common_data.encode("CGTA")
    alignment = overlap_align(seqa, seqb, common_data.score_matrix, common_data.gap_open, common_data.gap_extend)

    assert alignment.score == 4
    assert alignment.frag_count == 1

    expected_frags = [
        AlignFrag(frag_type=FragType.Match, sa_start=1, sb_start=0, len=4),
    ]
    assert alignment.align_frag == expected_frags


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
    assert alignment.score == 4  # AGCT match
    assert alignment.frag_count == 1
    expected_frags = [
        AlignFrag(
            frag_type=FragType.Match,
            sa_start=5,
            sb_start=3,
            len=4,
        ),  # AGCT in seqa starts at index 5, in seqb at index 3
    ]
    assert alignment.align_frag == expected_frags


def test_local_align_multi_fragment(multi_fragment_data: AlignmentInput) -> None:
    alignment = local_align(
        multi_fragment_data.seqa,
        multi_fragment_data.seqb,
        multi_fragment_data.score_matrix,
        multi_fragment_data.gap_open,
        multi_fragment_data.gap_extend,
    )
    assert alignment.score == 4
    assert alignment.frag_count == 5
    expected_frags = [
        AlignFrag(frag_type=FragType.Match, sa_start=4, sb_start=0, len=2),
        AlignFrag(frag_type=FragType.AGap, sa_start=6, sb_start=2, len=1),
        AlignFrag(frag_type=FragType.Match, sa_start=6, sb_start=3, len=2),
        AlignFrag(frag_type=FragType.AGap, sa_start=8, sb_start=5, len=1),
        AlignFrag(frag_type=FragType.Match, sa_start=8, sb_start=6, len=2),
    ]
    assert alignment.align_frag == expected_frags


def test_global_align_multi_fragment(multi_fragment_data: AlignmentInput) -> None:
    alignment = global_align(
        multi_fragment_data.seqa,
        multi_fragment_data.seqb,
        multi_fragment_data.score_matrix,
        multi_fragment_data.gap_open,
        multi_fragment_data.gap_extend,
    )
    assert alignment.score == 2
    assert alignment.frag_count == 8
    expected_frags = [
        AlignFrag(frag_type=FragType.Match, sa_start=0, sb_start=0, len=2),
        AlignFrag(frag_type=FragType.AGap, sa_start=2, sb_start=2, len=1),
        AlignFrag(frag_type=FragType.Match, sa_start=2, sb_start=3, len=2),
        AlignFrag(frag_type=FragType.AGap, sa_start=4, sb_start=5, len=1),
        AlignFrag(frag_type=FragType.Match, sa_start=4, sb_start=6, len=2),
        AlignFrag(frag_type=FragType.BGap, sa_start=6, sb_start=8, len=1),
        AlignFrag(frag_type=FragType.Match, sa_start=7, sb_start=8, len=2),
        AlignFrag(frag_type=FragType.BGap, sa_start=9, sb_start=10, len=1),
    ]
    assert alignment.align_frag == expected_frags


def test_glocal_align_multi_fragment(multi_fragment_data: AlignmentInput) -> None:
    alignment = glocal_align(
        multi_fragment_data.seqa,
        multi_fragment_data.seqb,
        multi_fragment_data.score_matrix,
        multi_fragment_data.gap_open,
        multi_fragment_data.gap_extend,
    )
    assert alignment.score == 4
    assert alignment.frag_count == 7
    expected_frags = [
        AlignFrag(frag_type=FragType.Match, sa_start=2, sb_start=0, len=2),
        AlignFrag(frag_type=FragType.AGap, sa_start=4, sb_start=2, len=1),
        AlignFrag(frag_type=FragType.Match, sa_start=4, sb_start=3, len=2),
        AlignFrag(frag_type=FragType.AGap, sa_start=6, sb_start=5, len=1),
        AlignFrag(frag_type=FragType.Match, sa_start=6, sb_start=6, len=2),
        AlignFrag(frag_type=FragType.AGap, sa_start=8, sb_start=8, len=1),
        AlignFrag(frag_type=FragType.Match, sa_start=8, sb_start=9, len=1),
    ]
    assert alignment.align_frag == expected_frags


def test_overlap_align_multi_fragment(multi_fragment_data: AlignmentInput) -> None:
    alignment = overlap_align(
        multi_fragment_data.seqa,
        multi_fragment_data.seqb,
        multi_fragment_data.score_matrix,
        multi_fragment_data.gap_open,
        multi_fragment_data.gap_extend,
    )
    assert alignment.score == 4
    assert alignment.frag_count == 5
    expected_frags = [
        AlignFrag(frag_type=FragType.Match, sa_start=4, sb_start=0, len=2),
        AlignFrag(frag_type=FragType.AGap, sa_start=6, sb_start=2, len=1),
        AlignFrag(frag_type=FragType.Match, sa_start=6, sb_start=3, len=2),
        AlignFrag(frag_type=FragType.AGap, sa_start=8, sb_start=5, len=1),
        AlignFrag(frag_type=FragType.Match, sa_start=8, sb_start=6, len=2),
    ]
    assert alignment.align_frag == expected_frags


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


def test_glocal_align_empty_seqa(common_data: AlignmentInput) -> None:
    with pytest.raises(ValueError, match=r"Input sequences cannot be empty."):
        glocal_align(b"", common_data.seqb, common_data.score_matrix, common_data.gap_open, common_data.gap_extend)


def test_glocal_align_empty_seqb(common_data: AlignmentInput) -> None:
    with pytest.raises(ValueError, match=r"Input sequences cannot be empty."):
        glocal_align(common_data.seqa, b"", common_data.score_matrix, common_data.gap_open, common_data.gap_extend)


def test_overlap_align_empty_seqa(common_data: AlignmentInput) -> None:
    with pytest.raises(ValueError, match=r"Input sequences cannot be empty."):
        overlap_align(b"", common_data.seqb, common_data.score_matrix, common_data.gap_open, common_data.gap_extend)


def test_overlap_align_empty_seqb(common_data: AlignmentInput) -> None:
    with pytest.raises(ValueError, match=r"Input sequences cannot be empty."):
        overlap_align(common_data.seqa, b"", common_data.score_matrix, common_data.gap_open, common_data.gap_extend)
