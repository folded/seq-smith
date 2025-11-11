import numpy as np
import pytest
from c_align import (
    AlignParams,
    CAlignmentPtr,
    FragType,
    global_align,
    glocal_align,
    lib,
    local_align,
    overlap_align,
)


def _encode(seq: str, alphabet: str) -> bytes:
    """Encode a sequence using the provided alphabet."""
    mapping = {char: i for i, char in enumerate(alphabet)}
    return bytes(mapping[char] for char in seq)


# Python representation of AlignFrag and Alignment
class AlignFrag:
    def __init__(
        self,
        frag_type: FragType,
        sa_start: int,
        sb_start: int,
        hsp_len: int,
    ) -> None:
        self.frag_type = frag_type
        self.sa_start = sa_start
        self.sb_start = sb_start
        self.hsp_len = hsp_len

    def __eq__(self, other: object) -> bool:
        return (
            isinstance(other, AlignFrag)
            and self.frag_type == other.frag_type
            and self.sa_start == other.sa_start
            and self.sb_start == other.sb_start
            and self.hsp_len == other.hsp_len
        )

    def __hash__(self) -> int:
        return hash((self.frag_type, self.sa_start, self.sb_start, self.hsp_len))

    def __repr__(self) -> str:
        return (
            f"AlignFrag(frag_type={self.frag_type}, sa_start={self.sa_start}, "
            f"sb_start={self.sb_start}, hsp_len={self.hsp_len})"
        )


class Alignment:
    def __init__(self, score: int, frag_count: int, align_frag: list[AlignFrag]) -> None:
        self.score = score
        self.frag_count = frag_count
        self.align_frag = align_frag

    def __eq__(self, other: object) -> bool:
        return (
            isinstance(other, Alignment)
            and self.score == other.score
            and self.frag_count == other.frag_count
            and self.align_frag == other.align_frag
        )

    def __hash__(self) -> int:
        return hash((self.score, self.frag_count, tuple(self.align_frag)))

    def __repr__(self) -> str:
        return f"Alignment(score={self.score}, frag_count={self.frag_count}, align_frag={self.align_frag})"


def convert_c_alignment_to_py_alignment(
    c_alignment_ptr: CAlignmentPtr | None,
) -> Alignment | None:
    if not c_alignment_ptr:
        return None

    c_alignment = c_alignment_ptr.contents  # type: ignore[attr-defined]
    score = c_alignment.score
    frag_count = c_alignment.frag_count

    py_align_frags = []
    current_c_frag_ptr = c_alignment.align_frag
    while current_c_frag_ptr:
        current_c_frag = current_c_frag_ptr.contents
        py_align_frags.append(
            AlignFrag(
                frag_type=FragType(current_c_frag.type),
                sa_start=current_c_frag.sa_start,
                sb_start=current_c_frag.sb_start,
                hsp_len=current_c_frag.hsp_len,
            ),
        )
        current_c_frag_ptr = current_c_frag.next

    # Free the C allocated memory
    lib.alignment_free(c_alignment_ptr)

    return Alignment(score, frag_count, py_align_frags)


# Fixtures from test_seq_align.py
@pytest.fixture
def common_data() -> AlignParams:
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
    return AlignParams(seqa, seqb, alpha_len, score_matrix, gap_open, gap_extend)


@pytest.fixture
def multi_fragment_data() -> AlignParams:
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
    return AlignParams(seqa, seqb, alpha_len, score_matrix, gap_open, gap_extend)


# Test functions
def test_c_local_align_simple(common_data: AlignParams) -> None:
    c_alignment_ptr = local_align(common_data)
    alignment = convert_c_alignment_to_py_alignment(c_alignment_ptr)

    assert alignment.score == 1
    assert alignment.frag_count == 1
    expected_frags = [
        AlignFrag(frag_type=FragType.MATCH, sa_start=3, sb_start=3, hsp_len=1),
    ]
    assert alignment.align_frag == expected_frags


def test_c_global_align_simple(common_data: AlignParams) -> None:
    c_alignment_ptr = global_align(common_data)
    alignment = convert_c_alignment_to_py_alignment(c_alignment_ptr)

    assert alignment.score == 0
    assert alignment.frag_count == 1
    expected_frags = [
        AlignFrag(frag_type=FragType.MATCH, sa_start=0, sb_start=0, hsp_len=4),
    ]
    assert alignment.align_frag == expected_frags


def test_c_glocal_align_simple(common_data: AlignParams) -> None:
    c_alignment_ptr = glocal_align(common_data)
    alignment = convert_c_alignment_to_py_alignment(c_alignment_ptr)

    assert alignment.score == 0
    assert alignment.frag_count == 1
    expected_frags = [
        AlignFrag(frag_type=FragType.MATCH, sa_start=0, sb_start=0, hsp_len=4),
    ]
    assert alignment.align_frag == expected_frags


def test_c_align_simple(common_data: AlignParams) -> None:
    c_alignment_ptr = overlap_align(common_data)
    alignment = convert_c_alignment_to_py_alignment(c_alignment_ptr)

    assert alignment.score == 0
    assert alignment.frag_count == 1
    expected_frags = [
        AlignFrag(frag_type=FragType.MATCH, sa_start=0, sb_start=0, hsp_len=4),
    ]
    assert alignment.align_frag == expected_frags


def test_c_local_align_multi_fragment(multi_fragment_data: AlignParams) -> None:
    c_alignment_ptr = local_align(multi_fragment_data)
    alignment = convert_c_alignment_to_py_alignment(c_alignment_ptr)

    assert alignment.score == 4
    assert alignment.frag_count == 5
    expected_frags = [
        AlignFrag(frag_type=FragType.MATCH, sa_start=4, sb_start=0, hsp_len=2),
        AlignFrag(frag_type=FragType.A_GAP, sa_start=6, sb_start=2, hsp_len=1),
        AlignFrag(frag_type=FragType.MATCH, sa_start=6, sb_start=3, hsp_len=2),
        AlignFrag(frag_type=FragType.A_GAP, sa_start=8, sb_start=5, hsp_len=1),
        AlignFrag(frag_type=FragType.MATCH, sa_start=8, sb_start=6, hsp_len=2),
    ]
    assert alignment.align_frag == expected_frags


def test_c_global_align_multi_fragment(multi_fragment_data: AlignParams) -> None:
    c_alignment_ptr = global_align(multi_fragment_data)
    alignment = convert_c_alignment_to_py_alignment(c_alignment_ptr)

    assert alignment.score == 2
    assert alignment.frag_count == 8
    expected_frags = [
        AlignFrag(frag_type=FragType.MATCH, sa_start=0, sb_start=0, hsp_len=2),
        AlignFrag(frag_type=FragType.A_GAP, sa_start=2, sb_start=2, hsp_len=1),
        AlignFrag(frag_type=FragType.MATCH, sa_start=2, sb_start=3, hsp_len=2),
        AlignFrag(frag_type=FragType.A_GAP, sa_start=4, sb_start=5, hsp_len=1),
        AlignFrag(frag_type=FragType.MATCH, sa_start=4, sb_start=6, hsp_len=2),
        AlignFrag(frag_type=FragType.B_GAP, sa_start=6, sb_start=8, hsp_len=1),
        AlignFrag(frag_type=FragType.MATCH, sa_start=7, sb_start=8, hsp_len=2),
        AlignFrag(frag_type=FragType.B_GAP, sa_start=9, sb_start=10, hsp_len=1),
    ]
    assert alignment.align_frag == expected_frags


def test_c_glocal_align_multi_fragment(multi_fragment_data: AlignParams) -> None:
    c_alignment_ptr = glocal_align(multi_fragment_data)
    alignment = convert_c_alignment_to_py_alignment(c_alignment_ptr)

    assert alignment.score == 4
    assert alignment.frag_count == 7
    expected_frags = [
        AlignFrag(frag_type=FragType.MATCH, sa_start=2, sb_start=0, hsp_len=2),
        AlignFrag(frag_type=FragType.A_GAP, sa_start=4, sb_start=2, hsp_len=1),
        AlignFrag(frag_type=FragType.MATCH, sa_start=4, sb_start=3, hsp_len=2),
        AlignFrag(frag_type=FragType.A_GAP, sa_start=6, sb_start=5, hsp_len=1),
        AlignFrag(frag_type=FragType.MATCH, sa_start=6, sb_start=6, hsp_len=2),
        AlignFrag(frag_type=FragType.A_GAP, sa_start=8, sb_start=8, hsp_len=1),
        AlignFrag(frag_type=FragType.MATCH, sa_start=8, sb_start=9, hsp_len=1),
    ]
    assert alignment.align_frag == expected_frags


def test_c_align_multi_fragment(multi_fragment_data: AlignParams) -> None:
    c_alignment_ptr = overlap_align(multi_fragment_data)
    alignment = convert_c_alignment_to_py_alignment(c_alignment_ptr)

    assert alignment.score == 4
    assert alignment.frag_count == 5
    expected_frags = [
        AlignFrag(frag_type=FragType.MATCH, sa_start=4, sb_start=0, hsp_len=2),
        AlignFrag(frag_type=FragType.A_GAP, sa_start=6, sb_start=2, hsp_len=1),
        AlignFrag(frag_type=FragType.MATCH, sa_start=6, sb_start=3, hsp_len=2),
        AlignFrag(frag_type=FragType.A_GAP, sa_start=8, sb_start=5, hsp_len=1),
        AlignFrag(frag_type=FragType.MATCH, sa_start=8, sb_start=6, hsp_len=2),
    ]
    assert alignment.align_frag == expected_frags


# Test with empty sequences (expecting ValueError from Python wrapper, not C)
def test_c_local_align_empty_seqa(common_data: AlignParams) -> None:
    common_data.seqa = b""
    c_alignment_ptr = local_align(common_data)
    alignment = convert_c_alignment_to_py_alignment(c_alignment_ptr)
    assert alignment is None or (alignment.score == 0 and alignment.frag_count == 0)


def test_c_local_align_empty_seqb(common_data: AlignParams) -> None:
    common_data.seqb = b""
    c_alignment_ptr = local_align(common_data)
    alignment = convert_c_alignment_to_py_alignment(c_alignment_ptr)
    assert alignment is None or (alignment.score == 0 and alignment.frag_count == 0)


def test_c_global_align_empty_seqa(common_data: AlignParams) -> None:
    common_data.seqa = b""
    c_alignment_ptr = global_align(common_data)
    alignment = convert_c_alignment_to_py_alignment(c_alignment_ptr)
    assert alignment is None or (alignment.score == 0 and alignment.frag_count == 0)


def test_c_global_align_empty_seqb(common_data: AlignParams) -> None:
    common_data.seqb = b""
    c_alignment_ptr = global_align(common_data)
    alignment = convert_c_alignment_to_py_alignment(c_alignment_ptr)
    assert alignment is None or (alignment.score == 0 and alignment.frag_count == 0)


def test_c_glocal_align_empty_seqa(common_data: AlignParams) -> None:
    common_data.seqa = b""
    c_alignment_ptr = glocal_align(common_data)
    alignment = convert_c_alignment_to_py_alignment(c_alignment_ptr)
    assert alignment is None or (alignment.score == 0 and alignment.frag_count == 0)


def test_c_glocal_align_empty_seqb(common_data: AlignParams) -> None:
    common_data.seqb = b""
    c_alignment_ptr = glocal_align(common_data)
    alignment = convert_c_alignment_to_py_alignment(c_alignment_ptr)
    assert alignment is None or (alignment.score == 0 and alignment.frag_count == 0)


def test_c_align_empty_seqa(common_data: AlignParams) -> None:
    common_data.seqa = b""
    c_alignment_ptr = overlap_align(common_data)
    alignment = convert_c_alignment_to_py_alignment(c_alignment_ptr)
    assert alignment is None or (alignment.score == 0 and alignment.frag_count == 0)


def test_c_align_empty_seqb(common_data: AlignParams) -> None:
    common_data.seqb = b""
    c_alignment_ptr = overlap_align(common_data)
    alignment = convert_c_alignment_to_py_alignment(c_alignment_ptr)
    assert alignment is None or (alignment.score == 0 and alignment.frag_count == 0)


def test_local_align_perfect_match_subsegment() -> None:
    alphabet = "ACGTXYZW"
    seqa = _encode("XXXXXAGCTYYYYY", alphabet)
    seqb = _encode("ZZZAGCTWWW", alphabet)
    score_matrix = np.eye(8, 8, dtype=np.int32) * 2 - 1
    gap_open = -2
    gap_extend = -1
    params = AlignParams(seqa, seqb, 8, score_matrix, gap_open, gap_extend)
    alignment = local_align(params)
    alignment = convert_c_alignment_to_py_alignment(alignment)
    assert alignment.score == 4  # AGCT match
    assert alignment.frag_count == 1
    expected_frags = [
        AlignFrag(
            frag_type=FragType.MATCH,
            sa_start=5,
            sb_start=3,
            hsp_len=4,
        ),  # AGCT in seqa starts at index 5, in seqb at index 3
    ]
    assert alignment.align_frag == expected_frags
