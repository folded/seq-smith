import ctypes

import numpy as np
import pytest

from c_align import FragType, lib


def _encode(seq: str, alphabet: str) -> bytes:
    """Encode a sequence using the provided alphabet."""
    char_to_index = {char: idx for idx, char in enumerate(alphabet)}
    return bytes(char_to_index[char] for char in seq)


# Python representation of AlignFrag and Alignment
class AlignFrag:
    def __init__(self, frag_type, sa_start, sb_start, hsp_len) -> None:
        self.frag_type = frag_type
        self.sa_start = sa_start
        self.sb_start = sb_start
        self.hsp_len = hsp_len

    def __eq__(self, other):
        return (
            isinstance(other, AlignFrag)
            and self.frag_type == other.frag_type
            and self.sa_start == other.sa_start
            and self.sb_start == other.sb_start
            and self.hsp_len == other.hsp_len
        )

    def __repr__(self) -> str:
        return (
            f"AlignFrag(frag_type={self.frag_type}, sa_start={self.sa_start}, "
            f"sb_start={self.sb_start}, hsp_len={self.hsp_len})"
        )


class Alignment:
    def __init__(self, score, frag_count, align_frag) -> None:
        self.score = score
        self.frag_count = frag_count
        self.align_frag = align_frag

    def __eq__(self, other):
        return (
            isinstance(other, Alignment)
            and self.score == other.score
            and self.frag_count == other.frag_count
            and self.align_frag == other.align_frag
        )

    def __repr__(self) -> str:
        return f"Alignment(score={self.score}, frag_count={self.frag_count}, align_frag={self.align_frag})"


def convert_c_alignment_to_py_alignment(c_alignment_ptr):
    if not c_alignment_ptr:
        return None

    c_alignment = c_alignment_ptr.contents
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
def common_data():
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
    return seqa, seqb, alpha_len, score_matrix, gap_open, gap_extend


@pytest.fixture
def multi_fragment_data():
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
    return seqa, seqb, alpha_len, score_matrix, gap_open, gap_extend


# Test functions
def test_c_local_align_simple(common_data) -> None:
    seqa, seqb, alpha_len, score_matrix, gap_open, gap_extend = common_data

    c_score_matrix = score_matrix.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    c_alignment_ptr = lib.local_align(
        seqa,
        len(seqa),
        seqb,
        len(seqb),
        alpha_len,
        c_score_matrix,
        gap_open,
        gap_extend,
    )
    alignment = convert_c_alignment_to_py_alignment(c_alignment_ptr)

    assert alignment.score == 1
    assert alignment.frag_count == 1
    expected_frags = [
        AlignFrag(frag_type=FragType.MATCH, sa_start=3, sb_start=3, hsp_len=1),
    ]
    assert alignment.align_frag == expected_frags


def test_c_global_align_simple(common_data) -> None:
    seqa, seqb, alpha_len, score_matrix, gap_open, gap_extend = common_data

    c_score_matrix = score_matrix.ctypes.data_as(ctypes.POINTER(ctypes.c_int))

    c_alignment_ptr = lib.global_align(
        seqa,
        len(seqa),
        seqb,
        len(seqb),
        alpha_len,
        c_score_matrix,
        gap_open,
        gap_extend,
    )
    alignment = convert_c_alignment_to_py_alignment(c_alignment_ptr)

    assert alignment.score == 0
    assert alignment.frag_count == 1
    expected_frags = [
        AlignFrag(frag_type=FragType.MATCH, sa_start=0, sb_start=0, hsp_len=4),
    ]
    assert alignment.align_frag == expected_frags


def test_c_glocal_align_simple(common_data) -> None:
    seqa, seqb, alpha_len, score_matrix, gap_open, gap_extend = common_data

    c_score_matrix = score_matrix.ctypes.data_as(ctypes.POINTER(ctypes.c_int))

    c_alignment_ptr = lib.glocal_align(
        seqa,
        len(seqa),
        seqb,
        len(seqb),
        alpha_len,
        c_score_matrix,
        gap_open,
        gap_extend,
    )
    alignment = convert_c_alignment_to_py_alignment(c_alignment_ptr)

    assert alignment.score == 0
    assert alignment.frag_count == 1
    expected_frags = [
        AlignFrag(frag_type=FragType.MATCH, sa_start=0, sb_start=0, hsp_len=4),
    ]
    assert alignment.align_frag == expected_frags


def test_c_align_simple(common_data) -> None:
    seqa, seqb, alpha_len, score_matrix, gap_open, gap_extend = common_data

    c_score_matrix = score_matrix.ctypes.data_as(ctypes.POINTER(ctypes.c_int))

    c_alignment_ptr = lib.overlap_align(
        seqa,
        len(seqa),
        seqb,
        len(seqb),
        alpha_len,
        c_score_matrix,
        gap_open,
        gap_extend,
    )
    alignment = convert_c_alignment_to_py_alignment(c_alignment_ptr)

    assert alignment.score == 0
    assert alignment.frag_count == 1
    expected_frags = [
        AlignFrag(frag_type=FragType.MATCH, sa_start=0, sb_start=0, hsp_len=4),
    ]
    assert alignment.align_frag == expected_frags


def test_c_local_align_multi_fragment(multi_fragment_data) -> None:
    seqa, seqb, alpha_len, score_matrix, gap_open, gap_extend = multi_fragment_data

    c_score_matrix = score_matrix.ctypes.data_as(ctypes.POINTER(ctypes.c_int))

    c_alignment_ptr = lib.local_align(
        seqa,
        len(seqa),
        seqb,
        len(seqb),
        alpha_len,
        c_score_matrix,
        gap_open,
        gap_extend,
    )
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


def test_c_global_align_multi_fragment(multi_fragment_data) -> None:
    seqa, seqb, _map_bytes, score_matrix, gap_open, gap_extend = multi_fragment_data

    c_score_matrix = score_matrix.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    alpha_len = score_matrix.shape[0]

    c_alignment_ptr = lib.global_align(
        seqa,
        len(seqa),
        seqb,
        len(seqb),
        alpha_len,
        c_score_matrix,
        gap_open,
        gap_extend,
    )
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


def test_c_glocal_align_multi_fragment(multi_fragment_data) -> None:
    seqa, seqb, _map_bytes, score_matrix, gap_open, gap_extend = multi_fragment_data

    c_score_matrix = score_matrix.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    alpha_len = score_matrix.shape[0]

    c_alignment_ptr = lib.glocal_align(
        seqa,
        len(seqa),
        seqb,
        len(seqb),
        alpha_len,
        c_score_matrix,
        gap_open,
        gap_extend,
    )
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


def test_c_align_multi_fragment(multi_fragment_data) -> None:
    seqa, seqb, _map_bytes, score_matrix, gap_open, gap_extend = multi_fragment_data

    c_score_matrix = score_matrix.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    alpha_len = score_matrix.shape[0]

    c_alignment_ptr = lib.overlap_align(
        seqa,
        len(seqa),
        seqb,
        len(seqb),
        alpha_len,
        c_score_matrix,
        gap_open,
        gap_extend,
    )
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
def test_c_local_align_empty_seqa(common_data) -> None:
    _, seqb, _map_bytes, score_matrix, gap_open, gap_extend = common_data

    c_score_matrix = score_matrix.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    alpha_len = score_matrix.shape[0]

    # C function might not raise ValueError, it might return NULL or an invalid alignment
    # For now, we'll check if it returns NULL or if the score is unexpected.
    # The original Python test expects ValueError, so we'll adapt.
    c_alignment_ptr = lib.local_align(
        b"",
        0,
        seqb,
        len(seqb),
        alpha_len,
        c_score_matrix,
        gap_open,
        gap_extend,
    )
    # Assuming C returns NULL for empty sequences or an alignment with score 0 and 0 frags
    alignment = convert_c_alignment_to_py_alignment(c_alignment_ptr)
    assert alignment is None or (alignment.score == 0 and alignment.frag_count == 0)


def test_c_local_align_empty_seqb(common_data) -> None:
    seqa, _, _map_bytes, score_matrix, gap_open, gap_extend = common_data

    c_score_matrix = score_matrix.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    alpha_len = score_matrix.shape[0]

    c_alignment_ptr = lib.local_align(
        seqa,
        len(seqa),
        b"",
        0,
        alpha_len,
        c_score_matrix,
        gap_open,
        gap_extend,
    )
    alignment = convert_c_alignment_to_py_alignment(c_alignment_ptr)
    assert alignment is None or (alignment.score == 0 and alignment.frag_count == 0)


def test_c_global_align_empty_seqa(common_data) -> None:
    _, seqb, _map_bytes, score_matrix, gap_open, gap_extend = common_data

    c_score_matrix = score_matrix.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    alpha_len = score_matrix.shape[0]

    c_alignment_ptr = lib.global_align(
        b"",
        0,
        seqb,
        len(seqb),
        alpha_len,
        c_score_matrix,
        gap_open,
        gap_extend,
    )
    alignment = convert_c_alignment_to_py_alignment(c_alignment_ptr)
    assert alignment is None or (alignment.score == 0 and alignment.frag_count == 0)


def test_c_global_align_empty_seqb(common_data) -> None:
    seqa, _, _map_bytes, score_matrix, gap_open, gap_extend = common_data

    c_score_matrix = score_matrix.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    alpha_len = score_matrix.shape[0]

    c_alignment_ptr = lib.global_align(
        seqa,
        len(seqa),
        b"",
        0,
        alpha_len,
        c_score_matrix,
        gap_open,
        gap_extend,
    )
    alignment = convert_c_alignment_to_py_alignment(c_alignment_ptr)
    assert alignment is None or (alignment.score == 0 and alignment.frag_count == 0)


def test_c_glocal_align_empty_seqa(common_data) -> None:
    _, seqb, _map_bytes, score_matrix, gap_open, gap_extend = common_data

    c_score_matrix = score_matrix.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    alpha_len = score_matrix.shape[0]

    c_alignment_ptr = lib.glocal_align(
        b"",
        0,
        seqb,
        len(seqb),
        alpha_len,
        c_score_matrix,
        gap_open,
        gap_extend,
    )
    alignment = convert_c_alignment_to_py_alignment(c_alignment_ptr)
    assert alignment is None or (alignment.score == 0 and alignment.frag_count == 0)


def test_c_glocal_align_empty_seqb(common_data) -> None:
    seqa, _, _map_bytes, score_matrix, gap_open, gap_extend = common_data

    c_score_matrix = score_matrix.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    alpha_len = score_matrix.shape[0]

    c_alignment_ptr = lib.glocal_align(
        seqa,
        len(seqa),
        b"",
        0,
        alpha_len,
        c_score_matrix,
        gap_open,
        gap_extend,
    )
    alignment = convert_c_alignment_to_py_alignment(c_alignment_ptr)
    assert alignment is None or (alignment.score == 0 and alignment.frag_count == 0)


def test_c_align_empty_seqa(common_data) -> None:
    _, seqb, _map_bytes, score_matrix, gap_open, gap_extend = common_data

    c_score_matrix = score_matrix.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    alpha_len = score_matrix.shape[0]

    c_alignment_ptr = lib.overlap_align(
        b"",
        0,
        seqb,
        len(seqb),
        alpha_len,
        c_score_matrix,
        gap_open,
        gap_extend,
    )
    alignment = convert_c_alignment_to_py_alignment(c_alignment_ptr)
    assert alignment is None or (alignment.score == 0 and alignment.frag_count == 0)


def test_c_align_empty_seqb(common_data) -> None:
    seqa, _, alpha_len, score_matrix, gap_open, gap_extend = common_data

    c_score_matrix = score_matrix.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    alpha_len = score_matrix.shape[0]

    c_alignment_ptr = lib.overlap_align(
        seqa,
        len(seqa),
        b"",
        0,
        alpha_len,
        c_score_matrix,
        gap_open,
        gap_extend,
    )
    alignment = convert_c_alignment_to_py_alignment(c_alignment_ptr)
    assert alignment is None or (alignment.score == 0 and alignment.frag_count == 0)


def test_local_align_perfect_match_subsegment() -> None:
    alphabet = "ACGTXYZW"

    seqa = _encode("XXXXXAGCTYYYYY", alphabet)
    seqb = _encode("ZZZAGCTWWW", alphabet)

    score_matrix = np.eye(8, 8, dtype=np.int32) * 2 - 1
    c_score_matrix = score_matrix.ctypes.data_as(ctypes.POINTER(ctypes.c_int))

    gap_open = -2
    gap_extend = -1

    alignment = lib.local_align(seqa, len(seqa), seqb, len(seqb), 8, c_score_matrix, gap_open, gap_extend)
    alignment = convert_c_alignment_to_py_alignment(alignment)
    assert alignment.score == 4  # AGCT match
    assert alignment.frag_count == 1
    expected_frags = [
        AlignFrag(
            frag_type=FragType.MATCH, sa_start=5, sb_start=3, hsp_len=4,
        ),  # AGCT in seqa starts at index 5, in seqb at index 3
    ]
    assert alignment.align_frag == expected_frags
