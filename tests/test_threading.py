from collections.abc import Callable

import pytest
from conftest import AlignmentDataMany

from seq_smith import (
    global_align,
    global_align_many,
    local_align,
    local_align_many,
    local_global_align,
    local_global_align_many,
    overlap_align,
    overlap_align_many,
)


@pytest.mark.parametrize(
    ("align_func", "align_many_func"),
    [
        (global_align, global_align_many),
        (local_align, local_align_many),
        (local_global_align, local_global_align_many),
        (overlap_align, overlap_align_many),
    ],
)
def test_align_many_vs_single(
    threading_data: AlignmentDataMany,
    align_func: Callable,
    align_many_func: Callable,
) -> None:
    seqa = threading_data.seqa
    seqbs = threading_data.seqbs
    score_matrix = threading_data.score_matrix
    gap_open = threading_data.gap_open
    gap_extend = threading_data.gap_extend

    # Single threaded loop
    expected = []
    for seqb in seqbs:
        expected.append(align_func(seqa, seqb, score_matrix, gap_open, gap_extend))

    # Multi-threaded
    actual = align_many_func(seqa, seqbs, score_matrix, gap_open, gap_extend)

    assert len(actual) == len(expected)
    for a, e in zip(actual, expected, strict=True):
        assert a.score == e.score
        assert a.fragments == e.fragments


def test_threading_num_threads(threading_data: AlignmentDataMany) -> None:
    actual = global_align_many(
        threading_data.seqa,
        threading_data.seqbs,
        threading_data.score_matrix,
        threading_data.gap_open,
        threading_data.gap_extend,
        num_threads=2,
    )
    assert len(actual) == len(threading_data.seqbs)


def test_empty_input(threading_data: AlignmentDataMany) -> None:
    actual = global_align_many(
        threading_data.seqa,
        [],
        threading_data.score_matrix,
        threading_data.gap_open,
        threading_data.gap_extend,
    )
    assert actual == []
