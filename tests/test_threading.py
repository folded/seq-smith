import numpy as np
import pytest
from seq_smith import (
    encode,
    global_align,
    global_align_many,
    local_align,
    local_align_many,
    local_global_align,
    local_global_align_many,
    make_score_matrix,
    overlap_align,
    overlap_align_many,
)


@pytest.fixture
def alignment_data():
    alphabet = "ACGT"
    seqa = encode("ACGTACGT", alphabet)
    seqbs = [
        encode("ACGTACGT", alphabet),
        encode("ACGT", alphabet),
        encode("CGTA", alphabet),
        encode("AAAA", alphabet),
        encode("GGGG", alphabet),
    ]
    score_matrix = make_score_matrix(alphabet, 1, -1)
    gap_open = -2
    gap_extend = -1
    return seqa, seqbs, score_matrix, gap_open, gap_extend


@pytest.mark.parametrize(
    "align_func, align_many_func",
    [
        (global_align, global_align_many),
        (local_align, local_align_many),
        (local_global_align, local_global_align_many),
        (overlap_align, overlap_align_many),
    ],
)
def test_align_many_vs_single(alignment_data, align_func, align_many_func):
    seqa, seqbs, score_matrix, gap_open, gap_extend = alignment_data

    # Single threaded loop
    expected = []
    for seqb in seqbs:
        expected.append(align_func(seqa, seqb, score_matrix, gap_open, gap_extend))

    # Multi-threaded
    actual = align_many_func(seqa, seqbs, score_matrix, gap_open, gap_extend)

    assert len(actual) == len(expected)
    for a, e in zip(actual, expected):
        assert a.score == e.score
        assert a.fragments == e.fragments


def test_threading_num_threads(alignment_data):
    seqa, seqbs, score_matrix, gap_open, gap_extend = alignment_data
    
    # Test with specific number of threads
    actual = global_align_many(seqa, seqbs, score_matrix, gap_open, gap_extend, num_threads=2)
    assert len(actual) == len(seqbs)


def test_empty_input(alignment_data):
    seqa, _, score_matrix, gap_open, gap_extend = alignment_data
    
    actual = global_align_many(seqa, [], score_matrix, gap_open, gap_extend)
    assert actual == []
