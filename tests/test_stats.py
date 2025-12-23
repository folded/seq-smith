import numpy as np
from conftest import AlignmentData

from seq_smith import encode, global_align


def test_stats_exact_match(common_data: AlignmentData) -> None:
    seqa = encode("ACGT", common_data.alphabet)
    seqb = encode("ACGT", common_data.alphabet)
    alignment = global_align(
        seqa,
        seqb,
        common_data.score_matrix,
        common_data.gap_open,
        common_data.gap_extend,
    )

    assert alignment.stats.num_exact_matches == 4
    assert alignment.stats.num_positive_mismatches == 0
    assert alignment.stats.num_negative_mismatches == 0
    assert alignment.stats.num_a_gaps == 0
    assert alignment.stats.num_b_gaps == 0
    assert alignment.stats.len == 4


def test_stats_mismatches() -> None:
    alphabet = "ABC"
    score_matrix = np.array([[2, 1, -1], [1, 2, -1], [-1, -1, 2]], dtype=np.int32)

    seqa = encode("AAA", alphabet)
    seqb = encode("ABC", alphabet)

    alignment = global_align(seqa, seqb, score_matrix, -10, -10)  # High gap penalty to force alignment

    assert alignment.stats.num_exact_matches == 1  # A-A
    assert alignment.stats.num_positive_mismatches == 1  # A-B
    assert alignment.stats.num_negative_mismatches == 1  # A-C
    assert alignment.stats.num_a_gaps == 0
    assert alignment.stats.num_b_gaps == 0
    assert alignment.stats.len == 3


def test_stats_gaps(common_data: AlignmentData) -> None:
    seqa = encode("A", common_data.alphabet)
    seqb = encode("AAA", common_data.alphabet)

    # A-- vs AAA (or --A or -A-)
    alignment = global_align(
        seqa,
        seqb,
        common_data.score_matrix,
        common_data.gap_open,
        common_data.gap_extend,
    )

    assert alignment.stats.num_exact_matches == 1
    assert alignment.stats.num_a_gaps == 2
    assert alignment.stats.num_b_gaps == 0
    assert alignment.stats.len == 3


def test_stats_gap_b(common_data: AlignmentData) -> None:
    seqa = encode("AAA", common_data.alphabet)
    seqb = encode("A", common_data.alphabet)

    alignment = global_align(
        seqa,
        seqb,
        common_data.score_matrix,
        common_data.gap_open,
        common_data.gap_extend,
    )

    assert alignment.stats.num_exact_matches == 1
    assert alignment.stats.num_b_gaps == 2
    assert alignment.stats.num_a_gaps == 0
    assert alignment.stats.len == 3
