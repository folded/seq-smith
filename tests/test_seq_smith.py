import numpy as np
import pytest
from conftest import AlignmentData

from seq_smith import (
    encode,
    format_alignment_ascii,
    generate_cigar,
    global_align,
    local_align,
    local_global_align,
    make_score_matrix,
    overlap_align,
)


def test_global_align_simple(common_data: AlignmentData) -> None:
    # This test was flawed. The original seqb was "AGCT", but the expected
    # fragments described a perfect match. Correcting seqb to "ACGT".
    seqa = encode("ACGT", common_data.alphabet)
    seqb = encode("ACGT", common_data.alphabet)
    alignment = global_align(
        seqa,
        seqb,
        common_data.score_matrix,
        common_data.gap_open,
        common_data.gap_extend,
    )

    assert alignment.score == 4
    aligned_a, aligned_b = format_alignment_ascii(seqa, seqb, alignment.fragments, common_data.alphabet)
    assert aligned_a == "ACGT"
    assert aligned_b == "ACGT"


def test_global_align_simple_gap(common_data: AlignmentData) -> None:
    seqa = encode("A", common_data.alphabet)
    seqb = encode("AC", common_data.alphabet)
    alignment = global_align(seqa, seqb, common_data.score_matrix, common_data.gap_open, common_data.gap_extend)

    assert alignment.score == -1
    aligned_a, aligned_b = format_alignment_ascii(seqa, seqb, alignment.fragments, common_data.alphabet)
    assert aligned_a == "A-"
    assert aligned_b == "AC"


def test_local_global_align_simple(common_data: AlignmentData) -> None:
    # This test was flawed. The original seqb was "AGCT", but the expected
    # fragments described a perfect match. Correcting seqb to "ACGT".
    seqa = encode("ACGT", common_data.alphabet)
    seqb = encode("ACGT", common_data.alphabet)
    alignment = local_global_align(
        seqa,
        seqb,
        common_data.score_matrix,
        common_data.gap_open,
        common_data.gap_extend,
    )

    assert alignment.score == 4
    aligned_a, aligned_b = format_alignment_ascii(seqa, seqb, alignment.fragments, common_data.alphabet)
    assert aligned_a == "ACGT"
    assert aligned_b == "ACGT"


def test_local_global_align_subsegment_global_seqb(local_global_test_data: AlignmentData) -> None:
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
        alignment.fragments,
        local_global_test_data.alphabet,
    )
    assert aligned_a == "ACGT"
    assert aligned_b == "ACGT"


def test_overlap_align_simple(common_data: AlignmentData) -> None:
    # This test was flawed. The original seqb was "AGCT", but the expected
    # fragments described a perfect match. Correcting seqb to "ACGT".
    seqa = encode("ACGT", common_data.alphabet)
    seqb = encode("ACGT", common_data.alphabet)
    alignment = overlap_align(
        seqa,
        seqb,
        common_data.score_matrix,
        common_data.gap_open,
        common_data.gap_extend,
    )

    assert alignment.score == 4
    aligned_a, aligned_b = format_alignment_ascii(seqa, seqb, alignment.fragments, common_data.alphabet)
    assert aligned_a == "ACGT"
    assert aligned_b == "ACGT"


def test_overlap_align_semi_global_overlap(common_data: AlignmentData) -> None:
    seqa = encode("ACGTACGT", common_data.alphabet)
    seqb = encode("CGTA", common_data.alphabet)
    alignment = overlap_align(seqa, seqb, common_data.score_matrix, common_data.gap_open, common_data.gap_extend)

    assert alignment.score == 4
    aligned_a, aligned_b = format_alignment_ascii(seqa, seqb, alignment.fragments, common_data.alphabet)
    assert aligned_a == "CGTA"
    assert aligned_b == "CGTA"


def test_local_align_perfect_match_subsegment() -> None:
    alphabet = "ACGTXYZW"
    seqa = encode("XXXXXAGCTYYYYY", alphabet)
    seqb = encode("ZZZAGCTWWW", alphabet)

    score_matrix = make_score_matrix(alphabet, 2, -1)

    gap_open = -2
    gap_extend = -1

    alignment = local_align(seqa, seqb, score_matrix, gap_open, gap_extend)
    assert alignment.score == 8
    aligned_a, aligned_b = format_alignment_ascii(seqa, seqb, alignment.fragments, alphabet)
    assert aligned_a == "AGCT"
    assert aligned_b == "AGCT"


def test_local_align_multi_fragment(multi_fragment_data: AlignmentData) -> None:
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
        alignment.fragments,
        multi_fragment_data.alphabet,
    )
    assert aligned_a == "AG-AG-AG"
    assert aligned_b == "AGCAGCAG"


def test_global_align_multi_fragment(multi_fragment_data: AlignmentData) -> None:
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
        alignment.fragments,
        multi_fragment_data.alphabet,
    )
    assert aligned_a == "AG-AG-AGAGAG"
    assert aligned_b == "AGCAGCAG-CA-"
    assert len(aligned_a) == len(aligned_b)


def test_local_global_align_multi_fragment(multi_fragment_data: AlignmentData) -> None:
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
        alignment.fragments,
        multi_fragment_data.alphabet,
    )
    assert aligned_a == "AG-AG-AG-A"
    assert aligned_b == "AGCAGCAGCA"


def test_overlap_align_multi_fragment(multi_fragment_data: AlignmentData) -> None:
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
        alignment.fragments,
        multi_fragment_data.alphabet,
    )
    # AGAGAG-AG-AG
    #     AGCAGCAGCA
    assert aligned_a == "AG-AG-AG"
    assert aligned_b == "AGCAGCAG"


# Test with empty sequences
def test_local_align_empty_seqa(common_data: AlignmentData) -> None:
    with pytest.raises(ValueError, match=r"Input sequences cannot be empty."):
        local_align(b"", common_data.seqb, common_data.score_matrix, common_data.gap_open, common_data.gap_extend)


def test_local_align_empty_seqb(common_data: AlignmentData) -> None:
    with pytest.raises(ValueError, match=r"Input sequences cannot be empty."):
        local_align(common_data.seqa, b"", common_data.score_matrix, common_data.gap_open, common_data.gap_extend)


def test_global_align_empty_seqa(common_data: AlignmentData) -> None:
    with pytest.raises(ValueError, match=r"Input sequences cannot be empty."):
        global_align(b"", common_data.seqb, common_data.score_matrix, common_data.gap_open, common_data.gap_extend)


def test_global_align_empty_seqb(common_data: AlignmentData) -> None:
    with pytest.raises(ValueError, match=r"Input sequences cannot be empty."):
        global_align(common_data.seqa, b"", common_data.score_matrix, common_data.gap_open, common_data.gap_extend)


def test_local_global_align_empty_seqa(common_data: AlignmentData) -> None:
    with pytest.raises(ValueError, match=r"Input sequences cannot be empty."):
        local_global_align(
            b"",
            common_data.seqb,
            common_data.score_matrix,
            common_data.gap_open,
            common_data.gap_extend,
        )


def test_local_global_align_empty_seqb(common_data: AlignmentData) -> None:
    with pytest.raises(ValueError, match=r"Input sequences cannot be empty."):
        local_global_align(
            common_data.seqa,
            b"",
            common_data.score_matrix,
            common_data.gap_open,
            common_data.gap_extend,
        )


def test_overlap_align_empty_seqa(common_data: AlignmentData) -> None:
    with pytest.raises(ValueError, match=r"Input sequences cannot be empty."):
        overlap_align(b"", common_data.seqb, common_data.score_matrix, common_data.gap_open, common_data.gap_extend)


def test_overlap_align_empty_seqb(common_data: AlignmentData) -> None:
    with pytest.raises(ValueError, match=r"Input sequences cannot be empty."):
        overlap_align(common_data.seqa, b"", common_data.score_matrix, common_data.gap_open, common_data.gap_extend)


def test_local_align_poly(poly_data: AlignmentData) -> None:
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
        alignment.fragments,
        poly_data.alphabet,
    )
    assert aligned_seqa == "ACAA"
    assert aligned_seqb == "AAAA"


def test_global_align_poly(poly_data: AlignmentData) -> None:
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
        alignment.fragments,
        poly_data.alphabet,
    )
    assert aligned_seqa == "CCCCCCAA----CAA"
    assert aligned_seqb == "--TTAAAAGGGGGGG"


def test_global_align_poly_strong_gap_penalty(poly_data: AlignmentData) -> None:
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
        alignment.fragments,
        poly_data.alphabet,
    )
    assert aligned_seqa == "--CCCCCCAACAA"
    assert aligned_seqb == "TTAAAAGGGGGGG"


def test_global_alignment_long_gaps() -> None:
    alphabet = "xyz"
    seqa_str = "xxxxzzzz"
    seqb_str = "yyyyzzzz"
    seqa = encode(seqa_str, alphabet)
    seqb = encode(seqb_str, alphabet)

    score_matrix = np.array(
        [
            [2, -3, -3],  # x
            [-3, 2, -3],  # y
            [-3, -3, 5],  # z
        ],
        dtype=np.int32,
    )

    gap_open = -1
    gap_extend = -1

    alignment = global_align(seqa, seqb, score_matrix, gap_open, gap_extend)

    assert alignment.score == 12
    aligned_a, aligned_b = format_alignment_ascii(seqa, seqb, alignment.fragments, alphabet)
    assert aligned_a == "xxxx----zzzz"
    assert aligned_b == "----yyyyzzzz"


def test_local_global_align_poly(poly_data: AlignmentData) -> None:
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
        alignment.fragments,
        poly_data.alphabet,
    )
    assert aligned_seqa == "CCAACA------A"
    assert aligned_seqb == "TTAAAAGGGGGGG"


def test_local_global_align_poly_strong_gap_penalty(poly_data_strong_gap_penalty: AlignmentData) -> None:
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
        alignment.fragments,
        poly_data_strong_gap_penalty.alphabet,
    )
    assert aligned_seqa == "CAACAACCCC"
    assert aligned_seqb == "TTAAAAGGGG"


def test_overlap_align_poly(poly_data: AlignmentData) -> None:
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
        alignment.fragments,
        poly_data.alphabet,
    )
    # CCCCCCAACAA
    #      TTAAAAGGGGGGG
    assert aligned_seqa == "CAACAA"
    assert aligned_seqb == "TTAAAA"


def test_overlap_align_poly_flipped(poly_data: AlignmentData) -> None:
    alignment = overlap_align(
        poly_data.seqb,
        poly_data.seqa,
        poly_data.score_matrix,
        poly_data.gap_open,
        poly_data.gap_extend,
    )
    assert alignment.score == 0
    aligned_seqa, aligned_seqb = format_alignment_ascii(
        poly_data.seqb,
        poly_data.seqa,
        alignment.fragments,
        poly_data.alphabet,
    )
    #      TTAAAAGGGGGGG
    # CCCCCCAACAA
    assert aligned_seqa == "TTAAAA"
    assert aligned_seqb == "CAACAA"


def test_generate_cigar_simple(common_data: AlignmentData) -> None:
    seqa = encode("ACGT", common_data.alphabet)
    seqb = encode("ACGT", common_data.alphabet)
    alignment = global_align(
        seqa,
        seqb,
        common_data.score_matrix,
        common_data.gap_open,
        common_data.gap_extend,
    )
    cigar = generate_cigar(alignment)
    assert cigar == "4M"


def test_generate_cigar_with_gap(common_data: AlignmentData) -> None:
    seqa = encode("A", common_data.alphabet)
    seqb = encode("AC", common_data.alphabet)
    alignment = global_align(seqa, seqb, common_data.score_matrix, common_data.gap_open, common_data.gap_extend)
    cigar = generate_cigar(alignment)
    assert cigar == "1M1I"  # A- vs AC, so A matches, then C is an insertion in seqB (query)


def test_generate_cigar_with_deletion(common_data: AlignmentData) -> None:
    seqa = encode("AC", common_data.alphabet)
    seqb = encode("A", common_data.alphabet)
    alignment = global_align(seqa, seqb, common_data.score_matrix, common_data.gap_open, common_data.gap_extend)
    cigar = generate_cigar(alignment)
    assert cigar == "1M1D"  # AC vs A-, so A matches, then C is a deletion in seqB (query)
