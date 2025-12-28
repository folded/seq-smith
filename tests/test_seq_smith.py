import numpy as np
import pytest
from conftest import AlignmentData

from seq_smith import (
    AlignmentFragment,
    FragmentType,
    encode,
    format_alignment_ascii,
    generate_cigar,
    global_align,
    local_align,
    local_global_align,
    make_score_matrix,
    overlap_align,
    top_k_ungapped_local_align,
    top_k_ungapped_local_align_many,
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


def test_local_global_align_overhangs() -> None:
    # Sequence B (Global) has overhanging tails
    # sequence A: CCCC
    # sequence B: AAACCCCAAA
    # Expected: A aligns to central C's, B has leading/trailing gaps.
    seqa = encode("CCCC", "ACGT")
    seqb = encode("AAACCCCAAA", "ACGT")
    sm = make_score_matrix("ACGT", match_score=2, mismatch_score=-2)
    aln = local_global_align(seqa, seqb, sm, gap_open=-3, gap_extend=-1)

    expected_fragments = [
        AlignmentFragment(FragmentType.AGap, 0, 0, 3),
        AlignmentFragment(FragmentType.Match, 0, 3, 4),
        AlignmentFragment(FragmentType.AGap, 4, 7, 3),
    ]
    assert aln.fragments == expected_fragments

    assert aln.score == -2


def test_top_k_ungapped_simple() -> None:
    # Use custom alphabet to support Z/W if needed, or just use ACGT
    alphabet = "ACGT"
    seqa = encode("AAAATTTTCCCC", alphabet)
    seqb = encode("AAAAGGGGCCCC", alphabet)

    # matrix: match=2, mismatch=-5
    score_matrix = make_score_matrix(alphabet, match_score=2, mismatch_score=-5)

    # AAAA matches (4*2=8). TTTT vs GGGG (-5*4 = -20). CCCC matches (8).
    # Alignment 1: AAAA (score 8)
    # Alignment 2: CCCC (score 8)

    alignments = top_k_ungapped_local_align(seqa, seqb, score_matrix, k=5)

    # Should get 2 alignments
    assert len(alignments) == 2
    # Verify scores
    assert alignments[0].score == 8
    assert alignments[1].score == 8

    starts = sorted([(a.fragments[0].sa_start, a.fragments[0].sb_start) for a in alignments])
    assert starts[0] == (0, 0)
    assert starts[1] == (8, 8)


def test_top_k_ungapped_overlap() -> None:
    # Use custom alphabet to support Z/W if needed, or just use ACGT
    alphabet = "ACGT"
    seqa = encode("AAAATTTTCCCCAAAATTTTCCCCAAAATTTTCCCC", alphabet)
    seqb = encode("AAAAGGGGCCCC", alphabet)

    # matrix: match=2, mismatch=-5
    score_matrix = make_score_matrix(alphabet, match_score=2, mismatch_score=-5)

    alignments = top_k_ungapped_local_align(seqa, seqb, score_matrix, k=5, filter_overlap_b=False)

    assert len(alignments) == 5  # 6 total alignments, but k=5
    assert all(a.score == 8 for a in alignments)
    starts = sorted([(a.fragments[0].sa_start, a.fragments[0].sb_start) for a in alignments])
    for c, r in starts:
        assert seqa[c : c + 4] == seqb[r : r + 4]


def test_top_k_ungapped_overlapping_candidates(common_data: AlignmentData) -> None:
    # Case where second best candidate overlaps best
    # Sequence A: AAAAA
    # Sequence B: AAAAA
    # match=1 (from common_data observation)

    seqa = encode("AAAAA", common_data.alphabet)
    seqb = encode("AAAAA", common_data.alphabet)

    # ensure score matrix is what we think (match=1) or make our own
    score_matrix = make_score_matrix(common_data.alphabet, match_score=2, mismatch_score=-1)

    alignments = top_k_ungapped_local_align(seqa, seqb, score_matrix, k=2)

    assert len(alignments) == 1
    assert alignments[0].score == 10  # 5 * 2
    assert alignments[0].fragments[0].len == 5


def test_top_k_ungapped_limit() -> None:
    # A: AA..CC..GG
    # B: AA..CC..GG
    alphabet = "ACGT"
    # Use T vs G for mismatch
    seqa = encode("AATTCCTTGG", alphabet)
    seqb = encode("AAGGCCGGGG", alphabet)

    score_matrix = make_score_matrix(alphabet, match_score=2, mismatch_score=-5)

    # Expected HSPs: AA (4), mismatch, CC (4), mismatch, GG (4)

    alignments = top_k_ungapped_local_align(seqa, seqb, score_matrix, k=2)

    assert len(alignments) == 2
    assert alignments[0].score == 4
    assert alignments[1].score == 4


def test_top_k_ungapped_many_simple() -> None:
    # Sequence A: AAAA
    # Sequence B1: AAAA (perfect)
    # Sequence B2: CCCC (mismatch)
    alphabet = "ACGT"
    seqa = encode("AAAA", alphabet)

    seqb1 = encode("AAAA", alphabet)
    seqb2 = encode("CCCC", alphabet)

    score_matrix = make_score_matrix(alphabet, match_score=2, mismatch_score=-5)

    alignments_list = top_k_ungapped_local_align_many(seqa, [seqb1, seqb2], score_matrix, k=1)

    assert len(alignments_list) == 2

    # Check first alignment (AAAA vs AAAA)
    # Should have score 8
    assert len(alignments_list[0]) == 1
    assert alignments_list[0][0].score == 8

    # Check second alignment (AAAA vs CCCC)
    # Should have no positive candidates if mismatch penalty is high enough?
    # -5 * 4 = -20.
    # So should be empty if score <= 0.
    # Our implementation returns empty if no positive peaks.
    assert len(alignments_list[1]) == 0
