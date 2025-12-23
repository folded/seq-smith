import dataclasses

import numpy as np
import pytest

from seq_smith import encode, make_score_matrix


@dataclasses.dataclass
class AlignmentData:
    alphabet: str
    seqa: bytes
    seqb: bytes
    score_matrix: np.ndarray
    gap_open: int
    gap_extend: int


@pytest.fixture
def common_data() -> AlignmentData:
    alphabet = "ACGT"
    seqa = encode("ACGT", alphabet)
    seqb = encode("ACGT", alphabet)

    score_matrix = make_score_matrix(alphabet, 1, -1)
    gap_open = -2
    gap_extend = -1

    return AlignmentData(
        alphabet=alphabet,
        seqa=seqa,
        seqb=seqb,
        score_matrix=score_matrix,
        gap_open=gap_open,
        gap_extend=gap_extend,
    )


@pytest.fixture
def complex_data() -> AlignmentData:
    alphabet = "ACGT"
    seqa = encode("GAATTCAGTTA", alphabet)
    seqb = encode("GGATCGA", alphabet)

    score_matrix = make_score_matrix(alphabet, 2, -1)
    gap_open = -3
    gap_extend = -1
    return AlignmentData(
        alphabet=alphabet,
        seqa=seqa,
        seqb=seqb,
        score_matrix=score_matrix,
        gap_open=gap_open,
        gap_extend=gap_extend,
    )


@pytest.fixture
def local_global_test_data() -> AlignmentData:
    alphabet = "ACGTX"
    seqa = encode("XACGTX", alphabet)
    seqb = encode("ACGT", alphabet)

    score_matrix = make_score_matrix(alphabet, 1, -1)
    score_matrix[4, 4] = 0  # X vs X
    gap_open = -2
    gap_extend = -1
    return AlignmentData(
        alphabet=alphabet,
        seqa=seqa,
        seqb=seqb,
        score_matrix=score_matrix,
        gap_open=gap_open,
        gap_extend=gap_extend,
    )


@pytest.fixture
def multi_fragment_data() -> AlignmentData:
    alphabet = "ACGT"
    seqa = encode("AGAGAGAGAG", alphabet)
    seqb = encode("AGCAGCAGCA", alphabet)

    score_matrix = make_score_matrix(alphabet, 1, -1)

    gap_open = -1
    gap_extend = -1
    return AlignmentData(
        alphabet=alphabet,
        seqa=seqa,
        seqb=seqb,
        score_matrix=score_matrix,
        gap_open=gap_open,
        gap_extend=gap_extend,
    )


@pytest.fixture
def poly_data() -> AlignmentData:
    alphabet = "ACGT"
    seqa = encode("CCCCCCAACAA", alphabet)
    seqb = encode("TTAAAAGGGGGGG", alphabet)

    score_matrix = make_score_matrix(alphabet, 1, -1)

    gap_open = -2
    gap_extend = -1
    return AlignmentData(
        alphabet=alphabet,
        seqa=seqa,
        seqb=seqb,
        score_matrix=score_matrix,
        gap_open=gap_open,
        gap_extend=gap_extend,
    )


@pytest.fixture
def poly_data_strong_gap_penalty() -> AlignmentData:
    alphabet = "ACGT"
    seqa = encode("CCCCCCAACAACCCCCCC", alphabet)
    seqb = encode("TTAAAAGGGG", alphabet)

    score_matrix = make_score_matrix(alphabet, 1, -1)

    gap_open = -100
    gap_extend = -100
    return AlignmentData(
        alphabet=alphabet,
        seqa=seqa,
        seqb=seqb,
        score_matrix=score_matrix,
        gap_open=gap_open,
        gap_extend=gap_extend,
    )


@dataclasses.dataclass
class AlignmentDataMany:
    alphabet: str
    seqa: bytes
    seqbs: list[bytes]
    score_matrix: np.ndarray
    gap_open: int
    gap_extend: int


@pytest.fixture
def threading_data() -> AlignmentDataMany:
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
    return AlignmentDataMany(
        alphabet=alphabet,
        seqa=seqa,
        seqbs=seqbs,
        score_matrix=score_matrix,
        gap_open=gap_open,
        gap_extend=gap_extend,
    )
