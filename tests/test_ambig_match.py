import numpy as np
import pytest

import seq_smith


def test_make_score_matrix_basic() -> None:
    alphabet = "ACGT"
    match_score = 10
    mismatch_score = -5
    expected_matrix = np.array([[10, -5, -5, -5], [-5, 10, -5, -5], [-5, -5, 10, -5], [-5, -5, -5, 10]], dtype=np.int32)
    assert np.array_equal(seq_smith.make_score_matrix(alphabet, match_score, mismatch_score), expected_matrix)


def test_make_score_matrix_with_ambiguous_char() -> None:
    alphabet = "ACGTN"
    match_score = 10
    mismatch_score = -5
    ambiguous = "N"
    ambiguous_match_score = 0
    expected_matrix = np.array(
        [[10, -5, -5, -5, 0], [-5, 10, -5, -5, 0], [-5, -5, 10, -5, 0], [-5, -5, -5, 10, 0], [0, 0, 0, 0, 0]],
        dtype=np.int32,
    )
    result_matrix = seq_smith.make_score_matrix(
        alphabet,
        match_score,
        mismatch_score,
        ambiguous=ambiguous,
        ambiguous_match_score=ambiguous_match_score,
    )
    assert np.array_equal(result_matrix, expected_matrix)


def test_make_score_matrix_with_multiple_ambiguous_chars() -> None:
    alphabet = "ACGTXY"
    match_score = 5
    mismatch_score = -2
    ambiguous = "XY"
    ambiguous_match_score = 1
    expected_matrix = np.array(
        [
            [5, -2, -2, -2, 1, 1],
            [-2, 5, -2, -2, 1, 1],
            [-2, -2, 5, -2, 1, 1],
            [-2, -2, -2, 5, 1, 1],
            [1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 1],
        ],
        dtype=np.int32,
    )
    result_matrix = seq_smith.make_score_matrix(
        alphabet,
        match_score,
        mismatch_score,
        ambiguous=ambiguous,
        ambiguous_match_score=ambiguous_match_score,
    )
    assert np.array_equal(result_matrix, expected_matrix)


def test_make_score_matrix_ambiguous_no_ambiguous_match_score_raises_error() -> None:
    alphabet = "ACGTN"
    match_score = 10
    mismatch_score = -5
    ambiguous = "N"
    with pytest.raises(ValueError, match="ambiguous_match_score and ambiguous must both be provided"):
        seq_smith.make_score_matrix(alphabet, match_score, mismatch_score, ambiguous=ambiguous)


def test_make_score_matrix_ambiguous_char_not_in_alphabet_raises_error() -> None:
    alphabet = "ACGT"
    match_score = 10
    mismatch_score = -5
    ambiguous = "N"
    ambiguous_match_score = 0
    with pytest.raises(ValueError, match="all ambiguous characters must be included in the alphabet"):
        seq_smith.make_score_matrix(
            alphabet,
            match_score,
            mismatch_score,
            ambiguous=ambiguous,
            ambiguous_match_score=ambiguous_match_score,
        )


def test_make_score_matrix_empty_ambiguous_string() -> None:
    alphabet = "ACGT"
    match_score = 10
    mismatch_score = -5
    ambiguous = ""
    ambiguous_match_score = 0  # This value should be ignored
    expected_matrix = np.array([[10, -5, -5, -5], [-5, 10, -5, -5], [-5, -5, 10, -5], [-5, -5, -5, 10]], dtype=np.int32)
    assert np.array_equal(
        seq_smith.make_score_matrix(
            alphabet,
            match_score,
            mismatch_score,
            ambiguous=ambiguous,
            ambiguous_match_score=ambiguous_match_score,
        ),
        expected_matrix,
    )


def test_make_score_matrix_ambiguous_is_entire_alphabet() -> None:
    alphabet = "ACGT"
    match_score = 10
    mismatch_score = -5
    ambiguous = "ACGT"
    ambiguous_match_score = 5
    expected_matrix = np.array([[5, 5, 5, 5], [5, 5, 5, 5], [5, 5, 5, 5], [5, 5, 5, 5]], dtype=np.int32)
    result_matrix = seq_smith.make_score_matrix(
        alphabet,
        match_score,
        mismatch_score,
        ambiguous=ambiguous,
        ambiguous_match_score=ambiguous_match_score,
    )
    assert np.array_equal(result_matrix, expected_matrix)
