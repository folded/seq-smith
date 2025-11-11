import ctypes
from collections.abc import Callable
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import TYPE_CHECKING, TypeAlias

import numpy as np


# Define C structures
class FragType(Enum):
    A_GAP = 0
    B_GAP = 1
    MATCH = 2


class CAlignFrag(ctypes.Structure):
    pass  # Forward declaration for self-referential pointer


CAlignFrag._fields_ = [
    ("next", ctypes.POINTER(CAlignFrag)),
    ("type", ctypes.c_int),
    ("sa_start", ctypes.c_int),
    ("sb_start", ctypes.c_int),
    ("hsp_len", ctypes.c_int),
]


class CAlignment(ctypes.Structure):
    pass


CAlignment._fields_ = [
    ("align_frag", ctypes.POINTER(CAlignFrag)),
    ("frag_count", ctypes.c_int),
    ("score", ctypes.c_int),
]


if TYPE_CHECKING:
    CAlignmentPtr: TypeAlias = ctypes.POINTER[CAlignment]  # type: ignore[valid-type]
else:
    CAlignmentPtr = ctypes.POINTER(CAlignment)


UBytePtr = ctypes.POINTER(ctypes.c_char)


@dataclass
class AlignParams:
    seqa: bytes
    seqb: bytes
    alpha_len: int
    score_matrix: np.ndarray
    gap_open: int
    gap_extend: int


# Load the shared library
lib = ctypes.CDLL(str(Path(__file__).parent / "libalign.so"))

# Define argtypes and restype for C functions
lib.local_align.argtypes = [
    UBytePtr,
    ctypes.c_int,
    UBytePtr,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.POINTER(ctypes.c_int),
    ctypes.c_int,
    ctypes.c_int,
]
lib.local_align.restype = ctypes.POINTER(CAlignment)

lib.global_align.argtypes = [
    UBytePtr,
    ctypes.c_int,
    UBytePtr,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.POINTER(ctypes.c_int),
    ctypes.c_int,
    ctypes.c_int,
]
lib.global_align.restype = ctypes.POINTER(CAlignment)

lib.glocal_align.argtypes = [
    UBytePtr,
    ctypes.c_int,
    UBytePtr,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.POINTER(ctypes.c_int),
    ctypes.c_int,
    ctypes.c_int,
]
lib.glocal_align.restype = ctypes.POINTER(CAlignment)

lib.overlap_align.argtypes = [
    UBytePtr,
    ctypes.c_int,
    UBytePtr,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.POINTER(ctypes.c_int),
    ctypes.c_int,
    ctypes.c_int,
]
lib.overlap_align.restype = ctypes.POINTER(CAlignment)

# Helper to free C Alignment structure
lib.alignment_free.argtypes = [ctypes.POINTER(CAlignment)]
lib.alignment_free.restype = None


def _align_wrapper(
    func: Callable[[bytes, int, bytes, int, int, np.ndarray, int, int], CAlignmentPtr],
    params: AlignParams,
) -> CAlignmentPtr:
    c_score_matrix = (ctypes.c_int * (params.alpha_len * params.alpha_len))(*params.score_matrix.flatten())

    return func(
        params.seqa,
        len(params.seqa),
        params.seqb,
        len(params.seqb),
        params.alpha_len,
        c_score_matrix,
        params.gap_open,
        params.gap_extend,
    )


def local_align(params: AlignParams) -> CAlignmentPtr:
    return _align_wrapper(lib.local_align, params)


def global_align(params: AlignParams) -> CAlignmentPtr:
    return _align_wrapper(lib.global_align, params)


def glocal_align(params: AlignParams) -> CAlignmentPtr:
    return _align_wrapper(lib.glocal_align, params)


def overlap_align(params: AlignParams) -> CAlignmentPtr:
    return _align_wrapper(lib.overlap_align, params)
