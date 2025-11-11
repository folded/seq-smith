import ctypes
from enum import Enum


# Define C structures
class FragType(Enum):
    A_GAP = 0
    B_GAP = 1
    MATCH = 2


class C_AlignFrag(ctypes.Structure):
    pass  # Forward declaration for self-referential pointer


C_AlignFrag._fields_ = [
    ("next", ctypes.POINTER(C_AlignFrag)),
    ("type", ctypes.c_int),
    ("sa_start", ctypes.c_int),
    ("sb_start", ctypes.c_int),
    ("hsp_len", ctypes.c_int),
]


class C_Alignment(ctypes.Structure):
    _fields_ = [
        ("align_frag", ctypes.POINTER(C_AlignFrag)),
        ("frag_count", ctypes.c_int),
        ("score", ctypes.c_int),
    ]


# Load the shared library
lib = ctypes.CDLL("/Users/tjs/repo/seq-align/libalign.so")

# Define argtypes and restype for C functions
# Alignment *local_align_raw(const unsigned char *sa, int sa_len, const unsigned char *sb, int sb_len, int alpha_len, int *score_matrix, int gap_open, int gap_extend);
lib.local_align.argtypes = [
    ctypes.c_char_p,
    ctypes.c_int,
    ctypes.c_char_p,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.POINTER(ctypes.c_int),
    ctypes.c_int,
    ctypes.c_int,
]
lib.local_align.restype = ctypes.POINTER(C_Alignment)

# Alignment *global_align_raw(const unsigned char *sa, int sa_len, const unsigned char *sb, int sb_len, int alpha_len, int *score_matrix, int gap_open, int gap_extend);
lib.global_align.argtypes = [
    ctypes.c_char_p,
    ctypes.c_int,
    ctypes.c_char_p,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.POINTER(ctypes.c_int),
    ctypes.c_int,
    ctypes.c_int,
]
lib.global_align.restype = ctypes.POINTER(C_Alignment)

# Alignment *glocal_align_raw(const unsigned char *sa, int sa_len, const unsigned char *sb, int sb_len, int alpha_len, int *score_matrix, int gap_open, int gap_extend);
lib.glocal_align.argtypes = [
    ctypes.c_char_p,
    ctypes.c_int,
    ctypes.c_char_p,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.POINTER(ctypes.c_int),
    ctypes.c_int,
    ctypes.c_int,
]
lib.glocal_align.restype = ctypes.POINTER(C_Alignment)

# Alignment *overlap_align_raw(const unsigned char *sa, int sa_len, const unsigned char *sb, int sb_len, int alpha_len, int *score_matrix, int gap_open, int gap_extend);
lib.overlap_align.argtypes = [
    ctypes.c_char_p,
    ctypes.c_int,
    ctypes.c_char_p,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.POINTER(ctypes.c_int),
    ctypes.c_int,
    ctypes.c_int,
]
lib.overlap_align.restype = ctypes.POINTER(C_Alignment)

# Helper to free C Alignment structure
lib.alignment_free.argtypes = [ctypes.POINTER(C_Alignment)]
lib.alignment_free.restype = None


def _align_wrapper(func, seqa: str, seqb: str, alpha_len: int, score_matrix, gap_open: int, gap_extend: int):
    sa_bytes = seqa.encode("ascii")
    sb_bytes = seqb.encode("ascii")

    c_score_matrix = (ctypes.c_int * (alpha_len * alpha_len))(*score_matrix.flatten())

    return func(
        sa_bytes,
        len(sa_bytes),
        sb_bytes,
        len(sb_bytes),
        alpha_len,
        c_score_matrix,
        gap_open,
        gap_extend,
    )


def local_align(seqa: str, seqb: str, alpha_len: int, score_matrix, gap_open: int, gap_extend: int):
    return _align_wrapper(lib.local_align, seqa, seqb, alpha_len, score_matrix, gap_open, gap_extend)


def global_align(seqa: str, seqb: str, alpha_len: int, score_matrix, gap_open: int, gap_extend: int):
    return _align_wrapper(lib.global_align, seqa, seqb, alpha_len, score_matrix, gap_open, gap_extend)


def glocal_align(seqa: str, seqb: str, alpha_len: int, score_matrix, gap_open: int, gap_extend: int):
    return _align_wrapper(lib.glocal_align, seqa, seqb, alpha_len, score_matrix, gap_open, gap_extend)


def overlap_align(seqa: str, seqb: str, alpha_len: int, score_matrix, gap_open: int, gap_extend: int):
    return _align_wrapper(lib.overlap_align, seqa, seqb, alpha_len, score_matrix, gap_open, gap_extend)
