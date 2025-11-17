# seq-smith

A Rust-based sequence alignment library for Python.

## Installation

You can install `seq-smith` using pip:

```bash
pip install seq-smith
```

## Usage

`seq-smith` provides several alignment functions and helper functions to make sequence alignment easy. Here's a basic example of how to perform a global alignment:

```python
from seq_smith import global_align, make_score_matrix, encode

# Define your alphabet
alphabet = "ACGT"

# Create a scoring matrix
score_matrix = make_score_matrix(alphabet, match_score=1, mismatch_score=-1)

# Encode sequences
seqa = encode("ACGT", alphabet)
seqb = encode("AGCT", alphabet)

# Define gap penalties
gap_open = -2
gap_extend = -1

# Perform the alignment
alignment = global_align(seqa, seqb, score_matrix, gap_open, gap_extend)

# Print the alignment score
print(f"Alignment score: {alignment.score}")

# Print the alignment fragments
for frag in alignment.align_frag:
    print(frag)
```

## Helper Functions

### `make_score_matrix(alphabet, match_score, mismatch_score)`

Creates a scoring matrix for a given alphabet.

**Example:**

```python
from seq_smith import make_score_matrix

alphabet = "ACGT"
score_matrix = make_score_matrix(alphabet, match_score=2, mismatch_score=-1)
print(score_matrix)
```

### `encode(seq, alphabet)`

Encodes a sequence into a byte array using the provided alphabet.

**Example:**

```python
from seq_smith import encode

alphabet = "ACGT"
encoded_seq = encode("AGCT", alphabet)
print(encoded_seq)
# Output: b'\x00\x02\x01\x03'
```

### `decode(encoded_seq, alphabet)`

Decodes a byte-encoded sequence back to a string using the provided alphabet.

**Example:**

```python
from seq_smith import decode

alphabet = "ACGT"
decoded_seq = decode(b'\x00\x02\x01\x03', alphabet)
print(decoded_seq)
# Output: AGCT
```

## Alignment Types

### Global Alignment (`global_align`)

Performs a global alignment between two sequences using the Needleman-Wunsch algorithm. This alignment type attempts to align every residue in both sequences.

**Example:**

```python
from seq_smith import global_align, make_score_matrix, encode

alphabet = "GATTACA"
score_matrix = make_score_matrix(alphabet, 2, -1)
seqa = encode("GATTACA", alphabet)
seqb = encode("GCATGCA", alphabet)

alignment = global_align(seqa, seqb, score_matrix, -2, -1)
# Expected score: 0
```

### Local Alignment (`local_align`)

Performs a local alignment between two sequences using the Smith-Waterman algorithm. This alignment type finds the best-scoring local region of similarity between the two sequences.

**Example:**

```python
from seq_smith import local_align, make_score_matrix, encode

alphabet = "ACGTXYZW"
score_matrix = make_score_matrix(alphabet, 2, -1)
seqa = encode("XXXXXAGCTYYYYY", alphabet)
seqb = encode("ZZZAGCTWWW", alphabet)

alignment = local_align(seqa, seqb, score_matrix, -2, -1)
# Expected score: 8 (for "AGCT")
```

### Local-Global Alignment (`local_global_align`)

This alignment finds the best local alignment of `seqa` within `seqb`, but `seqb` must be aligned globally.

**Example:**

```python
from seq_smith import local_global_align, make_score_matrix, encode

alphabet = "ACGTX"
score_matrix = make_score_matrix(alphabet, 2, -1)
seqa = encode("XACGTX", alphabet)
seqb = encode("ACGT", alphabet)

alignment = local_global_align(seqa, seqb, score_matrix, -2, -1)
# Expected score: 8
```

### Overlap Alignment (`overlap_align`)

Performs an overlap alignment between two sequences. This alignment type does not penalize gaps at the start or end of either sequence, making it suitable for finding overlaps between sequences, such as in sequence assembly.

**Example:**

```python
from seq_smith import overlap_align, make_score_matrix, encode

alphabet = "ACGT"
score_matrix = make_score_matrix(alphabet, 2, -1)
seqa = encode("ACGTACGT", alphabet)
seqb = encode("CGTA", alphabet)

alignment = overlap_align(seqa, seqb, score_matrix, -2, -1)
# Expected score: 8
```
