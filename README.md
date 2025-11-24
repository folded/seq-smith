# seq-smith

<img src="https://raw.githubusercontent.com/folded/seq-smith/main/docs/source/_static/seq-smith.svg" alt="seq-smith" width=200>

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
for frag in alignment.fragments:
    print(frag)
```

## Alignment Types

`seq-smith` supports the following alignment strategies:

- **Global Alignment (`global_align`):** Uses the Needleman-Wunsch algorithm to align every residue in both sequences.
- **Local Alignment (`local_align`):** Uses the Smith-Waterman algorithm to find the best-scoring local region of similarity.
- **Local-Global Alignment (`local_global_align`):** Finds the best local alignment of the first sequence within the second, requiring the second sequence to be aligned globally.
- **Overlap Alignment (`overlap_align`):** Does not penalize gaps at the start or end of either sequence, making it ideal for finding overlaps between sequences (e.g., in sequence assembly).

## Multi-threaded Alignment

`seq-smith` supports multi-threaded alignment for 1-vs-many scenarios. This is useful when you have a query sequence and want to align it against a large database of target sequences.

Available functions:

- `global_align_many`
- `local_align_many`
- `local_global_align_many`
- `overlap_align_many`

Example:

```python
from seq_smith import global_align_many

# ... setup alphabet, score_matrix, gap penalties ...

seqa = encode("ACGT", alphabet)
seqbs = [encode("ACGT", alphabet), encode("AGCT", alphabet), encode("AAAA", alphabet)]

# Align seqa against all sequences in seqbs in parallel
alignments = global_align_many(seqa, seqbs, score_matrix, gap_open, gap_extend, num_threads=4)

for aln in alignments:
    print(aln.score)
```

[full documentation](https://seq-smith.readthedocs.io).
