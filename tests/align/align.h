/*
align.h

Created by Toby Sargeant.
Copyright (c) 2013-2015  Toby Sargeant and The University of Melbourne. All
rights reserved.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

__author__ = "Toby Sargeant"
__copyright__ = "Copyright 2013-2015, Toby Sargeant and The University of
Melbourne"
__credits__ = ["Toby Sargeant","Matthew Wakefield",]
__license__ = "GPLv3"
__version__ = "0.5.1"
__maintainer__ = "Matthew Wakefield"
__email__ = "matthew.wakefield@unimelb.edu.au"
__status__ = "Development"
*/
#ifndef ALIGN_H_INCLUDED
#define ALIGN_H_INCLUDED

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

#include "helpers.h"

typedef enum { A_GAP = 0, B_GAP = 1, MATCH = 2 } FragType;

typedef struct AlignFrag {
  struct AlignFrag *next;
  FragType type;
  int sa_start, sb_start, hsp_len;
} AlignFrag;

typedef struct Alignment {
  AlignFrag *align_frag;
  int frag_count;
  int score;
} Alignment;

int align_frag_count(AlignFrag *f);
void align_frag_free(AlignFrag *f);
void alignment_free(Alignment *alignment);
Alignment *alignment_new(AlignFrag *align_frag, int score);

typedef Alignment *(RawAlignFunc)(const unsigned char *, int,
                                  const unsigned char *, int, int, int *, int,
                                  int);
typedef Alignment *(AlignFunc)(const unsigned char *, int,
                               const unsigned char *, int, int,
                               const unsigned char *, int *_matrix, int, int);

Alignment *align_raw(const unsigned char *sa, int sa_len,
                     const unsigned char *sb, int sb_len, int alpha_len,
                     int *score_matrix, int gap_open, int gap_extend);

Alignment *overlap_align(const unsigned char *seqa, int sa_len,
                         const unsigned char *seqb, int sb_len, int alpha_len,
                         int *score_matrix, int gap_open, int gap_extend);

Alignment *local_align_raw(const unsigned char *sa, int sa_len,
                           const unsigned char *sb, int sb_len, int alpha_len,
                           int *score_matrix, int gap_open, int gap_extend);

Alignment *local_align(const unsigned char *seqa, int sa_len,
                       const unsigned char *seqb, int sb_len, int alpha_len,
                       int *score_matrix, int gap_open, int gap_extend);

Alignment *global_align_raw(const unsigned char *sa, int sa_len,
                            const unsigned char *sb, int sb_len, int alpha_len,
                            int *score_matrix, int gap_open, int gap_extend);

Alignment *global_align(const unsigned char *seqa, int sa_len,
                        const unsigned char *seqb, int sb_len, int alpha_len,
                        int *score_matrix, int gap_open, int gap_extend);

Alignment *glocal_align_raw(const unsigned char *sa, int sa_len,
                            const unsigned char *sb, int sb_len, int alpha_len,
                            int *score_matrix, int gap_open, int gap_extend);

Alignment *glocal_align(const unsigned char *seqa, int sa_len,
                        const unsigned char *seqb, int sb_len, int alpha_len,
                        int *score_matrix, int gap_open, int gap_extend);

#endif
