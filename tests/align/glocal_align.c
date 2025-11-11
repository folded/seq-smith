/*
glocal_align.c

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
#include "align.h"

#define SWAP(x, y)                                                             \
  do {                                                                         \
    typeof(x) _v = x;                                                          \
    x = y;                                                                     \
    y = _v;                                                                    \
  } while (0)

// ============================================================================
static AlignFrag *traceback(const unsigned char *sa, int sa_len,
                            const unsigned char *sb, int sb_len,
                            int *dir_matrix, int s_col, int s_row) {
  AlignFrag *result = NULL;

  while (s_col >= 0 && s_row >= 0) {
    AlignFrag *temp = malloc(sizeof(AlignFrag));
    int d;
    if (!temp) {
      align_frag_free(result);
      return NULL;
    }

    d = dir_matrix[s_col * sb_len + s_row];
    if (d < 0) {
      s_row -= -d;
      temp->type = A_GAP;
      temp->hsp_len = -d;
    } else if (d > 0) {
      s_col -= +d;
      temp->type = B_GAP;
      temp->hsp_len = +d;
    } else {
      int count = 0;
      do {
        --s_col;
        --s_row;
        count++;
      } while (s_col >= 0 && s_row >= 0 &&
               (d = dir_matrix[s_col * sb_len + s_row]) == 0);
      temp->type = MATCH;
      temp->hsp_len = count;
    }
    temp->sa_start = s_col + 1;
    temp->sb_start = s_row + 1;

    temp->next = result;
    result = temp;
  }

  return result;
}

// ============================================================================
Alignment *glocal_align_raw(const unsigned char *sa, int sa_len,
                            const unsigned char *sb, int sb_len, int alpha_len,
                            int *score_matrix, int gap_open, int gap_extend) {
  Alignment *result = NULL;
#define SCORE(a, b) score_matrix[(a) * alpha_len + (b)]
  int score, dir;
  int row, col;
  int vgap_pos;
  int vgap_score;
  int max_score, max_row, max_col;

  int *curr_score = NULL;
  int *prev_score = NULL;
  int *dir_matrix = NULL;
  int *hgap_pos = NULL;
  int *hgap_score = NULL;

  int *dirp = NULL;

  if (sa_len <= 0 || sb_len <= 0 || !sa || !sb || !score_matrix ||
      gap_open > 0 || gap_extend > 0) {
    goto fail;
  }

  curr_score = malloc(sb_len * sizeof(int));
  prev_score = malloc(sb_len * sizeof(int));
  dir_matrix = malloc(sa_len * sb_len * sizeof(int));
  hgap_pos = malloc(sb_len * sizeof(int));
  hgap_score = malloc(sb_len * sizeof(int));

  if (!curr_score || !prev_score || !dir_matrix || !hgap_pos || !hgap_score) {
    goto fail;
  }

  dirp = dir_matrix;

  max_score = score = 0;
  max_row = 0;
  max_col = 0;

  for (row = 0; row < sb_len; row++) {
    prev_score[row] = INT_MIN / 2;
    hgap_pos[row] = -1;
    hgap_score[row] = INT_MIN / 2;
  }
  // fprintf(stderr, "\n");

  for (col = 0; col < sa_len; col++) {
    // first row has to be a match.
    score = SCORE(sa[col], sb[0]);
    dir = 0;

    curr_score[0] = score;
    dirp[0] = dir;
    // fprintf(stderr, "%2d,%2d %2d/%2d:%2d  ", sa[col], sb[0], SCORE(sa[col],
    // sb[0]), curr_score[0], dirp[0]);

    vgap_score = score + gap_open;
    vgap_pos = 0;

    for (row = 1; row < sb_len - 1; row++) {
      score = prev_score[row - 1] + SCORE(sa[col], sb[row]);
      dir = 0;

      if (score < vgap_score) {
        score = vgap_score;
        dir = -(row - vgap_pos);
      }

      if (score < hgap_score[row]) {
        score = hgap_score[row];
        dir = col - hgap_pos[row];
      }

      curr_score[row] = score;
      dirp[row] = dir;
      // fprintf(stderr, "%2d,%2d %2d/%2d:%2d  ", sa[col], sb[row],
      // SCORE(sa[col], sb[row]), curr_score[row], dirp[row]);

      if (score + gap_open >= vgap_score + gap_extend) {
        vgap_score = score + gap_open;
        vgap_pos = row;
      } else {
        vgap_score += gap_extend;
      }

      if (score + gap_open >= hgap_score[row] + gap_extend) {
        hgap_score[row] = score + gap_open;
        hgap_pos[row] = col;
      } else {
        hgap_score[row] += gap_extend;
      }
    }

    // last row has to be a match.
    score = prev_score[row - 1] + SCORE(sa[col], sb[row]);
    dir = 0;

    curr_score[row] = score;
    dirp[row] = dir;
    // fprintf(stderr, "%2d,%2d %2d/%2d:%2d  ", sa[col], sb[row], SCORE(sa[col],
    // sb[row]), curr_score[row], dirp[row]);

    if (score >= max_score) {
      max_row = row;
      max_col = col;
      max_score = score;
    }

    // fprintf(stderr, "\n");

    dirp = dirp + sb_len;
    SWAP(curr_score, prev_score);
  }

  result = alignment_new(
      traceback(sa, sa_len, sb, sb_len, dir_matrix, max_col, max_row),
      max_score);

fail:
  if (curr_score)
    free(curr_score);
  if (prev_score)
    free(prev_score);
  if (dir_matrix)
    free(dir_matrix);
  if (hgap_pos)
    free(hgap_pos);
  if (hgap_score)
    free(hgap_score);

  return result;
#undef SCORE
}

// ============================================================================
Alignment *glocal_align(const unsigned char *seqa, int sa_len,
                        const unsigned char *seqb, int sb_len, int alpha_len,
                        int *score_matrix, int gap_open, int gap_extend) {
  return glocal_align_raw(seqa, sa_len, seqb, sb_len, alpha_len, score_matrix,
                          gap_open, gap_extend);
}
