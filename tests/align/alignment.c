/*
align.c

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

int align_frag_count(AlignFrag *f) {
  int i = 0;
  while (f) {
    i++;
    f = f->next;
  }
  return i;
}

void align_frag_free(AlignFrag *f) {
  while (f) {
    AlignFrag *n = f->next;
    free(f);
    f = n;
  }
}

void alignment_free(Alignment *alignment) {
  align_frag_free(alignment->align_frag);
  free(alignment);
}

Alignment *alignment_new(AlignFrag *align_frag, int score) {
  Alignment *alignment = malloc(sizeof(Alignment));
  if (alignment) {
    alignment->align_frag = align_frag;
    alignment->frag_count = align_frag_count(align_frag);
    alignment->score = score;
  }
  return alignment;
}
