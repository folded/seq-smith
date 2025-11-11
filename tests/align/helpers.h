/*
helpers.h

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
#ifndef HELPERS_H_INCLUDED
#define HELPERS_H_INCLUDED

#define SWAP(x, y)                                                             \
  do {                                                                         \
    typeof(x) _v = x;                                                          \
    x = y;                                                                     \
    y = _v;                                                                    \
  } while (0)

// ============================================================================
static inline void to_raw(const char *in_seq, unsigned char *out_seq,
                          int seq_len, const unsigned char *map) {
  int i;
  for (i = 0; i < seq_len; i++)
    out_seq[i] = (unsigned char)map[(unsigned)in_seq[i]];
}

#endif
