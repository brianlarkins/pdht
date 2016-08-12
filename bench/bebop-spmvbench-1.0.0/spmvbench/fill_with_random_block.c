/**
 * @file fill_with_random_block.c
 * @author mfh
 * @date 2004 Sep 28
 *
 * Renamed from fill_with_random_block.cc and converted back to C++ 
 * (from C) by mfh on 2004 Jan 14
 ********************************************************************/
#include "config.h"
#include "fill_with_random_block.h"
#include <random_number.h>

#ifdef DEBUG
#  include <stdio.h>
#endif /* DEBUG */


void 
fill_with_random_block (double values[], const int row_size, 
                        const int col_size, const int start_index,
                        const int reproducible,
                        const int row_major_block)
{
  int i, j;
  int block_size = row_size * col_size;

  if (row_major_block)
    {
      for (i = start_index; i < block_size + start_index; i++)
        values[i] = smvm_random_double (-1,1);
    }
  else /* Block in column-major order, but iterate through in row-major order. */
    {
      for (i = 0; i < row_size; i++)
        for (j = 0; j < col_size; j++)
          values[row_size*j + i] = smvm_random_double (-1,1);
    }
}


