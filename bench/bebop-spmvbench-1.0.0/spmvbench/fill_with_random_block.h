#ifndef _fill_with_random_block_h
#define _fill_with_random_block_h
/**
 * @file fill_with_random_block.h
 * Maintainer: mfh
 * Last modified: 2003 Dec 27
 * Source language:  C (converted back from C++ by mfh on 2004 Jan 14)
 * For internal use.
 ****************************************************************************/

/**
 * @brief Fills in a block at a given place in values[] with random nonzeros.
 *   
 * @param values         Array used to store nonzeros in block CSR matrix format.
 * @param row_size       Size of row in a block.
 * @param col_size       Size of column in a block.
 * @param start_index    Index of the start of the block in the values array.
 * @param reproducible   true (nonzero) if we want reproducible "random" numbers
 *                       (same random seed with each invocation), false (zero) 
 *                       if we want a different seed each invocation.
 * @param row_major_block   If true (nonzero), blocks will be filled in with 
 *                          unit stride, thus preserving the caller's row-major 
 *                          orientation.  If false (zero), blocks will be filled 
 *                          in with stride r, which takes a column-major oriented 
 *                          block and fills it in a row-oriented way.
 * @return Nothing.
 *
 * Fills in an r x c block at a given position in the values[] array with
 * random nonzeros in [-1,1].  The "row_major_block" option is useful if we
 * want to compare row-major block SMVM with column-major SMVM.  If
 * reproducible is set to true, and row_major_block is set to true for
 * row-major and false for column-major, their test data sets and results of
 * the matrix-vector multiplication should be exactly the same (as interpreted
 * in their respective data structures).  If reproducible is set to false, then
 * there is no need to set row_major_block to false, since in that case the
 * specific values in the block are irrelevant and one can fill it any way one
 * pleases.
 */
void 
fill_with_random_block (double values[], const int row_size, 
                        const int col_size, const int start_index,
                        const int reproducible, const int row_major_block);

#endif /* NOT _fill_with_random_block_h */

