#ifndef _stencil_h
#define _stencil_h
/**
 * @file stencil.h
 * @author mfh
 * @since 05 Jan 2005 
 * @date 06 Sep 2005
 *
 * Functions for creating (B)CSR format sparse matrices corresponding to various 
 * stencils.
 */

struct bcsr_matrix_t;
struct coo_matrix_t;


/**
 * Creates a COO format sparse matrix which has the nonzero pattern of a 
 * 9-point stencil matrix, and is symmetric positive definite and diagonally 
 * dominant.  Only the rows [start_row, end_row] are included.  If start_row
 * == 0 and end_row == n^2 - 1, the resulting matrix will be n^2 x n^2.
 *
 * @param n [IN]  Dimension of one side of the (square) grid 
 *
 * @return NULL if error in creation, else pointer to the CSR matrix.
 */
struct coo_matrix_t*
stencil_9pt_coo_matrix (const int n, const int start_row, const int end_row);


/**
 * Creates a COO format sparse matrix which has the nonzero pattern of a 
 * 27-point stencil matrix, and is symmetric positive definite and diagonally 
 * dominant.  Only the rows [start_row, end_row] are included.  If start_row
 * == 0 and end_row == n^3 - 1, the resulting matrix will be n^3 x n^3.
 *
 * @param n [IN]  Dimension of one side of the (cubic) grid 
 *
 * @return NULL if error in creation, else pointer to the CSR matrix.
 */
struct coo_matrix_t*
stencil_27pt_coo_matrix (const int n, const int start_row, const int end_row);

/**
 * Creates a COO format sparse matrix which has the nonzero pattern of a 
 * 3-point stencil matrix, and is symmetric positive definite and diagonally 
 * dominant.  Only the rows [start_row, end_row] are included.  If start_row
 * == 0 and end_row == n - 1, the resulting matrix will be n x n.
 *
 * @param n [IN]  Dimension of one side of the (linear) grid 
 *
 * @return NULL if error in creation, else pointer to the CSR matrix.
 */
struct coo_matrix_t*
stencil_3pt_coo_matrix (const int n, const int start_row, const int end_row);


/**
 * Creates a CSR format sparse matrix which has the nonzero structure of the
 * 9-point stencil (2-D grid, points ordered consecutively column-wise) and 
 * nonzero values that are chosen randomly according to a uniform [-1,1] 
 * distribution.  The resulting matrix will be n^2 x n^2.
 *
 * @param n [IN]  Dimension of one side of the (square) grid 
 *
 * @return NULL if error in creation, else pointer to the CSR matrix.
 */
struct bcsr_matrix_t*
stencil_9pt_random_csr_matrix (const int n);

/**
 * Creates a BCSR format sparse matrix with r x c blocks.  The blocks are 
 * arranged in the pattern of a 9-point stencil (2-D grid, points ordered 
 * consecutively column-wise).  The nonzero values are chosen randomly 
 * according to a uniform [-1,1] distribution.  The resulting matrix will
 * be rn^2 x cn^2.
 *
 * @param n [IN]  Dimension of one side of the (square) grid
 * @param r [IN]  Number of rows in each dense register block
 * @param c [IN]  Number of columns in each dense register block
 *
 * @return NULL if error in creation, else pointer to the CSR matrix.
 */
struct bcsr_matrix_t*
stencil_9pt_random_bcsr_matrix (const int n, const int r, const int c);

/**
 * Creates a BCSR format sparse matrix with r x c blocks.  The blocks are 
 * arranged in the pattern of a 3-point stencil (1-D grid, points ordered 
 * consecutively).  The nonzero values are chosen randomly according to a 
 * uniform [-1,1] distribution.  The resulting matrix will be rn x cn.
 *
 * @param n [IN]  Dimension of the grid used to generate the matrix
 * @param r [IN]  Number of rows in each dense register block
 * @param c [IN]  Number of columns in each dense register block
 *
 * @return NULL if error in creation, else pointer to the CSR matrix.
 */
struct bcsr_matrix_t*
stencil_3pt_random_bcsr_matrix (const int n, const int r, const int c);

/**
 * Creates a BCSR format sparse matrix with r x c blocks.  The blocks are 
 * arranged in the pattern of a 27-point stencil (3-D grid, points in
 * the usual consecutive order).  The nonzero values are chosen randomly 
 * according to a uniform [-1,1] distribution.  The resulting matrix will 
 * be rN^3 x cN^3.
 *
 * @param N [IN]  Number of points on one side of the 3-D grid used to generate 
 *                the matrix
 * @param r [IN]  Number of rows in each dense register block
 * @param c [IN]  Number of columns in each dense register block
 *
 * @return NULL if error in creation, else pointer to the CSR matrix.
 */
struct bcsr_matrix_t*
stencil_27pt_random_bcsr_matrix (const int N, const int r, const int c);

/**
 * Creates a BCSR format sparse matrix structure with r x c blocks.  The 
 * blocks are arranged in the pattern of a 3-point stencil (1-D grid, points 
 * in the natural order).  There are no nonzero values -- only the structure
 * is stored.  Only rows [start_block_row, end_block_row] will be generated.
 * The resulting rowptr array will be relative to the specified chunk of the
 * matrix.  This means that, for example, if one matrix is generated with 
 * rows 0 to 10, and another with rows 11 to 20, conjoining the two rowptr 
 * arrays doesn't produce a usable BCSR matrix.  However, this feature makes
 * the output of this function suitable for input into the ParMETIS graph 
 * partitioning algorithms, if the row ranges are distributed across processes.
 *
 * The resulting matrix has c*n columns.
 *
 * @param n [IN]  Number of points in the 1-D grid used to generate the matrix
 * @param r [IN]  Number of rows in each dense register block
 * @param c [IN]  Number of columns in each dense register block
 * @param pnnzb [OUT]  Pointer to nnzb (number of nonzero blocks)
 * @param pptr [OUT]   Pointer to the rowptr array in the BCSR data structure
 * @param pidx [OUT]   Pointer to the colind array in the BCSR data structure
 * @param start_block_row
 * @param end_block_row  
 *
 * @return errcode (nonzero if error, else zero)
 *
 * @warning NOT TESTED!!!
 */
int 
stencil_3pt_dist_bcsr_matrix_structure (const int n, const int r, const int c, 
					int *pnnzb, int **pptr, int **pidx, 
					const int start_block_row, 
					const int end_block_row);

/**
 * Creates a BCSR format sparse matrix structure with r x c blocks.  The 
 * blocks are arranged in the pattern of a 9-point stencil (2-D grid, points 
 * in the natural order).  There are no nonzero values -- only the structure
 * is stored.  Only rows [start_block_row, end_block_row] will be generated.
 * The resulting rowptr array will be relative to the specified chunk of the
 * matrix.  This means that, for example, if one matrix is generated with 
 * rows 0 to 10, and another with rows 11 to 20, conjoining the two rowptr 
 * arrays doesn't produce a usable BCSR matrix.  However, this feature makes
 * the output of this function suitable for input into the ParMETIS graph 
 * partitioning algorithms, if the row ranges are distributed across processes.
 *
 * The resulting matrix has c*n*n columns.
 *
 * @param n [IN]  Number of points in the 2-D grid used to generate the matrix
 * @param r [IN]  Number of rows in each dense register block
 * @param c [IN]  Number of columns in each dense register block
 * @param pnnzb [OUT]  Pointer to nnzb (number of nonzero blocks)
 * @param pptr [OUT]   Pointer to the rowptr array in the BCSR data structure
 * @param pidx [OUT]   Pointer to the colind array in the BCSR data structure
 * @param start_block_row
 * @param end_block_row  
 *
 * @return errcode (nonzero if error, else zero)
 *
 * @warning NOT TESTED!!!
 */
int 
stencil_9pt_dist_bcsr_matrix_structure (const int n, const int r, const int c, 
					int *pnnzb, int **pptr, int **pidx, 
					const int start_block_row, 
					const int end_block_row);

/**
 * Creates a BCSR format sparse matrix structure with r x c blocks.  The 
 * blocks are arranged in the pattern of a 27-point stencil (3-D grid, points 
 * in the natural order).  There are no nonzero values -- only the structure
 * is stored.  Only rows [start_block_row, end_block_row] will be generated.
 * The resulting rowptr array will be relative to the specified chunk of the
 * matrix.  This means that, for example, if one matrix is generated with 
 * rows 0 to 10, and another with rows 11 to 20, conjoining the two rowptr 
 * arrays doesn't produce a usable BCSR matrix.  However, this feature makes
 * the output of this function suitable for input into the ParMETIS graph 
 * partitioning algorithms, if the row ranges are distributed across processes.
 *
 * The resulting matrix has c*n*n*n columns.
 *
 * @param n [IN]  Number of points in the 3-D grid used to generate the matrix
 * @param r [IN]  Number of rows in each dense register block
 * @param c [IN]  Number of columns in each dense register block
 * @param pnnzb [OUT]  Pointer to nnzb (number of nonzero blocks)
 * @param pptr [OUT]   Pointer to the rowptr array in the BCSR data structure
 * @param pidx [OUT]   Pointer to the colind array in the BCSR data structure
 * @param start_block_row
 * @param end_block_row  
 *
 * @return errcode (nonzero if error, else zero)
 *
 * @warning NOT TESTED!!!
 */
int 
stencil_27pt_dist_bcsr_matrix_structure (const int n, const int r, const int c, 
					 int *pnnzb, int **pptr, int **pidx, 
					 const int start_block_row, 
					 const int end_block_row);


#endif /* _stencil_h */
