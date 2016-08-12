/**
 * @file stencil.c
 * @author mfh
 * @since 05 Jan 2005 
 * @date 12 Sep 2005
 */
#include <config.h>
#include <stencil.h>
#include <bcsr_matrix.h>
#include <coo_matrix.h>
#include <random_number.h>
#include <smvm_malloc.h>
#include <smvm_util.h>

#include <assert.h>



/**
 * Creates a BCSR format sparse matrix with r x c blocks.  The blocks are 
 * arranged in the pattern of a 27-point stencil (3-D grid, points in
 * the usual consecutive order).  The nonzero values are chosen randomly 
 * according to a uniform [-1,1] distribution.  The resulting matrix will 
 * be rN^3 x cN^3.
 *
 * @param N [IN]  Dimension of one side of the 3-D grid used to generate 
 *                the matrix
 * @param r [IN]  Number of rows in each dense register block
 * @param c [IN]  Number of columns in each dense register block
 * @param pnnzb [OUT]  Pointer to the number of nonzero blocks
 * @param pptr [OUT]  Pointer to the row_start array of the sparse matrix
 * @param pidx [OUT]  Pointer to the col_idx array of the sparse matrix
 * @param pval [OUT]  Pointer to the values array of the sparse matrix
 *
 * @return Nonzero if error in creation, else zero.
 */
static int 
stencil_27pt_regblocked_random (const int N, const int r, const int c, 
				int *pnnzb, int **pptr, int **pidx, 
				double **pval)
{
  int *rowptr = NULL;
  int *colidx = NULL;
  double *values = NULL;
  int i = 0; /* Index of the current row */
  int j = 0; /* Index within each register block */
  int k = 0; /* Number of nonzero blocks added thus far */
  const int block_size = r * c;
  const int N_cubed = N * N * N;
  const int m = r * N_cubed;
  const int n = c * N_cubed;
  const int num_block_cols = n / c;
  int offset_array[27];  /* column index offsets */

  /** 
   * \bug Need to test N_cubed, m and n for overflowing int!!!  
   */

  WITH_DEBUG( fprintf(stderr, "=== stencil_27pt_regblocked_random ===\n") );
  if (N <= 0) return 1;
  if (r <= 0) return 2;
  if (c <= 0) return 2;

  rowptr = smvm_malloc ((N_cubed + 1) * sizeof (int));
  /* We overestimate nnz here; we will later revise the overestimate and scale
   * down the vectors. */
  colidx = smvm_malloc (27 * N_cubed * sizeof (int));
  values = smvm_malloc (block_size * 27 * m * sizeof (double));

  if (N == 1)
    {
      WITH_DEBUG2(fprintf(stderr, "Special case N == 1\n"));

      rowptr[0] = 0;
      for (j = 0; j < block_size; j++)
	values[j] = smvm_random_double(-1,1);
      colidx[0] = 0;
      rowptr[1] = 1;

      WITH_DEBUG(fprintf(stderr, "Reallocating colidx and values arrays..."));
      colidx = smvm_realloc (colidx, 1 * sizeof (int));
      values = smvm_realloc (values, block_size * sizeof (double));
      WITH_DEBUG(fprintf(stderr, "done."));

      *pnnzb = 1;
      *pptr = rowptr;
      *pidx = colidx;
      *pval = values;

      WITH_DEBUG(fprintf(stderr, "=== Done with stencil_27_pt_regblocked_random ===\n"));
      return 0;
    }
  else if (N == 2)
    {
      const int num_offset_array_entries = 15;
      offset_array[0] = -7;
      offset_array[1] = -6;
      offset_array[2] = -5;
      offset_array[3] = -4;
      offset_array[4] = -3;
      offset_array[5] = -2;
      offset_array[6] = -1;
      offset_array[7] = 0;
      offset_array[8] = 1;
      offset_array[9] = 2;
      offset_array[10] = 3;
      offset_array[11] = 4;
      offset_array[12] = 5;
      offset_array[13] = 6;
      offset_array[14] = 7;

      rowptr[0] = 0;
      for (i = 0; i < N_cubed; i++)
	{
	  int ia;
	  for (ia = 0; ia < num_offset_array_entries; ia++)
	    {
	      int t = i + offset_array[ia];

	      if (t >= 0 && t < num_block_cols)
		{
		  for (j = 0; j < block_size; j++) 
		    values[k*block_size+j] = smvm_random_double(-1,1);
		  colidx[k] = c*t;
		  k++;
		}
	    }
	  rowptr[i+1] = k;
	}

      WITH_DEBUG(fprintf(stderr, "Reallocating colidx and values arrays..."));
      colidx = smvm_realloc (colidx, k * sizeof (int));
      values = smvm_realloc (values, k * block_size * sizeof (double));
      WITH_DEBUG(fprintf(stderr, "done."));

      *pnnzb = k;
      *pptr = rowptr;
      *pidx = colidx;
      *pval = values;

      WITH_DEBUG(fprintf(stderr, "=== Done with stencil_27_pt_regblocked_random ===\n"));
      return 0;
    }

  WITH_DEBUG(fprintf(stderr,"nnz upper bound:  %d\n", block_size * 27 * m));

  offset_array[0] = -N*N - N - 1;
  offset_array[1] = -N*N - N;
  offset_array[2] = -N*N - N + 1;
  offset_array[3] = -N*N - 1;
  offset_array[4] = -N*N;
  offset_array[5] = -N*N + 1;
  offset_array[6] = -N*N + N - 1;
  offset_array[7] = -N*N + N;
  offset_array[8] = -N*N + N + 1;
  offset_array[9] = -N - 1;
  offset_array[10] = -N;
  offset_array[11] = -N + 1;
  offset_array[12] = -1;
  offset_array[13] = 0;
  offset_array[14] = 1;
  offset_array[15] = N - 1;
  offset_array[16] = N;
  offset_array[17] = N + 1;
  offset_array[18] = N*N - N - 1;
  offset_array[19] = N*N - N;
  offset_array[20] = N*N - N + 1;
  offset_array[21] = N*N - 1;
  offset_array[22] = N*N;
  offset_array[23] = N*N + 1;
  offset_array[24] = N*N + N - 1;
  offset_array[25] = N*N + N;
  offset_array[26] = N*N + N + 1;

  /* Check offset array for duplicates */
  for (j = 1; j < 27; j++)
    {
      int jj;
      for (jj = 0; jj < j; jj++)
	{
	  if (offset_array[jj] == offset_array[j])
	    fprintf(stderr, "*** stencil_27pt_regblocked_random: WARNING: "
		    "offset_array[%d] == offset_array[%d] == %d ***\n", 
		    jj, j, offset_array[j]);
	}
    }

  rowptr[0] = 0;
  for (i = 0; i < N_cubed; i++)
    {
      int ia;
      for (ia = 0; ia < 27; ia++)
	{
	  int t = i + offset_array[ia];

	  if (t >= 0 && t < num_block_cols)
	    {
	      for (j = 0; j < block_size; j++) 
		values[k*block_size+j] = smvm_random_double(-1,1);
	      colidx[k] = c*t;
	      k++;
	    }
	}
      rowptr[i+1] = k;
    }

  /* When we first allocated colidx and values, we overestimated their size.
   * Now use the actual nnz value to cut back the sizes of the arrays. */
  WITH_DEBUG(fprintf(stderr, "Reallocating colidx and values arrays..."));
  colidx = smvm_realloc (colidx, k * sizeof (int));
  values = smvm_realloc (values, block_size * k * sizeof (double));
  WITH_DEBUG(fprintf(stderr, "done.\n"));

  *pnnzb = k;
  *pptr = rowptr;
  *pidx = colidx;
  *pval = values;
  WITH_DEBUG(fprintf(stderr, "=== Done with stencil_27_pt_regblocked_random ===\n"));
  return 0;
}


/**
 * Creates a BCSR format sparse matrix with r x c blocks.  The blocks are 
 * arranged in the pattern of a 3-point stencil (1-D grid, points ordered 
 * consecutively).  The nonzero values are chosen randomly according to a 
 * uniform [-1,1] distribution.  The resulting matrix will be rn x cn.
 *
 * @param n [IN]  Dimension of the grid used to generate the matrix
 * @param r [IN]  Number of rows in each dense register block
 * @param c [IN]  Number of columns in each dense register block
 * @param pnnzb [OUT]  Pointer to the number of nonzero blocks
 * @param pptr [OUT]  Pointer to the row_start array of the sparse matrix
 * @param pidx [OUT]  Pointer to the col_idx array of the sparse matrix
 * @param pval [OUT]  Pointer to the values array of the sparse matrix
 *
 * @return Nonzero if error in creation, else zero.
 */
static int
stencil_3pt_regblocked_random (const int n, const int r, const int c, 
			       int *pnnzb, int **pptr, int **pidx, 
			       double **pval)
{
  int *rowptr = NULL;
  int *colidx = NULL;
  double *values = NULL;
  int nnz = 0; 
  int nnzb = 0;
  int i = 0;
  int j = 0;
  int k = 0; /* Number of nonzero blocks added thus far */
  const int block_size = r * c;

  WITH_DEBUG( fprintf(stderr, "=== stencil_3pt_regblocked_random ===\n") );
  if (n <= 0) return 1;
  if (r <= 0) return 2;
  if (c <= 0) return 2;

  /*
   * 2 nonzero blocks for block row 0, 3 nonzero blocks for block rows 1:n-2, 
   * 2 nonzero blocks for block row n-1.
   */
  nnzb = 3 * n - 2;
  nnz = block_size * nnzb;
  WITH_DEBUG( fprintf(stderr, "\tn = %d, nnzb = %d, nnz = %d\n", n, nnzb, nnz) );
  if (nnzb <= 0) return 3; 

  rowptr = smvm_malloc ((r*n + 1) * sizeof (int));
  colidx = smvm_malloc (nnzb * sizeof (int));
  values = smvm_malloc (nnz * sizeof (double));

  /* Special cases */
  if (n == 1)
    {
      rowptr[0] = 0;
      WITH_DEBUG(fprintf(stderr, "Special case: n == 1\n"));
      for (j = 0; j < block_size; j++)
	values[j] = 4.0;

      colidx[k] = 0;
      rowptr[1] = 1;

      *pnnzb = 1;
      *pptr  = rowptr; 
      *pidx  = colidx;
      *pval  = values;
      WITH_DEBUG( fprintf(stderr, "=== Done with stencil_3pt_regblocked_random ===\n") );
      return 0;
    }



  /* Block row 0: 2 nonzero blocks */
  WITH_DEBUG(fprintf(stderr, "Block row 0\n"));
  rowptr[0] = 0;
  for (j = 0; j < block_size; j++) values[k*block_size+j] = 4.0;
  colidx[k] = c * i;
  k++;
  for (j = 0; j < block_size; j++) values[k*block_size+j] = -1.0;
  colidx[k] = c * (i + 1);
  k++;
  rowptr[1] = k;
  i++;

  /* Block rows 1 to n-2: 3 nonzero blocks */
  WITH_DEBUG(fprintf(stderr, "Block rows 1 to %d\n", n - 2));
  for (; i < n - 1; i++)
    {
      WITH_DEBUG(fprintf(stderr, "\tblock row %d\n", i));
      for (j = 0; j < block_size; j++) values[k*block_size+j] = -1.0;
      colidx[k] = c * (i - 1);
      k++;
      for (j = 0; j < block_size; j++) values[k*block_size+j] = 4.0;
      colidx[k] = c * i;
      k++;
      for (j = 0; j < block_size; j++) values[k*block_size+j] = -1.0;
      colidx[k] = c * (i + 1);
      k++;
      rowptr[i+1] = k;
    }

  /* Block row n-1: 2 nonzero blocks */
  WITH_DEBUG(fprintf(stderr, "Block row %d\n", n - 1));
  for (j = 0; j < block_size; j++) values[k*block_size+j] = -1.0;
  colidx[k] = c * (i - 1);
  k++;
  for (j = 0; j < block_size; j++) values[k*block_size+j] = 4.0;
  colidx[k] = c * i;
  k++;
  rowptr[i+1] = k;

  *pnnzb = k;
  *pptr = rowptr;
  *pidx = colidx;
  *pval = values;
  WITH_DEBUG( fprintf(stderr, "=== Done with stencil_3pt_regblocked_random ===\n") );
  return 0;
}


static int
__stencil_9pt_coo_matrix (const int n, int *pnnz, int **pI, int **pJ, 
			  double **pval, const int start_row, 
			  const int end_row)
{
  int nnz = 0; 
  int i = 0;
  int k = 0; /* Number of nonzeros added thus far */
  int *II = NULL;
  int *JJ = NULL;
  double *val = NULL;

  WITH_DEBUG2(fprintf(stderr, "=== __stencil_9pt_coo_matrix ===\n"));
  if (n <= 0) return 1;
  if (start_row < 0) return 1;
  if (end_row >= n*n) return 1;
  if (end_row < start_row) return 1;

  /*
   * If we wanted the exact nnz count, it would be the following:
   *
   * 5 nonzeros for row 0, 
   * 6 nonzeros for rows 1 to n-2, 
   * 7 nonzeros for row n-1, 
   * 8 nonzeros for row n, 
   * 9 nonzeros for rows n+1 to n^2 - n - 2, 
   * 8 nonzeros for row n^2 - n - 1, 
   * 7 nonzeros for row n^2 - n, 
   * 6 nonzeros for rows n^2 - n + 1 to n^2 - 2, 
   * 5 nonzeros for row n^2 - 1.
   *
   * Instead, we just overestimate nnz and reallocate the arrays -- a lazy 
   * approach but less error-prone.
   */
  nnz = 9 * (end_row - start_row + 1);
  assert (nnz > 0);

  II = smvm_calloc (nnz, sizeof (int));
  JJ = smvm_calloc (nnz, sizeof (int));
  val = smvm_calloc (nnz, sizeof (double));

#ifdef add_elt
#  undef add_elt
#endif
#define add_elt( value, row, col ) do {  if(i>=start_row && i<=end_row) {val[k] = value; II[k] = row; JJ[k] = col; k++;} } while(0)

  /* Special case: n == 1 */
  if (n == 1)
    {
      assert (start_row == end_row);

      add_elt( 16.0, 0, 0 );
      II = smvm_realloc (II, k * sizeof (int));
      JJ = smvm_realloc (JJ, k * sizeof (int));
      val = smvm_realloc (val, k * sizeof (double));

      *pnnz = k;
      *pI = II;
      *pJ = JJ;
      *pval = val;
      return 0;
    }

  if (start_row >= 0)
    {
      /* Row 0: 5 nonzeros */
      add_elt( 16.0, i, i );
      add_elt( -1.0, i, i+1 ); 
      add_elt( -1.0, i, i+n-1 );
      add_elt( -1.0, i, i+n );
      add_elt( -1.0, i, i+n+1 );
      i++;
    }

  /* Rows 1 to n-2: 6 nonzero blocks */
  for (i = 1; i < MIN(n - 1, end_row+1); i++)
    {
      add_elt( -1.0, i, i-1 );
      add_elt( 16.0, i, i );
      add_elt( -1.0, i, i+1 );
      add_elt( -1.0, i, i+n-1 );
      add_elt( -1.0, i, i+n );
      add_elt( -1.0, i, i+n+1 );
    }

  if (n-1 <= end_row)
    {
      i = n-1; /* just to make sure */
      /* Row n-1: 7 nonzeros */
      add_elt( -1.0, i, i-n+1 );
      add_elt( -1.0, i, i-1 );
      add_elt( 16.0, i, i );
      add_elt( -1.0, i, i+1 );
      add_elt( -1.0, i, i+n-1 );
      add_elt( -1.0, i, i+n );
      add_elt( -1.0, i, i+n+1 );
      i++;
    }

  if (n <= end_row)
    {
      /* Row n: 8 nonzero blocks */
      add_elt( -1.0, i, i-n );
      add_elt( -1.0, i, i-n+1 );
      add_elt( -1.0, i, i-1 );
      add_elt( 16.0, i, i );
      add_elt( -1.0, i, i+1 );
      add_elt( -1.0, i, i+n-1 );
      add_elt( -1.0, i, i+n );
      add_elt( -1.0, i, i+n+1 );
      i++;
    }

   /* Rows n+1 to n^2 - n - 2: 9 nonzeros */
   for (i = n+1; i < MIN(n*n - n - 1, end_row+1); i++)
     {
       add_elt( -1.0, i, i-n-1 );
       add_elt( -1.0, i, i-n );
       add_elt( -1.0, i, i-n+1 );
       add_elt( -1.0, i, i-1 );
       add_elt( 16.0, i, i );
       add_elt( -1.0, i, i+1 );
       add_elt( -1.0, i, i+n-1 );
       add_elt( -1.0, i, i+n );
       add_elt( -1.0, i, i+n+1 );
     }

   if (n*n - n - 1 <= end_row)
     {
       /* Row n^2 - n - 1: 8 nonzeros */
       add_elt( -1.0, i, i-n-1 );
       add_elt( -1.0, i, i-n );
       add_elt( -1.0, i, i-n+1 );
       add_elt( -1.0, i, i-1 );
       add_elt( 16.0, i, i );
       add_elt( -1.0, i, i+1 );
       add_elt( -1.0, i, i+n-1 );
       add_elt( -1.0, i, i+n );
       i++;
     }

   if (n*n - n <= end_row)
     {
       /* Row n^2 - n: 7 nonzeros */
       add_elt( -1.0, i, i-n-1 );
       add_elt( -1.0, i, i-n );
       add_elt( -1.0, i, i-n+1 );
       add_elt( -1.0, i, i-1 );
       add_elt( 16.0, i, i );
       add_elt( -1.0, i, i+1 );
       add_elt( -1.0, i, i+n-1 );
       i++;
     }

  /* Rows n^2 - n + 1 to n^2 - 2: 6 nonzeros */
  for (; i < MIN(n*n - 1, end_row+1); i++)
    {
      add_elt( -1.0, i, i-n-1 );
      add_elt( -1.0, i, i-n );
      add_elt( -1.0, i, i-n+1 );
      add_elt( -1.0, i, i-1 );
      add_elt( 16.0, i, i );
      add_elt( -1.0, i, i+1 );
    }

  if (n*n - 1 <= end_row)
    {
      /* Row n^2 - 1: 5 nonzero blocks */
      add_elt( -1.0, i, i-n-1 );
      add_elt( -1.0, i, i-n );
      add_elt( -1.0, i, i-n+1 );
      add_elt( -1.0, i, i-1 );
      add_elt( 16.0, i, i );
    }

  WITH_DEBUG2(fprintf(stderr, "Reallocating arrays..."));
  II = smvm_realloc (II, k * sizeof (int));
  JJ = smvm_realloc (JJ, k * sizeof (int));
  val = smvm_realloc (val, k * sizeof (double));
  WITH_DEBUG2(fprintf(stderr, "done."));

  *pnnz = k;
  *pI = II;
  *pJ = JJ;
  *pval = val;
  return 0;
}


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
 * @param pnnzb [OUT]  Pointer to the number of nonzero blocks
 * @param pptr [OUT]  Pointer to the row_start array of the sparse matrix
 * @param pidx [OUT]  Pointer to the col_idx array of the sparse matrix
 * @param pval [OUT]  Pointer to the values array of the sparse matrix
 *
 * @return Nonzero if error in creation, else zero.
 */
int
stencil_9pt_regblocked_random (const int n, const int r, const int c, 
			       int *pnnzb, int **pptr, int **pidx, 
			       double **pval)
{
  int *rowptr = NULL;
  int *colidx = NULL;
  double *values = NULL;
  int nnz = 0; 
  int nnzb = 0;
  int i = 0;
  int j = 0;
  int k = 0; /* Number of nonzero blocks added thus far */
  const int block_size = r * c;

  WITH_DEBUG( fprintf(stderr, "=== stencil_9pt_regblocked_random ===\n") );
  if (n <= 0) return 1;
  if (r <= 0) return 2;
  if (c <= 0) return 2;

  /*
   * 5 nonzero blocks for block row 0, 
   * 6 nonzero blocks for block rows 1 to n-2, 
   * 7 nonzero blocks for block row n-1, 
   * 8 nonzero blocks for block row n, 
   * 9 nonzero blocks for block rows n+1 to n^2 - n - 2, 
   * 8 nonzero blocks for block row n^2 - n - 1, 
   * 7 nonzero blocks for block row n^2 - n, 
   * 6 nonzero blocks for block rows n^2 - n + 1 to n^2 - 2, 
   * 5 nonzero blocks for block row n^2 - 1.
   */
  nnzb = 2*(5 + 6*(n-2) + 7 + 8) + 9 * (n*n - 2*n - 2);
  nnz = block_size * nnzb;
  WITH_DEBUG( fprintf(stderr, "\tnnzb = %d, nnz = %d\n", nnzb, nnz) );
  if (nnzb <= 0) return 3; 

  rowptr = smvm_malloc ((n*n+1) * sizeof (int));
  colidx = smvm_malloc (nnzb * sizeof (int));
  values = smvm_malloc (nnz * sizeof (double));

  if (n == 1)
    {
      WITH_DEBUG(fprintf(stderr, "Special case: n == 1\n"));

      rowptr[0] = 0;
      for (j = 0; j < block_size; j++) 
	values[j] = smvm_random_double (-1,1);

      colidx[0] = 0;
      rowptr[1] = 1;

      *pnnzb = 1;
      *pptr = rowptr;
      *pidx = colidx;
      *pval = values;

      return 0;
    }
  else if (n == 2)
    {
      WITH_DEBUG(fprintf (stderr, "Special case: n == 2\n"));

      /* Block row 0 */

      rowptr[0] = 0;
      for (j = 0; j < block_size; j++) 
	values[j] = smvm_random_double (-1,1);
      colidx[0] = 0;

      for (j = 0; j < block_size; j++)
	values[1*block_size + j] = smvm_random_double (-1, 1);
      colidx[1] = c;

      /* Block row 1 */

      rowptr[1] = 2;
      for (j = 0; j < block_size; j++) 
	values[2*block_size + j] = smvm_random_double (-1,1);
      colidx[2] = 0;

      for (j = 0; j < block_size; j++)
	values[3*block_size + j] = smvm_random_double (-1, 1);
      colidx[3] = c;

      for (j = 0; j < block_size; j++)
	values[4*block_size + j] = smvm_random_double (-1, 1);
      colidx[4] = 2 * c;

      /* Block row 2 */

      rowptr[2] = 5;
      for (j = 0; j < block_size; j++) 
	values[5*block_size + j] = smvm_random_double (-1,1);
      colidx[5] = c;

      for (j = 0; j < block_size; j++)
	values[6*block_size + j] = smvm_random_double (-1, 1);
      colidx[6] = 2*c;

      for (j = 0; j < block_size; j++)
	values[7*block_size + j] = smvm_random_double (-1, 1);
      colidx[7] = 3 * c;

      /* Block row 3 */

      rowptr[3] = 8;
      for (j = 0; j < block_size; j++) 
	values[8*block_size + j] = smvm_random_double (-1,1);
      colidx[8] = 2*c;

      for (j = 0; j < block_size; j++)
	values[9*block_size + j] = smvm_random_double (-1, 1);
      colidx[9] = 3*c;

      rowptr[4] = 10;

      *pnnzb = 10;
      *pptr = rowptr;
      *pidx = colidx;
      *pval = values;

      return 0;
    }

  /* Block row 0: 5 nonzero blocks */
  rowptr[0] = 0;
  for (j = 0; j < block_size; j++) 
    values[k*block_size + j] = smvm_random_double (-1,1);
  colidx[k] = c * i;
  k++;

  for (j = 0; j < block_size; j++) 
    values[k*block_size + j] = smvm_random_double (-1,1);
  colidx[k] = c * (i + 1);
  k++;

  for (j = 0; j < block_size; j++) 
    values[k*block_size + j] = smvm_random_double (-1,1);
  colidx[k] = c * (i + n - 1);
  k++;

  for (j = 0; j < block_size; j++) 
    values[k*block_size + j] = smvm_random_double (-1,1);
  colidx[k] = c * (i + n);
  k++;

  for (j = 0; j < block_size; j++) 
    values[k*block_size + j] = smvm_random_double (-1,1);
  colidx[k] = c * (i + n + 1);
  k++;

  rowptr[1] = k;
  i++;

  /* Block rows 1 to n-2: 6 nonzero blocks */
  for (; i < n - 1; i++)
    {
      for (j = 0; j < block_size; j++) 
	values[k*block_size + j] = smvm_random_double (-1,1);
      colidx[k] = c * (i - 1);
      k++;

      for (j = 0; j < block_size; j++) 
	values[k*block_size + j] = smvm_random_double (-1,1);
      colidx[k] = c * i;
      k++;

      for (j = 0; j < block_size; j++) 
	values[k*block_size + j] = smvm_random_double (-1,1);
      colidx[k] = c * (i + 1);
      k++;

      for (j = 0; j < block_size; j++) 
	values[k*block_size + j] = smvm_random_double (-1,1);
      colidx[k] = c * (i + n - 1);
      k++;

      for (j = 0; j < block_size; j++) 
	values[k*block_size + j] = smvm_random_double (-1,1);
      colidx[k] = c * (i + n);
      k++;

      for (j = 0; j < block_size; j++) 
	values[k*block_size + j] = smvm_random_double (-1,1);
      colidx[k] = c * (i + n + 1);
      k++;

      rowptr[i+1] = k;
    }

  /* Block row n-1: 7 nonzero blocks */
  for (j = 0; j < block_size; j++) 
    values[k*block_size + j] = smvm_random_double (-1,1);
  colidx[k] = c * (i - n + 1);
  k++;

  for (j = 0; j < block_size; j++) 
    values[k*block_size + j] = smvm_random_double (-1,1);
  colidx[k] = c * (i - 1);
  k++;

  for (j = 0; j < block_size; j++) 
    values[k*block_size + j] = smvm_random_double (-1,1);
  colidx[k] = c * i;
  k++;

  for (j = 0; j < block_size; j++) 
    values[k*block_size + j] = smvm_random_double (-1,1);
  colidx[k] = c * (i + 1);
  k++;

  for (j = 0; j < block_size; j++) 
    values[k*block_size + j] = smvm_random_double (-1,1);
  colidx[k] = c * (i + n - 1);
  k++;

  for (j = 0; j < block_size; j++) 
    values[k*block_size + j] = smvm_random_double (-1,1);
  colidx[k] = c * (i + n);
  k++;

  for (j = 0; j < block_size; j++) 
    values[k*block_size + j] = smvm_random_double (-1,1);
  colidx[k] = c * (i + n + 1);
  k++;

  rowptr[i+1] = k;
  i++;

  /* Block row n: 8 nonzero blocks */
  for (j = 0; j < block_size; j++) 
    values[k*block_size + j] = smvm_random_double (-1,1);
  colidx[k] = c * (i - n);
  k++;

  for (j = 0; j < block_size; j++) 
    values[k*block_size + j] = smvm_random_double (-1,1);
  colidx[k] = c * (i - n + 1);
  k++;

  for (j = 0; j < block_size; j++) 
    values[k*block_size + j] = smvm_random_double (-1,1);
  colidx[k] = c * (i - 1);
  k++;

  for (j = 0; j < block_size; j++) 
    values[k*block_size + j] = smvm_random_double (-1,1);
  colidx[k] = c * i;
  k++;

  for (j = 0; j < block_size; j++) 
    values[k*block_size + j] = smvm_random_double (-1,1);
  colidx[k] = c * (i + 1);
  k++;

  for (j = 0; j < block_size; j++) 
    values[k*block_size + j] = smvm_random_double (-1,1);
  colidx[k] = c * (i + n - 1);
  k++;

  for (j = 0; j < block_size; j++) 
    values[k*block_size + j] = smvm_random_double (-1,1);
  colidx[k] = c * (i + n);
  k++;

  for (j = 0; j < block_size; j++) 
    values[k*block_size + j] = smvm_random_double (-1,1);
  colidx[k] = c * (i + n + 1);
  k++;

  rowptr[i+1] = k;
  i++;

  /* Block rows n+1 to n^2 - n - 2: 9 nonzero blocks */
  for (; i < n*n - n - 1; i++)
    {
      for (j = 0; j < block_size; j++) 
	values[k*block_size + j] = smvm_random_double (-1,1);
      colidx[k] = c * (i - n - 1);
      k++;

      for (j = 0; j < block_size; j++) 
	values[k*block_size + j] = smvm_random_double (-1,1);
      colidx[k] = c * (i - n);
      k++;

      for (j = 0; j < block_size; j++) 
	values[k*block_size + j] = smvm_random_double (-1,1);
      colidx[k] = c * (i - n + 1);
      k++;

      for (j = 0; j < block_size; j++) 
	values[k*block_size + j] = smvm_random_double (-1,1);
      colidx[k] = c * (i - 1);
      k++;

      for (j = 0; j < block_size; j++) 
	values[k*block_size + j] = smvm_random_double (-1,1);
      colidx[k] = c * i;
      k++;

      for (j = 0; j < block_size; j++) 
	values[k*block_size + j] = smvm_random_double (-1,1);
      colidx[k] = c * (i + 1);
      k++;

      for (j = 0; j < block_size; j++) 
	values[k*block_size + j] = smvm_random_double (-1,1);
      colidx[k] = c * (i + n - 1);
      k++;

      for (j = 0; j < block_size; j++) 
	values[k*block_size + j] = smvm_random_double (-1,1);
      colidx[k] = c * (i + n);
      k++;

      for (j = 0; j < block_size; j++) 
	values[k*block_size + j] = smvm_random_double (-1,1);
      colidx[k] = c * (i + n + 1);
      k++;

      rowptr[i+1] = k;
    }

  /* Block row n^2 - n - 1: 8 nonzero blocks */
  for (j = 0; j < block_size; j++) 
    values[k*block_size + j] = smvm_random_double (-1,1);
  colidx[k] = c * (i - n - 1);
  k++;

  for (j = 0; j < block_size; j++) 
    values[k*block_size + j] = smvm_random_double (-1,1);
  colidx[k] = c * (i - n);
  k++;

  for (j = 0; j < block_size; j++) 
    values[k*block_size + j] = smvm_random_double (-1,1);
  colidx[k] = c * (i - n + 1);
  k++;

  for (j = 0; j < block_size; j++) 
    values[k*block_size + j] = smvm_random_double (-1,1);
  colidx[k] = c * (i - 1);
  k++;

  for (j = 0; j < block_size; j++) 
    values[k*block_size + j] = smvm_random_double (-1,1);
  colidx[k] = c * i;
  k++;

  for (j = 0; j < block_size; j++) 
    values[k*block_size + j] = smvm_random_double (-1,1);
  colidx[k] = c * (i + 1);
  k++;

  for (j = 0; j < block_size; j++) 
    values[k*block_size + j] = smvm_random_double (-1,1);
  colidx[k] = c * (i + n - 1);
  k++;

  for (j = 0; j < block_size; j++) 
    values[k*block_size + j] = smvm_random_double (-1,1);
  colidx[k] = c * (i + n);
  k++;

  rowptr[i+1] = k;
  i++;


  /* Block row n^2 - n: 7 nonzero blocks */
  for (j = 0; j < block_size; j++) 
    values[k*block_size + j] = smvm_random_double (-1,1);
  colidx[k] = c * (i - n - 1);
  k++;

  for (j = 0; j < block_size; j++) 
    values[k*block_size + j] = smvm_random_double (-1,1);
  colidx[k] = c * (i - n);
  k++;

  for (j = 0; j < block_size; j++) 
    values[k*block_size + j] = smvm_random_double (-1,1);
  colidx[k] = c * (i - n + 1);
  k++;

  for (j = 0; j < block_size; j++) 
    values[k*block_size + j] = smvm_random_double (-1,1);
  colidx[k] = c * (i - 1);
  k++;

  for (j = 0; j < block_size; j++) 
    values[k*block_size + j] = smvm_random_double (-1,1);
  colidx[k] = c * i;
  k++;

  for (j = 0; j < block_size; j++) 
    values[k*block_size + j] = smvm_random_double (-1,1);
  colidx[k] = c * (i + 1);
  k++;

  for (j = 0; j < block_size; j++) 
    values[k*block_size + j] = smvm_random_double (-1,1);
  colidx[k] = c * (i + n - 1);
  k++;

  rowptr[i+1] = k;
  i++;


  /* Rows n^2 - n + 1 to n^2 - 2: 6 nonzero blocks */
  for (; i < n*n - 1; i++)
    {
      for (j = 0; j < block_size; j++) 
	values[k*block_size + j] = smvm_random_double (-1,1);
      colidx[k] = c * (i - n - 1);
      k++;
      for (j = 0; j < block_size; j++) 
	values[k*block_size + j] = smvm_random_double (-1,1);
      colidx[k] = c * (i - n);
      k++;
      for (j = 0; j < block_size; j++) 
	values[k*block_size + j] = smvm_random_double (-1,1);
      colidx[k] = c * (i - n + 1);
      k++;
      for (j = 0; j < block_size; j++) 
	values[k*block_size + j] = smvm_random_double (-1,1);
      colidx[k] = c * (i - 1);
      k++;
      for (j = 0; j < block_size; j++) 
	values[k*block_size + j] = smvm_random_double (-1,1);
      colidx[k] = c * i;
      k++;
      for (j = 0; j < block_size; j++) 
	values[k*block_size + j] = smvm_random_double (-1,1);
      colidx[k] = c * (i + 1);
      k++;
      rowptr[i+1] = k;
    }

  /* Row n^2 - 1: 5 nonzero blocks */
  for (j = 0; j < block_size; j++) 
    values[k*block_size + j] = smvm_random_double (-1,1);
  colidx[k] = c * (i - n - 1);
  k++;
  for (j = 0; j < block_size; j++) 
    values[k*block_size + j] = smvm_random_double (-1,1);
  colidx[k] = c * (i - n);
  k++;
  for (j = 0; j < block_size; j++) 
    values[k*block_size + j] = smvm_random_double (-1,1);
  colidx[k] = c * (i - n + 1);
  k++;
  for (j = 0; j < block_size; j++) 
    values[k*block_size + j] = smvm_random_double (-1,1);
  colidx[k] = c * (i - 1);
  k++;
  for (j = 0; j < block_size; j++) 
    values[k*block_size + j] = smvm_random_double (-1,1);
  colidx[k] = c * i;
  k++;
  rowptr[i+1] = k;

  *pnnzb = k;
  *pptr = rowptr;
  *pidx = colidx;
  *pval = values;
  return 0;
}


/**
 * Creates a CSR format sparse matrix which has the nonzero structure of the
 * 9-point stencil (2-D grid, points ordered consecutively column-wise) and 
 * nonzero values that are chosen randomly according to a uniform [-1,1] 
 * distribution.  The resulting matrix will be n^2 x n^2.
 *
 * @param n [IN]  Dimension of one side of the (square) grid 
 * @param pnnz [OUT]  Pointer to the number of nonzeros
 * @param pptr [OUT]  Pointer to the row_start array of the sparse matrix
 * @param pidx [OUT]  Pointer to the col_idx array of the sparse matrix
 * @param pval [OUT]  Pointer to the values array of the sparse matrix
 *
 * @return Nonzero if error in creation, else zero.
 */
static int
stencil_9pt_random (const int n, 
		    int *pnnz, int **pptr, int **pidx, double **pval)
{
  return stencil_9pt_regblocked_random (n, 1, 1, pnnz, pptr, pidx, pval);
}




struct bcsr_matrix_t*
stencil_9pt_random_csr_matrix (const int n)
{
  int nnz = 0;
  int* ptr = NULL;
  int* idx = NULL;
  double *val = NULL;
  int errcode = 0;

  errcode = stencil_9pt_random (n, &nnz, &ptr, &idx, &val);
  if (errcode != 0)
    return NULL;
  return create_bcsr_matrix (n*n, n*n, 1, 1, nnz, val, idx, ptr, UNSYMMETRIC, 0, REAL, 0);
}


struct bcsr_matrix_t*
stencil_27pt_random_bcsr_matrix (const int N, const int r, const int c)
{
  int nnzb = 0;
  int *ptr = NULL;
  int *idx = NULL;
  double *val = NULL;
  int errcode = 0;
  int N_cubed = N*N*N;

  errcode = stencil_27pt_regblocked_random (N, r, c, &nnzb, &ptr, &idx, &val);
  if (errcode != 0)
    return NULL;

  return create_bcsr_matrix (N_cubed, N_cubed, r, c, nnzb, val, idx, ptr, UNSYMMETRIC, 0, REAL, 0);
}


struct bcsr_matrix_t*
stencil_3pt_random_bcsr_matrix (const int n, const int r, const int c)
{
  int nnzb = 0;
  int *ptr = NULL;
  int *idx = NULL;
  double *val = NULL;
  int errcode = 0;

  errcode = stencil_3pt_regblocked_random (n, r, c, &nnzb, &ptr, &idx, &val);
  if (errcode != 0)
    return NULL;

  return create_bcsr_matrix (n, n, r, c, nnzb, val, idx, ptr, UNSYMMETRIC, 0, REAL, 0);
}

struct bcsr_matrix_t*
stencil_9pt_random_bcsr_matrix (const int n, const int r, const int c)
{
  int nnzb = 0;
  int *ptr = NULL;
  int *idx = NULL;
  double *val = NULL;
  int errcode = 0;
  int n_squared = n * n;

  errcode = stencil_9pt_regblocked_random (n, r, c, &nnzb, &ptr, &idx, &val);
  if (errcode != 0)
    return NULL;

  return create_bcsr_matrix (n_squared, n_squared, r, c, nnzb, val, idx, ptr, UNSYMMETRIC, 0, REAL, 0);
}



static int
__stencil_27pt_coo_matrix (const int n, int *pnnz, int **pI, int **pJ, double **pval, const int start_row, const int end_row)
{
  int nnz = 0; 
  int i = 0;
  int k = 0; /* Number of nonzeros added thus far */
  int *II = NULL;
  int *JJ = NULL;
  double *val = NULL;
  const int n_cubed = n*n*n;
  int offset_array[27];  /* column index offsets */

  WITH_DEBUG2(fprintf(stderr, "=== __stencil_27pt_coo_matrix ===\n"));
  if (n <= 0) return 1;
  if (end_row >= n_cubed) return 1;
  if (start_row < 0) return 1;
  if (end_row < start_row) return 1;

  /* At first we overestimate nnz */
  nnz = 27 * (end_row - start_row + 1);
  if (nnz <= 0) return 3;

  II = smvm_calloc (nnz, sizeof (int));
  JJ = smvm_calloc (nnz, sizeof (int));
  val = smvm_calloc (nnz, sizeof (double));

#ifdef add_elt
#  undef add_elt
#endif 
#define add_elt( value, row, col ) do { val[k] = value; II[k] = row; JJ[k] = col; k++; } while(0)

  if (n == 1)
    {
      WITH_DEBUG2(fprintf(stderr, "Special case n == 1\n"));
      add_elt( 32.0, 1, 1 );
    }
  else if (n == 2)
    {
      const int num_offset_array_entries = 15;
      offset_array[0] = -7;
      offset_array[1] = -6;
      offset_array[2] = -5;
      offset_array[3] = -4;
      offset_array[4] = -3;
      offset_array[5] = -2;
      offset_array[6] = -1;
      offset_array[7] = 0;
      offset_array[8] = 1;
      offset_array[9] = 2;
      offset_array[10] = 3;
      offset_array[11] = 4;
      offset_array[12] = 5;
      offset_array[13] = 6;
      offset_array[14] = 7;

      for (i = start_row; i <= end_row; i++)
	{
	  int ia;
	  for (ia = 0; ia < num_offset_array_entries; ia++)
	    {
	      int t = i + offset_array[ia];

	      if (t >= 0 && t < n_cubed)
		{
		  if (t == i)
		    add_elt( 32.0, i, t );
		  else
		    add_elt( -1.0, i, t );
		}
	    }
	}
    }
  else
    {
      const int num_offset_array_entries = 27;
      int j;

      offset_array[0] = -n*n - n - 1;
      offset_array[1] = -n*n - n;
      offset_array[2] = -n*n - n + 1;
      offset_array[3] = -n*n - 1;
      offset_array[4] = -n*n;
      offset_array[5] = -n*n + 1;
      offset_array[6] = -n*n + n - 1;
      offset_array[7] = -n*n + n;
      offset_array[8] = -n*n + n + 1;
      offset_array[9] = -n - 1;
      offset_array[10] = -n;
      offset_array[11] = -n + 1;
      offset_array[12] = -1;
      offset_array[13] = 0;
      offset_array[14] = 1;
      offset_array[15] = n - 1;
      offset_array[16] = n;
      offset_array[17] = n + 1;
      offset_array[18] = n*n - n - 1;
      offset_array[19] = n*n - n;
      offset_array[20] = n*n - n + 1;
      offset_array[21] = n*n - 1;
      offset_array[22] = n*n;
      offset_array[23] = n*n + 1;
      offset_array[24] = n*n + n - 1;
      offset_array[25] = n*n + n;
      offset_array[26] = n*n + n + 1;

      /* Check offset array for duplicates */
      for (j = 1; j < num_offset_array_entries; j++)
	{
	  int jj;
	  for (jj = 0; jj < j; jj++)
	    {
	      if (offset_array[jj] == offset_array[j])
		fprintf(stderr, "*** stencil_27pt_regblocked_random: WARNING: "
			"offset_array[%d] == offset_array[%d] == %d ***\n", 
			jj, j, offset_array[j]);
	    }
	}

      /* Add the nonzeros to the matrix */
      for (i = start_row; i <= end_row; i++)
	{
	  int ia;
	  for (ia = 0; ia < num_offset_array_entries; ia++)
	    {
	      int t = i + offset_array[ia];

	      if (t >= 0 && t < n_cubed)
		{
		  if (t == i)
		    add_elt( 32.0, i, t );
		  else
		    add_elt( -1.0, i, t );
		}
	    }
	}
    }

  WITH_DEBUG(fprintf(stderr, "Reallocating arrays..."));
  II = smvm_realloc (II, k * sizeof (int));
  JJ = smvm_realloc (JJ, k * sizeof (int));
  val = smvm_realloc (val, k * sizeof (double));
  WITH_DEBUG(fprintf(stderr, "done."));

  *pnnz = k;
  *pI = II;
  *pJ = JJ;
  *pval = val;
  return 0;
}


static int
__stencil_3pt_coo_matrix (const int n, int *pnnz, int **pI, int **pJ, double **pval, const int start_row, const int end_row)
{
  int nnz = 0; 
  int i = 0;
  int k = 0; /* Number of nonzeros added thus far */
  int *II = NULL;
  int *JJ = NULL;
  double *val = NULL;
  int offset_array[3];  /* column index offsets */

  WITH_DEBUG2(fprintf(stderr, "=== __stencil_3pt_coo_matrix ===\n"));
  if (n <= 0) return 1;
  if (start_row < 0) return 1;
  if (end_row >= n) return 1;
  if (end_row < start_row) return 1;

  /* At first we overestimate nnz */
  nnz = 3 * (end_row - start_row + 1);
  if (nnz <= 0) return 3;

  II = smvm_calloc (nnz, sizeof (int));
  JJ = smvm_calloc (nnz, sizeof (int));
  val = smvm_calloc (nnz, sizeof (double));

#ifdef add_elt
#  undef add_elt
#endif 
#define add_elt( value, row, col ) do { val[k] = value; II[k] = row; JJ[k] = col; k++; } while(0)

  if (n == 1)
    {
      WITH_DEBUG2(fprintf(stderr, "Special case n == 1\n"));
      add_elt( 4.0, 1, 1 );
    }
  else 
    {
      const int num_offset_array_entries = 3;
      offset_array[0] = -1;
      offset_array[1] = 0;
      offset_array[2] = +1;

      for (i = start_row; i <= end_row; i++)
	{
	  int ia;
	  for (ia = 0; ia < num_offset_array_entries; ia++)
	    {
	      int t = i + offset_array[ia];

	      if (t >= 0 && t < n)
		{
		  if (t == i)
		    add_elt( 4.0, i, t );
		  else
		    add_elt( -1.0, i, t );
		}
	    }
	}
    }

  WITH_DEBUG(fprintf(stderr, "Reallocating arrays..."));
  II = smvm_realloc (II, k * sizeof (int));
  JJ = smvm_realloc (JJ, k * sizeof (int));
  val = smvm_realloc (val, k * sizeof (double));
  WITH_DEBUG(fprintf(stderr, "done."));

  *pnnz = k;
  *pI = II;
  *pJ = JJ;
  *pval = val;
  return 0;
}


struct coo_matrix_t*
stencil_3pt_coo_matrix (const int n, const int start_row, const int end_row)
{
  int errcode = 0;
  int nnz = 0;
  int* II = NULL;
  int* JJ = NULL;
  double *val = NULL;

  errcode = __stencil_3pt_coo_matrix (n, &nnz, &II, &JJ, &val, start_row, end_row);
  if (errcode != 0)
    {
      fprintf (stderr, "*** stencil_3pt_coo_matrix: failed to generate "
	       "3-point stencil matrix! ***\n");
      return NULL;
    }

  return create_coo_matrix (end_row - start_row + 1, n, nnz, II, JJ, val, 
			    ZERO, UNSYMMETRIC, 0, REAL);
}


struct coo_matrix_t*
stencil_9pt_coo_matrix (const int n, const int start_row, const int end_row)
{
  int errcode = 0;
  int nnz = 0;
  int* II = NULL;
  int* JJ = NULL;
  double *val = NULL;

  errcode = __stencil_9pt_coo_matrix (n, &nnz, &II, &JJ, &val, start_row, end_row);
  if (errcode != 0)
    {
      fprintf (stderr, "*** stencil_9pt_coo_matrix: failed to generate "
	       "9-point stencil matrix! ***\n");
      return NULL;
    }

  return create_coo_matrix (end_row - start_row + 1, n*n, nnz, II, JJ, val, 
			    ZERO, UNSYMMETRIC, 0, REAL);
}


struct coo_matrix_t*
stencil_27pt_coo_matrix (const int n, const int start_row, const int end_row)
{
  int errcode = 0;
  int nnz = 0;
  int* II = NULL;
  int* JJ = NULL;
  double *val = NULL;
  const int n_cubed = n * n * n;

  errcode = __stencil_27pt_coo_matrix (n, &nnz, &II, &JJ, &val, start_row, end_row);
  if (errcode != 0)
    {
      fprintf (stderr, "*** stencil_27pt_coo_matrix: failed to generate "
	       "27-point stencil matrix! ***\n");
      return NULL;
    }

  return create_coo_matrix (end_row - start_row + 1, n_cubed, nnz, II, JJ, val, 
			    ZERO, UNSYMMETRIC, 0, REAL);
}




int 
stencil_27pt_dist_bcsr_matrix_structure (const int n, const int r, const int c, 
					 int *pnnzb, int **pptr, int **pidx, 
					 const int start_block_row, 
					 const int end_block_row)
{
  int *rowptr = NULL;
  int *colidx = NULL;
  int i = 0; /* Index of the current row */
  int j = 0; /* Index within each register block */
  int k = 0; /* Number of nonzero blocks added thus far */
  const int block_size = r * c;
  const int n_cubed = n * n * n;
  const int M = r * n_cubed;
  const int N = c * n_cubed;
  const int num_block_cols = N / c;
  int num_offset_array_entries = 27;
  int offset_array[27];  /* column index offsets */

  /** 
   * \bug Need to test n_cubed, M and N for overflowing int!!!  
   */

  WITH_DEBUG( fprintf(stderr, "=== stencil_27pt_dist_bcsr_matrix_structure ===\n") );
  if (n <= 0) return 1;
  if (r <= 0) return 2;
  if (c <= 0) return 2;

  rowptr = smvm_malloc ((end_block_row - start_block_row + 2) * sizeof (int));
  /* We overestimate nnz here; we will later revise the overestimate and scale
   * down the vectors. */
  colidx = smvm_malloc (27 * n_cubed * sizeof (int));

  if (n == 1)
    {
      WITH_DEBUG2(fprintf(stderr, "Special case n == 1\n"));

      rowptr[0] = 0;
      colidx[0] = 0;
      rowptr[1] = 1;

      WITH_DEBUG(fprintf(stderr, "Reallocating colidx array..."));
      colidx = smvm_realloc (colidx, 1 * sizeof (int));
      WITH_DEBUG(fprintf(stderr, "done."));

      *pnnzb = 1;
      *pptr = rowptr;
      *pidx = colidx;

      WITH_DEBUG(fprintf(stderr, "=== Done with stencil_27pt_dist_bcsr_matrix_structure ===\n"));
      return 0;
    }
  else if (n == 2)
    {
      num_offset_array_entries = 15;
      offset_array[0] = -7;
      offset_array[1] = -6;
      offset_array[2] = -5;
      offset_array[3] = -4;
      offset_array[4] = -3;
      offset_array[5] = -2;
      offset_array[6] = -1;
      offset_array[7] = 0;
      offset_array[8] = 1;
      offset_array[9] = 2;
      offset_array[10] = 3;
      offset_array[11] = 4;
      offset_array[12] = 5;
      offset_array[13] = 6;
      offset_array[14] = 7;

      rowptr[0] = 0;
      for (i = start_block_row; i < MIN( end_block_row+1, n_cubed ); i++)
	{
	  int ia;
	  for (ia = 0; ia < num_offset_array_entries; ia++)
	    {
	      int t = i + offset_array[ia];

	      if (t >= 0 && t < num_block_cols)
		{
		  colidx[k] = c*t;
		  k++;
		}
	    }
	  rowptr[i+1 - start_block_row] = k;
	}

      WITH_DEBUG(fprintf(stderr, "Reallocating colidx array..."));
      colidx = smvm_realloc (colidx, k * sizeof (int));
      WITH_DEBUG(fprintf(stderr, "done."));

      *pnnzb = k;
      *pptr = rowptr;
      *pidx = colidx;

      WITH_DEBUG(fprintf(stderr, "=== Done with stencil_27pt_dist_bcsr_matrix_structure ===\n"));
      return 0;
    }

  WITH_DEBUG(fprintf(stderr,"nnz upper bound:  %d\n", block_size * 27 * M));

  offset_array[0] = -n*n - n - 1;
  offset_array[1] = -n*n - n;
  offset_array[2] = -n*n - n + 1;
  offset_array[3] = -n*n - 1;
  offset_array[4] = -n*n;
  offset_array[5] = -n*n + 1;
  offset_array[6] = -n*n + n - 1;
  offset_array[7] = -n*n + n;
  offset_array[8] = -n*n + n + 1;
  offset_array[9] = -n - 1;
  offset_array[10] = -n;
  offset_array[11] = -n + 1;
  offset_array[12] = -1;
  offset_array[13] = 0;
  offset_array[14] = 1;
  offset_array[15] = n - 1;
  offset_array[16] = n;
  offset_array[17] = n + 1;
  offset_array[18] = n*n - n - 1;
  offset_array[19] = n*n - n;
  offset_array[20] = n*n - n + 1;
  offset_array[21] = n*n - 1;
  offset_array[22] = n*n;
  offset_array[23] = n*n + 1;
  offset_array[24] = n*n + n - 1;
  offset_array[25] = n*n + n;
  offset_array[26] = n*n + n + 1;

  /* Check offset array for duplicates */
  for (j = 1; j < 27; j++)
    {
      int jj;
      for (jj = 0; jj < j; jj++)
	{
	  if (offset_array[jj] == offset_array[j])
	    fprintf(stderr, "*** stencil_27pt_dist_bcsr_matrix_structure: WARNING: "
		    "offset_array[%d] == offset_array[%d] == %d ***\n", 
		    jj, j, offset_array[j]);
	}
    }

  rowptr[0] = 0;
  for (i = start_block_row; i < MIN( end_block_row+1, n_cubed ); i++)
    {
      int ia;
      for (ia = 0; ia < num_offset_array_entries; ia++)
	{
	  int t = i + offset_array[ia];

	  if (t >= 0 && t < num_block_cols)
	    {
	      colidx[k] = c*t;
	      k++;
	    }
	}
      rowptr[i+1 - start_block_row] = k;
    }

  /* When we first allocated colidx, we overestimated its size.
   * Now use the actual nnz value to cut back the size of the array. */
  WITH_DEBUG(fprintf(stderr, "Reallocating colidx array..."));
  colidx = smvm_realloc (colidx, k * sizeof (int));
  WITH_DEBUG(fprintf(stderr, "done.\n"));

  *pnnzb = k;
  *pptr = rowptr;
  *pidx = colidx;
  WITH_DEBUG(fprintf(stderr, "=== Done with stencil_27pt_dist_bcsr_matrix_structure ===\n"));
  return 0;
}


int
stencil_9pt_dist_bcsr_matrix_structure (const int n, const int r, const int c, 
					int *pnnzb, int **pptr, int **pidx, 
					const int start_block_row, 
					const int end_block_row)
{
  int *rowptr = NULL;
  int *colidx = NULL;
  int nnzb = 0;
  int i = 0;
  int k = 0; /* Number of nonzero blocks added thus far */
  int num_block_cols = n*n*c;
  int offset_array[9];
  int num_offset_array_entries = 9;

  WITH_DEBUG( fprintf(stderr, "=== stencil_9pt_dist_bcsr_matrix_structure ===\n") );
  if (n <= 0) return 1;
  if (end_block_row < start_block_row) return 1;
  if (r <= 0) return 2;
  if (c <= 0) return 2;

  offset_array[0] = -n - 1;
  offset_array[1] = -n;
  offset_array[2] = -n + 1;
  offset_array[3] = -1;
  offset_array[4] = 0;
  offset_array[5] = +1;
  offset_array[6] = n - 1;
  offset_array[7] = n;
  offset_array[8] = n + 1;

  rowptr = smvm_malloc ((end_block_row - start_block_row + 2) * sizeof (int));
  /* Overestimate nnzb; reallocate it at the end */
  nnzb = 9 * (end_block_row - start_block_row + 1);
  colidx = smvm_malloc (nnzb * sizeof (int));

#if 0
  /*
   * 5 nonzero blocks for block row 0, 
   * 6 nonzero blocks for block rows 1 to n-2, 
   * 7 nonzero blocks for block row n-1, 
   * 8 nonzero blocks for block row n, 
   * 9 nonzero blocks for block rows n+1 to n^2 - n - 2, 
   * 8 nonzero blocks for block row n^2 - n - 1, 
   * 7 nonzero blocks for block row n^2 - n, 
   * 6 nonzero blocks for block rows n^2 - n + 1 to n^2 - 2, 
   * 5 nonzero blocks for block row n^2 - 1.
   */
  nnzb = 2*(5 + 6*(n-2) + 7 + 8) + 9 * (n*n - 2*n - 2);
  nnz = block_size * nnzb;
  WITH_DEBUG( fprintf(stderr, "\tnnzb = %d, nnz = %d\n", nnzb, nnz) );
  if (nnzb <= 0) return 3; 
#endif 

  if (n == 1)
    {
      WITH_DEBUG(fprintf(stderr, "Special case: n == 1\n"));

      rowptr[0] = 0;
      colidx[0] = 0;
      rowptr[1] = 1;

      colidx = smvm_realloc (colidx, 1 * sizeof (int));
      *pnnzb = 1;
      *pptr = rowptr;
      *pidx = colidx;
      return 0;
    }
  else if (n == 2)
    {
      int num_offset_array_entries = 3;
      WITH_DEBUG(fprintf (stderr, "Special case: n == 2\n"));

      offset_array[0] = -1;
      offset_array[1] = 0;
      offset_array[2] = +1;
      num_offset_array_entries = 3;

      /* Block row 0 */
      rowptr[0] = 0;
      for (i = start_block_row; i < MIN( end_block_row+1, n*n ); i++)
	{
	  int ia;
	  for (ia = 0; ia < num_offset_array_entries; ia++)
	    {
	      int t = i + offset_array[ia];

	      if (t >= 0 && t < num_block_cols)
		{
		  colidx[k] = c*t;
		  k++;
		}
	    }
	  rowptr[i+1 - start_block_row] = k;
	}

      *pnnzb = 10;
      *pptr = rowptr;
      *pidx = colidx;
      return 0;
    }


  rowptr[0] = 0;
  for (i = start_block_row; i < MIN( end_block_row+1, n*n ); i++)
    {
      int ia;
      for (ia = 0; ia < num_offset_array_entries; ia++)
	{
	  int t = i + offset_array[ia];

	  if (t >= 0 && t < num_block_cols)
	    {
	      colidx[k] = c*t;
	      k++;
	    }
	}
      rowptr[i+1 - start_block_row] = k;
    }

  /* Reallocate colidx to use only exactly as much space as we need */
  colidx = smvm_realloc (colidx, k * sizeof (int));

  *pnnzb = k;
  *pptr = rowptr;
  *pidx = colidx;
  return 0;
}



int
stencil_3pt_dist_bcsr_matrix_structure (const int n, const int r, const int c, 
					int *pnnzb, int **pptr, int **pidx, 
					const int start_block_row, 
					const int end_block_row)
{
  int *rowptr = NULL;
  int *colidx = NULL;
  int nnzb = 0;
  int i = 0; /* Current row index */
  int k = 0; /* Number of nonzero blocks added thus far */
  int num_block_cols = n * c;
  int offset_array[3];

  WITH_DEBUG( fprintf(stderr, "=== stencil_3pt_dist_bcsr_matrix_structure ===\n") );
  if (n <= 0) return 1;
  if (end_block_row < start_block_row) return 1;
  if (r <= 0) return 2;
  if (c <= 0) return 2;

  offset_array[0] = -1;
  offset_array[1] = 0;
  offset_array[2] = +1;

  /* Overestimate nnzb and readjust later. */
  nnzb = 3 * (end_block_row - start_block_row + 1); 
  if (nnzb <= 0) return 3; 

  /* Special cases */
  if (n == 1)
    {
      rowptr = smvm_malloc (2 * sizeof (int));
      colidx = smvm_malloc (1 * sizeof (int));

      rowptr[0] = 0;
      k = 0;
      colidx[k] = 0;
      rowptr[1] = 1;

      *pnnzb = 1;
      *pptr  = rowptr; 
      *pidx  = colidx;
      WITH_DEBUG2( fprintf(stderr, "=== Done with stencil_3pt_dist_bcsr_matri"
			   "x_structure ===\n") );
      return 0;
    }

  rowptr = smvm_malloc ((end_block_row - start_block_row + 2) * sizeof (int));
  colidx = smvm_malloc (nnzb * sizeof (int));

  rowptr[0] = 0;
  for (i = start_block_row; i < MIN( end_block_row + 1, n ); i++)
    {
      int ia;
      for (ia = 0; ia < 3; ia++)
	{
	  int t = i + offset_array[ia];

	  if (t >= 0 && t < num_block_cols)
	    {
	      colidx[k] = c*t;
	      k++;
	    }
	}
      rowptr[i+1 - start_block_row] = k;
    }

  /* Reallocate colidx to use only exactly as much space as we need */
  colidx = smvm_realloc (colidx, k * sizeof (int));

  *pnnzb = k;
  *pptr = rowptr;
  *pidx = colidx;
  WITH_DEBUG2( fprintf(stderr, "=== Done with stencil_3pt_dist_bcsr_matrix_structure ===\n") );
  return 0;
}


