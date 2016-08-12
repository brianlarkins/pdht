/**
 * @file csr_matrix.c
 * @author Mark Hoemmen
 * @date 02 Mar 2006
 * @since 31 May 2005
 * @version 1.0
 *
 * Copyright (c) 2006, Regents of the University of California 
 * All rights reserved.
 * Redistribution and use in source and binary forms, with or
 * without modification, are permitted provided that the
 * following conditions are met:
 * 
 * * Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 *
 * * Redistributions in binary form must reproduce the above copyright 
 *   notice, this list of conditions and the following disclaimer in 
 *   the documentation and/or other materials provided with the 
 *   distribution.
 *
 * * Neither the name of the University of California, Berkeley, nor
 *   the names of its contributors may be used to endorse or promote
 *   products derived from this software without specific prior
 *   written permission.  
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#include <config.h>
#include <csr_matrix.h>
#include <csr_matmatmult.h>
#include <coo_matrix.h>
#include <csc_matrix.h>
#include <iohb.h>
#include <read_mm.h>
#include <sparse_vector.h>

#include <__complex.h>
#include <smvm_malloc.h>
#include <smvm_util.h>

#include <assert.h>
#include <math.h>   /* fabs */
#include <string.h> /* memcpy */



/**
 * An (index,value) pair.  In our context of sparse matrices in CSR format, 
 * it represents a single row index and its associated nonzero value.
 */
struct index_real_value_pair_t
{
  int index;
  double value;
};

struct index_complex_value_pair_t
{
  int index;
  double_Complex value;
};


/**
 * Compares two (index,value) pairs by their indices.
 */
static int
compare_index_real_value_pairs (const void *pa, const void *pb)
{
  const struct index_real_value_pair_t* a = (const struct index_real_value_pair_t*) pa;
  const struct index_real_value_pair_t* b = (const struct index_real_value_pair_t*) pb;

  if (a->index > b->index)
    return 1;
  else if (a->index < b->index)
    return -1;
  else
    return 0;
}

/**
 * Compares two (index,value) pairs by their indices.
 */
static int
compare_index_complex_value_pairs (const void *pa, const void *pb)
{
  const struct index_complex_value_pair_t* a = (const struct index_complex_value_pair_t*) pa;
  const struct index_complex_value_pair_t* b = (const struct index_complex_value_pair_t*) pb;

  if (a->index > b->index)
    return 1;
  else if (a->index < b->index)
    return -1;
  else
    return 0;
}

static int
compare_ints (const void *pa, const void *pb)
{
  int a = *((int*) pa);
  int b = *((int*) pb);

  if (a > b)
    return 1;
  else if (a < b)
    return -1;
  else 
    return 0;
}




/** 
 * Takes row number i of the given sparse matrix A, and copies all the 
 * indices (and values, if there are any values) of that column into "col".
 *
 * @param A [IN]    Sparse matrix in CSR format
 * @param i [IN]    Index of the row to copy (zero-based indices)
 * @param col [OUT] An array of either index_real_value_pair_t, 
 *                  index_complex_value_pair_t or int, depending on whether 
 *                  the value type of the matrix is REAL, COMPLEX, or PATTERN,
 *                  respectively.
 * @param max_nnz [IN]  The max number of nonzeros in any column of A
 */
static void
copy_row2pairs (const struct csr_matrix_t* A, int i, void* row, int max_nnz)
{
  int a = A->rowptr[i];
  int b = A->rowptr[i+1];
  int nnz = b - a;
  int k;

  assert (nnz <= max_nnz);

  if (A->value_type == REAL)
    {
      struct index_real_value_pair_t* _row = (struct index_real_value_pair_t*) row;
      const double* const values = (const double* const) (A->values);

      for (k = 0; k < nnz; k++)
	{
	  assert ((a+k) < b);
	  _row[k].index = A->colidx[a+k];
	  _row[k].value = values[a+k];
	}
    }
  else if (A->value_type == COMPLEX)
    {
      struct index_complex_value_pair_t* _row = (struct index_complex_value_pair_t*) row;
      const double_Complex* const values = (const double_Complex* const) (A->values);

      for (k = 0; k < nnz; k++)
	{
	  assert ((a+k) < b);
	  _row[k].index = A->colidx[a+k];
	  _row[k].value = values[a+k];
	}
    }
  else if (A->value_type == PATTERN)
    {
      int* _row = (int*) row;

      for (k = 0; k < nnz; k++)
	{
	  assert ((a+k) < b);
	  _row[k] = A->colidx[a+k];
	}
    }
}


/** 
 * Given a sparse matrix A in CSR format, a row number i, and a list of
 * (index,value) pairs, copies the list of (index,value) pairs back into that
 * row of the matrix.
 */
static void
copy_pairs2row (const void* row, int max_nnz, struct csr_matrix_t* A, int i)
{
  int a = A->rowptr[i];
  int b = A->rowptr[i+1];
  int nnz = b - a;

  int k;

  assert (nnz <= max_nnz);

  if (A->value_type == REAL)
    {
      struct index_real_value_pair_t* _row = (struct index_real_value_pair_t*) row;
      double* values = (double*) (A->values);

      for (k = 0; k < nnz; k++)
	{
	  assert ((a+k) < b);
	  A->colidx[a+k] = _row[k].index;
	  values[a+k] = _row[k].value;
	}
    }
  else if (A->value_type == COMPLEX)
    {
      struct index_complex_value_pair_t* _row = (struct index_complex_value_pair_t*) row;
      double_Complex* values = (double_Complex*) (A->values);

      for (k = 0; k < nnz; k++)
	{
	  assert ((a+k) < b);
	  A->colidx[a+k] = _row[k].index;
	  values[a+k] = _row[k].value;
	}
    }
  else if (A->value_type == PATTERN)
    {
      int* _row = (int*) row;

      for (k = 0; k < nnz; k++)
	{
	  assert ((a+k) < b);
	  A->colidx[a+k] = _row[k];
	}
    }
}


/**
 * Sorts the column indices within each column of the sparse matrix A.
 */
static void
csr_matrix_sort_colidx (struct csr_matrix_t* A)
{
  int i;
  int max_nnz;  /* Will be the max # of nonzeros in each column */

  if (A->m <= 0)
    return;

  /* Find the max # of nonzeros in each row.  This lets us allocate 
   * workspace large enough to accomodate any row (otherwise we would
   * have to allocate it separately for _each_ row, which would waste
   * a lot of malloc calls). */
  max_nnz = A->rowptr[1] - A->rowptr[0];
  for (i = 1; i < A->m; i++)
    {
      int nnz = A->rowptr[i+1] - A->rowptr[i];
      max_nnz = MAX(nnz, max_nnz);
    }

  /* Sort the column indices in each column, reordering the corresponding values accordingly. */
  if (A->value_type == REAL)
    {
      /* row: Workspace for sorting (index,value) pairs */
      struct index_real_value_pair_t* row = smvm_malloc (max_nnz * sizeof (struct index_real_value_pair_t));

      for (i = 0; i < A->m; i++)
	{
	  int nnz = A->rowptr[i+1] - A->rowptr[i];
	  /* Copy the (index,value) pairs into the temp location "row". */
	  copy_row2pairs (A, i, row, max_nnz);
	  /* Do the sorting in "row". */
	  qsort (row, nnz, sizeof (struct index_real_value_pair_t), compare_index_real_value_pairs);
	  /* Copy the (index,value) pairs back from "row" into the matrix. */
	  copy_pairs2row (row, max_nnz, A, i);
	}

      /* Free up the temp storage. */
      smvm_free (row);
    }
  else if (A->value_type == COMPLEX)
    {
      struct index_complex_value_pair_t* row = smvm_malloc (max_nnz * sizeof (struct index_complex_value_pair_t));
      for (i = 0; i < A->m; i++)
	{
	  int nnz = A->rowptr[i+1] - A->rowptr[i];
	  copy_row2pairs (A, i, row, max_nnz);
	  qsort (row, nnz, sizeof (struct index_real_value_pair_t), compare_index_complex_value_pairs);
	  copy_pairs2row (row, max_nnz, A, i);
	}

      smvm_free (row);
    }
  else if (A->value_type == PATTERN)
    {
      int* row = smvm_malloc (max_nnz * sizeof (int));

      for (i = 0; i < A->m; i++)
	{
	  int nnz = A->rowptr[i+1] - A->rowptr[i];
	  copy_row2pairs (A, i, row, max_nnz);
	  qsort (row, nnz, sizeof (int), compare_ints);
	  copy_pairs2row (row, max_nnz, A, i);
	}

      smvm_free (row);
    }
}




/*=====================================================================*/
void
pack_csr_matrix (struct csr_matrix_t* A, 
		 const int m, const int n, const int nnz, 
		 void* values, int* colidx, int* rowptr,
		 enum symmetry_type_t symmetry_type,
		 enum symmetric_storage_location_t symmetric_storage_location,
		 enum value_type_t value_type)
{
  WITH_DEBUG2(fprintf(stderr, "=== pack_csr_matrix ===\n")); 

  A->m = m;
  A->n = n;
  A->nnz = nnz;
  A->values = values;
  A->colidx = colidx;
  A->rowptr = rowptr;
  A->symmetry_type = symmetry_type;
  A->symmetric_storage_location = symmetric_storage_location;
  A->value_type = value_type;

  WITH_DEBUG2(fprintf(stderr, "=== Done with pack_csr_matrix ===\n")); 
}


/*=====================================================================*/
void
init_csr_matrix (struct csr_matrix_t* A, 
		 const int m, const int n, const int nnz, 
		 void* values, int* colidx, int* rowptr,
		 enum symmetry_type_t symmetry_type,
		 enum symmetric_storage_location_t symmetric_storage_location,
		 enum value_type_t value_type)
{
  pack_csr_matrix (A, m, n, nnz, values, colidx, rowptr, symmetry_type, 
		   symmetric_storage_location, value_type);
}



/*======================================================================*/
void
dealloc_csr_matrix (struct csr_matrix_t* A)
{
  WITH_DEBUG2(fprintf(stderr, "=== dealloc_csr_matrix ===\n"));

  if (A->values != NULL && A->nnz > 0) 
    {
      WITH_DEBUG2(fprintf(stderr,"\tFreeing values pointer...\n"));
      smvm_free (A->values);
      WITH_DEBUG2(fprintf(stderr,"...done.\n"));
      A->values = NULL;
    }

  if (A->colidx != NULL && A->nnz > 0)
    {
      WITH_DEBUG2(fprintf(stderr,"\tFreeing colidx pointer...\n"));
      smvm_free (A->colidx);
      WITH_DEBUG2(fprintf(stderr,"...done.\n"));
      A->colidx = NULL;
    }

  if (A->rowptr != NULL && A->n > 0)
    {
      WITH_DEBUG2(fprintf(stderr,"\tFreeing rowptr pointer...\n"));
      smvm_free (A->rowptr);
      WITH_DEBUG2(fprintf(stderr,"...done.\n"));
      A->rowptr = NULL;
    }

  A->m = 0;
  A->n = 0;
  A->nnz = 0;

  WITH_DEBUG2(fprintf(stderr, "=== Done with dealloc_csr_matrix ===\n"));
}


/*======================================================================*/
void
copy_csr_matrix (struct csr_matrix_t* dest, const struct csr_matrix_t* src)
{
  const int m = src->m;
  const int n = src->n;
  const int nnz = src->nnz;

  WITH_DEBUG2(fprintf(stderr, "=== copy_csr_matrix ===\n"));
  dest->m   = m;
  dest->n   = n;
  dest->nnz = nnz;

  if (src->value_type == REAL)
    {
      dest->values = smvm_malloc (nnz * sizeof (double));
      memcpy ((double*) (dest->values), (double*) (src->values), 
	      nnz * sizeof (double));
    }
  else if (src->value_type == COMPLEX)
    {
      dest->values = smvm_malloc (nnz * sizeof (double_Complex));
      memcpy ((double_Complex*) (dest->values), 
	      (double_Complex*) (src->values), 
	      nnz * sizeof (double_Complex));
    }
  else if (src->value_type == PATTERN)
    dest->values = NULL;

  dest->colidx = smvm_malloc (nnz * sizeof (int));
  dest->rowptr = smvm_malloc ((m+1) * sizeof (int));

  memcpy (dest->colidx, src->colidx, nnz * sizeof (int));
  memcpy (dest->rowptr, src->rowptr, (m+1) * sizeof (int));

  dest->symmetry_type = src->symmetry_type;
  dest->symmetric_storage_location = src->symmetric_storage_location;
  dest->value_type = src->value_type;

  WITH_DEBUG2(fprintf(stderr, "=== Done with copy_csr_matrix ===\n"));
}


/*======================================================================*/
void
unpack_csr_matrix (const struct csr_matrix_t* A,
		   int* m, int* n, int* nnz,
		   void** values, int** colidx, int** rowptr,
		   enum symmetry_type_t* symmetry_type,
		   enum symmetric_storage_location_t* symmetric_storage_location,
		   enum value_type_t* value_type)
{
  WITH_DEBUG2(fprintf(stderr, "=== unpack_csr_matrix ===\n"));

  *m = A->m;
  *n = A->n;
  *nnz = A->nnz;

  *values = A->values;
  *colidx = A->colidx;
  *rowptr = A->rowptr;

  *symmetry_type = A->symmetry_type;
  *symmetric_storage_location = A->symmetric_storage_location;
  *value_type = A->value_type;

  WITH_DEBUG2(fprintf(stderr, "=== Done with unpack_csr_matrix ===\n"));
}



/*======================================================================*/
struct csr_matrix_t*
create_csr_matrix (const int m, const int n, const int nnz, 
		   void* values, int* colidx, int* rowptr,
		   enum symmetry_type_t symmetry_type,
		   enum symmetric_storage_location_t symmetric_storage_location,
		   enum value_type_t value_type)
{
  struct csr_matrix_t *A = smvm_malloc (sizeof (struct csr_matrix_t));

  WITH_DEBUG2(fprintf(stderr, "=== create_csr_matrix ===\n")); 
  init_csr_matrix (A, m, n, nnz, values, colidx, rowptr, symmetry_type, 
		   symmetric_storage_location, value_type);
  WITH_DEBUG2(fprintf(stderr, "=== Done with create_csr_matrix ===\n")); 
  return A;
}


/*======================================================================*/
void
destroy_csr_matrix (struct csr_matrix_t* A)
{
  WITH_DEBUG2(fprintf(stderr, "=== destroy_csr_matrix ===\n")); 
  if (A != NULL)
    {
      dealloc_csr_matrix (A);
      smvm_free (A);
    }
  WITH_DEBUG2(fprintf(stderr, "=== Done with destroy_csr_matrix ===\n")); 
}

struct csc_matrix_t*
csr_to_csc (const struct csr_matrix_t* A)
{
  int *col_nnz = NULL; /* # of nonzeros in each column */
  int i;
  struct csc_matrix_t* B = smvm_calloc (1, sizeof (struct csc_matrix_t));

  B->m = A->m;
  B->n = A->n;
  B->nnz = A->nnz;
  B->symmetry_type = A->symmetry_type;
  B->symmetric_storage_location = A->symmetric_storage_location;
  B->value_type = A->value_type;

  B->colptr = smvm_malloc ((A->n + 1) * sizeof (int));
  B->rowidx = smvm_malloc (A->nnz * sizeof (int));
  if (A->value_type == REAL)
    B->values = smvm_malloc (A->nnz * sizeof (double));
  else if (A->value_type == COMPLEX)
    B->values = smvm_malloc (A->nnz * sizeof (double_Complex));
  else if (A->value_type == PATTERN)
    B->values = NULL;

  if (A->symmetry_type == SYMMETRIC)
    {
      /* Just give it the transpose, since $A^T = A$ for a symmetric matrix. 
       * TODO:  Think about a clever way to do this for a skew-symmetric or
       * a Hermitian matrix!  */

      assert (A->m == A->n);

      memcpy (B->colptr, A->rowptr, (A->n + 1) * sizeof (int));
      memcpy (B->rowidx, A->colidx, A->nnz * sizeof (int));

      if (A->value_type == REAL)
	memcpy ((double*) (B->values), (double*) (A->values), 
		A->nnz * sizeof (double));
      else if (A->value_type == COMPLEX)
	memcpy ((double_Complex*) (B->values), 
		(double_Complex*) (A->values), 
		A->nnz * sizeof (double));

      return B;
    }

  /* count # of non-zeros per column */

  col_nnz = smvm_calloc (A->n, sizeof(int));

  for (i = 0; i < A->nnz; i++)
    {
      int k = A->colidx[i];
      col_nnz[k]++;
    }

  /*
   *  initialize CSC's column pointers.
   *  reset col_nnz to zero.
   */
  B->colptr[0] = 0;
  for (i = 1; i <= A->n; i++)
    {
      B->colptr[i] = B->colptr[i-1] + col_nnz[i-1];
      col_nnz[i-1] = 0;
    }

  /*
   *  convert from CSR to CSC.
   *  use col_nnz to keep track of the number of non-zeros
   *  added to each column.
   */
  if (A->value_type == REAL)
    {
      double* in_values  = (double*) (A->values);
      double* out_values = (double*) (B->values);

      for (i = 0; i < A->m; i++)
	{
	  int k;
	  int nnz_row;    /* # of non-zeros in row i of A */

	  nnz_row = A->rowptr[i+1] - A->rowptr[i];

	  for (k = A->rowptr[i]; k < (A->rowptr[i]+nnz_row); k++)
	    {
	      int j = A->colidx[ k ];             /* col index */
	      double a = in_values[ k ];          /* non-zero value */
	      int h = B->colptr[j] + col_nnz[j];  /* non-zero position */

	      /* add the non-zero A(i,j) to B */

	      B->rowidx[ h ] = i;
	      out_values[ h ] = a;

	      col_nnz[j]++;
	    }
	}
    }
  else if (A->value_type == COMPLEX)
    {
      double_Complex* in_values  = (double_Complex*) (A->values);
      double_Complex* out_values = (double_Complex*) (B->values);

      for (i = 0; i < A->m; i++)
	{
	  int k;
	  int nnz_row;    /* # of non-zeros in row i of A */

	  nnz_row = A->rowptr[i+1] - A->rowptr[i];

	  for (k = A->rowptr[i]; k < (A->rowptr[i]+nnz_row); k++)
	    {
	      int j = A->colidx[ k ];               /* col index */
	      double_Complex a = in_values[ k ];   /* non-zero value */
	      int h = B->colptr[j] + col_nnz[j];    /* non-zero position */

	      /* add the non-zero A(i,j) to B */

	      B->rowidx[ h ] = i;
	      out_values[ h ] = a;

	      col_nnz[j]++;
	    }
	}
    }
  else if (A->value_type == PATTERN)
    {
      for (i = 0; i < A->m; i++)
	{
	  int k;
	  int nnz_row;    /* # of non-zeros in row i of A */

	  nnz_row = A->rowptr[i+1] - A->rowptr[i];

	  for (k = A->rowptr[i]; k < (A->rowptr[i]+nnz_row); k++)
	    {
	      int j = A->colidx[ k ];             /* col index */
	      int h = B->colptr[j] + col_nnz[j];  /* non-zero position */

	      /* add the non-zero A(i,j) to B */
	      B->rowidx[ h ] = i;
	      col_nnz[j]++;
	    }
	}
    }

  /* clean-up */
  smvm_free (col_nnz);

  return B;
}


struct csr_matrix_t*
csc_to_csr (struct csc_matrix_t* A)
{
  int *row_nnz = NULL; /* # of nonzeros in each row */
  int i, j;
  struct csr_matrix_t* B = smvm_calloc (1, sizeof (struct csr_matrix_t));

  B->m = A->m;
  B->n = A->n;
  B->nnz = A->nnz;
  B->symmetry_type = A->symmetry_type;
  B->symmetric_storage_location = A->symmetric_storage_location;
  B->value_type = A->value_type;

  B->rowptr = smvm_malloc ((A->m + 1) * sizeof (int));
  B->colidx = smvm_malloc (A->nnz * sizeof (int));
  if (A->value_type == REAL)
    B->values = smvm_malloc (A->nnz * sizeof (double));
  else if (A->value_type == COMPLEX)
    B->values = smvm_malloc (A->nnz * sizeof (double_Complex));
  else if (A->value_type == PATTERN)
    B->values = NULL;

  if (A->symmetry_type == SYMMETRIC)
    {
      /* Just give it the transpose, since $A^T = A$ for a symmetric matrix. 
       * TODO:  Think about a clever way to do this for a skew-symmetric or
       * a Hermitian matrix!  */

      assert (A->m == A->n);

      memcpy (B->rowptr, A->colptr, (A->m + 1) * sizeof (int));
      memcpy (B->colidx, A->rowidx, A->nnz * sizeof (int));

      if (A->value_type == REAL)
	memcpy ((double*) (B->values), (double*) (A->values), 
		A->nnz * sizeof (double));
      else if (A->value_type == COMPLEX)
	memcpy ((double_Complex*) (B->values), 
		(double_Complex*) (A->values), 
		A->nnz * sizeof (double));

      return B;
    }

  /* count # of non-zeros per row */

  row_nnz = smvm_calloc (A->m, sizeof(int));

  for (i = 0; i < A->nnz; i++)
    {
      int k = A->rowidx[i];
      row_nnz[k]++;
    }

  /*
   *  initialize CSR's row pointers.
   *  reset row_nnz to zero.
   */
  B->rowptr[0] = 0;
  for (i = 1; i <= A->m; i++)
    {
      B->rowptr[i] = B->rowptr[i-1] + row_nnz[i-1];
      row_nnz[i-1] = 0;
    }

  /*
   *  convert from CSC to CSR.
   *  use row_nnz to keep track of the number of non-zeros
   *  added to each row.
   */
  if (A->value_type == REAL)
    {
      double* in_values  = (double*) (A->values);
      double* out_values = (double*) (B->values);

      for (j = 0; j < A->n; j++)
	{
	  int k;
	  int nnz_col;    /* # of non-zeros in column j of A */

	  nnz_col = A->colptr[j+1] - A->colptr[j];

	  for (k = A->colptr[j]; k < (A->colptr[j]+nnz_col); k++)
	    {
	      int i = A->rowidx[ k ];             /* row index */
	      double a = in_values[ k ];          /* non-zero value */
	      int h = B->rowptr[i] + row_nnz[i];  /* non-zero position */

	      /* add the non-zero A(i,j) to B */

	      B->colidx[ h ] = j;
	      out_values[ h ] = a;

	      row_nnz[i]++;
	    }
	}
    }
  else if (A->value_type == COMPLEX)
    {
      double_Complex* in_values  = (double_Complex*) (A->values);
      double_Complex* out_values = (double_Complex*) (B->values);

      for (j = 0; j < A->n; j++)
	{
	  int k;
	  int nnz_col;    /* # of non-zeros in column j of A */

	  nnz_col = A->colptr[j+1] - A->colptr[j];

	  for (k = A->colptr[j]; k < (A->colptr[j]+nnz_col); k++)
	    {
	      int i = A->rowidx[ k ];               /* row index */
	      double_Complex a = in_values[ k ];   /* non-zero value */
	      int h = B->rowptr[i] + row_nnz[i];    /* non-zero position */

	      /* add the non-zero A(i,j) to B */

	      B->colidx[ h ] = j;
	      out_values[ h ] = a;

	      row_nnz[i]++;
	    }
	}
    }
  else if (A->value_type == PATTERN)
    {
      for (j = 0; j < A->n; j++)
	{
	  int k;
	  int nnz_col;    /* # of non-zeros in column j of A */

	  nnz_col = A->colptr[j+1] - A->colptr[j];

	  for (k = A->colptr[j]; k < (A->colptr[j]+nnz_col); k++)
	    {
	      int i = A->rowidx[ k ];  /* row index */
	      int h = B->rowptr[i] + row_nnz[i];  /* non-zero position */

	      /* add the non-zero A(i,j) to B */
	      B->colidx[ h ] = j;
	      row_nnz[i]++;
	    }
	}
    }

  /* clean-up */
  smvm_free (row_nnz);

  return B;
}


int
save_csr_matrix_in_harwell_boeing_format (const char* const filename, 
					  const struct csr_matrix_t* A)
{
  struct csc_matrix_t *B = csr_to_csc (A);

  if (B == NULL)
    return -1;
  else
    {
      const int errcode = save_csc_matrix_in_harwell_boeing_format (filename, B);
      destroy_csc_matrix (B);
      return errcode;
    }
}


int
print_csr_matrix_in_matrix_market_format (FILE* out, const struct csr_matrix_t* A)
{
  int start, end;
  int i, k;
  char value_type_label[20];
  char symmetry_type_label[20];

  WITH_DEBUG2(fprintf(stderr, "=== print_csr_matrix_in_matrix_market_format ===\n")); 

  if (A->value_type == REAL)
    strncpy (value_type_label, "real", 19);
  else if (A->value_type == COMPLEX)
    strncpy (value_type_label, "complex", 19);
  else if (A->value_type == PATTERN)
    strncpy (value_type_label, "pattern", 19);
  else 
    {
      fprintf (stderr, "*** print_csr_matrix_in_matrix_market_format: Unsupported value type! ***\n");
      WITH_DEBUG2(fprintf(stderr, "=== Done with print_csr_matrix_in_matrix_market_format ===\n")); 
      return -1;
    }

  if (A->symmetry_type == UNSYMMETRIC)
    strncpy (symmetry_type_label, "general", 19);
  else if (A->symmetry_type == SYMMETRIC)
    strncpy (symmetry_type_label, "symmetric", 19);
  else if (A->symmetry_type == SKEW_SYMMETRIC)
    strncpy (symmetry_type_label, "skew-symmetric", 19);
  else if (A->symmetry_type == HERMITIAN)
    strncpy (symmetry_type_label, "hermitian", 19);
  else
    {
      fprintf (stderr, "*** print_csr_matrix_in_matrix_market_format: Unsupported symmetry type! ***\n");
      WITH_DEBUG2(fprintf(stderr, "=== Done with print_csr_matrix_in_matrix_market_format ===\n")); 
      return -1;
    }

  fprintf (out, "%%%%MatrixMarket matrix coordinate %s %s\n", 
	   value_type_label, symmetry_type_label);

  if (smvm_debug_level() > 1)
    {
      fprintf (out, "%% rowptr[%d]: ", A->m + 1);
      for (i = 0; i <= A->m; i++)
	fprintf (out, " %d", A->rowptr[i]);

      fprintf (out, "\n%% colidx[%d]: ", A->nnz);
      for (i = 0; i < A->nnz; i++)
	fprintf (out, " %d", A->colidx[i]);

      if (A->value_type != PATTERN)
	{
	  fprintf (out, "\n%% values[%d]: ", A->nnz);

	  if (A->value_type == REAL)
	    {
	      const double* const values = (const double* const) (A->values);

	      for (i = 0; i < A->nnz; i++)
		fprintf (out, " %g", values[i]);
	    }
	  else if (A->value_type == COMPLEX)
	    {
	      const double_Complex* const values = (const double_Complex* const) (A->values);

	      for (i = 0; i < A->nnz; i++)
		fprintf (out, " %g+I*%g", double_Complex_real_part(values[i]), double_Complex_imag_part(values[i]));
	    }
	}

      fprintf (out, "\n");
    }

  fprintf (out, "%d %d %d\n", A->m, A->n, A->nnz);

  if (A->value_type == REAL)
    {
      const double* const values = (const double* const) (A->values);

      for (i = 0; i < A->m; i++)
	{
	  start = A->rowptr[i];
	  end   = A->rowptr[i+1];

	  /* MatrixMarket files use 1-based indices. */
	  for (k = start; k < end; k++)
	    fprintf (out, "%d %d %.13e\n", i+1, A->colidx[k] + 1, values[k]);
	}
    }
  else if (A->value_type == COMPLEX)
    {
      const double_Complex* const values = (const double_Complex* const) (A->values);

      for (i = 0; i < A->m; i++)
	{
	  start = A->rowptr[i];
	  end   = A->rowptr[i+1];

	  /* MatrixMarket files use 1-based indices. */
	  for (k = start; k < end; k++)
	    fprintf (out, "%d %d %.13e %.13e\n", i+1, A->colidx[k] + 1, 
		     double_Complex_real_part(values[k]), 
		     double_Complex_imag_part(values[k]));
	}
    }
  else if (A->value_type == PATTERN)
    {
      for (i = 0; i < A->m; i++)
	{
	  start = A->rowptr[i];
	  end   = A->rowptr[i+1];

	  /* MatrixMarket files use 1-based indices. */
	  for (k = start; k < end; k++)
	    fprintf (out, "%d %d\n", i+1, A->colidx[k] + 1);
	}
    }

  WITH_DEBUG2(fprintf(stderr, "=== Done with print_csc_matrix_in_matrix_mar"
		      "ket_format ===\n")); 
  return 0;
}



int 
save_csr_matrix_in_matrix_market_format (const char* const filename, 
					 const struct csr_matrix_t* A)
{
  FILE* out = NULL;
  int errcode = 0;

  WITH_DEBUG2(fprintf(stderr, "=== save_csr_matrix_in_matrix_market_format"
		      " ===\n")); 
  out = fopen (filename, "w");
  if (out == NULL)
    {
      fprintf (stderr, "*** save_csr_matrix_in_matrix_market_format: failed "
	       "to open output file \"%s\" ***\n", filename);
      return -1;
    }

  errcode = print_csr_matrix_in_matrix_market_format (out, A);
  if (0 != fclose (out))
    {
      fprintf (stderr, "*** save_csr_matrix_in_matrix_market_format: failed "
	       "to close output file \"%s\" ***\n", filename);
      return -1;
    }
  WITH_DEBUG2(fprintf(stderr, "=== Done with save_csr_matrix_in_matrix_marke"
		      "t_format ===\n")); 
  return errcode;
}


int
csr_matrix_expand_symmetric_storage (struct csr_matrix_t* A)
{
  /* # of non-zeros in each row of original symmetric matrix */
  int* cur_row_nnz;   
  /* # of non-zeros in each row of final expanded matrix */
  int* new_row_nnz;   
  /* total # of non-zeros in final expanded matrix */
  int new_nnz;        
  int i;
  int* new_rowptr = NULL;
  int* new_colidx = NULL;
  void* new_values = NULL;

  WITH_DEBUG2(fprintf(stderr, "=== csr_matrix_expand_symmetric_storage ===\n"));

  if (A->m != A->n)
    {
      fprintf (stderr, "*** csr_matrix_expand_symmetric_storage: A is not "
	       "square! ***\n");
      return -1;
    }
  else if (A->symmetry_type == UNSYMMETRIC)
    {
      WITH_DEBUG(fprintf (stderr, "*** csr_matrix_expand_symmetric_storage: A"
			  " is already stored in an unsymmetric format! ***\n"));
      return 0; /* Not an error -- just means that we don't have any work to do */
    }

  /* set-up */
  cur_row_nnz = smvm_calloc (A->m, sizeof (int));
  new_row_nnz = smvm_calloc (A->m, sizeof(int));

  /*
   * Scan A and count how many new non-zeros we'll need to create.
   *
   * Post:
   *   cur_row_nnz[i] == # of non-zeros in row i of the original symmetric 
   *                     matrix.
   *   new_row_nnz[i] == # of non-zeros to be stored in row i of the final 
   *                     expanded matrix.
   *   new_nnz == total # of non-zeros to be stored in the final expanded 
   *              matrix.
   */
  new_nnz = 0;
  for (i = 0; i < A->m; i++)
    {
      cur_row_nnz[i] = A->rowptr[i+1] - A->rowptr[i];
      new_row_nnz[i] = cur_row_nnz[i];
      new_nnz += new_row_nnz[i];
    }
  for (i = 0; i < A->m; i++)
    {
      int k;
      for (k = A->rowptr[i]; k < A->rowptr[i+1]; k++)
	{
	  int j = A->colidx[k];
	  if (j != i)
	    {
	      new_row_nnz[j]++;
	      new_nnz++;
	    }
	}
    }

  /*
   *  Initialize row pointers in expanded matrix.
   *
   *  Post:
   *    new_rowptr initialized to the correct, final values.
   *    new_row_nnz[i] reset to be equal to cur_row_nnz[i].
   */
  new_rowptr = smvm_calloc (A->m + 1, sizeof(int));
  for (i = 1; i <= A->m; i++)
    {
      new_rowptr[i] = new_rowptr[i-1] + new_row_nnz[i-1];
      new_row_nnz[i-1] = cur_row_nnz[i-1];
    }
  new_rowptr[A->m] = new_nnz;

  new_colidx = smvm_calloc (new_nnz, sizeof (int));
  if (A->value_type == REAL)
    new_values = smvm_malloc (new_nnz * sizeof (double));
  else if (A->value_type == COMPLEX)
    new_values = smvm_malloc (new_nnz * sizeof (double_Complex));
  else if (A->value_type == PATTERN)
    new_values = NULL;
 
  /*
   *  Complete expansion of A to full storage.
   *
   *  Post:
   *    (new_rowptr, new_colidx, new_values) is the full-storage equivalent of A.
   *    new_row_nnz[i] == # of non-zeros in row i of A.
   */
  for (i = 0; i < A->m; i++)
    {
      int cur_nnz = cur_row_nnz[i];
      int k_cur   = A->rowptr[i];
      int k_new   = new_rowptr[i];

      /* copy current non-zeros from old matrix to new matrix */
      memcpy (new_colidx + k_new, A->colidx + k_cur, cur_nnz * sizeof(int));
      /* mfh 24 Feb 2006:  fixed bug (had forgotten to check for value type) */
      if (A->value_type == REAL)
	memcpy ((double*)new_values + k_new, 
		(double*)new_values + k_cur, 
		cur_nnz * sizeof(double));
      else if (A->value_type == COMPLEX)
	memcpy ((double_Complex*)new_values + k_new, 
		(double_Complex*)new_values + k_cur, 
		cur_nnz * sizeof(double_Complex));

      /* fill in the symmetric "missing" values */
      while (k_cur < A->rowptr[i+1])
	{
	  /* non-zero of original matrix */
	  int j = A->colidx[k_cur];

	  if (j != i)  /* if not a non-diagonal element ... */
	    {
	      /* position of this transposed element in new matrix */
	      k_new = new_rowptr[j] + new_row_nnz[j];

	      /* store */
	      new_colidx[k_new] = i;
	      if (A->value_type == REAL)
		{
		  double* old_values = (double*) (A->values);
		  double* __new_values = (double*) new_values;
		  double a = old_values[k_cur];

		  if (A->symmetry_type == SYMMETRIC)
		    __new_values[k_new] = a;
		  else if (A->symmetry_type == SKEW_SYMMETRIC)
		    __new_values[k_new] = -a;
		  else if (A->symmetry_type == HERMITIAN)
		    __new_values[k_new] = a;
		}
	      else if (A->value_type == COMPLEX)
		{
		  double_Complex* old_values = (double_Complex*) (A->values);
		  double_Complex* __new_values = (double_Complex*) new_values;
		  double_Complex a = old_values[k_cur];

		  if (A->symmetry_type == SYMMETRIC)
		    __new_values[k_new] = a;
		  else if (A->symmetry_type == SKEW_SYMMETRIC)
		    __new_values[k_new] = double_Complex_negate (a);
		  else if (A->symmetry_type == HERMITIAN)
		    __new_values[k_new] = double_Complex_conj (a);
		}
		    
	      /*  update so next element stored at row j will appear
	       *  at the right place.  */
	      new_row_nnz[j]++;
	    }

	  k_cur++;
	}
    }

  /* clean-up */
  free (cur_row_nnz);
  free (new_row_nnz);

  A->nnz = new_nnz;
  smvm_free (A->values);
  A->values = new_values;
  smvm_free (A->colidx);
  A->colidx = new_colidx;
  smvm_free (A->rowptr);
  A->rowptr = new_rowptr;

  csr_matrix_sort_colidx (A);
  return 0;
}


int
print_csr_matrix_in_matlab_format (FILE* out, struct csr_matrix_t* A)
{
  int errcode = 0;
  struct coo_matrix_t* B = csr_to_coo (A);
  if (B == NULL)
    {
      fprintf (stderr, "*** print_csr_matrix_in_matlab_format: Failed to con"
	       "vert CSR-format sparse matrix into COO format for printing! "
	       "***\n");
      return -1;
    }
  errcode = print_coo_matrix_in_matlab_format (out, B);
  destroy_coo_matrix (B);
  return errcode;
}


int
save_csr_matrix_in_matlab_format (const char* const filename, struct csr_matrix_t* A)
{
  FILE* out = NULL;
  int errcode = 0;

  WITH_DEBUG2(fprintf(stderr, "=== save_csr_matrix_in_matlab_format ===\n"));
  out = fopen (filename, "w");
  if (out == NULL)
    {
      fprintf (stderr, "*** save_csr_matrix_in_matlab_format: failed "
	       "to open output file \"%s\" ***\n", filename);
      WITH_DEBUG2(fprintf(stderr, "=== Done with save_csr_matrix_in_matlab_format ===\n"));
      return -1;
    }

  errcode = print_csr_matrix_in_matlab_format (out, A);
  if (0 != fclose (out))
    {
      fprintf (stderr, "*** save_csr_matrix_in_matlab_format: failed "
	       "to close output file \"%s\" ***\n", filename);
      WITH_DEBUG2(fprintf(stderr, "=== Done with save_csr_matrix_in_matlab_format ===\n"));
      return -1;
    }
  WITH_DEBUG2(fprintf(stderr, "=== Done with save_csr_matrix_in_matlab_format ===\n"));
  return errcode;
}


int
valid_csr_matrix_p (const struct csr_matrix_t* A)
{
  const int m = A->m;
  const int n = A->n;
  const int nnz = A->nnz;
  const int* rowptr = A->rowptr;
  const int* colidx = A->colidx;
  int reached_nnz_in_rowptr_p = 0;
  int i, j;

  WITH_DEBUG2(fprintf(stderr, "=== valid_csr_matrix_p ===\n")); 

  if (m < 1)
    {
      fprintf(stderr,"*** valid_csr_matrix_p: Matrix has m = %d < 1 ***\n", m);
      WITH_DEBUG2(fprintf(stderr, "=== Done with valid_csr_matrix_p ===\n")); 
      return 0;
    }
  else if (n < 1)
    {
      fprintf(stderr,"*** valid_csr_matrix_p: Matrix has n = %d < 1 ***\n", nnz);
      WITH_DEBUG2(fprintf(stderr, "=== Done with valid_csr_matrix_p ===\n")); 
      return 0;
    }
  else if (nnz < 0)
    {
      fprintf(stderr,"*** valid_csr_matrix_p: Matrix has nnz = %d < 0 ***\n", nnz);
      WITH_DEBUG2(fprintf(stderr, "=== Done with valid_csr_matrix_p ===\n")); 
      return 0;
    }

  for (i = 0; i < m; i++)
    {
      const int cur = rowptr[i];
      const int next = rowptr[i+1];

      if (cur == nnz)
	{
	  if (! reached_nnz_in_rowptr_p)
	    reached_nnz_in_rowptr_p = 1;
	  if (next != cur)
	    {
	      fprintf (stderr, "*** valid_csr_matrix_p: Although rowptr[%d]==nnz==%d, "
		       "rowptr[%d+1]==%d != rowptr[%d] ***\n", 
		       i, nnz, i, next, i);
	      WITH_DEBUG2(fprintf(stderr, "=== Done with valid_csr_matrix_p ===\n")); 
	      return 0;
	    }
	}

      /* Verify that rowptr entries are in order */
      if (cur > next)
	{
	  fprintf (stderr, "*** valid_csr_matrix_p: Matrix: rowptr[%d] = %d > "
		   "rowptr[%d] = %d ***\n", i, cur, i+1, next);
	  WITH_DEBUG2(fprintf(stderr, "=== Done with valid_csr_matrix_p ===\n")); 
	  return 0;
	}

      /* Verify that current rowptr entry doesn't exceed nnz */
      if (cur >= nnz) /** \bug Should be cur >= nnz, right???? */
	{
	  fprintf (stderr, "*** valid_csr_matrix_p: Matrix: At col %d, rowptr[i] = "
		   "%d >= nnz = %d ***\n", i, rowptr[i], nnz);
	  WITH_DEBUG2(fprintf(stderr, "=== Done with valid_csr_matrix_p ===\n")); 
	  return 0;
	}

      /* Verify that colidx entries are in range */
      for (j = cur; j < next; j++)
	{
	  if (colidx[j] < 0 || colidx[j] >= n)
	    {
	      fprintf (stderr, "*** valid_csr_matrix_p: Matrix: at col %d, "
		       "colidx[%d]=%d out of range [%d,%d) ***\n", 
		       i, j, colidx[j], 0, n);
	      WITH_DEBUG2(fprintf(stderr, "=== Done with valid_csr_matrix_p ===\n")); 
	      return 0;
	    }
	}
    }

  /* Verify that the last rowptr entry is nnz */
  if (rowptr[m] != nnz)
    {
      fprintf (stderr, "*** valid_csr_matrix_p: rowptr[m=%d] = %d != nnz=%d ***\n", 
	       m, rowptr[m], nnz);
      WITH_DEBUG2(fprintf(stderr, "=== Done with valid_csr_matrix_p ===\n")); 
      return 0;
    }

  WITH_DEBUG2(fprintf(stderr, "=== Done with valid_csr_matrix_p ===\n"));
  return 1;
}


struct csr_matrix_t*
csr_matrix_matmatmult (struct csr_matrix_t* B, struct csr_matrix_t* A)
{
  /* B is m x p and A is p x n. */
  const int m = B->m;
  const int p = B->n;
  const int n = A->n;
  const enum value_type_t value_type = B->value_type;
  int errcode = 0;

  if (p != A->m)
    {
      fprintf (stderr, "*** csr_matrix_matmatmult: incompatible dimensions ***\n");
      return NULL;
    }
  else if (B->symmetry_type != UNSYMMETRIC || A->symmetry_type != UNSYMMETRIC)
    {
      fprintf (stderr, "*** csr_matrix_matmatmult: currently sparse matrix-matrix multiplication is only supported for matrices using unsymmetric storage ***\n");
      return NULL;
    }
  else if (B->value_type != A->value_type)
    {
      fprintf (stderr, "*** csr_matrix_matmatmult: currently sparse matrix-matrix multiplication is only supported for matrices of the same value type ***\n");
      return NULL;
    }

  if (value_type == REAL)
    {
      struct csr_matrix_t* C = smvm_malloc (sizeof (struct csr_matrix_t));
      double* Cvalues = (double*) C->values;
      double* Bvalues = (double*) B->values;
      double* Avalues = (double*) A->values;

      errcode = csr_matmatmult_double (&(C->rowptr), &(C->colidx), &Cvalues, &(C->nnz), 
				       1.0, 
				       B->rowptr, B->colidx, Bvalues,
				       A->rowptr, A->colidx, Avalues,
				       m, p, n);
      if (errcode != 0)
	{
	  fprintf (stderr, "*** csr_matrix_matmatmult: error code %d ***\n", errcode);
	  return NULL;
	}
      C->values = (void*) Cvalues;
      C->m = m;
      C->n = n;
      C->symmetry_type = UNSYMMETRIC;
      C->symmetric_storage_location = 0;
      C->value_type = REAL;
      return C;
    }
  else if (value_type == COMPLEX)
    {
      struct csr_matrix_t* C = smvm_malloc (sizeof (struct csr_matrix_t));
      const double_Complex one = new_double_Complex (1.0, 0.0);
      double_Complex* Cvalues = (double_Complex*) C->values;
      double_Complex* Bvalues = (double_Complex*) B->values;
      double_Complex* Avalues = (double_Complex*) A->values;

      errcode = csr_matmatmult_complex (&(C->rowptr), &(C->colidx), &Cvalues,
					&(C->nnz), one, B->rowptr, B->colidx, 
					Bvalues, A->rowptr, A->colidx, Avalues,
					m, p, n);
      if (errcode != 0)
	{
	  fprintf (stderr, "*** csr_matrix_matmatmult: error code %d ***\n", errcode);
	  return NULL;
	}
      C->values = (void*) Cvalues;
      C->m = m;
      C->n = n;
      C->symmetry_type = UNSYMMETRIC;
      C->symmetric_storage_location = 0;
      C->value_type = COMPLEX;
      return C;
    }
  else /* value_type == PATTERN */
    {
      struct csr_matrix_t* C = smvm_malloc (sizeof (struct csr_matrix_t));

      errcode = csr_matmatmult_pattern (&(C->rowptr), &(C->colidx), &(C->nnz), 
					B->rowptr, B->colidx, 
					A->rowptr, A->colidx, 
					m, p, n);
      if (errcode != 0)
	{
	  fprintf (stderr, "*** csr_matrix_matmatmult: error code %d ***\n", errcode);
	  return NULL;
	}
      C->m = m;
      C->n = n;
      C->symmetry_type = UNSYMMETRIC;
      C->symmetric_storage_location = 0;
      C->value_type = PATTERN;
      return C;
    }
}
