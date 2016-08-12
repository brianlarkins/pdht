/**
 * @file coo_matrix.c
 * @author Mark Hoemmen
 * @date 24 Feb 2006
 * @since Summer 2004
 *
 * Coordinate format sparse matrix data structure.
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
#include <coo_matrix.h>
#include <csc_matrix.h>
#include <read_mm.h>

#include <avltree_intpair.h>
#include <__complex.h>
#include <random_number.h>
#include <smvm_malloc.h>
#include <smvm_util.h>

#include <assert.h>
#include <stdio.h>
#include <string.h> /* memcpy */

/**
 * Works like (void*) (&array[idx]), in which array is treated as an
 * array of objects, each of which satisfies sizeof(object) == size.
 * "size" must be a variable or constant in scope ("deliberate variable
 * capture").
 */
#ifdef VOIDAREF
#  undef VOIDAREF
#endif /* VOIDAREF */
#define VOIDAREF( array, idx )  ((void*) ((char*) (array) + size*(idx)))                  


static int 
my_compare_lexicographically (const int x1, const int y1, const int x2, const int y2)
{
  if (x1 < y1)
    return -1;
  else if (x1 > y1)
    return +1;
  else
    {
      if (x2 < y2)
	return -1;
      else if (x2 > y2)
	return +1;
      else 
	return 0;
    }

  return 0; /* to pacify the compiler */
}


/**
 * If A has value_type REAL, workspace should be 
 * 2*max(sizeof(double), sizeof(int)) * r * c;
 * if A has value_type COMPLEX, workspace should be
 * 2*max(sizeof(double_Complex), sizeof(int)) * r * c.
 */
static int
partition (struct coo_matrix_t* A,
	   int low, 
	   int high, 
	   void *workspace)
{
  int left, right, pivot;
  int int_pivot_item1 = 0;
  int int_pivot_item2 = 0;
  int int_swapspace = 0;
  void *pivot_item = NULL; 
  void *swapspace = NULL;   
  size_t size = 0;
  const int block_size = 1;

  WITH_DEBUG2(fprintf(stderr, "=== partition ===\n"));
  WITH_DEBUG2(fprintf(stderr, "[low,high] = [%d,%d]\n", low, high));

  if (A->value_type == REAL)
    size = sizeof (double);
  else if (A->value_type == COMPLEX)
    size = sizeof (double_Complex);
  else 
    size = 0;

  pivot_item = workspace;
  swapspace = VOIDAREF(workspace, block_size);

#ifdef SWAP_INTS
#  undef SWAP_INTS
#endif /* SWAP_INTS */
#define SWAP_INTS( array, i, j )  do { \
  int_swapspace = array[i]; \
  array[i] = array[j]; \
  array[j] = int_swapspace; \
} while(0)

#ifdef SWAP
#  undef SWAP
#endif /* SWAP */
#define SWAP( array, i, j )  do { \
  if (A->value_type == REAL)  \
    { \
      double* __v = (double*) (array);                     \
      memcpy (swapspace, &__v[block_size * i], block_size * sizeof(double));            \
      memcpy (&__v[block_size * i], &__v[block_size * j], block_size * sizeof(double)); \
      memcpy (&__v[block_size * j], swapspace, block_size * sizeof(double));         \
    } \
  else if (A->value_type == COMPLEX) \
    { \
      double_Complex* __v = (double_Complex*) (array);     \
      memcpy (swapspace, &__v[block_size * i], block_size * sizeof(double_Complex)); \
      memcpy (&__v[block_size * i], &__v[block_size * j], block_size * sizeof(double_Complex));   \
      memcpy (&__v[block_size * j], swapspace, block_size * sizeof(double_Complex)); \
    } \
} while(0)

#ifdef ASSIGN
#  undef ASSIGN
#endif /* ASSIGN */
#define ASSIGN( voidp1, voidp2 )  do { \
  if (A->value_type == REAL)        \
    memcpy ((voidp1), (voidp2), block_size * sizeof(double)); \
  else if (A->value_type == COMPLEX) \
    memcpy ((voidp1), (voidp2), block_size * sizeof(double_Complex)); \
} while(0)

  /* mfh 21 Oct 2005: For randomized quicksort, pick a random integer i 
   * between low and high inclusive, and swap A[low] and A[i]. */
  {
    int i = smvm_random_integer (low, high);

    SWAP_INTS( (A->II), low, i );
    SWAP_INTS( (A->JJ), low, i );
    SWAP( (A->val), low, i );
  }         

  int_pivot_item1 = A->II[low];
  int_pivot_item2 = A->JJ[low];
  ASSIGN (pivot_item, VOIDAREF(A->val, block_size*low));

  pivot = left = low;
  right = high;

  while (left < right)
    {
      /* Move left while item < pivot */
      while (my_compare_lexicographically (A->II[left], int_pivot_item1, 
					   A->JJ[left], int_pivot_item2) <= 0 && 
	     left <= high)
	left++;

      /* Move right while item > pivot */
      while (my_compare_lexicographically (A->II[right], int_pivot_item1, 
					   A->JJ[right], int_pivot_item2) > 0 && 
	     right >= low)
	right--;

      if (left < right)
	{
	  SWAP_INTS( (A->II), left, right );
	  SWAP_INTS( (A->JJ), left, right );
	  SWAP( (A->val), left, right );
	}
    }

  /* right is final position for the pivot */
  A->II[low] = A->II[right];
  A->JJ[low] = A->JJ[right];
  ASSIGN( VOIDAREF(A->val, block_size*low), VOIDAREF(A->val, block_size*right) );

  A->II[right] = int_pivot_item1;
  A->JJ[right] = int_pivot_item2;
  ASSIGN( VOIDAREF(A->val, block_size*right), pivot_item );

  return right;
}

static void
quicksort (struct coo_matrix_t* A,
	   int low, 
	   int high, 
	   void* workspace)
{
  int pivot;

  WITH_DEBUG2(fprintf(stderr, "=== quicksort ===\n"));
  WITH_DEBUG2(fprintf(stderr, "[low,high] = [%d,%d]\n", low, high));

  if (high > low)
    {
      pivot = partition (A, low, high, workspace);
      quicksort (A, low, pivot - 1, workspace);
      quicksort (A, pivot + 1, high, workspace);
    }
}
	
static int
sort_coo_matrix_by_rows_then_columns (struct coo_matrix_t* A)
{
  void* workspace = NULL;

  WITH_DEBUG2(fprintf(stderr,"=== sort_coo_matrix_by_rows_then_columns ===\n"));

  if (A->value_type == REAL)
    workspace = smvm_calloc (3, MAX(sizeof(double), sizeof(int))); 
  else if (A->value_type == COMPLEX)
    workspace = smvm_calloc (3, MAX(sizeof(double), sizeof(int)));
  else if (A->value_type == PATTERN)
    workspace = NULL;
  else 
    {
      fprintf (stderr, "*** sort_bcoo_matrix_by_rows_then_columns: A has an invalid value_type! ***\n");
      return -1;
    }
  quicksort (A, 0, A->nnz - 1, workspace);
  if (workspace != NULL)
    {
      WITH_DEBUG2(fprintf(stderr, "\tsort_bcoo_matrix_by_rows_then_columns: freeing workspace\n"));
      smvm_free (workspace);
    }
  WITH_DEBUG2(fprintf(stderr,"=== Done with sort_bcoo_matrix_by_rows_then_columns ===\n"));

  return 0;
}



struct coo_matrix_t*
csc_to_coo_matrix (const struct csc_matrix_t* A, int index_base)
{
  struct coo_matrix_t *B = NULL;
  const int m = A->m;
  const int n = A->n;
  const int nnz = A->nnz;
  const int *colptr = A->colptr;
  const int *rowidx = A->rowidx;
  const double *values = A->values;
  int count = 0;
  int j, k;

  WITH_DEBUG2(fprintf(stderr, "=== csc_to_coo_matrix ==="));
  B = smvm_calloc (1, sizeof (struct coo_matrix_t));
  B->m = m;
  B->n = n;
  B->nnz = nnz;
  B->II = smvm_malloc (nnz * sizeof (int));
  B->JJ = smvm_malloc (nnz * sizeof (int));
  B->val = smvm_malloc (nnz * sizeof (double)); /* csc_matrix_t contains real values only */

  B->index_base = index_base;
  B->symmetry_type = A->symmetry_type;
  B->symmetric_storage_location = A->symmetric_storage_location;
  B->value_type = REAL; /* csc_matrix_t contains real values */

  if (index_base != 0 && index_base != 1)
    {
      fprintf (stderr,"*** csc_to_coo_matrix: invalid index base %d ***\n", index_base);
      smvm_exit (EXIT_FAILURE);
    }

  /* 18 oct 2004:  For some reason, in the j loop, the colptr array
   * starts being overwritten after a certain point !!!   Problem
   * solved:  you were allocating B->val with sizeof(int) instead of
   * sizeof(double), so writing to B->val was eventually overwriting the
   * colptr array. */

  for (j = 0; j < n; j++)
    {
      const int start = colptr[j];
      const int end   = colptr[j+1];

      if (start < 0 || end > nnz + 1 || start > end)
	{
	  fprintf(stderr,"*** ERROR: csc_to_coo_matrix: start,end wrong ***\n");
	  fprintf(stderr, "start,end = %d,%d\n", start, end);
	  smvm_exit (EXIT_FAILURE);
	}
                                                                                                                             
      for (k = start; k < end; k++)
        {
	  if (count >= nnz)
	    {
	      fprintf(stderr,"*** ERROR: csc_to_coo_matrix: count >= nnz ***\n");
	      fprintf(stderr, "j = %d, count = %d, nnz = %d\n", j, count, nnz);
	      smvm_exit (EXIT_FAILURE);
	    }
	  else if (count < 0)
	    {
	      fprintf(stderr,"*** ERROR: count < 0 ***\n");
	      smvm_exit (EXIT_FAILURE);
	    }
          B->II[count]   = rowidx[k] + index_base;
          B->JJ[count]   = j + index_base;
          ((double*) (B->val))[count] = values[k];
          count++;
        }
    }
                                                                                                                             
  WITH_DEBUG2(fprintf(stderr, "=== Done with csc_to_coo_matrix ==="));
  return B;
}


void
alloc_coo_matrix (struct coo_matrix_t* A, 
		  const int m, const int n, const int nnz,
		  enum index_base_t index_base, 
		  enum symmetry_type_t symmetry_type, 
		  enum symmetric_storage_location_t symmetric_storage_location,
		  enum value_type_t value_type)
{
  WITH_DEBUG2(fprintf(stderr, "=== alloc_coo_matrix ===\n"));

  A->m = m;
  A->n = n;
  A->nnz = nnz;

  A->II = smvm_calloc (nnz, sizeof (int));
  A->JJ = smvm_calloc (nnz, sizeof (int));
  if (value_type == REAL)
    A->val = smvm_calloc (nnz, sizeof (double));
  else if (value_type == COMPLEX)
    A->val = smvm_calloc (nnz, sizeof (double_Complex));
  else if (value_type == PATTERN)
    A->val = NULL;

  A->index_base = index_base;
  A->symmetry_type = symmetry_type;
  A->symmetric_storage_location = symmetric_storage_location;
  A->value_type = value_type;

  WITH_DEBUG2(fprintf(stderr, "=== Done with alloc_coo_matrix ===\n"));
}


void 
dealloc_coo_matrix (struct coo_matrix_t* A)
{
  WITH_DEBUG2(fprintf(stderr, "=== dealloc_coo_matrix ===\n"));
  if (A != NULL)
    {
      if (A->II != NULL) smvm_free (A->II);
      if (A->JJ != NULL) smvm_free (A->JJ);
      if (A->val != NULL) smvm_free (A->val);

      A->II = NULL;
      A->JJ = NULL;
      A->val = NULL;
      A->m = 0;
      A->n = 0;
      A->nnz = 0;
    }
  WITH_DEBUG2(fprintf(stderr, "=== Done with dealloc_coo_matrix ===\n"));
}


int
print_coo_matrix_in_matrix_market_format (FILE* out, const struct coo_matrix_t* A)
{
  int i;

  char symmetry_type_label[20];
  char value_type_label[20];

  WITH_DEBUG2(fprintf(stderr, "=== print_coo_matrix_in_matrix_market_format ===\n"));

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
      fprintf (stderr, "*** print_coo_matrix_in_matrix_market_format: "
	       "Invalid symmetry type %d of the given coo_matrix_t sparse "
	       "matrix! ***\n", A->symmetry_type);
      WITH_DEBUG2(fprintf(stderr, "=== Done with print_coo_matrix_in_matrix_market_format ===\n"));
      return -1;
    }

  if (A->value_type == REAL)
    strncpy (value_type_label, "real", 19);
  else if (A->value_type == COMPLEX)
    strncpy (value_type_label, "complex", 19);
  else if (A->value_type == PATTERN)
    strncpy (value_type_label, "pattern", 19);
  else
    {
      fprintf (stderr, "*** print_coo_matrix_in_matrix_market_format: "
	       "Unsupported value type! ***\n");
      WITH_DEBUG2(fprintf(stderr, "=== Done with print_coo_matrix_in_matrix_market_format ===\n"));
      return -1;
    }


  fprintf (out, "%%%%MatrixMarket matrix coordinate %s %s\n", value_type_label, symmetry_type_label);
  fprintf (out, "%d %d %d\n", A->m, A->n, A->nnz);

  if (A->value_type == REAL)
    {
      double* val = (double*) (A->val);

      /* MatrixMarket format uses 1-based indices, so we have to convert
       * to 1-based when we print out the matrix.  */
      for (i = 0; i < A->nnz; i++)
	fprintf (out, "%d %d %.13e\n", A->II[i] + (1 - A->index_base), A->JJ[i] + (1 - A->index_base), val[i]);
    }
  else if (A->value_type == COMPLEX)
    {
      double_Complex* val = (double_Complex*) (A->val);

      for (i = 0; i < A->nnz; i++)
	fprintf (out, "%d %d %.13e %.13e\n", 
		 A->II[i] + (1 - A->index_base), 
		 A->JJ[i] + (1 - A->index_base), 
		 double_Complex_real_part(val[i]), 
		 double_Complex_imag_part(val[i]));
    }
  else if (A->value_type == PATTERN)
    {
      for (i = 0; i < A->nnz; i++)
	fprintf (out, "%d %d\n", A->II[i] + (1 - A->index_base), A->JJ[i] + (1 - A->index_base));
    }

  WITH_DEBUG2(fprintf(stderr, "=== Done with print_coo_matrix_in_matrix_market_format ===\n"));
  return 0;
}


void
destroy_coo_matrix (struct coo_matrix_t* A)
{
  WITH_DEBUG2(fprintf(stderr, "=== destroy_coo_matrix ===\n"));
  if (A != NULL)
    {
      if (A->II != NULL) smvm_free (A->II);
      if (A->JJ != NULL) smvm_free (A->JJ);
      if (A->val != NULL) smvm_free (A->val);

      smvm_free (A);
    } 
  WITH_DEBUG2(fprintf(stderr, "=== Done with destroy_coo_matrix ===\n"));
}


void
coo_c_to_fortran (struct coo_matrix_t* A)
{
  const int nnz = A->nnz;
  int* II = A->II;
  int* JJ = A->JJ;
  int k;

  WITH_DEBUG2(fprintf(stderr, "=== coo_c_to_fortran ===\n"));
  if (A->index_base == ONE)
    {
      fprintf(stderr,"*** coo_c_to_fortran: matrix already is one-based ***\n");
      return;
    }

  for (k = 0; k < nnz; k++)
    {
      II[k]++;
      JJ[k]++;
    }

  A->index_base = ONE;
  WITH_DEBUG2(fprintf(stderr, "=== Done with coo_c_to_fortran ===\n"));
}


void
coo_fortran_to_c (struct coo_matrix_t* A)
{
  const int nnz = A->nnz;
  int* II = A->II;
  int* JJ = A->JJ;
  int k;

  WITH_DEBUG2(fprintf(stderr, "=== coo_fortran_to_c ===\n"));
  if (A->index_base == ZERO)
    {
      fprintf(stderr,"*** coo_fortran_to_c: matrix already is zero-based ***\n");
      return;
    }

  for (k = 0; k < nnz; k++)
    {
      II[k]--;
      JJ[k]--;
    }

  A->index_base = ZERO;
  WITH_DEBUG2(fprintf(stderr, "=== Done with coo_fortran_to_c ===\n"));
}


int
valid_coo_matrix_p (struct coo_matrix_t *A)
{
  /* pairset: Set of pairs of indices already seen in the matrix */
  avl_tree_intpair pairset = NULL; 
  int m, n, nnz, k;
  int *II;
  int *JJ;
  void *val;
  enum index_base_t index_base;

  WITH_DEBUG2(fprintf(stderr, "=== valid_coo_matrix_p ===\n"));
  if (A == NULL) 
    {
      WITH_DEBUG2(fprintf(stderr, "=== Done with valid_coo_matrix_p ===\n"));
      return 0;
    }

  m = A->m;
  if (m <= 0) 
    {
      WITH_DEBUG2(fprintf(stderr, "=== Done with valid_coo_matrix_p ===\n"));
      return 0;
    }

  n = A->n;
  if (n <= 0) 
    {
      WITH_DEBUG2(fprintf(stderr, "=== Done with valid_coo_matrix_p ===\n"));
      return 0;
    }

  nnz = A->nnz;
  if (nnz < 0) 
    {
      WITH_DEBUG2(fprintf(stderr, "=== Done with valid_coo_matrix_p ===\n"));
      return 0;
    }
  if (nnz == 0) 
    {
      WITH_DEBUG2(fprintf(stderr, "=== Done with valid_coo_matrix_p ===\n"));
      return 1; /* If the matrix is empty, we've checked enough */
    }

  /* Ensure that the arrays are non-NULL */
  II = A->II;
  if (II == NULL) 
    {
      WITH_DEBUG2(fprintf(stderr, "=== Done with valid_coo_matrix_p ===\n"));
      return 0;
    }
  JJ = A->JJ;
  if (JJ == NULL) 
    { 
      WITH_DEBUG2(fprintf(stderr, "=== Done with valid_coo_matrix_p ===\n"));
      return 0;
    }
  val = A->val;
  if (val == NULL) 
    {
      WITH_DEBUG2(fprintf(stderr, "=== Done with valid_coo_matrix_p ===\n"));
      return 0;
    }

  /* Ensure that the index base is valid */
  index_base = A->index_base;
  if (index_base != ZERO && index_base != ONE)
    {
      fprintf(stderr,"*** valid_coo_matrix_p: invalid index base ***\n");
      WITH_DEBUG2(fprintf(stderr, "=== Done with valid_coo_matrix_p ===\n"));
      return 0;
    }

  /* Initialize the set of pairs of indices */
  pairset = make_empty_avl_tree_intpair (NULL);

  /* 
   * Ensure that all the indices are within range, and check for duplicate 
   * indices.  Account for the index base.
   */
  for (k = 0; k < nnz; k++)
    {
      struct int_pair pair;

      pair.first = II[k] - index_base;
      pair.second = JJ[k] - index_base;

      if (pair.first < 0 || pair.first >= m)
	{
	  fprintf(stderr, "*** valid_coo_matrix_p: at k = %d, row index %d "
		 "out of valid range %d:%d ***", 
		 k, pair.first + index_base, index_base, m-1 + index_base);

	  pairset = make_empty_avl_tree_intpair (pairset);
	  WITH_DEBUG2(fprintf(stderr, "=== Done with valid_coo_matrix_p ===\n"));
	  return 0;
	}
      else if (pair.second < 0 || pair.second >= n)
	{
	  fprintf(stderr, "*** valid_coo_matrix_p: at k = %d, column index %d "
		 "out of valid range %d:%d ***", 
		 k, pair.second + index_base, index_base, n-1 + index_base);

	  pairset = make_empty_avl_tree_intpair (pairset);
	  WITH_DEBUG2(fprintf(stderr, "=== Done with valid_coo_matrix_p ===\n"));
	  return 0;
	}

      /* Check if the pair of indices has already been seen */
      if (NULL != find_intpair (pair, pairset))
	{
	  fprintf(stderr,"*** valid_coo_matrix_p: at k = %d, index pair (%d,%d) "
		 "has already been seen in the matrix! ***", 
		 k, pair.first + index_base, pair.second + index_base);
	  pairset = make_empty_avl_tree_intpair (pairset);
	  WITH_DEBUG2(fprintf(stderr, "=== Done with valid_coo_matrix_p ===\n"));
	  return 0;
	}
      else
	pairset = insert_intpair (pair, pairset);
    }

  pairset = make_empty_avl_tree_intpair (pairset);
  WITH_DEBUG2(fprintf(stderr, "=== Done with valid_coo_matrix_p ===\n"));
  return 1;  /* We've made it through the gauntlet of tests */
}


struct coo_matrix_t*
copy_coo_matrix (const struct coo_matrix_t* A)
{
  const int nnz = A->nnz;
  const int m   = A->m;
  const int n   = A->n;
  const enum index_base_t index_base = A->index_base;
  const enum symmetry_type_t symmetry_type = A->symmetry_type;
  const enum symmetric_storage_location_t symmetric_storage_location = A->symmetric_storage_location;
  const enum value_type_t value_type = A->value_type;
  struct coo_matrix_t *B = smvm_calloc (1, sizeof (struct coo_matrix_t));

  WITH_DEBUG2(fprintf(stderr, "=== copy_coo_matrix ===\n"));

  B->m   = m;
  B->n   = n;
  B->nnz = nnz;
  B->index_base = index_base;
  B->symmetry_type = symmetry_type;
  B->symmetric_storage_location = symmetric_storage_location;
  B->value_type = value_type;
  B->II = smvm_malloc (nnz * sizeof (int));
  B->JJ = smvm_malloc (nnz * sizeof (int));

  if (value_type == REAL)
    B->val = smvm_malloc (nnz * sizeof (double));
  else if (value_type == COMPLEX)
    B->val = smvm_malloc (nnz * sizeof (double_Complex));
  else if (value_type == PATTERN)
    B->val = NULL;

  memcpy (B->II, A->II, nnz * sizeof (int));
  memcpy (B->JJ, A->JJ, nnz * sizeof (int));

  if (value_type == REAL)
    memcpy (B->val, A->val, nnz * sizeof (double));
  else if (value_type == COMPLEX)
    memcpy (B->val, A->val, nnz * sizeof (double_Complex));

  WITH_DEBUG2(fprintf(stderr, "=== Done with copy_coo_matrix ===\n"));
  return B;
}


int
coo_matrix_equal_p (struct coo_matrix_t *A, struct coo_matrix_t *B)
{
  const int nnz = A->nnz;
  int errcode = 0;
  int k;

  WITH_DEBUG2(fprintf(stderr, "=== coo_matrix_equal_p ===\n"));

  assert (A != NULL);
  assert (B != NULL);

  if (A->m != B->m) return 0;
  if (A->n != B->n) return 0;
  if (A->nnz != B->nnz) return 0;
  if (A->index_base != B->index_base) return 0;
  if (A->symmetry_type != B->symmetry_type) return 0;
  if (A->symmetric_storage_location != B->symmetric_storage_location) return 0;
  if (A->value_type != B->value_type) return 0;

  /*
   * Put the rows and columns of A and B in a canonical order so we can 
   * easily test if all the entries are equal.
   */

  errcode = sort_coo_matrix_by_rows_then_columns (A);
  if (errcode != 0)
    {
      fprintf (stderr, "*** coo_matrix_equal_p: Failed to sort COO format matrix A ***\n");
      exit (EXIT_FAILURE);
    }
  errcode = sort_coo_matrix_by_rows_then_columns (B);
  if (errcode != 0)
    {
      fprintf (stderr, "*** coo_matrix_equal_p: Failed to sort COO format matrix B ***\n");
      exit (EXIT_FAILURE);
    }
  
  if (A->value_type == REAL)
    {
      const double* const aval = (double*) (A->val);
      const double* const bval = (double*) (B->val);

      for (k = 0; k < nnz; k++)
	{
	  if (A->II[k] != B->II[k]) return 0;
	  if (A->JJ[k] != B->JJ[k]) return 0;
	  if (aval[k] != bval[k]) return 0;
	}
    }
  else if (A->value_type == COMPLEX)
    {
      const double_Complex* const aval = (double_Complex*) (A->val);
      const double_Complex* const bval = (double_Complex*) (B->val);

      for (k = 0; k < nnz; k++)
	{
	  if (A->II[k] != B->II[k]) return 0;
	  if (A->JJ[k] != B->JJ[k]) return 0;
	  if (double_Complex_not_equal(aval[k], bval[k])) return 0;
	}
    }
  else if (A->value_type == PATTERN)
    {
      for (k = 0; k < nnz; k++)
	{
	  if (A->II[k] != B->II[k]) return 0;
	  if (A->JJ[k] != B->JJ[k]) return 0;
	}
    }

  return 1;
}


void
init_coo_matrix (struct coo_matrix_t *A, int m, int n, int nnz, int *II, 
		 int *JJ, void *val, enum index_base_t index_base,
		 enum symmetry_type_t symmetry_type, 
		 enum symmetric_storage_location_t symmetric_storage_location,
		 enum value_type_t value_type)
{
  WITH_DEBUG2(fprintf(stderr, "=== init_coo_matrix ===\n"));
  A->m = m;
  A->n = n;
  A->nnz = nnz;
  A->II   = II;
  A->JJ   = JJ;
  A->val = val;
  A->index_base = index_base;
  A->symmetry_type = symmetry_type;
  A->symmetric_storage_location = symmetric_storage_location;
  A->value_type = value_type;
  WITH_DEBUG2(fprintf(stderr, "=== Done with init_coo_matrix ===\n"));
}


/** 
* Checks whether the indices in the given COO-format sparse matrix are 
* within their valid ranges, for Fortran-style (1-based) indices.
*
* @param m [IN]    Number of rows in the matrix
* @param n [IN]    Number of columns in the matrix
* @param nnz [IN]  Number of structural nonzeros (stored elements)
* @param II [IN]    Array of row indices
* @param JJ [IN]    Array of column indices
*
* @return Nonzero if indices in range, 0 if an index is out of range.
*/
int
coo_matrix_in_fortran_format_p (const int m, const int n, const int nnz, 
				const int II[], const int JJ[])
{
  int k;
  int bad = 0;

  WITH_DEBUG2(fprintf(stderr, "=== coo_matrix_in_fortran_format_p ===\n"));
  for (k = 0; k < nnz; k++)
    {
      int i = II[k];
      int j = JJ[k];

      if (i < 1 || i > m || j < 1 || j > n)
	{
	  fprintf (stderr, "\n*** Element k=%d of array is out of valid index range, "
		 "using Fortran-style indexing: m,n = %d,%d, (i,j) = (%d,%d)"
		 " ***\n", k, m, n, i, j);
	  bad = 1;
	}
    }
  WITH_DEBUG2(fprintf(stderr, "=== Done with coo_matrix_in_fortran_format_p ===\n"));
  return (bad == 0);
}


struct coo_matrix_t*
reserve_coo_matrix (const int m, const int n, const int nnz, 
		    enum index_base_t index_base,
		    enum symmetry_type_t symmetry_type, 
		    enum symmetric_storage_location_t symmetric_storage_location,
		    enum value_type_t value_type)
{
  struct coo_matrix_t* A = smvm_malloc (sizeof (struct coo_matrix_t));
  int *II = smvm_malloc (nnz * sizeof (int));
  int *JJ = smvm_malloc (nnz * sizeof (int));
  void *val = NULL;
  
  WITH_DEBUG2(fprintf(stderr, "=== reserve_coo_matrix ===\n"));
  if (value_type == REAL)
    val = smvm_malloc (nnz * sizeof (double));
  else if (value_type == COMPLEX)
    val = smvm_malloc (nnz * sizeof (double_Complex));
  else if (value_type == PATTERN)
    val = NULL;

  init_coo_matrix (A, m, n, nnz, II, JJ, val, index_base, symmetry_type, 
		   symmetric_storage_location, value_type);

  WITH_DEBUG2(fprintf(stderr, "=== Done with reserve_coo_matrix ===\n"));
  return A;
}


int
save_coo_matrix_in_harwell_boeing_format (const char* const filename, 
					  const struct coo_matrix_t* A)
{
  int errcode = 0;
  struct csc_matrix_t* B = smvm_malloc (sizeof (struct csc_matrix_t));

  WITH_DEBUG2(fprintf(stderr, "=== save_coo_matrix_in_harwell_boeing_format ===\n"));

  /* We already have a function for saving CSC-format sparse matrices in
   * Harwell-Boeing format, so just create a CSC-format version of A and call
   * the appropriate save function on it. */

  /* coo_to_csc_matrix is a bit hackish:  The struct B has to be allocated, but
   * its internal storage should not previously be allocated (this function does
   * that for us). */
  coo_to_csc_matrix (B, A);

  errcode = save_csc_matrix_in_harwell_boeing_format (filename, B);

  /* Get rid of the CSC-format copy. */
  destroy_csc_matrix (B);

  WITH_DEBUG2(fprintf(stderr, "=== Done with save_coo_matrix_in_harwell_boeing_format ===\n"));
  return errcode;
}

int
save_coo_matrix_in_matrix_market_format (const char* const filename, const struct coo_matrix_t* A)
{
  FILE* out = NULL;
  int errcode = 0;

  WITH_DEBUG2(fprintf(stderr, "=== save_coo_matrix_in_matrix_market_format ===\n"));
  out = fopen (filename, "w");
  if (out == NULL)
    {
      fprintf (stderr, "*** save_coo_matrix_in_matrix_market_format: failed "
	       "to open output file \"%s\" ***\n", filename);
      WITH_DEBUG2(fprintf(stderr, "=== Done with save_coo_matrix_in_matrix_market_format ===\n"));
      return -1;
    }

  errcode = print_coo_matrix_in_matrix_market_format (out, A);
  if (0 != fclose (out))
    {
      fprintf (stderr, "*** save_coo_matrix_in_matrix_market_format: failed "
	       "to close output file \"%s\" ***\n", filename);
      WITH_DEBUG2(fprintf(stderr, "=== Done with save_coo_matrix_in_matrix_market_format ===\n"));
      return -1;
    }
  WITH_DEBUG2(fprintf(stderr, "=== Done with save_coo_matrix_in_matrix_market_format ===\n"));
  return errcode;
}


int
coo_matrix_expand_symmetric_storage (struct coo_matrix_t* A)
{
  int old_nnz = A->nnz;
  int new_nnz = A->nnz; /* To be changed */
  int num_diag_entries = 0;
  int k;
  int current_idx = A->nnz;

  WITH_DEBUG2(fprintf(stderr, "=== coo_matrix_expand_symmetric_storage ===\n"));

  if (A->symmetry_type == UNSYMMETRIC)
    {
      WITH_DEBUG(fprintf(stderr, "*** coo_matrix_expand_symmetric_storage: ma"
			 "trix is already stored in an unsymmetric format! ***\n"));
      return 0; /* Not an error -- just means that we don't have any work to do */
    }

  /* Count the number of nonzeros on the diagonal.  These entries do not need to be replicated. */
  for (k = 0; k < old_nnz; k++)
    {
      if (A->II[k] == A->JJ[k])
	num_diag_entries++;
    }

  new_nnz = 2*old_nnz - num_diag_entries;
  assert (new_nnz >= old_nnz);
  assert (new_nnz <= 2*old_nnz);

  /* Expand the arrays */
  A->II = smvm_realloc (A->II, new_nnz * sizeof (int));
  A->JJ = smvm_realloc (A->JJ, new_nnz * sizeof (int));
  if (A->value_type == REAL)
    A->val = smvm_realloc (A->val, new_nnz * sizeof (double));
  else if (A->value_type == COMPLEX)
    A->val = smvm_realloc (A->val, new_nnz * sizeof (double_Complex));

  /* Expand the off-diagonal elements, so that if (i,j) is in the matrix, (j,i)
   * is there also.  This algorithm assumes correctness of the symmetric
   * representation.  Depending upon whether the matrix is symmetric or
   * skew-symmetric, we reflect the element resp. the negative of the element
   * across the diagonal.  (The symmetric/skew-symmetric choice is outside the 
   * for-loop for efficiency (to avoid an extra branch in the inner loop).) */

  current_idx = old_nnz;
  if (A->symmetry_type == SYMMETRIC)
    {
      for (k = 0; k < old_nnz; k++)
	{
	  if (A->II[k] != A->JJ[k]) /* If not a diagonal element */
	    {
	      const int i = A->II[k];
	      const int j = A->JJ[k];

	      A->II[current_idx] = j;
	      A->JJ[current_idx] = i;
	      if (A->value_type == REAL)
		{
		  double* val = (double*) (A->val);
		  double temp = val[k];
		  val[current_idx] = temp;
		}
	      else if (A->value_type == COMPLEX)
		{
		  double_Complex* val = (double_Complex*) (A->val);
		  double_Complex temp = val[k];
		  val[current_idx] = temp;
		}

	      current_idx++;
	    }
	}
    }
  else if (A->symmetry_type == SKEW_SYMMETRIC)
    {
      for (k = 0; k < old_nnz; k++)
	{
	  if (A->II[k] != A->JJ[k]) /* If not a diagonal element */
	    {
	      const int i = A->II[k];
	      const int j = A->JJ[k];

	      A->II[current_idx] = j;
	      A->JJ[current_idx] = i;
	      if (A->value_type == REAL)
		{
		  double* val = (double*) (A->val);
		  double temp = val[k];
		  val[current_idx] = -temp;
		}
	      else if (A->value_type == COMPLEX)
		{
		  double_Complex* val = (double_Complex*) (A->val);
		  double_Complex temp = val[k];
		  val[current_idx] = double_Complex_negate(temp);
		}

	      current_idx++;
	    }
	}
    }
  else if (A->symmetry_type == HERMITIAN)
    {
      for (k = 0; k < old_nnz; k++)
	{
	  if (A->II[k] != A->JJ[k]) /* If not a diagonal element */
	    {
	      const int i = A->II[k];
	      const int j = A->JJ[k];

	      A->II[current_idx] = j;
	      A->JJ[current_idx] = i;
	      if (A->value_type == REAL)
		{
		  /* Real matrices can be labeled as "hermitian," though they really 
		   * should be called "symmetric." */
		  double* val = (double*) (A->val);
		  double temp = val[k];
		  val[current_idx] = temp;
		}
	      else if (A->value_type == COMPLEX)
		{
		  double_Complex* val = (double_Complex*) (A->val);
		  double_Complex temp = val[k];
		  val[current_idx] = double_Complex_conj (temp);
		}

	      current_idx++;
	    }
	}

    }

  /* Change the symmetry type to unsymmetric */
  A->symmetry_type = UNSYMMETRIC;
  A->nnz = new_nnz;

  WITH_DEBUG2(fprintf(stderr, "=== coo_matrix_expand_symmetric_storage ===\n"));
  return 0;
}


int
print_coo_matrix_in_matlab_format (FILE* out, struct coo_matrix_t* A)
{
  const int m = A->m;
  const int n = A->n;
  const int index_base = A->index_base;
  int found_mn_elt_p = 0;
  struct coo_matrix_t* B = NULL;
  int i;

  WITH_DEBUG2(fprintf(stderr, "=== print_coo_matrix_in_matlab_format ===\n"));

  if (A->symmetry_type != UNSYMMETRIC)
    {
      B = coo_matrix_expand_symmetric_storage_copy (A);
      if (B == NULL)
	{
	  fprintf (stderr, "=== print_coo_matrix_in_matlab_format: Failed to "
		   "expand matrix A into unsymmetric storage! ***\n");
	  return -1;
	}
    }

  /* 
   * In a Matlab format file, the matrix dimensions m and n are not explicitly
   * stored, so if the last row or column of the matrix contains all zeros, we
   * have to output an explicit zero at (m,n) in order that the matrix
   * dimensions are correctly stored.  So we need to search the matrix for an
   * (m,n) entry (rather, (m - index_base, n - index_base), and if there is not
   * already one, then we add one with an explicit zero value.
   */
  if (B != NULL)
    {
      for (i = 0; i < B->nnz; i++)
	{
	  if (B->II[i] == m - index_base && B->JJ[i] == n - index_base)
	    found_mn_elt_p = 1;
	}

      if (! found_mn_elt_p)
	{
	  B->II = smvm_realloc (B->II, (B->nnz + 1) * sizeof (int));
	  B->JJ = smvm_realloc (B->JJ, (B->nnz + 1) * sizeof (int));

	  B->II[B->nnz] = m - index_base;
	  B->JJ[B->nnz] = n - index_base;

	  if (B->value_type == REAL)
	    {
	      double* val;

	      B->val = smvm_realloc (B->val, (B->nnz + 1) * sizeof (double));
	      val = (double*) (B->val);
	      val[B->nnz] = 0.0;
	    }
	  else if (B->value_type == COMPLEX)
	    {
	      double_Complex* val;

	      B->val = smvm_realloc (B->val, (B->nnz + 1) * sizeof (double_Complex));
	      val = (double_Complex*) (B->val);
	      val[B->nnz] = double_Complex_ZERO;
	    }

	  B->nnz = B->nnz + 1;
	}
    }
  else
    {
      for (i = 0; i < A->nnz; i++)
	{
	  if (A->II[i] == m - index_base && A->JJ[i] == n - index_base)
	    found_mn_elt_p = 1;
	}

      if (! found_mn_elt_p)
	{
	  A->II = smvm_realloc (A->II, (A->nnz + 1) * sizeof (int));
	  A->JJ = smvm_realloc (A->JJ, (A->nnz + 1) * sizeof (int));

	  A->II[A->nnz] = m - index_base;
	  A->JJ[A->nnz] = n - index_base;

	  if (A->value_type == REAL)
	    {
	      double* val;

	      A->val = smvm_realloc (A->val, (A->nnz + 1) * sizeof (double));
	      val = (double*) (A->val);
	      val[A->nnz] = 0.0;
	    }
	  else if (A->value_type == COMPLEX)
	    {
	      double_Complex* val;

	      A->val = smvm_realloc (A->val, (A->nnz + 1) * sizeof (double_Complex));
	      val = (double_Complex*) (A->val);
	      val[A->nnz] = double_Complex_ZERO;
	    }

	  A->nnz = A->nnz + 1;
	}

    }

  /*
   * Output the (i,j) indices and values.
   */

  if (A->value_type == REAL)
    {
      if (B != NULL)
	{
	  const double* const val = (const double* const) (B->val);

	  /* Matlab format uses 1-based indices, so we have to convert
	   * to 1-based when we print out the matrix.  */
	  for (i = 0; i < B->nnz; i++)
	    fprintf (out, "%d %d %.13e\n", B->II[i] + (1 - B->index_base), B->JJ[i] + (1 - B->index_base), val[i]);
	}
      else 
	{
	  const double* const val = (const double* const) (A->val);

	  for (i = 0; i < A->nnz; i++)
	    fprintf (out, "%d %d %.13e\n", A->II[i] + (1 - A->index_base), A->JJ[i] + (1 - A->index_base), val[i]);
	}
    }
  else if (A->value_type == COMPLEX)
    {
      if (B != NULL)
	{
	  const double_Complex* const val = (const double_Complex* const) (B->val);

	  for (i = 0; i < B->nnz; i++)
	    fprintf (out, "%d %d %.13e %.13e\n", 
		     B->II[i] + (1 - B->index_base), 
		     B->JJ[i] + (1 - B->index_base), 
		     double_Complex_real_part(val[i]), 
		     double_Complex_imag_part(val[i]));
	}
      else
	{
	  const double_Complex* const val = (const double_Complex* const) (A->val);

	  for (i = 0; i < A->nnz; i++)
	    fprintf (out, "%d %d %.13e %.13e\n", 
		     A->II[i] + (1 - A->index_base), 
		     A->JJ[i] + (1 - A->index_base), 
		     double_Complex_real_part(val[i]), 
		     double_Complex_imag_part(val[i]));
	}
    }
  else if (A->value_type == PATTERN)
    {
      if (B != NULL)
	{
	  for (i = 0; i < B->nnz; i++)
	    fprintf (out, "%d %d\n", B->II[i] + (1 - B->index_base), B->JJ[i] + (1 - B->index_base));
	}
      else
	{
	  for (i = 0; i < A->nnz; i++)
	    fprintf (out, "%d %d\n", A->II[i] + (1 - A->index_base), A->JJ[i] + (1 - A->index_base));
	}
    }

  if (B != NULL)
    destroy_coo_matrix (B);

  WITH_DEBUG2(fprintf(stderr, "=== Done with print_coo_matrix_in_matlab_format ===\n"));
  return 0;
}


int
save_coo_matrix_in_matlab_format (const char* const filename, struct coo_matrix_t* A)
{
  FILE* out = NULL;
  int errcode = 0;

  WITH_DEBUG2(fprintf(stderr, "=== save_coo_matrix_in_matlab_format ===\n"));
  out = fopen (filename, "w");
  if (out == NULL)
    {
      fprintf (stderr, "*** save_coo_matrix_in_matlab_format: failed "
	       "to open output file \"%s\" ***\n", filename);
      WITH_DEBUG2(fprintf(stderr, "=== Done with save_coo_matrix_in_matlab_format ===\n"));
      return -1;
    }

  errcode = print_coo_matrix_in_matlab_format (out, A);
  if (0 != fclose (out))
    {
      fprintf (stderr, "*** save_coo_matrix_in_matlab_format: failed "
	       "to close output file \"%s\" ***\n", filename);
      WITH_DEBUG2(fprintf(stderr, "=== Done with save_coo_matrix_in_matlab_format ===\n"));
      return -1;
    }
  WITH_DEBUG2(fprintf(stderr, "=== Done with save_coo_matrix_in_matlab_format ===\n"));
  return errcode;
}



struct coo_matrix_t*
coo_matrix_expand_symmetric_storage_copy (struct coo_matrix_t* A)
{
  struct coo_matrix_t* B = smvm_calloc (1, sizeof (struct coo_matrix_t));
  int old_nnz = A->nnz;
  int new_nnz = A->nnz; /* To be changed */
  int num_diag_entries = 0;
  int k;
  int current_idx = A->nnz;

  WITH_DEBUG2(fprintf(stderr, "=== coo_matrix_expand_symmetric_storage_copy ===\n"));

  if (A->symmetry_type == UNSYMMETRIC)
    {
      WITH_DEBUG(fprintf(stderr, "*** coo_matrix_expand_symmetric_storage_cop"
			 "y: Matrix A is already in unsymmetric format -- mak"
			 "ing copy ***\n"));
      return copy_coo_matrix (A);
    }

  /* Count the number of nonzeros on the diagonal.  These entries do not need to be replicated. */
  for (k = 0; k < old_nnz; k++)
    {
      if (A->II[k] == A->JJ[k])
	num_diag_entries++;
    }

  new_nnz = 2*old_nnz - num_diag_entries;
  assert (new_nnz >= old_nnz);
  assert (new_nnz <= 2*old_nnz);

  init_coo_matrix (B, A->m, A->n, new_nnz, NULL, NULL, NULL, A->index_base, 
		   UNSYMMETRIC, 0, A->value_type);

  /* Allocate the arrays */
  B->II = smvm_calloc (new_nnz, sizeof (int));
  memcpy (B->II, A->II, new_nnz * sizeof (int));
  B->JJ = smvm_calloc (new_nnz, sizeof (int));
  memcpy (B->JJ, A->JJ, new_nnz * sizeof (int));

  if (A->value_type == REAL)
    {
      B->val = smvm_calloc (new_nnz, sizeof (double));
      memcpy (B->val, A->val, new_nnz * sizeof (double));
    }
  else if (A->value_type == COMPLEX)
    {
      B->val = smvm_calloc (new_nnz, sizeof (double_Complex));
      memcpy (B->val, A->val, new_nnz * sizeof (double_Complex));
    }

  /* Expand the off-diagonal elements, so that if (i,j) is in the matrix, (j,i)
   * is there also.  This algorithm assumes correctness of the symmetric
   * representation.  Depending upon whether the matrix is symmetric or
   * skew-symmetric, we reflect the element resp. the negative of the element
   * across the diagonal.  */

  current_idx = old_nnz;
  if (A->symmetry_type == SYMMETRIC)
    {
      for (k = 0; k < old_nnz; k++)
	{
	  if (A->II[k] != A->JJ[k]) /* If not a diagonal element */
	    {
	      const int i = A->II[k];
	      const int j = A->JJ[k];

	      B->II[current_idx] = j;
	      B->JJ[current_idx] = i;
	      if (A->value_type == REAL)
		{
		  double temp = ((double*) A->val)[k];
		  ((double*) B->val)[current_idx] = temp;
		}
	      else if (A->value_type == COMPLEX)
		{
		  double_Complex temp = ((double_Complex*) A->val)[k];
		  ((double_Complex*) B->val)[current_idx] = temp;
		}

	      current_idx++;
	    }
	}
    }
  else if (A->symmetry_type == SKEW_SYMMETRIC)
    {
      for (k = 0; k < old_nnz; k++)
	{
	  if (A->II[k] != A->JJ[k]) /* If not a diagonal element */
	    {
	      const int i = A->II[k];
	      const int j = A->JJ[k];

	      A->II[current_idx] = j;
	      A->JJ[current_idx] = i;
	      if (A->value_type == REAL)
		{
		  double temp = ((double*) A->val)[k];
		  ((double*) B->val)[current_idx] = -temp;
		}
	      else if (A->value_type == COMPLEX)
		{
		  double_Complex temp = ((double_Complex*) A->val)[k];
		  ((double_Complex*) B->val)[current_idx] = double_Complex_negate(temp);
		}

	      current_idx++;
	    }
	}
    }
  else if (A->symmetry_type == HERMITIAN)
    {
      for (k = 0; k < old_nnz; k++)
	{
	  if (A->II[k] != A->JJ[k]) /* If not a diagonal element */
	    {
	      const int i = A->II[k];
	      const int j = A->JJ[k];

	      A->II[current_idx] = j;
	      A->JJ[current_idx] = i;
	      if (A->value_type == REAL)
		{
		  /* Real matrices can be labeled as "hermitian," though they really 
		   * should be called "symmetric." */
		  double temp = ((double*) A->val)[k];
		  ((double*) B->val)[current_idx] = temp;
		}
	      else if (A->value_type == COMPLEX)
		{
		  double_Complex temp = ((double_Complex*) A->val)[k];
		  ((double_Complex*) B->val)[current_idx] = double_Complex_conj(temp);
		}

	      current_idx++;
	    }
	}

    }

  WITH_DEBUG2(fprintf(stderr, "=== coo_matrix_expand_symmetric_storage_copy ===\n"));
  return B;
}



struct coo_matrix_t* 
load_coo_matrix_in_matlab_format (const char* const filename)
{
  /* If there are only two columns of data, then we have a pattern matrix.  If
   * there are three columns, then the values are real.  If there are four
   * columns, then the values are complex (the third column has the real parts
   * and the fourth column the imaginary parts).  Any more columns than that 
   * should be ignored.
   */

  char line[400]; /* max line length is 399 */
  FILE *f = NULL;
  int m = 0; 
  int n = 0;
  int nnz = 0;
  int *II = NULL; 
  int *JJ = NULL;
  void *val = NULL;
  enum value_type_t value_type = -1; /* flag */
  int num_cols_read = 0;
  int num_cols_to_read = 0;
  int max_nnz = 0;
  int ii, jj;
  double xr, xi;

  WITH_DEBUG2(fprintf(stderr,"=== load_coo_matrix_in_matlab_format ===\n"));

  WITH_DEBUG2(fprintf(stderr, "\tOpening file\n" ));
  f = fopen (filename, "r");
  if (f == NULL)
    {
      fprintf(stderr, "*** load_coo_matrix_in_matlab_format: Failed to open "
	      "Matlab file %s ***\n", filename);
      return NULL;
    }

  WITH_DEBUG2(fprintf(stderr, "\tGetting and scanning the first data line\n" ));

  /* Get the first line.  
   *
   * The Linux Programmer's Manual says:  "Never use gets().  Because it is
   * impossible to tell without knowing the data in advance how many characters
   * gets() will read, and because gets() will continue to store characters
   * past the end of the buffer, it is extremely dangerous to use.  It has
   * been used to break computer security.  Use fgets() instead."
   *
   * fgets stores the newline, if it finds one.  */
  fgets (line, 399, f);

  /* Try to match a four-column line, and see how many columns were actually
   * read.  Note that whitespace in the scanf format string matches any amount
   * of whitespace, including none.  If there are two columns, then this is a
   * pattern matrix; if there are three columns, the values are real; if there
   * are four columns, the values are complex. */

  num_cols_read = sscanf (line, " %d %d %lg %lg", &ii, &jj, &xr, &xi);

  /* Skip whatever lines don't match until we find a line that does */
  while (num_cols_read <= 0 && ! feof (f))
    num_cols_read = sscanf (line, " %d %d %lg %lg", &ii, &jj, &xr, &xi);

  if (feof(f))
    {
      fprintf (stderr, "*** load_coo_matrix_in_matlab_format: File contains no data! ***\n");
      WITH_DEBUG2(fprintf(stderr, "\tClosing file\n" ));
      fclose (f);
      return NULL;
    }
  else if (num_cols_read < 2)
    {
      fprintf (stderr, "*** load_coo_matrix_in_matlab_format: File does not "
	       "contain at least two columns of data! ***\n");
      WITH_DEBUG2(fprintf(stderr, "\tClosing file\n" ));
      fclose (f);
      return NULL;
    }

  if (num_cols_read == 2)
    {
      value_type = PATTERN;
      num_cols_to_read = 2;
      WITH_DEBUG2(fprintf(stderr, "\t2 columns of data -- pattern matrix\n" ));
    }
  else if (num_cols_read == 3)
    {
      value_type = REAL;
      num_cols_to_read = 3;
      WITH_DEBUG2(fprintf(stderr, "\t3 columns of data -- real-valued matrix\n" ));
    }
  else if (num_cols_read >= 4)
    {
      value_type = COMPLEX;
      num_cols_to_read = 4;
      WITH_DEBUG2(fprintf(stderr, "\t4 columns of data -- complex-valued matrix\n" ));
    }

  WITH_DEBUG2(fprintf(stderr, "\tInitial allocation of II, JJ and val arrays\n" ));

  max_nnz = 16; /* Initial size of arrays */
  II = smvm_calloc (max_nnz, sizeof (int));
  JJ = smvm_calloc (max_nnz, sizeof (int));
  if (value_type == PATTERN)
    val = NULL;
  else if (value_type == REAL)
    val = smvm_calloc (max_nnz, sizeof (double));
  else if (value_type == COMPLEX)
    val = smvm_calloc (max_nnz, sizeof (double_Complex));

  /* Store the first read element */
  II[0] = ii;
  JJ[0] = jj;
  if (value_type == REAL)
    ((double*) (val))[0] = xr;
  else if (value_type == REAL)
    ((double_Complex*) (val))[0] = new_double_Complex(xr, xi);
  nnz++;

  /* Indices are one-based in the file.  The max i seen becomes m, and the max
   * j seen becomes n. */
  if (ii > m)
    m = ii;
  if (jj > n)
    n = jj;

  WITH_DEBUG2(fprintf(stderr, "\tReading the rest of the file\n" ));

  while (! feof(f))
    {
      fgets (line, 399, f);
      num_cols_read = sscanf (line, " %d %d %lg %lg", &ii, &jj, &xr, &xi);

      /* Skip whatever lines don't match until we find a line that does */
      while (num_cols_read < num_cols_to_read && ! feof (f))
	num_cols_read = sscanf (line, " %d %d %lg %lg", &ii, &jj, &xr, &xi);

      if (num_cols_read < num_cols_to_read)
	break; /* There is nothing left to read */

      /* Expand the arrays if necessary */
      if (nnz + 1 > max_nnz)
	{
	  II = smvm_realloc (II, 2 * (nnz+1) * sizeof (int));
	  JJ = smvm_realloc (JJ, 2 * (nnz+1) * sizeof (int));

	  if (value_type == REAL)
	    val = smvm_realloc (val, 2 * (nnz+1) * sizeof (double));
	  else if (value_type == REAL)
	    val = smvm_realloc (val, 2 * (nnz+1) * sizeof (double_Complex));
	}

      /* Extract the new nonzero element */
      II[nnz] = ii;
      JJ[nnz] = jj;
      if (value_type == REAL)
	((double*) (val))[nnz] = xr;
      else if (value_type == REAL)
	((double_Complex*) (val))[nnz] = new_double_Complex(xr, xi);
      nnz++;

      /* Indices are one-based in the file.  The max i seen becomes m, and the max
       * j seen becomes n. */
      if (ii > m)
	m = ii;
      if (jj > n)
	n = jj;
    }

  WITH_DEBUG2(fprintf(stderr, "\tClosing file\n" ));
  fclose (f);
  return create_coo_matrix (m, n, nnz, II, JJ, val, ONE, UNSYMMETRIC, -1, 
			    value_type);
}


struct coo_matrix_t*
create_coo_matrix (int m, int n, int nnz, int *II, 
		   int *JJ, void *val, enum index_base_t index_base, 
		   enum symmetry_type_t symmetry_type, 
		   enum symmetric_storage_location_t symmetric_storage_location,
		   enum value_type_t value_type)
{
  struct coo_matrix_t* A = smvm_calloc (1, sizeof (struct coo_matrix_t));

  init_coo_matrix (A, m, n, nnz, II, JJ, val, index_base, symmetry_type, 
		   symmetric_storage_location, value_type);
  return A;
}


