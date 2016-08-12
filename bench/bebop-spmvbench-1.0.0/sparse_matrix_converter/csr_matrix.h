/**
 * @file csr_matrix.h
 * @date 09 Jun 2006
 * @since 31 May 2005
 * @author Mark Hoemmen
 * @version 1.0
 * 
 * CSR (compressed sparse row) format sparse matrix member functions.
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
#ifndef _csr_matrix_h
#define _csr_matrix_h

#include <enumerations.h>
#include <stdio.h> /* FILE */


/** Evaluates to true iff A is a square matrix */
#define SQUARE_P( A )   (A->m == A->n)


/**
 * @struct csr_matrix_t
 * @author Mark Hoemmen
 * @since 31 May 2005
 *
 * Sparse matrix in Compressed Sparse Row format.
 */
struct
csr_matrix_t
{
  /** Number of rows in the matrix */
  int m;
  
  /** Number of columns in the matrix */
  int n;
  
  /** 
   * Number of stored (nonzero) entries.  If the matrix is stored in a 
   * symmetric (or skew, etc.) format, nnz only refers to the number of 
   * stored entries, not the actual number of nonzeros. 
   */
  int nnz;

  /** Array of stored (nonzero) entries of the matrix */
  void* values;

  /** Array of column indices of the stored (nonzero) entries of the matrix */
  int* colidx;

  /** Array of indices into the colidx and values arrays, for each column */
  int* rowptr;

  /**
   * Symmetry type of the matrix.
   */
  enum symmetry_type_t symmetry_type;

  /**
   * If the matrix has a kind of symmetry (or skew-symmetry): Where the actual 
   * elements of the matrix are stored: in the lower triangle or the upper 
   * triangle.
   */
  enum symmetric_storage_location_t symmetric_storage_location;

  /** 
   * Indicates the type of the entries in val:  REAL means "double", COMPLEX
   * means "double _Complex", PATTERN means the "val" array is NULL and contains
   * no entries.
   */ 
  enum value_type_t value_type;
};


/**
 * Fills in the data structure with shallow copies of the given arguments.
 *
 * @param A [OUT]
 * @param m [IN]
 * @param n [IN]
 * @param nnz [IN]
 * @param values [IN]
 * @param colidx [IN]
 * @param rowptr [IN]
 * @param symmetry_type [IN]
 * @param symmetric_storage_location [IN]
 * @param value_type [IN]
 */
void
pack_csr_matrix (struct csr_matrix_t* A, 
		 const int m, const int n, const int nnz, 
		 void* values, int* colidx, int* rowptr,
		 enum symmetry_type_t symmetry_type,
		 enum symmetric_storage_location_t symmetric_storage_location,
		 enum value_type_t value_type);

/**
 * Fills in the data structure with shallow copies of the given arguments.
 *
 * @param A [OUT]
 * @param m [IN]
 * @param n [IN]
 * @param nnz [IN]
 * @param values [IN]
 * @param colidx [IN]
 * @param rowptr [IN]
 * @param symmetry_type [IN]
 * @param symmetric_storage_location [IN]
 * @param value_type [IN]
 */
void
init_csr_matrix (struct csr_matrix_t* A, 
		 const int m, const int n, const int nnz, 
		 void* values, int* colidx, int* rowptr,
		 enum symmetry_type_t symmetry_type,
		 enum symmetric_storage_location_t symmetric_storage_location,
		 enum value_type_t value_type);

/**
 * Dynamically allocates a new csr_matrix_t object, and fills it in with
 * shallow copies of the given arguments.
 *
 * @param m [IN]
 * @param n [IN]
 * @param nnz [IN]
 * @param values [IN]
 * @param colidx [IN]
 * @param rowptr [IN]
 * @param symmetry_type [IN]
 * @param symmetric_storage_location [IN]
 * @param value_type [IN]
 *
 * @return Pointer to freshly allocated csr_matrix_t object
 */
struct csr_matrix_t*
create_csr_matrix (const int m, const int n, const int nnz, 
		   void* values, int* colidx, int* rowptr,
		   enum symmetry_type_t symmetry_type,
		   enum symmetric_storage_location_t symmetric_storage_location,
		   enum value_type_t value_type);

/**
 * Unpacks the CSC format sparse matrix into the given data structures.
 *
 * @param A       [IN]
 * @param m      [OUT]
 * @param n      [OUT]
 * @param nnz    [OUT]
 * @param values [OUT]
 * @param colidx [OUT]
 * @param rowptr [OUT]
 * @param symmetry_type [OUT]
 * @param symmetric_storage_location [OUT]
 * @param value_type [OUT]
 */
void
unpack_csr_matrix (const struct csr_matrix_t* A,
		   int* m, int* n, int* nnz,
		   void** values, int** colidx, int** rowptr,
		   enum symmetry_type_t* symmetry_type,
		   enum symmetric_storage_location_t* symmetric_storage_location,
		   enum value_type_t* value_type);

/**
 * Deep-copies src into dest.
 *
 * @param src [IN]
 * @param dest [OUT]
 */
void
copy_csr_matrix (struct csr_matrix_t* dest, const struct csr_matrix_t* src);


/**
 * Deallocates the data structures used by A.
 */
void
dealloc_csr_matrix (struct csr_matrix_t* A);


/**
 * Deallocates the data structures used by A, and the csr_matrix_t data
 * structure itself.
 */
void
destroy_csr_matrix (struct csr_matrix_t* A);


struct csc_matrix_t; /* forward declaration */

/**
 * Returns a copy of the given CSC-format sparse matrix in CSR format.
 *
 * @param A [IN]
 *
 * @return NULL if error, else a valid malloc'd CSR-format sparse matrix
 */
struct csr_matrix_t*
csc_to_csr (struct csc_matrix_t* A);

/**
 * Returns a copy of the given CSR-format sparse matrix in CSC format.
 *
 * @param A [IN]
 *
 * @return NULL if error, else a valid malloc'd CSC-format sparse matrix
 */
struct csc_matrix_t*
csr_to_csc (const struct csr_matrix_t* A);

/**
 * Writes the given matrix in Harwell-Boeing format to the given file.
 *
 * @param filename [IN] The file to write
 * @param A        [IN] The matrix (in Harwell-Boeing format)
 *
 * @return Nonzero if something went wrong, else zero.
 */
int
save_csr_matrix_in_harwell_boeing_format (const char* const filename, 
					  const struct csr_matrix_t* A);

/**
 * Prints the given CSR-format matrix to the given output stream in MatrixMarket format.
 *
 * @param out [OUT]  Valid (open) output stream
 * @param A [IN]     CSR-format sparse matrix
 *
 * @return Zero if no error, else nonzero.
 */
int
print_csr_matrix_in_matrix_market_format (FILE* out, const struct csr_matrix_t* A);


/**
 * Saves the given CSR-format sparse matrix A to the given filename, 
 * in MatrixMarket format.
 */
int 
save_csr_matrix_in_matrix_market_format (const char* const filename, 
					 const struct csr_matrix_t* A);

/**
 * Prints the given CSR-format matrix to the given output stream in Matlab format.
 *
 * @param out [OUT]  Valid (open) output stream
 * @param A [IN]     CSR-format sparse matrix
 *
 * @return Zero if no error, else nonzero.
 */
int
print_csr_matrix_in_matlab_format (FILE* out, struct csr_matrix_t* A);


/**
 * Saves the given CSR-format sparse matrix A to the given filename, 
 * in Matlab format.
 */
int 
save_csr_matrix_in_matlab_format (const char* const filename, 
				  struct csr_matrix_t* A);




/**
 * If the sparse matrix A is stored in a symmetric space-saving format,
 * expands out A so that all the nonzeros ((i,j) and (j,i)) are explicitly
 * stored.
 *
 * @param A [IN/OUT]  Sparse matrix in CSR format
 *
 * @return Nonzero if error, else zero.
 */
int
csr_matrix_expand_symmetric_storage (struct csr_matrix_t* A);

/**
 * Returns nonzero if the given CSR format sparse matrix has valid indices,
 * else returns zero.
 */
int
valid_csr_matrix_p (const struct csr_matrix_t* A);


/**
 * Returns C := B * A, a sparse matrix in CSR format. 
 * Returns NULL on error, e.g. if B and A don't have the same value type,
 * or if the data format isn't supported (e.g. symmetric storage may not
 * be implemented yet).
 */
struct csr_matrix_t*
csr_matrix_matmatmult (struct csr_matrix_t* B, struct csr_matrix_t* A);




#endif /* _csr_matrix_h */
