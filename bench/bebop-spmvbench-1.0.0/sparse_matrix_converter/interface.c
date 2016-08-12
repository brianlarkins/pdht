/**
 * @file interface.c
 * @author Mark Hoemmen
 * @since 09 Jun 2006
 * @date 09 Jun 2006
 */

#include "config.h"
#include "sparse_matrix.h"
#include "sparse_matrix_ops.h"

#include <get_options.h>
#include <smvm_malloc.h>
#include <smvm_util.h>
#include <timer.h>

#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>


struct sparse_matrix_t*
sp_load (const char* path, const char* fmt)
{
  struct sparse_matrix_t* A = NULL;
  double seconds = get_seconds ();
  A = load_sparse_matrix (sparse_matrix_file_format_string_to_enum (fmt), path);
  if (A == NULL)
    {
      fprintf (stderr, "*** sp_load: Failed to load sparse matrix! ***\n");
      return NULL;
    }
  seconds = get_seconds () - seconds;
  printf ("Loading the sparse matrix took %g seconds.\n", seconds); 
  return A;
}

int
sp_save (struct sparse_matrix_t* A, const char* path, 
      const char* fmt)
{
  double seconds = get_seconds ();
  save_sparse_matrix (path, A, sparse_matrix_file_format_string_to_enum (fmt));
  seconds = get_seconds () - seconds;
  printf ("Saving the sparse matrix took %g seconds.\n", seconds);
  return 0;
}

void
sp_format (struct sparse_matrix_t* A)
{
  if (A == NULL)
    printf ("Sparse matrix is NULL!\n");
  else
    printf ("%s\n", sparse_matrix_format_string (A));
}

int
sp_convert (struct sparse_matrix_t* A, const char* type)
{
  double seconds = get_seconds ();
  int errcode = sparse_matrix_convert (A, 
	sparse_matrix_storage_format_string_to_enum (type));
  seconds = get_seconds () - seconds;
  if (errcode != 0)
    {
      printf ("*** Failed to convert sparse matrix! ***\n");
      return errcode;
    }
  else
    {
      printf ("Converting the sparse matrix took %g seconds.\n", seconds);
      return 0;
    }
}


struct sparse_matrix_t* 
sp_mult (struct sparse_matrix_t* B, struct sparse_matrix_t* A)
{
  return sparse_matrix_matmatmult (B, A);
}

void
sp_destroy (struct sparse_matrix_t* A)
{
  destroy_sparse_matrix (A);
}
