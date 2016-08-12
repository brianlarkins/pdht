/**
 * @file test_sort_joint_arrays.c
 * @author Mark Hoemmen
 * @since 7 Jul 2005
 * @date 23 Feb 2006
 *
 * Driver program for testing functions defined in sort_joint_arrays.c.  
 * Returns EXIT_SUCCESS if all tests pass, else prints an error message 
 * and returns EXIT_FAILURE.
 */
#include "../sort_joint_arrays.h"
#include "../smvm_malloc.h"
#include "../smvm_util.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * Assumes that `array' is a void*, an array of elements, each of which has size
 * `size'.  `idx' is the index into the array.  This macro is needed because
 * `array' is a void pointer, so array[idx] is an invalid expression (since the
 * compiler doesn't know how big the elements of `array' are).  So we do the 
 * pointer math ourselves.  
 *
 * @note The usual `AREF' macro dereferences the pointer; this macro does NOT.
 *
 * @warn This won't work if sizeof(char) != 1.  The reason we cast to char* is 
 * to get a T* such that sizeof(T) == 1.
 */
#define VOIDAREF( array, idx )  ((void*) ((char*) (array) + size*(idx)))


void
test_sort_joint_arrays (void* A1, void* A2, const void* const sorted_A1, 
			const void* const sorted_A2, size_t nmemb, size_t size,
			int (*compar) (const void*, const void*, const void*, const void*))
{
  int i;

  sort_joint_arrays (A1, A2, nmemb, size, compar);

  for (i = 0; i < size; i++)
    {
      if (0 != compar (VOIDAREF(A1, i), VOIDAREF(A2, i), 
		       VOIDAREF(sorted_A1, i), VOIDAREF(sorted_A2, i)))
   	{
	  fprintf (stderr, "*** test_sort_joint_arrays: Element %d of %d "
		   "failed to compare correctly! ***\n", i, nmemb);
	  exit (EXIT_FAILURE);
	}
    }
}



int
main (int argc, char** argv)
{
  printf ("=== test_sort_joint_arrays ===\n");



  printf ("=== Passed test_sort_joint_arrays ===\n");
  return EXIT_SUCCESS;
}
