#ifndef _smvm_benchmark_results_h
#define _smvm_benchmark_results_h
/**
 * @file smvm_benchmark_results.h
 * @author mfh 
 * @date 2004 Nov 29
 *
 ******************************************************************************/
#include <stdio.h> /* FILE */


/**
 * A summary of results from the entire benchmark.  This is most useful for 
 * integration with the HPCC benchmark (which requires that results be reported
 * in the HPCC_Params struct, as well as output to a file).  `best_' refers
 * to a result or parameter of the (r x c) register block combination with the
 * fastest average results; `worst_' to the slowest average results for all 
 * (r x c) register block combinations tried.  
 */
struct SMVM_benchmark_results
{
  /** Set to nonzero if the BSR SpMV benchmark suffered a fatal error. */
  int failure;

  /** 
   * `Return value' of the BSR SpMV benchmark; nonzero if an error 
   * occurred. 
   */
  int errcode;

  /** Problem size (n x n matrix) */
  int n;

  /** 
   * "r" (number of rows in a register block) corresponding to the best
   * (mflop/s) performance out of all the runs with the various register 
   * block dimensions.
   */
  int best_r;

  /** 
   * "c" (number of columns in a register block) corresponding to the best
   * (mflop/s) performance out of all the runs with the various register block 
   * dimensions.
   */
  int best_c;

  int best_block_size;

  /** 
   * "r" (number of rows in a register block) corresponding to the worst
   * (mflop/s) performance out of all the runs with the various register 
   * block dimensions.
   */
  int worst_r;

  /** 
   * "c" (number of columns in a register block) corresponding to the worst
   * (mflop/s) performance out of all the runs with the various register 
   * block dimensions.
   */
  int worst_c;

  int worst_block_size;

  /** 
   * Best (megaflops per second) performance of all the runs with the various 
   * register block dimensions.
   */
  double best_mflops_rate;

  /** 
   * Fill of the run with the best (mflop/s) performance of all the runs 
   * with the various register block dimensions.
   */
  double best_fill;

  /** 
   * Worst (megaflops per second) performance of all the runs with the various 
   * register block dimensions.
   */
  double worst_mflops_rate;

  /** 
   * Fill of the run with the worst (mflop/s) performance of all the runs 
   * with the various register block dimensions.
   */
  double worst_fill;

  /**
   * Fill of the 1x1 run.
   */
  double one_by_one_fill;

  /**
   * Megaflops per second rate for the 1x1 run.
   */
  double one_by_one_mflops_rate;

  /**
   * Aggregate Megaflops per second rate for finite element matrices (matrices
   * 2-9 in Rich Vuduc's thesis).
   */
  double aggregate_fem_mflops_rate;
};


/**
 * Outputs the given set of benchmark results to the given output stream.
 * This function should only be called by proc 0 (in MPI terms).
 */
void
smvm_output_benchmark_results (struct SMVM_benchmark_results* p_results, 
			       FILE* out);

/**
 * Returns a deep copy (allocated using malloc) of the object pointed to by
 * "in".
 */
struct SMVM_benchmark_results* 
smvm_copy_benchmark_results (const struct SMVM_benchmark_results* in);


/**
 * Deallocates the given object, assuming that it was allocated dynamically.
 */
void
smvm_free_benchmark_results (struct SMVM_benchmark_results* results);


/**
 * Allocates memory for a new SMVM_benchmark_results object, and returns
 * a pointer to the object.
 */
struct SMVM_benchmark_results*
smvm_alloc_benchmark_results ();


/**
 * Represents an array of SMVM_benchmark_results objects.  The
 * individual objects are accessed as pointers through the "array" array
 * of pointers.  "array" is an array of length maxlen, but it holds only
 * curlen pointers.
 */
struct SMVM_benchmark_results_array
{
  struct SMVM_benchmark_results **array;
  int curlen;
  int maxlen;
};

/**
 * Initializes the given SMVM benchmark results array object to hold zero
 * objects.
 */
void
init_smvm_benchmark_results_array (struct SMVM_benchmark_results_array* p);

/**
 * Appends the given results struct to the given array of results.
 * Uses shallow copy of the results pointer.
 */
void
append_to_smvm_benchmark_results_array (struct SMVM_benchmark_results_array* p,
					struct SMVM_benchmark_results* results);

/**
 * Returns a pointer to the SMVM_benchmark_results object stored in the
 * given benchmark results array object.
 */
struct SMVM_benchmark_results* 
get_from_smvm_benchmark_results_array (struct SMVM_benchmark_results_array *p, 
				       const int index);

/**
 * Destroys all the internal storage used by the given
 * SMVM_benchmark_results_array object.  The pointer to the object
 * itself is not freed.
 */
void
destroy_smvm_benchmark_results_array_contents (struct SMVM_benchmark_results_array* p);

/**
 * Returns the current number of elements in the given array of
 * benchmark results.
 */
int 
smvm_benchmark_results_array_length (struct SMVM_benchmark_results_array *p);


#endif /* NOT _smvm_benchmark_results_h */
