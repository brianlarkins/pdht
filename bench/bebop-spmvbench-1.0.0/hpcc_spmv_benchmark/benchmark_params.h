#ifndef BENCHMARK_PARAMS
#define BENCHMARK_PARAMS

/********************************************
 * benchmark_params.h
 *
 * Contains parameters for the SpMV benchmark
 *
 * author: hormozd
 * last modified: Nov 2, 2005 21:48 PST
 ********************************************/

/* structure for input parameters:
 *   row_blockdims - array of row block dimensions
 *   col_blockdims - array of column block dimensions
 *   n_row_blockdims - number of row block dimensions
 *   n_col_blockdims - number of column block dimensions
 *   rmaxlen - maximum digit length of row blocksize
 *   cmaxlen - maximum digit length of column blocksize 
 *   nnz_per_row_min - minimum nnz/row to test
 *   nnz_per_row_max - maximum nnz/row to test
 *   mindim - minimum problem dimension
 *   tmax - time constraint (in minutes)
 *   memmax - memory constraint (in megabytes)
 *   maxdim - maximum problem dimension (set at runtime)
 *   numdims - number of dimensions to test (set at runtime)
 *   outfilebase - base name of output file for MFLOP rates
 *   est_runtime - estimated runtime of benchmark (in minutes) 
 *   actual_runtime - actual runtime of benchmark 
 *   runall - flag telling us whether to run every problem size
 *            in the test space instead of adhering to a time limit. 
 *   nnz_per_row_to_run - # of nnz/row values to test 
 *   threshold_dims - 2D array that contains "threshold dimensions"
 *                    (dimensions below which running the benchmark is
 *                     "cheap") for each blocksize 
 *   nnz_row_vals - array containing values of nnz/row that will
 *                  actually be tested 
 *   unblocked_mflop_max - max unblocked mflop rate
 *   blocked_mflop_max - max blocked mflop rate
 *   unblocked_mflop_med - median unblocked mflop rate
 *   blocked_mflop_med - median blocked mflop rate 
 *   interval_fracs - statistics for matrix generator. NEEDS to be of
 *                    length 10; no more, no less */
typedef struct {
    int *row_blockdims;
    int *col_blockdims;
    int n_row_blockdims;
    int n_col_blockdims;
    int rmaxlen;
    int cmaxlen;
    int nnz_per_row_min;
    int nnz_per_row_max;
    int mindim;
    int tmax;
    int memmax;
    int maxdim;
    int numdims;
    char *outfilebase;
    int est_runtime;
    int actual_runtime;
    int runall;
    int nnz_per_row_to_run;
    int **threshold_dims;
    int *nnz_row_vals;
    double unblocked_mflop_max;
    double blocked_mflop_max;
    double unblocked_mflop_med;
    double blocked_mflop_med;
    double *interval_fracs;
} hpcc_spmv_params;

#endif
