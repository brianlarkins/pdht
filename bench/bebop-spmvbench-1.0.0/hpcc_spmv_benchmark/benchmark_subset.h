#ifndef benchmark_subset_h
#define benchmark_subset_h

/**************************************************
 * benchmark_subset.h
 * 
 * Determines the "appropriate subset" of the SpMV
 * benchmark to run by estimating the benchmark's
 * run time.
 *
 * author: hormozd
 * last updated: 2 Nov 2005, 20:42 PST
 **************************************************/

#include "benchmark_params.h"

void compute_threshold_dims(hpcc_spmv_params *bench_params);

int max_problem_dimension(hpcc_spmv_params *bench_params);

int benchmark_subset(hpcc_spmv_params *bench_params);

#endif
