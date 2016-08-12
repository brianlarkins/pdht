/**************************************************
 * main.c
 *
 * does command-line version of HPCC spmv benchmark
 *
 * author: hormozd
 * last updated: Dec 18, 2005, 22:32 PST
 **************************************************/

#include "benchmark_params.h"
#include "hpcc_spmv_benchmark.h"
#include <smvm_util.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
  int benchcheck, t_max, mem_max;
  int rowdims[] = {1,2,3,4,6,8}, coldims[] = {1,2,3,4,6,8};
  double interval_fracs[] = {0.659,0.114,0.0584,0.0684,0.0285,0.0186,0.0144,0.0271,0.00774,0.00387};
  hpcc_spmv_params bench_params;
   
  /* obtain arguments. should be invoked like this from
     the command line:
  
       hpcc_benchmark [tmax] [memmax]
    
     where tmax is the maximum allowed time in minutes, and memmmax is
     the memory constraint in megabytes */
  if (argc < 3) {
    printf("\nIncorrect number of arguments. Usage is\n");
    printf("\n\thpcc_benchmark [tmax] [memmax] (runall)\n");  
    printf("\nwhere tmax is the time limit in minutes and\n");
    printf("memmax is the memory limit in megabytes.\n");
    printf("runall is an optional parameter that if nonzero\n");
    printf("says to run the benchmark with no time limit\n");
    printf("(default value is zero)\n\n");
    return 1;
  }
  
  /* process arguments and check them */
  t_max = atoi(argv[1]);
  mem_max = atoi(argv[2]);
  
  if (argc >= 4)
    bench_params.runall = atoi(argv[3]);
  else
    bench_params.runall = 0;
  
  if (t_max <= 0 || mem_max <= 0) {
    printf("Error: Both tmax and memmax must be positive!\n");
    return 2;
  }
    
  /* fill benchmark parameter struct */
  bench_params.row_blockdims = rowdims;   
  bench_params.col_blockdims = coldims;
  bench_params.n_row_blockdims = 6;
  bench_params.n_col_blockdims = 6;
  bench_params.nnz_per_row_min = 24;
  bench_params.nnz_per_row_max = 34;
  bench_params.mindim = 512;
  bench_params.outfilebase = "mflops";
  bench_params.tmax = t_max;
  bench_params.memmax = mem_max;
  bench_params.rmaxlen = 1;
  bench_params.cmaxlen = 1; 
  bench_params.interval_fracs = interval_fracs;
  
  /* call "driver routines" */
  benchcheck = run_benchmark(&bench_params, 1, 1);
   
  /* we're done! return value returned by driver routine */
  return benchcheck;
}
