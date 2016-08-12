#ifndef hpcc_spmv_benchmark_h
#define hpcc_spmv_benchmark_h

/**********************************************************
 * hpcc_spmv_benchmark.h
 *
 * Front-end for the SpMV Benchmark run HPC Challenge style
 *
 * author: hormozd
 * last updated: Oct 19 2005, 21:11 PDT
 **********************************************************/

double median(double *a, int nelems);

int hpcc_spmv_benchmark(hpcc_spmv_params* bench_params, int write_to_files, int text_output);

int run_benchmark(hpcc_spmv_params* bench_params, int write_to_files, int text_output);

#endif
