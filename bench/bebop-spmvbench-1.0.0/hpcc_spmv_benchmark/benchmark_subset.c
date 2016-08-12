/**************************************************
 * benchmark_subset.c
 * 
 * Determines the "appropriate subset" of the SpMV
 * benchmark to run by estimating the benchmark's
 * run time.
 *
 * author: hormozd
 * last updated: Nov 3, 2005 13:29 PST
 **************************************************/

#include "benchmark_params.h"
#include "benchmark_subset.h"
#include <timer.h>
#include <smvm_malloc.h>
#include <smvm_util.h>
#include <smvm_benchmark.h>
#include <smvm_timing_results.h>
#include <math.h>
#include <stdlib.h>

/*********************************************************************
 * compute_threshold_dims
 *
 * args:
 *
 * returns: pointer to 2D array of "timer thresholds", i.e. problem
 *          dimensions for each blocksize being tested at which
 *          generating the matrix becomes no longer a "free" operation
 *********************************************************************/
void compute_threshold_dims(hpcc_spmv_params *bench_params) {
  int dim, r, c, r_curr, c_curr;
  double time1, time2;
  struct SMVM_parameters spmv_params;
  struct SMVM_timing_results spmv_results;

  /* fill threshold_dims array by running SpMV timings */
  for (r = 0; r < bench_params->n_row_blockdims; r++) {
    for (c = 0; c < bench_params->n_col_blockdims; c++) {
      for (dim = bench_params->mindim; dim <= bench_params->maxdim; dim *= 2) {
        r_curr = bench_params->row_blockdims[r];
        c_curr = bench_params->col_blockdims[c];

        spmv_params.m = dim;
        spmv_params.n = dim;
        spmv_params.r = r_curr;
        spmv_params.c = c_curr;
        spmv_params.num_trials = 10;
        spmv_params.percent_fill = -bench_params->nnz_per_row_max;
        spmv_params.dataoutfile = stdout;
        spmv_params.interval_fracs = bench_params->interval_fracs;

        /* run benchmark 10 trials */
        time1 = get_seconds();
        smvm_benchmark_with_results(&spmv_params, &spmv_results);
        time2 = get_seconds();

        /* bug fix: need this here. if maxdim is less than what would
           have been the threshold dimension had we searched for it all
           the way, then the "computed" threshold dimension would be 0
           without this line */
        bench_params->threshold_dims[r][c] = dim;

        if (time2 - time1 >= .1) {
          bench_params->threshold_dims[r][c] /= 2;
          break;
        }
      }
    }
  }
}

/***********************************************************************
 * max_problem_dimension
 *
 * args: benchmark parameters
 *       threshold_dims - 2D array of "timer thresholds"
 *                        corresponding to different blocksizes
 *
 * returns: maximum problem dimension we can run
 ***********************************************************************/
int max_problem_dimension(hpcc_spmv_params *bench_params)
{
  int dim, logdim, r, c, allowed_time;
  int r_curr, c_curr, nnz_per_row_tested, nnz_curr;
  double ***accum_time_ests;
  double time1, time2, est_run_time = 0, time_est_prev, run_time_base;
  struct SMVM_parameters spmv_params;
  struct SMVM_timing_results spmv_results;

  allowed_time = 60*bench_params->tmax;
  nnz_per_row_tested = bench_params->nnz_per_row_max - 
                       bench_params->nnz_per_row_min + 1;

  /* allocate space to hold accumulated estimated running time for each
     problem dimension that could be tested */
  accum_time_ests = smvm_calloc(bench_params->numdims, sizeof(double**));

  for (logdim = 0; logdim < bench_params->numdims; logdim++) {
    accum_time_ests[logdim] = smvm_calloc(bench_params->n_row_blockdims, sizeof(double*));

    for (r = 0; r < bench_params->n_row_blockdims; r++) {
      accum_time_ests[logdim][r] = smvm_calloc(bench_params->n_col_blockdims, sizeof(double));
    }
  }

  /* fill in array that estimates benchmark runtime for each problem
     dimension. also use this to estimate runtime for running every
     trial in our test space */
  for (r = 0; r < bench_params->n_row_blockdims; r++) {
    for (c = 0; c < bench_params->n_col_blockdims; c++) {
      logdim = 0;

      for (dim = bench_params->mindim; dim < bench_params->threshold_dims[r][c]; dim *= 2) {
        logdim++;
      }

      r_curr = bench_params->row_blockdims[r];
      c_curr = bench_params->col_blockdims[c]; 

      spmv_params.m = dim;
      spmv_params.n = dim;
      spmv_params.r = r_curr;
      spmv_params.c = c_curr;
      spmv_params.num_trials = 10;
      spmv_params.percent_fill = -bench_params->nnz_per_row_max;
      spmv_params.dataoutfile = stdout;
      spmv_params.interval_fracs = bench_params->interval_fracs;

      /* run benchmark 10 trials */
      time1 = get_seconds();
      smvm_benchmark_with_results(&spmv_params, &spmv_results);
      time2 = get_seconds();

      accum_time_ests[logdim][r][c] = time2 - time1;
      est_run_time += nnz_per_row_tested*(time2 - time1);
      logdim++;

      /* now compute estimate of benchmark runtime for one nnz/row value
         for this particular blocksize by doubling this runtime until
         we reach maxdim */
      for (dim *= 2; dim <= bench_params->maxdim; dim *= 2) {
        accum_time_ests[logdim][r][c] = 2*accum_time_ests[logdim-1][r][c];
        est_run_time += nnz_per_row_tested*accum_time_ests[logdim][r][c];
        logdim++;
      }
    }
  }

  /* now compute largest problem dimension we can go to. first step is
     to cut down on the amount of nnz/row trials. if that doesn't work,
     then cut the dimension and try again */
  logdim--;
  dim = bench_params->maxdim;
  time_est_prev = est_run_time;

  /* bug fix: hormozd 2/6/2006 2:06 PST - need this initialization here
     in case bottom loop doesn't run (which happens when the memory
     constraint is small!) */
  nnz_curr = nnz_per_row_tested;

  while (dim >= bench_params->mindim && est_run_time >= 1.1*allowed_time) {
    run_time_base = time_est_prev;

    /* if we've gone through the loop before, cut off the largest
       problem dimension from the test space */
    if (dim < bench_params->maxdim) {
      for (r = 0; r < bench_params->n_row_blockdims; r++) {
        for (c = 0; c < bench_params->n_col_blockdims; c++) {
          run_time_base -= nnz_per_row_tested*accum_time_ests[logdim+1][r][c];
        }
      }
    }

    time_est_prev = run_time_base;

    /* are we within alotted time? if not, cut down on nnz/row */
    nnz_curr = nnz_per_row_tested;

    while (est_run_time > 1.1*allowed_time && nnz_curr > 1) {
      nnz_curr--;
      est_run_time = run_time_base*nnz_curr/nnz_per_row_tested;
    }

    logdim--;
    dim = dim/2;
  }

  /* collect final results */
  if (dim < bench_params->maxdim) {
    dim *= 2;
    logdim += 2;
  }

  /* major bug fix here: if we don't cut any from our maximum problem 
     dimension (i.e. memory constrains problem size rather than time,
     except when we cut nnz/row but not dimension), we need to increment
     logdim to compensate for the decrement before the above loop */

  else
    logdim++;

  bench_params->maxdim = dim;
  bench_params->numdims = logdim;
  bench_params->est_runtime = est_run_time/60 + 1;
  bench_params->nnz_per_row_to_run = nnz_curr;

  /* cleanup */
  for (logdim = 0; logdim < bench_params->numdims; logdim++) {
    for (r = 0; r < bench_params->n_row_blockdims; r++) {
      smvm_free(accum_time_ests[logdim][r]);
    }

    smvm_free(accum_time_ests[logdim]);
  }

  smvm_free(accum_time_ests);
  return dim;
}

/*******************************************************
 * round_to_int
 * 
 * rounds a floating-point number to the nearest integer
 * only works for positive arguments
 *******************************************************/
double round_to_int(double x) {
  return floor(x+0.5);
}

/*****************************************
 * pick_nnz_row_vals
 * 
 * selects which nnz/row trials to conduct
 *****************************************/
void pick_nnz_row_vals(int *vals, int nnz_row_run, int nnz_row_min, int nnz_row_max) {
  int idx;
  double spacing, nnz_to_test; 
  
  if (nnz_row_run == 1)
    vals[0] = (nnz_row_min + nnz_row_max)/2;
  else if (nnz_row_run >= 2) {
    vals[0] = nnz_row_min;
    vals[nnz_row_run-1] = nnz_row_max;
 
    if (nnz_row_run > 2) {
      spacing = (double)(nnz_row_max - nnz_row_min)/(nnz_row_run-1);
      nnz_to_test = nnz_row_min + spacing;
  
      for (idx = 1; idx <= nnz_row_run-2; idx ++) {
        vals[idx] = (int)round_to_int(nnz_to_test);
        nnz_to_test += spacing;
      }
    }
  }
}

/******************************************
 * benchmark_subset
 *
 * determines subset of benchmark to be run
 ******************************************/
int benchmark_subset(hpcc_spmv_params *bench_params) {
  long int maxdim_intermediate;
  int nnz_per_row_med, dim, r, dimcheck;

  /* determine the absolute largest problem dimension allowed
   * us by the memory restriction */

  /* BUG NOTE (hormozd, 23:10 PST, 11/5/2005): for some reason, this
     calculation overflows, even though I'm using a long. the end result
     for, to give an example, memmax = 3000 should fit in an int, even
     though the numerator would overflow in an int but not a long (but
     shouldn't using a long to handle this numerator solve it?) */
    
#if 0
  maxdim_intermediate = 1048576*((long)bench_params->memmax);
  maxdim_intermediate = maxdim_intermediate/
                           (12*bench_params->nnz_per_row_max+20);
#endif
  /* BUG FIX (mhoemmen, 10:20 PST, 11/6/2005): the usual way to fix
   * int overflow bugs in intermediate calculations is to perform the
   * calculations using double-precision floating point (with its much
   * larger exponent range) and then cast the results back to int. */ 
  {
    double mi = 1048576.0 * (double) (bench_params->memmax);
    mi = mi / (12.0 * (double) (bench_params->nnz_per_row_max) + 20.0);
    maxdim_intermediate = mi;
  }
                           
  bench_params->maxdim = maxdim_intermediate;
  
  /* now find largest problem dimension we can test given that we're  
     doing problem doubling */
  bench_params->numdims = 0;
    
  for (dim = bench_params->mindim; dim <= bench_params->maxdim; dim *= 2) {
    bench_params->numdims++;
  }                        
  
  bench_params->maxdim = dim/2;
   
  /* if the number of dimension sizes that we're going to test is zero, we
     have a problem... */
  if (bench_params->numdims == 0)
    return 1;

  /* allocate space for threshold_dims array */
  bench_params->threshold_dims = smvm_calloc(bench_params->n_row_blockdims,sizeof(int*));
  
  for (r = 0; r < bench_params->n_row_blockdims; r++) {
    bench_params->threshold_dims[r] = smvm_calloc(bench_params->n_col_blockdims,sizeof(int));
  }
  
  /* these steps are only needed if we're sticking to the time limit */
  if (bench_params->runall == 0) {
    /* compute "threshold dimensions", for which the estimated runtime of
       the benchmark is on the order of seconds */
    compute_threshold_dims (bench_params);
  
    /* now find maximum allowed problem dimension */
    dimcheck = max_problem_dimension (bench_params);
  
    if (dimcheck < 0) {
      for (r = 0; r < bench_params->n_row_blockdims; r++)
        smvm_free(bench_params->threshold_dims[r]);

      smvm_free(bench_params->threshold_dims);
      return 1;
    }
  }
  else {
    bench_params->nnz_per_row_to_run = bench_params->nnz_per_row_max -
                                       bench_params->nnz_per_row_min + 1;
  }

  /* allocate array that contains nnz/row values to run and fill it */
  nnz_per_row_med = (bench_params->nnz_per_row_min + bench_params->nnz_per_row_max)/2;
  bench_params->nnz_row_vals = smvm_malloc(bench_params->nnz_per_row_to_run*sizeof(int));
  pick_nnz_row_vals(bench_params->nnz_row_vals,bench_params->nnz_per_row_to_run,bench_params->nnz_per_row_min,bench_params->nnz_per_row_max);

  /* done */
  return 0;
}
