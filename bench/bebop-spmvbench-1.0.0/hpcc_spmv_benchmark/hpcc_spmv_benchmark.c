/**********************************************************
 * hpcc_spmv_benchmark.c
 *
 * Front-end for the SpMV Benchmark run HPC Challenge style
 *
 * author: hormozd
 * last updated: Nov 2 2005, 20:43 PST
 **********************************************************/

#include "benchmark_params.h"
#include "benchmark_subset.h"
#include <bcsr_matrix.h>
#include <create_rand.h>
#include <smvm_benchmark.h>
#include <smvm_timing_results.h>
#include <smvm_timing_run.h>
#include <smvm_verify_result.h>
#include <smvm_malloc.h>
#include <smvm_util.h>
#include <math.h>
#include <timer.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/****************************
 * median
 * 
 * returns median of an array
 ****************************/
double median(double *a, int nelems) {
  /* sort input array */
  qsort(a, nelems, sizeof(double), smvm_compare_double);

  /* what we do depends on whether nelems is even or odd */
  if (nelems % 2 == 1)
    return a[nelems/2];
  else
    return (a[nelems/2 - 1]+a[nelems/2])/2;
}

/******************************************************
 * blockdim_idx
 *
 * find index for mflop data for a particular blocksize
 * assumes 1x1 starts the data array, and the blocksize
 * referenced by r_idx = 0, c_idx = 0 is 1x1
 ******************************************************/
int blockdim_idx(int r_idx, int c_idx, int numdims, int nnz_vals, int n_rdims) {
  return r_idx*n_rdims*numdims*nnz_vals + c_idx*numdims*nnz_vals;
}

/********************************************************
 * mflop_idx
 *
 * find index for mflop data based on blocksize, dim, and
 * nnz/row value
 ********************************************************/
int mflop_idx(int r_idx, int c_idx, int nnz_idx, int dim_idx, int numdims, int nnz_vals, int n_rdims) {
  return blockdim_idx(r_idx,c_idx,numdims,nnz_vals,n_rdims)+
         nnz_idx*numdims + dim_idx;
}

/****************************************
 * find_timer_granularity
 *
 * routine to determine timer granularity
 ****************************************/
double find_timer_granularity(double *dummyval) {
  double t = 0.0;
  int i, k = 1;
  double x = 0.0;
  double delta = 1e-20;

  init_timer();

  while (t <= 0) {
    t = -get_seconds();

    /* We need to do something in the loop so the compiler
       won't optimize it away. */
    for (i = 0; i < k; i++)
      x = x + delta;

    t = t + get_seconds();
    k *= 2;
  }

  /* refine estimate, since it was obtained with coarse binary search */
  k = k/2;

  if (k > 2) {
    /* keep dividing interval in half and trying again, starting
       with [k/2, k] */
    int kbig = k; 
    int ksmall = k/2;
    int ktest;
         
    while (kbig - ksmall > 1) {
      ktest = (kbig + ksmall)/2;
          
      /* test midpoint (rounded down) of [ksmall,kbig] */
      t = -get_seconds();
    
      for (i = 0; i < ktest; i++)
        x = x + delta;
  
      t = t + get_seconds();
         
      if (t > 0) {
        /* if ktest works, then look in [ktest,kbig] */
        ksmall = ktest;
      }
      else {
        /* otherwise, look in [ksmall,ktest] */
        kbig = ktest;
      }
    }
  }

  *dummyval = x;
  return t;
} 

/****************************************************************
 * smallest_ntrials
 *
 * specialty function to run SpMV trials for the smallest problem
 * size. determines how many iterations to run to compensate for 
 * timer granularity. as a bonus, it also runs the benchmark for
 * the inputs given to it.
 *
 * returns: the number of trials to run
 ****************************************************************/
int smallest_ntrials(hpcc_spmv_params *bench_params, struct SMVM_parameters *spmv_params, struct SMVM_timing_results *timing_results, double timer_granularity) {
  int ntrials = 1, nnzb_per_row, num_block_rows, num_block_columns;
  int m_eff, n_eff, dimmax, b_verify = 1;
  double *src, *dest, *tmp, tol;
  struct bcsr_matrix_t *A = NULL;

  nnzb_per_row = (int)floor((-spmv_params->percent_fill/((double)spmv_params->c))+.5);
  num_block_rows = (int)floor(((double)spmv_params->m/(double)spmv_params->r)+.5);
  num_block_columns = (int)floor(((double)spmv_params->n/(double)spmv_params->c)+.5);
  m_eff = num_block_rows*spmv_params->r;
  n_eff = num_block_columns*spmv_params->c;
  tol = m_eff*1.11022e-16;

  if (m_eff > n_eff)
    dimmax = m_eff;
  else
    dimmax = n_eff;

  if (nnzb_per_row == 0)
    nnzb_per_row = 1;

  /* generate matrix */
  A = create_random_matrix_banded_by_statistics(num_block_rows,
         num_block_columns, spmv_params->r, spmv_params->c, nnzb_per_row,
         bench_params->interval_fracs);

  /* allocate source and destination vectors */
  src = (double *)smvm_malloc(dimmax*sizeof(double));
  dest = (double *)smvm_calloc(dimmax,sizeof(double));

  /* fill source vector with random numbers */
  smvm_init_vector_rand(dimmax,src);

  /* commence binary search to find how many trials. goal is 99% accuracy. */
  timing_results->t_median = 0;

  while (timing_results->t_median*ntrials < timer_granularity*100) {
    /* perform SpMV run. verify on the first multiplication (which
       only consists of one trial) */
    smvm_timing_run_with_results2(timing_results, m_eff, n_eff,
	                          spmv_params->r, spmv_params->c, A->rowptr,
                                  A->colind, A->values, A->nnzb,
                                  src, dest, ntrials, 0, 0, b_verify, tol);

    if (b_verify)
      b_verify = 0;

    ntrials *= 2;
  }

  ntrials = ntrials/2;

  /* refine initially determined value of ntrials if necessary */
  if (ntrials > 2) {
    /* keep dividing interval in half and trying again, starting
       with [ntrials/2, ntrials] */
    int ntrials_big = ntrials; 
    int ntrials_small = ntrials/2;
    int ntrials_test;

    while (ntrials_big - ntrials_small > 1) {
      ntrials_test = (ntrials_big + ntrials_small)/2;

      /* test midpoint (rounded down) of [ntrials_small,ntrials_big] */
      smvm_timing_run_with_results2(timing_results,m_eff,
          n_eff,spmv_params->r,spmv_params->c,A->rowptr,
          A->colind,A->values,A->nnzb,src,dest,ntrials_test,0,0,0,tol);

      if (timing_results->t_median/ntrials_test > timer_granularity*100) {
        /* if ktest works, then look in [ktest,kbig] */
        ntrials_small = ntrials_test;
      }
      else {
        /* otherwise, look in [ksmall,ktest] */
        ntrials_big = ntrials_test;
      }
    }

    ntrials = ntrials_test;
  }

  /* want at least 10 trials */
  if (ntrials < 10) {
    ntrials = 10;

    /* run 10 matrix multiplications to get an accurate number */
    smvm_timing_run_with_results2(timing_results,m_eff,n_eff,
       spmv_params->r,spmv_params->c,A->rowptr,
       A->colind,A->values,A->nnzb,src,dest,ntrials,0,0,0,tol);
  }

  /* free source and destination vectors and matrix */
  smvm_free(src);
  smvm_free(dest);
  destroy_bcsr_matrix(A);

  return ntrials;
}


/**************************************************************
 * hpcc_spmv_benchmark
 *
 * driver function for benchmark
 * arguments: benchmark parameters + 
 *            write_to_files flag (= 0 to suppress file output)
 *            text_output flag (= 0 to suppress text output)
 * returns: error code; "true" output is printed on the screen
 **************************************************************/
int hpcc_spmv_benchmark(hpcc_spmv_params* bench_params, int write_to_files, int text_output) {
  int dim, r, c, num_unblocked_runs, max_blocked_runs;
  int num_runs = 0, nnz_row, r_curr, c_curr;
  int len, mflop_max_idx = 0, nnz_idx, logdim, bigidx;
  int lookahead_idx, nnz_per_row_med, next_nz_row, seeded_ntrials;
  int orig_ntrials, ntrials;
  int **threshold_dims;
  int *r_at_max_mflop, *c_at_max_mflop, *nnz_row_vals;
  double time1, time2, unblocked_mflop_max = 0, blocked_mflop_max = 0;
  double mflop_rate, unblocked_mflop_median, blocked_mflop_median;
  double slope, timergranularity, dummy, totaltime;
  double *mflop_rates;
  char *outfilename;
  struct SMVM_timing_results timing_results;
  struct SMVM_parameters spmv_params;
  FILE *outfp;

  /* set up benchmark */
  nnz_row_vals = bench_params->nnz_row_vals;
  threshold_dims = bench_params->threshold_dims;
  nnz_per_row_med = (bench_params->nnz_per_row_min +
                        bench_params->nnz_per_row_max)/2;
  num_unblocked_runs = bench_params->numdims*(bench_params->nnz_per_row_max-
                          bench_params->nnz_per_row_min+1);
  max_blocked_runs = (bench_params->n_row_blockdims*bench_params->n_col_blockdims-1)*
                        num_unblocked_runs;
  mflop_rates = smvm_calloc(num_unblocked_runs+max_blocked_runs,sizeof(double));
  outfilename = smvm_calloc(strlen(bench_params->outfilebase)+7+
                   bench_params->rmaxlen+bench_params->cmaxlen,sizeof(char));
  r_at_max_mflop = smvm_calloc(bench_params->n_row_blockdims*
                      bench_params->n_col_blockdims,sizeof(int)+1);
  c_at_max_mflop = smvm_calloc(bench_params->n_row_blockdims*
                      bench_params->n_col_blockdims,sizeof(int)+1);

  /* determine timer granularity */
  timergranularity = find_timer_granularity(&dummy);

  /* print runtime estimate */
  if (text_output) {
    printf("estimated runtime = %d minutes\n", bench_params->est_runtime);
    printf("max problem dimension = %d\n", bench_params->maxdim);
  }

  time1 = get_seconds();

  for (r = 0; r < bench_params->n_row_blockdims; r++) {
    for (c = 0; c < bench_params->n_col_blockdims; c++) {
      r_curr = bench_params->row_blockdims[r];
      c_curr = bench_params->col_blockdims[c];
      nnz_idx = 0;
      seeded_ntrials = 0;

      for (nnz_row = bench_params->nnz_per_row_min; nnz_row <= bench_params->nnz_per_row_max; nnz_row++) {
        for (dim = bench_params->mindim; dim <= bench_params->maxdim; dim *= 2) {
          /* only run benchmark trial if
               (1) problem size is small enough, or
               (2) we're doing an appropriate value of nnz/row */
          if (dim < threshold_dims[r][c] || (nnz_idx < bench_params->nnz_per_row_to_run && nnz_row == nnz_row_vals[nnz_idx])) {
            /* set up benchmark run */
            spmv_params.m = dim;
            spmv_params.n = dim;
            spmv_params.r = r_curr;
            spmv_params.c = c_curr;
            spmv_params.percent_fill = -nnz_row;
            spmv_params.interval_fracs = bench_params->interval_fracs;

            /* hormozd 13:26 PST, 1.17.2006 -- need to set this, even
               though it's not used. though the benchmark function being
               called creates no file output, not setting this pointer
               could cause the benchmark to abort unnecessarily. */
            spmv_params.dataoutfile = stdout;

            /* have we found an appropriate ntrials for this blocksize or
               not? if not, we find it and use the results of that as our
               SpMV trial. if we have, proceed normally. */

            if (seeded_ntrials == 0) {
              orig_ntrials = smallest_ntrials(bench_params, &spmv_params, &timing_results, timergranularity);
              spmv_params.num_trials = orig_ntrials;
              ntrials = orig_ntrials;
              seeded_ntrials = 1;
            }
            else {
              smvm_benchmark_with_results(&spmv_params, &timing_results);
            }

            /* how long does one less trial take? */
            totaltime = (ntrials-1)*timing_results.t_median;
            mflop_rate = timing_results.mflops/timing_results.t_median;

            /* check time per SpMV and adjust down if needed */ 
            while (ntrials > 10 && totaltime > 100*timergranularity) {
              totaltime -= timing_results.t_median;
              ntrials--;
            } 

            spmv_params.num_trials = ntrials;

            /* collect data */
            mflop_rates[num_runs] = mflop_rate;

            if (r == 0 && c == 0 && mflop_rate > unblocked_mflop_max)
              unblocked_mflop_max = mflop_rate;
            else if (r != 0 || c != 0) {
              if (mflop_rate > blocked_mflop_max) {
                blocked_mflop_max = mflop_rate;
                r_at_max_mflop[0] = r_curr;
	        c_at_max_mflop[0] = c_curr;
	        r_at_max_mflop[1] = 0;
	        c_at_max_mflop[1] = 0;
	        mflop_max_idx = 1;
	      }
	      else if (mflop_rate == blocked_mflop_max) {
		r_at_max_mflop[mflop_max_idx] = r_curr;
		c_at_max_mflop[mflop_max_idx] = c_curr;
		mflop_max_idx++;
	      }
            }
          }

          num_runs++;
        }

        /* reset ntrials to original value */
        ntrials = orig_ntrials;
        spmv_params.num_trials = orig_ntrials;

        if (nnz_idx < bench_params->nnz_per_row_to_run && nnz_row == nnz_row_vals[nnz_idx])
          nnz_idx++;
      }

      /* reset seeded_ntrials */
      seeded_ntrials = 0;
    }
  }

  time2 = get_seconds();
  bench_params->actual_runtime = time2 - time1;

  /* fill in tests that weren't done using interpolation and write
     resulting data to output file */
  bigidx = 0;

  for (r = 0; r < bench_params->n_row_blockdims; r++) {
    for (c = 0; c < bench_params->n_col_blockdims; c++) {
      r_curr = bench_params->row_blockdims[r];
      c_curr = bench_params->col_blockdims[c];

      if (write_to_files) {
        len = sprintf(outfilename,"%s_%dx%d.csv",
                      bench_params->outfilebase, r_curr, c_curr);
        outfilename[len] = '\0';

        /* open output file for mflop rate info */
        outfp = fopen(outfilename,"w");

        if (outfp == NULL) {
          fprintf(stderr, "ERROR: unable to open output file for %d x %dcase!", r_curr, c_curr);
          /* free what needs to be freed here */
          return 1;
        }

        /* write top row of output file */
        for (dim=bench_params->mindim;dim<=bench_params->maxdim;dim*=2) {
          fprintf(outfp,",%d",dim);
        }
        fprintf(outfp,"\n");
      }

      for (nnz_row = bench_params->nnz_per_row_min; nnz_row <= bench_params->nnz_per_row_max; nnz_row++) {
        /* write current nnz/row value to output file, if requested */
        if (write_to_files)
          fprintf(outfp,"%d",nnz_row);

        for (logdim = 0; logdim < bench_params->numdims; logdim++) {
          /* is interpolation needed to fill this data point? */
          if (mflop_rates[bigidx] == 0) {

            if (bench_params->nnz_per_row_to_run == 1) {
              lookahead_idx = bigidx + bench_params->numdims*(nnz_per_row_med - nnz_row);
              mflop_rates[bigidx] = mflop_rates[lookahead_idx];
            }
            else {
              lookahead_idx = bigidx+bench_params->numdims;
              next_nz_row = nnz_row+1;

              while (mflop_rates[lookahead_idx] == 0) {
                lookahead_idx += bench_params->numdims;
                next_nz_row++;
              }

              mflop_rates[bigidx] = mflop_rates[lookahead_idx];

              /* compute slope for linear interpolation. note that bigidx-numdims
                 will always be filled. */
              slope = (mflop_rates[lookahead_idx] - mflop_rates[bigidx-bench_params->numdims])/(next_nz_row-nnz_row+1);

              /* now fill in value */
              mflop_rates[bigidx] = mflop_rates[bigidx-bench_params->numdims]+slope;
            }
          }

          /* write data to output file, if requested */
          if (write_to_files)
            fprintf(outfp,",%f",mflop_rates[bigidx]);

          bigidx++;
        }

        if (write_to_files)
          fprintf(outfp,"\n");
      }

      if (write_to_files)
        fclose(outfp);
    }
  }

  /* later on: organize data to make "comprehensive" graph */

  /* get median mflop numbers for both unblocked and blocked cases
     to complete data analysis (max numbers were already collected
     during the benchmark runs) */
  unblocked_mflop_median = median(mflop_rates, num_unblocked_runs);
  blocked_mflop_median = median(mflop_rates + num_unblocked_runs, max_blocked_runs);

  /* output the magic four numbers */
  if (text_output) {
    printf("\t|Max\t\t|Median\n");
    printf("-------\t|---\t\t|------\n");
    printf("1x1\t|%.4g\t\t|%.4g\n", unblocked_mflop_max, unblocked_mflop_median);
    printf("Blocked\t|%.4g\t\t|%.4g\n", blocked_mflop_max, blocked_mflop_median);
    printf("\nblocked max occurs for blocksizes\n");
  }

  mflop_max_idx = 0;
  while (r_at_max_mflop[mflop_max_idx] != 0
         && c_at_max_mflop[mflop_max_idx] != 0)
  {
    if (text_output)
      printf("\t%dx%d\n",r_at_max_mflop[mflop_max_idx],c_at_max_mflop[mflop_max_idx]);

    mflop_max_idx++;
  }

  if (text_output)
    printf("actual runtime: %d minutes\n", (int)((time2-time1)/60));

  /* cleanup */
  bench_params->unblocked_mflop_max = unblocked_mflop_max;
  bench_params->blocked_mflop_max = blocked_mflop_max;
  bench_params->unblocked_mflop_med = unblocked_mflop_median;
  bench_params->blocked_mflop_med = blocked_mflop_median;

  for (r = 0; r < bench_params->n_row_blockdims; r++)
    smvm_free(threshold_dims[r]);

  smvm_free(threshold_dims);
  smvm_free(mflop_rates);
  smvm_free(r_at_max_mflop);
  smvm_free(c_at_max_mflop);
  smvm_free(outfilename);
  smvm_free(nnz_row_vals);
  return 0;
}

/****************************************
 * run_benchmark
 *
 * wrapper function to run SpMV benchmark
 ****************************************/
int run_benchmark(hpcc_spmv_params *bench_params, int write_to_files, int text_output) {
  int subsetcheck, returnval;

  subsetcheck = benchmark_subset(bench_params);

  if (subsetcheck == 0) {
    returnval = hpcc_spmv_benchmark(bench_params, write_to_files, text_output);
    return returnval;
  }
  else
    return subsetcheck;
}
