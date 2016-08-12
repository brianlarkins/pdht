/**
 * @file smvm_benchmark_results.c
 * @author mfh
 * @date 2004 Nov 29
 */

#include "smvm_benchmark_results.h"
#include "config.h"
#include <smvm_malloc.h>
#include <smvm_util.h>
#include <assert.h>

void
smvm_output_benchmark_results (struct SMVM_benchmark_results* p_results, 
			       FILE* out)
{
  WITH_DEBUG( fprintf(stderr, "=== smvm_output_benchmark results ===\n") );

  /* 
   * Note that fill differs only slightly for different block dimensions, when
   * n is large.
   */
  fprintf (out, "For %d x %d BSR sparse matrix:\n", 
	   p_results->n, p_results->n);
  fprintf (out, "Best median Mflop/s rate is %g Mflop/s, for (r,c) = (%d,%d) "
	   "(fill = %g)\n", p_results->best_mflops_rate, p_results->best_r, 
	   p_results->best_c, p_results->best_fill);
  fprintf (out, "Worst median Mflop/s rate is %g Mflop/s, for (r,c) = (%d,%d)"
	   " (fill = %g)\n", p_results->worst_mflops_rate, p_results->worst_r, 
	   p_results->worst_c, p_results->worst_fill);
  fprintf (out, "1x1 Mflop/s rate is %g Mflop/s (fill = %g)\n",
	   p_results->one_by_one_mflops_rate, p_results->one_by_one_fill);
  fprintf (out, "Aggregate FEM Mflop/s rate is %g Mflop/s\n",
	   p_results->aggregate_fem_mflops_rate);
}


void
smvm_free_benchmark_results (struct SMVM_benchmark_results* results)
{
  if (results != NULL)
    smvm_free (results);
}


struct SMVM_benchmark_results* 
smvm_copy_benchmark_results (const struct SMVM_benchmark_results* in)
{
  struct SMVM_benchmark_results *out = smvm_alloc_benchmark_results ();

  out->failure = in->failure;
  out->errcode = in->errcode;
  out->n       = in->n;

  out->best_r            = in->best_r;
  out->best_c            = in->best_r;
  out->best_block_size   = in->best_block_size;

  out->worst_r           = in->worst_r;
  out->worst_c           = in->worst_c;
  out->worst_block_size  = in->worst_block_size;

  out->best_mflops_rate  = in->best_mflops_rate;
  out->best_fill         = in->best_fill;
  out->worst_mflops_rate = in->worst_mflops_rate;
  out->worst_fill        = in->worst_fill;

  out->one_by_one_mflops_rate = in->one_by_one_mflops_rate;
  out->one_by_one_fill        = in->one_by_one_fill;

  return out;
}


struct SMVM_benchmark_results*
smvm_alloc_benchmark_results ()
{
  return smvm_malloc (sizeof (struct SMVM_benchmark_results));
}


void
init_smvm_benchmark_results_array (struct SMVM_benchmark_results_array* p)
{
  p->array = NULL;
  p->curlen = 0;
  p->maxlen = 0;
}


void
append_to_smvm_benchmark_results_array (struct SMVM_benchmark_results_array* p,
					struct SMVM_benchmark_results* results)
{
  assert (p->curlen >= 0);
  assert (p->maxlen >= 0);
  assert (p->curlen <= p->maxlen);

  if (p->curlen >= p->maxlen)
    {
      p->maxlen = 1 + 2 * p->maxlen;
      p->array = smvm_realloc (p->array, p->maxlen * 
			       sizeof (struct SMVM_benchmark_results*));
    }
  p->array[p->curlen] = results;
  p->curlen++;
}


void
destroy_smvm_benchmark_results_array_contents (struct SMVM_benchmark_results_array* p)
{
  int i;
  if (p != NULL)
    {
      for (i = 0; i < p->curlen; i++) smvm_free (p->array[i]);
      if (p->array != NULL) smvm_free (p->array);
    }
}


struct SMVM_benchmark_results* 
get_from_smvm_benchmark_results_array (struct SMVM_benchmark_results_array *p, 
				       const int index)
{
  assert (index >= 0 && index < p->curlen);
  return p->array[index];
}


int 
smvm_benchmark_results_array_length (struct SMVM_benchmark_results_array *p)
{
  return p->curlen;
}
