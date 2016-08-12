/**
 * @file test_random_number.c
 * @author Mark Hoemmen
 * @since 7 Jul 2005
 * @date 23 Feb 2006
 *
 * Driver program for testing functions defined in random_number.c.  Returns 
 * EXIT_SUCCESS if all tests pass, else prints an error message and returns 
 * EXIT_FAILURE.
 *
 * Copyright (c) 2006, Regents of the University of California 
 * All rights reserved.
 * Redistribution and use in source and binary forms, with or
 * without modification, are permitted provided that the
 * following conditions are met:
 * 
 * * Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 *
 * * Redistributions in binary form must reproduce the above copyright 
 *   notice, this list of conditions and the following disclaimer in 
 *   the documentation and/or other materials provided with the 
 *   distribution.
 *
 * * Neither the name of the University of California, Berkeley, nor
 *   the names of its contributors may be used to endorse or promote
 *   products derived from this software without specific prior
 *   written permission.  
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#include "../random_number.h"
#include "../smvm_malloc.h"
#include "../smvm_util.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h> /* clock() */
#include <string.h>


void
test_random_number (int low, int high)
{
  const int total = high - low + 1;
  struct random_integer_from_range_without_replacement_generator_t* gen = NULL;
  int* returned = NULL;
  int errcode = 0;
  int p = 0;
  int i;

  if (low > high)
    {
      gen = create_random_integer_from_range_without_replacement_generator (low, high);
      if (gen->num_remaining != 0)
	{
	  fprintf (stderr, "*** test_random_number: create_random_integer_from_ra"
		   "nge_without_replacement_generator (%d, %d) should return a gen"
		   "erator with num_remaining == 0, but num_remaining = %d instea"
		   "d ***\n", low, high, gen->num_remaining);
	  exit (EXIT_FAILURE);
	}

      errcode = return_random_integer_from_range_without_replacement (&p, gen);
      if (errcode == 0)
	{
	  fprintf (stderr, "*** test_random_number: return_random_integer_from_ra"
		   "nge_without_replacement should always fail when high < low, b"
		   "ut did not fail! ***\n");
	  exit (EXIT_FAILURE);
	}

      destroy_random_integer_from_range_without_replacement_generator (gen);
      return;
    }

  returned = smvm_calloc (sizeof (int), total);
  gen = create_random_integer_from_range_without_replacement_generator (low, high);
  assert (gen != NULL);
  assert (gen->low == low);
  assert (gen->high == high);

  /* Try to get high-low+1 distinct integers within the range [low,high]. */

  for (i = 0; i < total; i++)
    {
      int j;

      errcode = return_random_integer_from_range_without_replacement (&p, gen);
      if (errcode != 0)
	{
	  fprintf (stderr, "*** test_random_number: at iteration %d of %d, re"
		   "turn_random_integer_from_range_without_replacement should"
		   " have returned zero for success, but returned nonzero ins"
		   "tead ***\n", i+1, total);
	  exit (EXIT_FAILURE);
	}
      if (p < low || p > high)
	{
	  fprintf (stderr, "*** test_random_number: at iteration %d of %d, re"
		   "turn_random_integer_from_range_without_replacement return"
		   "ed a number %d outside the valid range [%d,%d] ***\n", 
		   i+1, total, p, gen->low, gen->high);
	  exit (EXIT_FAILURE);
	}
      returned[i] = p;

      /* Make sure that p is unique (has not been chosen already). */
      for (j = 0; i < i-1; i++)
	{
	  if (p == returned[j])
	    {
	      fprintf (stderr, "*** test_random_number: at iteration %d of %d"
		       ", return_random_integer_from_range_without_replacemen"
		       "t return ed a number %d that was returned previously "
		       "at iteration %d ***\n", i+1, total, p, j+1);
	      exit (EXIT_FAILURE);
	    }
	}

      /* Make sure that all the remaining numbers are unique. */
      for (j = 1; j < gen->num_remaining; j++)
	{
	  int k;
	  for (k = 0; k < j; k++)
	    {
	      if (gen->remaining[j] == gen->remaining[k])
		{
		  fprintf (stderr, "*** test_random_number: at iteration %d o"
			   "f %d, gen->remaining has a duplicate entry! ***\n",
			   i+1, total);
		  exit (EXIT_FAILURE);
		}
	    }
	}
    }

  if (gen->num_remaining != 0)
    {
      fprintf (stderr, "*** test_random_number: after removing all %d numbers, "
	       "%d are still remaining, though 0 should remain! ***\n", 
	       total, gen->num_remaining);
      exit (EXIT_FAILURE);
    }

  /* Make sure that returning a random integer fails.  Try this twice. */

  errcode = return_random_integer_from_range_without_replacement (&p, gen);
  if (errcode == 0)
    {
      fprintf (stderr, "*** test_random_number: return_random_integer_from_ra"
	       "nge_without_replacement should fail after returning all %d nu"
	       "mbers, but did not fail! ***\n", 6);
      exit (EXIT_FAILURE);
    }

  errcode = return_random_integer_from_range_without_replacement (&p, gen);
  if (errcode == 0)
    {
      fprintf (stderr, "*** test_random_number: return_random_integer_from_ra"
	       "nge_without_replacement should fail after returning all %d nu"
	       "mbers, but did not fail! ***\n", 6);
      exit (EXIT_FAILURE);
    }

  destroy_random_integer_from_range_without_replacement_generator (gen);
  smvm_free (returned);
}


void
test_random_double (double low, double high)
{
  double d = smvm_random_double (low, high);

  printf ("Random double in [%g,%g): %g\n", low, high, d);
  assert (d >= low);
  assert (d < high);
}


int 
main (int argc, char** argv)
{
  int low, high;
  smvm_set_debug_level_from_environment ();
  smvm_init_randomizer ((int) clock());

  printf ("=== test_random_number ===\n");

  if (argc > 1)
    {
      if (0 == strcmp (argv[1], "without-replacement"))
	{
	  low = 1;
	  high = 6;
	  test_random_number (low, high);

	  /* Now try a larger interval. */

	  low = -50;
	  high = 200;
	  test_random_number (low, high);

	  /* Now make a generator where high < low -- this one should never return any values. */
	  low = 0;
	  high = -1;
	  test_random_number (low, high);
	}
      else if (0 == strcmp (argv[1], "random_double"))
	{
	  int i;

	  for (i = 0; i < 20; i++)
	    {
	      test_random_double (0.0, 1.0);
	      test_random_double (0.0, 1.0);
	      test_random_double (0.0, 1.0);
	    }

	  test_random_double (-10.0, 10.0);
	  test_random_double (3.14159, 7.654321);
	}
      else if (0 == strcmp (argv[1], "random_integer"))
	{
	  int i, die;
	  for (i = 0; i < 100; i++)
	    {
	      die = smvm_random_integer (1, 6);
	      assert (die >= 1);
	      assert (die <= 6);
	      fprintf (stderr, "%d\n", die);
	    }
	}
    }
  else
    {
      low = 1;
      high = 6;
      test_random_number (low, high);

      /* Now try a larger interval. */

      low = -50;
      high = 200;
      test_random_number (low, high);

      /* Now make a generator where high < low -- this one should never return any values. */
      low = 0;
      high = -1;
      test_random_number (low, high);

      test_random_double (0.0, 1.0);
      test_random_double (0.0, 1.0);
      test_random_double (0.0, 1.0);

      test_random_double (-10.0, 10.0);
      test_random_double (3.14159, 7.654321);
    }

  printf ("=== Passed test_random_number ===\n");
  return 0;
}

