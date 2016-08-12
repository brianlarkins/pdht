/*****************************************************************************
 * @file main.c
 * @author mfh
 * @date 2004 Nov 21
 * @brief Driver for invoker of sparse matrix-vector multiply benchmarks.
 *
 * Driver for invoker of sparse matrix-vector multiply benchmarks.  It uses 
 * the interface exposed by smvm_invoker.h.  That interface is an intermediate 
 * layer, between the raw benchmark and the interface exposed (in 
 * smvm_hpcc_interface.h) to the HPC Challenge benchmark suite.  This file 
 * remains in the project for historical reasons; users seeking a driver that 
 * produces a standalone executable for the SpMV benchmark should refer to 
 * smvm_hpcc_interface_tester.c.
 ****************************************************************************/
#include "config.h"

#include <errno.h>    /* errno used for checking string -> int conversion */
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>   /* getopt (POSIX command-line args reader) */

#include "smvm_invoker.h"
#include <smvm_util.h>


/****************************************************************************/
static void
usage ()
{
  fprintf (stderr, 
	   "Driver for sparse matrix-vector multiplication benchmark\n\n"
	   "Options:  -o <output_file> -m <maxmem> -s -d\n"
	   "-o <output_file>  File where results should be stored\n"
	   "-m <maxmem>       Max KB main memory that benchmark can use\n"
	   "-n <num_trials>   Number of trials to run each (r,c) pair\n"
	   "-s                Run sparse benchmark\n"
	   "-d                Run dense-matrix-in-sparse-format benchmark\n"
	   "The -s and -d flags may be combined.\n"
	   "\n\n");
}


  
/**
 * Shows given error and the usage, and exits with given error code.
 * Caller is responsible for any endlines in the error message.  
 *
 * @see usage()
 *
 * @param message  Message to print before displaying usage.
 * @param errcode  Argument of exit().
 *
 * @return  Does not return (calls exit(int)).
 */
static void
die_with_usage (char* msg, int errcode)
{
  fprintf (stderr, msg);
  usage ();
  exit (errcode);
}


/**
 * Command-line arguments for test wrapper for SMVM benchmark invoker.  These
 * are used for testing, and not meant to be included when the benchmark is
 * integrated into HPCC (MFH 2004 Jan 31).  In particular, a data output
 * filename should not be passed directly to the benchmark, because HPCC
 * writes data to stdout.
 *
 * @note Memory sizes are in KB, because if they were in bytes, we might have 
 *       to worry about exceeding bounds of int.
 */
struct
SMVM_invoker_cmdline_params
{
  /** Name of the file where data output should go. */
  char* dataoutfile;
 
  /** Max amount of mem in KB the matrix and vector data structures can use. */
  int max_mem_KB;
  
  /** Number of trials to run, for each register block dimension pair. */
  int num_trials;

  /** Set != 0 to run dense-matrix-in-sparse-format benchmark, else set = 0.*/
  int b_run_dense;

  /** Set != 0 to run usual sparse benchmark, else set = 0.*/
  int b_run_sparse;
};


/** 
 * Initializes the command-line arguments struct.  Assumes that the struct has
 * been allocated already.
 */
static void
smvm_init_invoker_cmdline_params (struct SMVM_invoker_cmdline_params* p_params)
{
  p_params->b_run_dense           = 0;
  p_params->b_run_sparse          = 0;
  p_params->dataoutfile           = NULL;
  p_params->max_mem_KB            = 0;
  p_params->num_trials            = 0;
}


/**
 * Get the command-line arguments (@see usage() for details).  We use
 * errno and the C standard functions that modify errno, to check for
 * correct conversion of command-line arguments.  Aborts if the
 * arguments cannot be retrieved correctly.
 *
 * @param p_params  Holds the command-line arguments' values.
 * @param argc      Pass in the corresponding parameter of main().
 * @param argv      Pass in the corresponding parameter of main().
 *
 * @return Nothing.
 */
static void
smvm_get_cmdline_args (struct SMVM_invoker_cmdline_params* p_params, 
		       int argc, char* argv[])
{
  int c;
#ifdef DEBUG
  int dbg = 1;
#else /* DEBUG */
  int dbg = 0;
#endif /* DEBUG */

  /* For required parameters. */
  int b_seen_m = 0;
  int b_seen_n = 0;

  /* If abortworthy error encountered. */
  int b_abortworthy_error = 0;

  /* opt* (static variables) are part of getopt(), which is a POSIX standard
   * function for reading single-letter command-line arguments with parameters.
   */
  extern char *optarg;
  extern int optind, optopt, opterr;

  smvm_init_invoker_cmdline_params (p_params);
  
  while ((c = getopt(argc, argv, ":sdo:1:2:m:n:")) != -1) 
    {
      switch(c) 
	{
	case 's':
	  p_params->b_run_sparse = 1;
	  if (dbg) fprintf (stderr, "-s:  Sparse benchmark\n");
	  break;
	case 'd':
	  p_params->b_run_dense = 1;
	  if (dbg) 
	    fprintf (stderr, "-d:  Dense-in-sparse-format benchmark\n");
	  break;
	case 'o':
	  /* If `-o' option unspecified, output defaults to stdout.  Since `o' 
	   * is followed by `:' in getopt string, optarg must be non-NULL. */
	  p_params->dataoutfile = optarg;
	  if (dbg) fprintf (stderr, "-o:  Output file is %s\n", 
			    p_params->dataoutfile);
	  break;
	case 'm':
	  p_params->max_mem_KB = (int) strtol (optarg, NULL, 10);
	  if (errno != 0)
	    die_with_usage ("*** ERROR: -m parameter must be an integer! ***\n", -1);
	  if (dbg) 
	    fprintf (stderr, "-m:  Max allocable mem is %d KB\n", 
		     p_params->max_mem_KB);
	  b_seen_m = 1;
	  break;
	case 'n':
	  p_params->num_trials = (int) strtol (optarg, NULL, 10);
	  if (errno != 0)
	    die_with_usage ("*** ERROR: -n parameter must be an integer! ***\n", -1);
	  b_seen_n = 1;
	  break;
	case ':':
	  fprintf (stderr, "*** ERROR:  Option -%c is missing its parameter! ***\n", optopt);
	  die_with_usage ("", -1);
	  break;
	case '?':
	  fprintf(stderr, "*** ERROR:  Unknown option -%c ***\n", optopt);
	  die_with_usage ("", -1);
	  break;
	}
    }

  /* Check for required parameters. */
  if (! b_seen_m)
    {
      b_abortworthy_error = 1;
      fprintf (stderr, "*** `-m\' parameter not specified! ***\n");
    }
  if (! b_seen_n)
    {
      b_abortworthy_error = 1;
      fprintf (stderr, "*** `-n\' parameter not specified! ***\n");
    }

  if (! (p_params->b_run_sparse || p_params->b_run_dense))
    {
      b_abortworthy_error = 1;
      fprintf (stderr, "*** Must specify at least one of -s or -d "
	       "parameters ***\n");
    }

  if (b_abortworthy_error)
    {
      die_with_usage ("\n", -1);
    }
}


/****************************************************************************/
int 
main (int argc, char* argv[])
{
#ifdef DEBUG
  int dbg = 1;
#else 
  int dbg = 0;
#endif /* DEBUG */

#ifdef WARN
  int warn = 1;
#else
  int warn = 0;
#endif /* WARN */

  FILE* output = NULL;
  struct SMVM_invoker_cmdline_params  params;
  struct SMVM_invoker_traits          traits;
  struct SMVM_memory_hierarchy_traits memtraits;

  smvm_set_debug_level_from_environment ();

  smvm_get_cmdline_args (&params, argc, argv);

  /* Default to stdout if no filename given. */
  if (params.dataoutfile == NULL)
    output = stdout;
  else
    {
      output = fopen (params.dataoutfile, "a");
      if (output == NULL)
	{
	  fprintf (stderr, "*** ERROR: Failed to open output file %s! ***\n", 
		   params.dataoutfile);
	  die (-1);
	}
    }

  smvm_init_memory_hierarchy_traits (&memtraits, params.max_mem_KB);
  /* FIXME:  hardcoded fill! */
  smvm_construct_invocation (&traits, params.num_trials, 0.01, memtraits, output, warn, dbg);
  smvm_invoke (&traits, params.b_run_sparse, params.b_run_dense);

  /* Close the data output file. */
  if (fclose (output) != 0) 
    {
      fprintf (stderr, "*** ERROR: main (invoker):  Failed to close "
	       "data output file! ***\n");
      die (-1);
    }

  /** 
   * @bug Should return the appropriate error code!  (Doesn't matter for HPCC
   *      though, since some implementations of MPI always return 0.) 
   */
  return 0;
}
