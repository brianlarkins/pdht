/**
 * @file benchmark_tester.c
 * @author mfh
 * @date 2004 Dec 23
 * @since 2004 Sep 28
 * @brief Driver for low-level testing of SpMV benchmark.
 *
 * Driver for low-level testing of SpMV benchmark.  Generates a block CSR 
 * matrix with randomly distributed blocks, and runs a timing loop on it for 
 * the specified number of times.  For command-line arguments, see "usage" 
 * function.  Functions in this file are meant only for testing the SMVM 
 * benchmark, and not for actual invocation of it in practice.  Users are 
 * recommended to refer to smvm_hpcc_interface_tester.c for a driver that 
 * tests the SpMV benchmark at the same level of functionality as the HPC 
 * Challenge suite of benchmarks would use.
 *
 * Program return values:
 * <UL>
 * <LI> 0   Successful run
 * <LI> 1   Incorrect command-line parameters ("usage")
 * <LI> 2   One or more parameters out of bounds
 * <LI> -1  Some other fatal error
 * </UL>
 *****************************************************************************/
#include "config.h"
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include "smvm_benchmark.h"
#include <smvm_util.h>


/**
 * @brief Explains the command-line arguments.
 *
 * The classic "usage" message that explains the command-line arguments.
 * Output goes to stderr.  
 *
 * @param prog_name  Name of the executable.  Usually one passes in argv[0].
 *
 * @return Nothing.
 */
static void
usage (const char* prog_name)
{
  fprintf (stderr, "Sparse matrix-vector multiplication benchmark\n\n");
  fprintf (stderr, 
           "usage: %s <m> <n> <r> <c> <percent_fill> <num_trials> "
	   "<data_output_filename>\n\n", 
           prog_name);
  fprintf (stderr, 
           "Generates a block CSR matrix with randomly distributed blocks.\n"
           "<m>, <n>\n\t"
	   "Dimensions of the matrix\n"
           "<r>, <c>\n\t"
	   "Dimensions of the blocks\n"
           "<percent_fill>\n\t"
	   "Percentage of nonzeros in the matrix\n"
           "<num_trials>\n\t"
	   "Number of trials to run each test (if zero, does a dry run)\n"
           "<data_output_filename>\n\t"
	   "Name of (CSV format) data output file (created new each time)\n");
  fprintf (stderr, "\n\n");
}
  
  
/**
 * Shows an error message and the usage, and exits with the given error code.
 * Caller is responsible for any endlines in the error message.
 * Pass it argv[0] so the usage statement works right.
 *
 * @param message  Message to print before displaying usage.
 * @param argv0    Pass this argv[0], the invocation name.
 * @param errcode  What the executable returns (argument of exit(int)).
 *
 * @return  Does not return (calls exit(int)).
 */
static void
die_with_usage (const char* message, const char* argv0, const int errcode)
{
  fprintf (stderr, message);
  usage (argv0);
  exit (errcode);
}




/**
 * Get the command-line arguments (see `usage' for details).  We use
 * errno and the C standard functions that modify errno, to check for
 * correct conversion of command-line arguments.
 *
 * @param p_params  Holds the command-line arguments' values.
 * @param argc      Pass in the corresponding parameter of main().
 * @param argv      Pass in the corresponding parameter of main().
 */
static void
smvm_read_command_line_arguments (struct SMVM_parameters *p_params, 
				  int argc, char *argv[])
{ 
  const int min_num_args = 7;
  /* const int max_num_args = 7; */

  /* 
   * argc is the length of the argv array and thus includes 
   * the program invocation name.  That's why we check values
   * one larger than they seemingly should be.
   */
  if (argc < min_num_args + 1) 
    die_with_usage ("Too few arguments.\n", argv[0], 1);

  errno = 0;
  p_params->m = (int) strtol (argv[1], NULL, 10);
  if (errno != 0)
    die_with_usage ("Failed to extract <m> command-line argument!\n", argv[0], 1);

  errno = 0;
  p_params->n = (int) strtol (argv[2], NULL, 10);
  if (errno != 0)
    die_with_usage ("Failed to extract <n> command-line argument!\n", argv[0], 1);

  errno = 0;
  p_params->r = (int) strtol (argv[3], NULL, 10);
  if (errno != 0)
    die_with_usage ("Failed to extract <r> command-line argument!\n", argv[0], 1);

  errno = 0;
  p_params->c = (int) strtol (argv[4], NULL, 10);
  if (errno != 0)
    die_with_usage ("Failed to extract <c> command-line argument!\n", argv[0], 1);

  errno = 0;
  p_params->percent_fill = (double) strtod (argv[5], NULL);
  if (errno != 0)
    die_with_usage ("Failed to extract <percent_fill> command-line argument!\n", argv[0], 1);

  errno = 0;
  p_params->num_trials = (int) strtol (argv[6], NULL, 10);
  if (errno != 0)
    die_with_usage ("Failed to extract <num_trials> command-line argument!\n", argv[0], 1);

  if (argv[7] == NULL)
    die_with_usage ("<data_output_filename> parameter is empty!", argv[0], 1);
  else 
    {
      p_params->dataoutfile = fopen (argv[7], "a");
      if (p_params->dataoutfile == NULL)
	{
	  fprintf (stderr, "*** ERROR: Failed to open data output file %s "
		   "for appending! ***", argv[7]);
	  die (-1);
	}
    }
}



/****************************************************************************/
/*                                  MAIN                                    */
/****************************************************************************/
int
main (int argc, char *argv[])
{
  int retval;
  struct SMVM_parameters cl_args;

  smvm_set_debug_level_from_environment ();

  smvm_read_command_line_arguments (&cl_args, argc, argv);
  retval = smvm_benchmark (&cl_args);
  if (fclose (cl_args.dataoutfile) != 0)
    {
      fprintf (stderr, "*** ERROR:  Failed to close data output file! ***\n");
      die (-1);
    }
  return retval;
}

