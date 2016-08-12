/**
 * @file main.c
 * @author mfh
 * @since 14 Jun 2005
 * @date 11 Nov 2005
 *
 * Driver program for the sparse matrix generator.
 */


#include <config.h>
#include <bcsr_matrix.h>
#include <coo_matrix.h>
#include <get_options.h>
#include <random_number.h>
#include <smvm_malloc.h>
#include <smvm_util.h>
#include <timer.h>

#include <create_rand.h>
#include <dict.h>
#include <stencil.h>

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/** 
 * Map from function name to function pointer.
 */
static dict_t* fnname_to_fnptr = NULL;
/**
 * Map from function name to the number of arguments taken by that function.
 */
static dict_t* fnname_to_nargs = NULL;
/**
 * Map from "matrix type" as specified on the command line, to the corresponding function name.
 */
static dict_t* command_line_name_to_fnname = NULL;
/**
 * List of numbers of arguments taken by the various generator functions.
 */
static int nargs_array[] = {3, 3, 3, 5, 7, 8, 3, 3, 3};
/**
 * Number of matrix generator functions referenced in this file.
 */
static const int num_dict_entries = 9;
/**
 * "Matrix types" as specified on the command line.
 */
static char* mattypes[] = {
  "stencil-3pt-bcsr",
  "stencil-9pt-bcsr",
  "stencil-27pt-bcsr",
  "random",
  "random-banded",
  "random-banded-per-row", 
  "stencil-3pt-coo", 
  "stencil-9pt-coo", 
  "stencil-27pt-coo"
};
/**
 * Names of matrix generator functions.
 */
static char* fnnames[] = {
  "stencil_3pt_random_bcsr_matrix",
  "stencil_9pt_random_bcsr_matrix",
  "stencil_27pt_random_bcsr_matrix",
  "create_random_matrix",
  "create_random_matrix_banded",
  "create_random_matrix_banded_by_nnz_per_row_using_algorithm",
  "stencil_3pt_coo_matrix",
  "stencil_9pt_coo_matrix",
  "stencil_27pt_coo_matrix"
};
static char* matrix_type = NULL;


static void
usage (FILE* out, const struct arginfo* arglist, const struct arginfo* ext_args)
{
  int i;

  fprintf (out, "Usage:\n");
  fprintf (out, "sparse_matrix_converter [options] <matrix-type> [additional-params]\n");
  fprintf (out, "-v             If specified, use verbose mode\n");
  fprintf (out, "<matrix-type>: Type of the matrix to generate\n");
  fprintf (out, "[additional-params]: Parameters that go with the matrix type\n");
  fprintf (out, "\n");
  fprintf (out, "Here are the supported matrix types and the number of arguments required for each:\n");
  for (i = 0; i < num_dict_entries; i++)
    fprintf (out, "\t%s  %d\n", mattypes[i], nargs_array[i]);

  fprintf (out, "\n");
}



#if 0
static int
compare_ints (const void *pa, const void *pb)
{
  const int a = *((int*) pa);
  const int b = *((int*) pb);

  if (a < b)
    return -1;
  else if (a > b)
    return +1;
  else 
    return 0;
}
#endif

static void
deinit_dicts ()
{
  /* We have to free all the nodes before we can call dict_destroy. */
#define FREE_DICT( d ) do {if (d != NULL) { dict_free_nodes (d); dict_destroy (d); d = NULL; }} while(0)

  FREE_DICT( fnname_to_fnptr );
  FREE_DICT( fnname_to_nargs );
  FREE_DICT( command_line_name_to_fnname );
}

static void
init_dicts ()
{
  int i;

  fnname_to_fnptr = dict_create (num_dict_entries, (dict_comp_t) strcmp);
  fnname_to_nargs = dict_create (num_dict_entries, (dict_comp_t) strcmp);
  command_line_name_to_fnname = dict_create (num_dict_entries, (dict_comp_t) strcmp);

  dict_alloc_insert (fnname_to_fnptr, "create_random_matrix", 
		     &create_random_matrix);
  dict_alloc_insert (fnname_to_fnptr, "create_random_matrix_banded", 
		     &create_random_matrix_banded);
  dict_alloc_insert (fnname_to_fnptr, "create_random_matrix_banded_by_nnz_per"
		     "_row_using_algorithm", 
		     &create_random_matrix_banded_by_nnz_per_row_using_algorithm);
  dict_alloc_insert (fnname_to_fnptr, "stencil_3pt_random_bcsr_matrix", 
		     &stencil_3pt_random_bcsr_matrix);
  dict_alloc_insert (fnname_to_fnptr, "stencil_9pt_random_bcsr_matrix", 
		     &stencil_9pt_random_bcsr_matrix);
  dict_alloc_insert (fnname_to_fnptr, "stencil_27pt_random_bcsr_matrix", 
		     &stencil_27pt_random_bcsr_matrix);
  dict_alloc_insert (fnname_to_fnptr, "stencil_3pt_coo_matrix",
		     &stencil_3pt_coo_matrix);
  dict_alloc_insert (fnname_to_fnptr, "stencil_9pt_coo_matrix",
		     &stencil_9pt_coo_matrix);
  dict_alloc_insert (fnname_to_fnptr, "stencil_27pt_coo_matrix",
		     &stencil_27pt_coo_matrix);

  for (i = 0; i < num_dict_entries; i++)
    {
      dict_alloc_insert (fnname_to_nargs, fnnames[i], &nargs_array[i]);
      dict_alloc_insert (command_line_name_to_fnname, mattypes[i], fnnames[i]);
    }
}


/**
 * Calls the given function with the given arguments (first converting the 
 * arguments to integers), and returns a pointer to the resulting generated 
 * sparse matrix.  The caller is responsible for correctly identifying the 
 * type of the generated matrix.
 */
static void*
dispatch_matrix_generator (const char* const function_name, 
			   char *argument_list[], 
			   const int num_args_given)
{
  /* 
   * This declaration exploits the non-typesafeness of C.  In particular, it 
   * lets us call the function with any number of arguments that we want.  
   * In this context this is a feature and not a bug.
   */
  typedef void* (*fnptr_t) ();
  fnptr_t fnptr = NULL;
  dnode_t* dn = NULL;
  int nargs = 0;
  const void* temp;

  dn = dict_lookup (fnname_to_fnptr, function_name);
  if (dn == NULL)
    {
      fprintf (stderr, "*** dispatch_matrix_generator: failed to find "
	       "function name in map! ***\n");
      return NULL;
    }
  fnptr = (fnptr_t) (dnode_get (dn));
  if (fnptr == NULL)
    {
      fprintf (stderr, "*** dispatch_matrix_generator: failed to get "
	       "function from node! ***\n");
      return NULL;
    }
  dn = dict_lookup (fnname_to_nargs, function_name);
  if (dn == NULL)
    {
      fprintf (stderr, "*** dispatch_matrix_generator: failed to find "
	       "function name in nargs map! ***\n");
      return NULL;
    }
  temp = dnode_get (dn);
  if (temp == NULL)
    {
      fprintf (stderr, "*** dispatch_matrix_generator: failed to get nargs "
	       "from node! ***\n");
      return NULL;
    }
  nargs = *((int*) temp);

  if (nargs < 0)
    {
      fprintf (stderr, "*** dispatch_matrix_generator: nargs < 0 !!! ***\n");
      return NULL;
    }
  else if (nargs > num_args_given)
    {
      fprintf (stderr, "*** dispatch_matrix_generator: %d args required,"
	       " but only %d given! ***\n", nargs, num_args_given);
      return NULL;
    }
  else if (nargs < num_args_given)
    {
      fprintf (stderr, "*** dispatch_matrix_generator: WARNING: %d excess "
	       "arguments given! ***\n", num_args_given - nargs);
    }

  if (nargs == 0)
    return fnptr ();
  else if (nargs == 1)
    return fnptr (strtol (argument_list[0], NULL, 10));
  else if (nargs == 2)
    return fnptr (strtol (argument_list[0], NULL, 10),
		  strtol (argument_list[1], NULL, 10));
  else if (nargs == 3)
    return fnptr (strtol (argument_list[0], NULL, 10),
		  strtol (argument_list[1], NULL, 10),
		  strtol (argument_list[2], NULL, 10));
  else if (nargs == 4)
    return fnptr (strtol (argument_list[0], NULL, 10),
		  strtol (argument_list[1], NULL, 10),
		  strtol (argument_list[2], NULL, 10),
		  strtol (argument_list[3], NULL, 10));
  else if (nargs == 5)
    return fnptr (strtol (argument_list[0], NULL, 10),
		  strtol (argument_list[1], NULL, 10),
		  strtol (argument_list[2], NULL, 10),
		  strtol (argument_list[3], NULL, 10),
		  strtol (argument_list[4], NULL, 10));
  else if (nargs == 6)
    return fnptr (strtol (argument_list[0], NULL, 10),
		  strtol (argument_list[1], NULL, 10),
		  strtol (argument_list[2], NULL, 10),
		  strtol (argument_list[3], NULL, 10),
		  strtol (argument_list[4], NULL, 10),
		  strtol (argument_list[5], NULL, 10));
  else if (nargs == 7)
    return fnptr (strtol (argument_list[0], NULL, 10),
		  strtol (argument_list[1], NULL, 10),
		  strtol (argument_list[2], NULL, 10),
		  strtol (argument_list[3], NULL, 10),
		  strtol (argument_list[4], NULL, 10),
		  strtol (argument_list[5], NULL, 10),
		  strtol (argument_list[6], NULL, 10));
  else if (nargs == 8)
    return fnptr (strtol (argument_list[0], NULL, 10),
		  strtol (argument_list[1], NULL, 10),
		  strtol (argument_list[2], NULL, 10),
		  strtol (argument_list[3], NULL, 10),
		  strtol (argument_list[4], NULL, 10),
		  strtol (argument_list[5], NULL, 10),
		  strtol (argument_list[6], NULL, 10),
		  strtol (argument_list[7], NULL, 10));
  else
    {
      fprintf (stderr, "*** dispatch_matrix_generator: unsupported number of args! ***\n");
      return NULL;
    }
}


static void*
create_matrix (int argc, char *argv[], struct arginfo* arglist)
{
  extern int optind;
  int num_args = argc - optind - 1;
  dnode_t* dn = NULL;
  char* fnname = NULL;

  if (argc - optind < 1)
    {
      fprintf (stderr, "*** No matrix type specified ***\n");
      dump_usage (stderr, argv[0], arglist, NULL);
      return NULL;
    }
  matrix_type = argv[optind];

  dn = dict_lookup (command_line_name_to_fnname, matrix_type);
  if (dn == NULL)
    {
      fprintf (stderr, "*** create_matrix: Invalid matrix type %s ***\n", matrix_type);
      return NULL;
    }
  fnname = (char*) (dnode_get (dn));
  if (fnname == NULL)
    {
      fprintf (stderr, "*** create_matrix: failed to get function name from node! ***\n");
      return NULL;
    }

  return dispatch_matrix_generator (fnname, &argv[optind+1], num_args);
}


struct arginfo *arglist = NULL;


int
init (int argc, char** argv)
{
  smvm_set_debug_level_from_environment ();
  /* Set the get_options usage function to "usage", instead of using the default
   * usage function.  This is necessary because the command-line arguments include
   * things that are not "options" in the strict sense, because they do not follow
   * a "-[a-z]" pattern.  */
  register_usage_function (usage);
  arglist = register_arginfo (arglist, 'p', NULLARG, NULL, "If specified, pr"
			      "int the resulting matrix in MatrixMarket form"
			      "at to stdout");
  arglist = register_arginfo (arglist, 'v', NULLARG, NULL, "If specified, ac"
			      "tivate verbose mode");
  get_options (argc, argv, arglist, NULL);

  if (got_arg_p (arglist, 'v'))
    {
      printf ("Calling init_dicts...\n");
      fflush (stdout);
    }
  init_dicts ();

  /* Initialize the random number generator */
  smvm_init_randomizer (1);

  init_timer ();

  return EXIT_SUCCESS;
}

int
deinit ()
{
  if (got_arg_p (arglist, 'v'))
    {
      printf ("Calling deinit_dicts...\n");
      fflush (stdout);
    }
  deinit_dicts ();
  destroy_arginfo_list (arglist);
  deinit_timer ();

  return EXIT_SUCCESS;
}


int
main (int argc, char *argv[])
{
  void* A = NULL;
  int errcode = 0;
  double t = 0.0;

  errcode = init (argc, argv);
  if (errcode == EXIT_FAILURE)
    return errcode;

  if (got_arg_p (arglist, 'v'))
    {
      printf ("Calling create_matrix...\n");
      fflush (stdout);
    }
  t = get_seconds ();
  A = create_matrix (argc, argv, arglist);
  t = get_seconds () - t;
  printf ("Time to generate matrix:  %g seconds\n", t);
  if (A == NULL)
    {
      fprintf (stderr, "*** Failed to create sparse matrix! ***\n");
      exit (EXIT_FAILURE);
    }
  else if (got_arg_p (arglist, 'v'))
    {
      printf ("Successfully created sparse matrix.\n");
      fflush (stdout);
    }

  if (strcmp (matrix_type, "stencil-3pt-coo") == 0 ||
      strcmp (matrix_type, "stencil-9pt-coo") == 0 ||
      strcmp (matrix_type, "stencil-27pt-coo") == 0)
    {
      if (got_arg_p (arglist, 'p'))
	{
	  errcode = print_coo_matrix_in_matrix_market_format (stdout, (struct coo_matrix_t*) A);
	  fflush (stdout);

	  if (errcode != 0)
	    fprintf (stderr, "*** Error %d in printing COO matrix! ***\n", errcode);
	}

      if (got_arg_p (arglist, 'v'))
	{
	  printf ("Calling destroy_coo_matrix...\n");
	  fflush (stdout);
	}
      destroy_coo_matrix (A);
    }
  else
    {
      if (got_arg_p (arglist, 'p'))
	{
	  errcode = print_bcsr_matrix_in_matrix_market_format (stdout, (struct bcsr_matrix_t*) A);
	  fflush (stdout);

	  if (errcode != 0)
	    fprintf (stderr, "*** Error %d in printing BCSR matrix! ***\n", errcode);
	}

      if (got_arg_p (arglist, 'v'))
	{
	  printf ("Calling destroy_bcsr_matrix...\n");
	  fflush (stdout);
	}
      destroy_bcsr_matrix (A);
    }

  return deinit ();
}



