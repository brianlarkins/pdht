/**
 * @file get_options.c
 * @author mfh
 * @date 22 Jul 2005
 *
 * Implementation of a library for parsing command-line arguments.  
 *
 * @note Implementation uses the POSIX standard unistd.h getopt functions, 
 * but does not use any GNU extensions (getopt.h).
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
 ***************************************************************************/
#include <config.h>
#include <get_options.h>
#include <smvm_malloc.h>
#include <smvm_util.h>

#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h> 

#ifdef DMALLOC
#  define WITH_DMALLOC( OP )  do {OP;} while(0)
#else
#  define WITH_DMALLOC( OP )
#endif /* DMALLOC */

/**
 * A user-specified usage function.  If not specified, a default is used.
 */
void (*the_usage_function) (FILE*, const struct arginfo*, const struct arginfo*) = NULL;




/**
 * Describes and stores the value of a single command-line argument, and also 
 * defines a linked list of these command-line arguments.
 *
 * @note Clients never need to see the internals of this struct.
 */
struct arginfo 
{
  /**
   * [IN] Command-line option name (e.g. 't' for "-t <blah>").
   */
  char c;

  /**
   * [IN] Type of the command-line argument.
   */
  arg_type_t type;

  /**
   * [OUT] Pointer to the value of the command-line argument.  Its
   * initial value is the default value of the parameter.  
   */
  void* val;

  /**
   * (Pointer to the) default value of the command-line argument.
   * A deep copy of the default value specified by register_arginfo.
   */
  union {
    int i;
    double d;
    char *s;
  } default_val;

  /**
   * [IN] English description of the command-line argument.
   * Suitable for a `usage' dump if the user gets the arguments wrong.
   */
  char* desc;

  /**
   * [OUT] Nonzero if we got this argument, else zero.  Should be init'd to zero.
   */
  int got_arg;

  /**
   * [OUT] Nonzero if we got a value for this argument, else zero.
   * Should be init'd to zero.
   */
  int got_arg_value;

  /**
   * Pointer to the next command-line argument in the linked list (NULL if 
   * this is the last one).
   */
  struct arginfo* next;
};

/**
 * Returns the default value of the given command-line argument.
 */
static int 
get_int_default_argval (const struct arginfo* arg)
{
  if (arg->type != INTARG)
    {
      fprintf (stderr, "*** Attempt to interpret parameter of command-line "
	       "argument -%c as integer, though it is not ***\n", arg->c);
      smvm_exit (EXIT_FAILURE);
    }
  return arg->default_val.i;
}

/**
 * Returns the default value of the given command-line argument.
 */
double
get_double_default_argval (const struct arginfo* arg)
{
  if (arg->type != DOUBLEARG)
    {
      fprintf (stderr, "*** Attempt to interpret parameter of command-line "
	       "argument -%c as a double, though it is not ***\n", arg->c);
      smvm_exit (EXIT_FAILURE);
    }
  return arg->default_val.d;
}


/**
 * Returns the default value of the given command-line argument.
 */
char*
get_string_default_argval (const struct arginfo* arg)
{
  if (arg->type != STRINGARG)
    {
      fprintf (stderr, "*** Attempt to interpret parameter of command-line "
	       "argument -%c as a string, though it is not ***\n", arg->c);
      smvm_exit (EXIT_FAILURE);
    }
  return arg->default_val.s;
}

/**
 * Returns a pointer to the arginfo object in the given list with the given 
 * label, or NULL if none such object is found in the list.
 */
struct arginfo*
find_arginfo (struct arginfo *list, const char c)
{
  if (list == NULL)
    return NULL;

  while (list != NULL)
    {
      if (list->c == c)
	return list;

      list = list->next;
    }

  return NULL;
}


/**
 * Returns a dynamically allocated copy of the given string.
 */
static char*
create_string_copy (const char* str)
{
  const int len = strlen (str);
  char* out = smvm_malloc ((len + 1) * sizeof (char));
  int i;

  /* strcpy() is not trustworthy, because it could go on "forever" if str is not
   * zero-terminated.  Of course strlen() has the same problem... */
  for (i = 0; i < len; i++)
    out[i] = str[i];

  /* Ensure out is zero-terminated */
  out[len] = '\0';

  return out;
}



/**
 * Dynamically allocates an arginfo object, stores in it the given parameters,
 * and returns a pointer to the object.
 *
 * @param c [IN]     Command-line option name (e.g. 't' for "-t <blah>")
 *
 * @param type [IN]  Type of the command-line argument
 *
 * @param val [IN]   Pointer to the default value of the command-line argument
 *                   e.g. if type==DOUBLEARG or INTARG, val is a pointer to a 
 *                   variable holding the initial value.  If type==STRINGARG, 
 *                   then val is a char* holding the string.  In all of these 
 *                   cases, value will be copied (a deep copy), and val must 
 *                   be a valid pointer.  (You must always set a default value 
 *                   for any parameter.)  If type==NULLARG, val will not be 
 *                   referenced.  
 *
 * @param desc [IN]  English description of the command-line argument.  String
 *                   will be (deeply) copied.
 *
 * @param next [IN]  Pointer to the "next" arginfo object in the list 
 *
 * @return Pointer to freshly allocated arginfo object
 */
struct arginfo* 
create_arginfo (const char c, const arg_type_t type, const void *val, 
		const char *desc, struct arginfo *next);

/**
 * Deallocates the dynamic memory used by the given arginfo object, without 
 * affecting the objects in the tail of the list.
 */
void
deinit_arginfo (struct arginfo *arg);

/**
 * Deallocates the head of the given arginfo list and all the dynamic memory 
 * that the head uses, without affecting the objects in the tail of the list.
 */
void
destroy_arginfo (struct arginfo *head);




/*======================================================================*/
void 
set_int_argval (struct arginfo* arg, const int val)
{
  *((int*) (arg->val)) = val;
}

/*======================================================================*/
void
set_double_argval (struct arginfo* arg, const double val)
{
  *((double*) (arg->val)) = val;
}

/*======================================================================*/
void
set_string_argval (struct arginfo* arg, char* val)
{
  if (arg->val)
    smvm_free (arg->val);

  arg->val = create_string_copy (val);
}

/**
 * Dump usage for all the command-line arguments.  
 *
 * @param out [OUT]   Valid output stream
 * @param arg [IN]    List of command-line arguments
 */
void
dump_args (FILE* out, const struct arginfo* arg)
{
  for (; arg != NULL; arg = arg->next) 
    {
      fprintf (out, "\t-%c: %s", arg->c, arg->desc);
      switch (arg->type) 
	{
	case INTARG:
	  fprintf (out, " (default: %d)\n", get_int_default_argval (arg));
	  break;
	case DOUBLEARG:
	  fprintf (out, " (default: %g)\n", get_double_default_argval (arg));
	  break;
	case STRINGARG:
	    {
	      const char *const s = get_string_default_argval (arg);
	      if (s != NULL && strlen (s) > 0)
		fprintf (out, " (default: %s)\n", s);
	      else
		fprintf (out, " (default: none)\n");
	    }
	  break;
	case NULLARG:
	default:
	  fprintf (out, "\n");
	  break;
	}
    }
}


/**
 * Dump usage for a single command-line argument.
 */
void
dump_arg (FILE* out, const struct arginfo* arg)
{
  if (arg == NULL) return;

  fprintf (out, "\t-%c: %s", arg->c, arg->desc);
  switch (arg->type) 
    {
    case INTARG:
      fprintf (out, " (default: %d)\n", get_int_default_argval (arg));
      break;
    case DOUBLEARG:
      fprintf (out, " (default: %g)\n", get_double_default_argval (arg));
      break;
    case STRINGARG:
	{
	  const char *const s = get_string_default_argval (arg);
	  if (s != NULL && strlen (s) > 0)
	    fprintf (out, " (default: %s)\n", s);
	  else
	    fprintf (out, " (default: none)\n");
	}
      break;
    case NULLARG:
    default:
      fprintf (out, "\n");
      break;
    }
}


/*======================================================================*/
void
dump_usage (FILE* out, 
	    const char *const pathname, 
	    struct arginfo *core_args, 
	    struct arginfo *ext_args)
{
  const struct arginfo* arg = NULL;

  if (the_usage_function != NULL)
    {
      the_usage_function (out, core_args, ext_args);
      return;
    }

  fprintf(out, "%s [-h|-?]", pathname);
  if (ext_args)
    for (arg = ext_args; arg != NULL; arg = arg->next) 
      {
	fprintf(out, " [-%c%s]", arg->c, 
		(INTARG == arg->type ? " #" : 
		 (DOUBLEARG == arg->type? " #.##" : 
		  (STRINGARG == arg->type? " <string>" :
		   ""))));
      }
  if (core_args)
    for (arg = core_args; arg != NULL; arg = arg->next) 
      {
	fprintf(out, " [-%c%s]", arg->c,
		(INTARG == arg->type ? " #" :
		 (DOUBLEARG == arg->type ? " #.##" :
		  (STRINGARG == arg->type ? " <string>" : 
		   ""))));
      }
  fprintf (out, "\n");

  fprintf (out, "\t-h or -?: Display this.\n");

  if (ext_args) 
    dump_args (out, ext_args);

  if (core_args)
    dump_args (out, core_args);
}


/*======================================================================*/
static int
parse_arg (int c, struct arginfo* arg)
{
  extern char *optarg;

  for (; arg != NULL; arg = arg->next) 
    {
      if (c == arg->c) 
	{
	  if (arg->type == NULLARG)
	    {
	      arg->got_arg = 1;
	      arg->got_arg_value = 0;
	      return 1;
	    }
	  else if (arg->type == INTARG)
	    {
	      const int tmp = (int) strtol (optarg, NULL, 10);
	      errno = 0;
	      if (errno)
		{
		  fprintf (stderr, "*** Error in converting parameter %s of "
			   "command-line option -%c to integer ***\n", optarg, c);
		  smvm_exit (EXIT_FAILURE);
		}
	      set_int_argval (arg, tmp);
	      arg->got_arg = 1;
	      arg->got_arg_value = 1;
	      return 1;
	    }
	  else if (arg->type == DOUBLEARG)
	    {
	      const double tmp = strtod (optarg, NULL);
	      errno = 0;
	      if (errno)
		{
		  fprintf (stderr, "*** Error in converting parameter %s of "
			   "command-line option -%c to double ***\n", optarg, c);
		  smvm_exit (EXIT_FAILURE);
		}
	      set_double_argval (arg, tmp);
	      arg->got_arg = 1;
	      arg->got_arg_value = 1;
	      return 1;
	    }
	  else if (arg->type == STRINGARG) 
	    {
	      set_string_argval (arg, optarg);
	      arg->got_arg = 1;
	      arg->got_arg_value = 1;
	      return 1;
	    }
	  else 
	    return 0;
	}
    }
  return 0;
}

/*======================================================================*/
void
get_options (int argc, 
	     char **argv, 
	     struct arginfo *core_args, 
	     struct arginfo *ext_args)
{
  extern int getopt (int argc, char* const argv[], const char* optstring);
  int c;

  extern char *optarg;
  /* extern int optind, opterr, optopt; */

  char* optstring;
  char* cptr;
  int len_optstring = 2;
  const struct arginfo* arg;

  /* Construct the options string for getopt() */

  /* int dbg_count = 0; */
  for (arg = core_args; arg != NULL; arg = arg->next) 
    {
      ++len_optstring;
      if (NULLARG != arg->type) ++len_optstring;
    }
  if (ext_args)
    for (arg = ext_args; arg != NULL; arg = arg->next) 
      {
	++len_optstring;
	if (NULLARG != arg->type) ++len_optstring;
      }

  optstring = smvm_malloc ((len_optstring+1) * sizeof(char));
  optstring[0] = 'h';
  optstring[1] = '?';
  optstring[len_optstring] = '\000';

  cptr = optstring + 2;
  for (arg = core_args; arg != NULL; arg = arg->next) 
    {
      *cptr++ = arg->c;
      if (NULLARG != arg->type) *cptr++ = ':';
    }
  if (ext_args)
    for (arg = ext_args; arg != NULL; arg = arg->next) 
      {
	*cptr++ = arg->c;
	if (NULLARG != arg->type) *cptr++ = ':';
      }


  /* Read in the command-line options using getopt() */

  while (-1 != (c = getopt(argc, argv, optstring))) 
    {
      if ('h' == c || '?' == c) 
	{
	  dump_usage (stderr, argv[0], core_args, ext_args);
	  smvm_exit (EXIT_SUCCESS);
	}
      if (parse_arg (c, core_args)) continue;
      if (ext_args && parse_arg (c, ext_args)) continue;

      dump_usage (stderr, argv[0], core_args, ext_args);
      smvm_exit (EXIT_FAILURE);
    }

  smvm_free (optstring);
}



/*======================================================================*/
int 
get_int_argval (const struct arginfo* arg)
{
  if (arg->type != INTARG)
    {
      fprintf (stderr, "*** Attempt to interpret parameter of command-line "
	       "argument -%c as integer, though it is not ***\n", arg->c);
      smvm_exit (EXIT_FAILURE);
    }
  return *((int*) (arg->val));
}


/*======================================================================*/
double
get_double_argval (const struct arginfo* arg)
{
  if (arg->type != DOUBLEARG)
    {
      fprintf (stderr, "*** Attempt to interpret parameter of command-line "
	       "argument -%c as a double, though it is not ***\n", arg->c);
      smvm_exit (EXIT_FAILURE);
    }
  return *((double*) (arg->val));
}

/*======================================================================*/
char*
get_string_argval (const struct arginfo* arg)
{
  if (arg->type != STRINGARG)
    {
      fprintf (stderr, "*** Attempt to interpret parameter of command-line "
	       "argument -%c as a string, though it is not ***\n", arg->c);
      smvm_exit (EXIT_FAILURE);
    }
  return (char*) arg->val;
}


/*======================================================================*/
struct arginfo*
create_arginfo (const char c, const arg_type_t type, const void *val, 
		const char *desc, struct arginfo *next)
{
  struct arginfo *arg = smvm_malloc (sizeof (struct arginfo));

  arg->c = c;
  arg->type = type;

  if (type == INTARG)
    {
      arg->val = smvm_malloc (sizeof (int));
      *((int*)(arg->val)) = *((int*) val); /* Ewww pointer yuckiness */
      arg->default_val.i = *((int*) val);
    }
  else if (type == DOUBLEARG)
    {
      arg->val = smvm_malloc (sizeof (double));
      *((double*)(arg->val)) = *((double*) val);
      arg->default_val.d = *((double*) val);
    }
  else if (type == STRINGARG)
    {
      arg->val = create_string_copy (val);
      arg->default_val.s = create_string_copy (val);
    }
  else if (type == NULLARG)
    {
      arg->val = NULL;
      arg->default_val.s = (char*) NULL;
    }
  else 
    {
      fprintf (stderr, "*** init_arginfo: Invalid type %d for command-line "
	       "option %c ***\n", type, c);
      smvm_exit (EXIT_FAILURE);
    }

  arg->desc = create_string_copy (desc);
  arg->got_arg = 0;
  arg->got_arg_value = 0;
  arg->next = next;

  return arg;
}


/*======================================================================*/
void
destroy_arginfo_list (struct arginfo *list)
{
  struct arginfo *tail;

  if (list != NULL)
    {
      tail = list->next;
      destroy_arginfo (list);
      destroy_arginfo_list (tail);
    }
}


/*======================================================================*/
void
destroy_arginfo (struct arginfo *head)
{
  deinit_arginfo (head);
  smvm_free (head);
}

/*======================================================================*/
void
deinit_arginfo (struct arginfo *arg)
{
#ifdef DMALLOC
#  include <sys/types.h>
#  define DMALLOC_SIZE size_t
#  define DMALLOC_PNT  void *
#  define DMALLOC_VERIFY_ERROR	  0   /* function failed */
#  define DMALLOC_VERIFY_NOERROR  1   /* function succeeded */

  extern void dmalloc_message(const char *format, ...);
  extern int dmalloc_verify(const DMALLOC_PNT pnt);

  extern int dmalloc_examine(const DMALLOC_PNT pnt, DMALLOC_SIZE *size_p,
	 char **file_p, unsigned int *line_p, DMALLOC_PNT *ret_attr_p);
#endif


  WITH_DMALLOC( dmalloc_message (" "); dmalloc_message ("=== deinit_arginfo ===") );

  if (arg == NULL)
    {
      WITH_DMALLOC( dmalloc_message ("### arg == NULL, returning ###\n") );
      return;
    }

  if (arg->val != NULL && arg->type != NULLARG)
    {
#ifdef DMALLOC
      dmalloc_message ("FREEING argument %c, pointer val = %p\n", arg->c, arg->val);

      if (dmalloc_verify (arg->val) == DMALLOC_VERIFY_ERROR)
	dmalloc_message ("*** ERROR:  pointer is INVALID ***\n");
      else
	{
	  DMALLOC_SIZE size;
	  char *filename;
	  int linenum;
	  DMALLOC_PNT ret_addr;
	  unsigned long user_mark;
	  unsigned long seen;
	  int errcode;

	  errcode = dmalloc_examine (arg->val, &size, &filename, &linenum, &ret_addr);
	  dmalloc_message ("### Pointer status:  size %d, from file %s line %d ###\n", 
			   size, filename, linenum);
	}
#endif /* DMALLOC */

      smvm_free (arg->val);
      arg->val = NULL;
      if (arg->type == STRINGARG && arg->default_val.s != NULL)
	{
	  smvm_free (arg->default_val.s);
	  arg->default_val.s = NULL;
	}
    }

  if (arg->desc != NULL)
    {
#ifdef DMALLOC 
      dmalloc_message ("FREEING argument %c, pointer desc = %p\n", arg->c, arg->desc);

      if (dmalloc_verify (arg->desc) == DMALLOC_VERIFY_ERROR)
	dmalloc_message ("*** ERROR:  pointer is INVALID ***\n");
      else
	{
	  DMALLOC_SIZE size;
	  char *filename;
	  int linenum;
	  DMALLOC_PNT ret_addr;
	  int errcode;

	  errcode = dmalloc_examine (arg->desc, &size, &filename, &linenum, &ret_addr);
	  dmalloc_message ("### Pointer status:  size %d, from file %s line %d ###\n", 
			   size, filename, linenum);
	}
#endif /* DMALLOC */

      smvm_free (arg->desc);
      arg->desc = NULL;
    }
}


/*======================================================================*/
struct arginfo*
register_arginfo (struct arginfo *arginfo_list, 
		  const char c, 
		  const arg_type_t type, 
		  const void *val, 
		  const char *desc)
{
  struct arginfo *newhead = create_arginfo (c, type, val, desc, arginfo_list);

  return newhead;
}


/*======================================================================*/
int
got_arg_p (struct arginfo *list, const char c)
{
  struct arginfo *arg = find_arginfo (list, c); 

  if (arg == NULL)
    return 0;
  else if (arg->got_arg)
    return 1;
  /* else */
  return 0;
}


/*======================================================================*/
int
got_arg_value_p (struct arginfo *list, const char c)
{
  struct arginfo *arg = find_arginfo (list, c); 

  if (arg == NULL)
    return 0;
  else if (arg->got_arg_value)
    return 1;
  /* else */
  return 0;
}


void 
register_usage_function (void (usage_function) (FILE*, const struct arginfo*, const struct arginfo*))
{
  the_usage_function = usage_function;
}
