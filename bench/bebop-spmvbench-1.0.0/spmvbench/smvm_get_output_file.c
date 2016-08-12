/**
 * @file smvm_get_output_file.c
 * Maintainer:  mfh
 * Created:  2004 Jan 24
 * Last modified: 2004 Mar 17
 ************************************************************************/
#include "smvm_get_output_file.h"
#include <smvm_util.h>
#include <stdlib.h>


/**
 * Returns nonzero if the given file exists and is readable, zero if not.
 */
int
file_exists (const char filename[])
{
  FILE* f = NULL;
  int b_exists = 0;

  f = fopen (filename, "r");
  b_exists = (f != NULL);
  
  if (b_exists)
    {
      if (fclose (f) != 0)
	fprintf (stderr, "*** ERROR: file_exists: failed to close file %s "
		 "while testing if it exists! ", filename);
    }

  return b_exists;
}


/************************************************************************/
FILE*
smvm_get_output_file (const char outfilename[], const char col_headers[], 
		      const int b_dbg)
{
  FILE* file = NULL;
  int b_file_exists = 0;

  if (b_dbg) fprintf (stderr, ">>> smvm_get_output_file:\n");

  if (b_dbg) fprintf (stderr, "Checking if data output file exists...");
  b_file_exists = file_exists (outfilename);

  if (b_file_exists)
    {
      if (b_dbg) fprintf (stderr, "it does.\nReopening it for appending...");

      /*
       * The file already exists, so assume the column headers are in it 
       * already.  Reopen it for appending.
       */
      file = fopen (outfilename, "a");
      if (file == NULL)
	{
	  die_with_message ("*** ERROR: smvm_get_output_file: Failed to open "
			    "output file for appending ", -1);
	}
      if (b_dbg) fprintf (stderr, "done.\n");
    }
  else
    {
      /* The file doesn't exist, so write the column headers to it. */
      if (b_dbg)
	fprintf (stderr, "output file %s does not yet exist.\n"
		 "Creating it and writing column headers...",
		 outfilename);
      
      file = fopen (outfilename, "a");
      if (file == NULL)
	{
	  die_with_message ("*** ERROR: smvm_get_output_file: Failed to open "
			    "output file for appending ", -1);
	}
      fprintf (file, "%s", col_headers);

      if (b_dbg) fprintf (stderr, "done.\n");
    }

  return file;
}

