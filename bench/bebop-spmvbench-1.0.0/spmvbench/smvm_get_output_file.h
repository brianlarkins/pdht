#ifndef _smvm_get_output_file_h
#define _smvm_get_output_file_h
/**
 * @file smvm_get_output_file.h
 * @author mfh
 * @since 2004 Jan 24
 * @date 2004 Oct 02
 ***********************************************************************/
#include "config.h"
#include <stdio.h>  /* FILE */

/**
 * Returns an open file handle to the data file with the given name.  If the
 * data file does not yet exist, creates it and writes the given line of column
 * headers to it.  If unable to open the file, aborts (so the file handle
 * returned is guaranteed good).
 *
 * @param outfilename   Path (name) of the data output file to open.
 * @param col_headers   Line of headers to put at top of data output file. 
 *                      Endline is client-supplied.
 * @param b_dbg         Set to nonzero to display copious debug output.
 * @return  Open file handle to data output file.
 */
FILE* 
smvm_get_output_file (const char outfilename[], const char col_headers[], const int b_dbg);


#endif /* NOT _smvm_get_output_file_h */
