/**
 * @file test_split_pathname.c
 * @author Mark Hoemmen
 * @since 27 Jul 2005
 * @date 23 Feb 2006
 *
 * Driver program for testing split_pathname, in smvm_util.h.  Returns 
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
#include "../smvm_malloc.h"
#include "../smvm_util.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void
test_split_pathname (const char* const path, 
		     const char* const expected_parentdir, 
		     const char* const expected_namestem, 
		     const char* const expected_extn)
{
  char* parentdir = NULL;
  char* namestem = NULL;
  char* extn = NULL;

  split_pathname (&parentdir, &namestem, &extn, path);
  assert (parentdir != NULL);
  assert (namestem != NULL);
  assert (extn != NULL);
  assert (0 == strcmp (parentdir, expected_parentdir));
  assert (0 == strcmp (namestem, expected_namestem));
  assert (0 == strcmp (extn, expected_extn));

  smvm_free (parentdir);
  smvm_free (namestem);
  smvm_free (extn);
}


int 
main (int argc, char** argv)
{
  printf ("=== test_split_pathname ===\n");

  if (argc >= 5)
    test_split_pathname (argv[1], argv[2], argv[3], argv[4]);
  else
    {
      test_split_pathname ("/usr/lib/libgen.h", "/usr/lib", "libgen", "h");
      test_split_pathname ("foo/bar123_/ugabuga", "foo/bar123_", "ugabuga", "");
    }

  printf ("=== Passed test_split_pathname ===\n");
  return EXIT_SUCCESS;
}

