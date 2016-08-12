# Makefile
# Author: Mark Hoemmen
# Since: 26 June 2005
# Date:  09 June 2006
#
# Utility functions for the BeBOP group.
#
# Copyright (c) 2006, Regents of the University of California 
# All rights reserved.
# Redistribution and use in source and binary forms, with or
# without modification, are permitted provided that the
# following conditions are met:
# 
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright 
#   notice, this list of conditions and the following disclaimer in 
#   the documentation and/or other materials provided with the 
#   distribution.
#
# * Neither the name of the University of California, Berkeley, nor
#   the names of its contributors may be used to endorse or promote
#   products derived from this software without specific prior
#   written permission.  
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
# STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
# OF THE POSSIBILITY OF SUCH DAMAGE.
##########################################################################

# Platform-specific configuration options.
ARCHIVER = ar
ARFLAGS = r
CC = gcc
CFLAGS = -fPIC -O3 -Wall -std=c99 
#CFLAGS = -fpic -O3 -Wall -std=c99 
CPPFLAGS = -I.
LIBS = 
LINKER = gcc
RANLIB = echo
RM = rm -f

#
# For the "utility library," to be shared between projects.
# 
LIBUTILOBJ = __complex.o avltree_intpair.o get_options.o merge_sort.o \
	     mt19937ar.o random_number.o smvm_malloc.o smvm_util.o \
	     sort_joint_arrays.o timer.o 
LIBUTIL_STATIC = libbebop_util.a
LIBUTIL_DYNAMIC = libbebop_util.dylib
#LIBUTIL_DYNAMIC = libbebop_util.so

#########################################################################
## OUTPUT RULES
#########################################################################

all: $(LIBUTIL_STATIC) tests/test_random_number tests/test_split_pathname

dynamic: $(LIBUTIL_DYNAMIC)

# FIXME: only works for gcc!
#$(LIBUTIL_DYNAMIC): $(LIBUTILOBJ)
#	$(CC) -shared $(CFLAGS) $(LIBUTILOBJ) -o $(LIBUTIL_DYNAMIC)
$(LIBUTIL_DYNAMIC): $(LIBUTILOBJ)
	$(CC) -dynamiclib $(CFLAGS) $(LIBUTILOBJ) -o $(LIBUTIL_DYNAMIC)

$(LIBUTIL_STATIC): $(LIBUTILOBJ)
	$(ARCHIVER) $(ARFLAGS) $(LIBUTIL_STATIC) $(LIBUTILOBJ)
	$(RANLIB) $(LIBUTIL_STATIC)

tests/test_random_number: $(LIBUTIL_STATIC) tests/test_random_number.o
	$(LINKER) $(CPPFLAGS) $(CFLAGS) -o tests/test_random_number \
		tests/test_random_number.o $(LIBUTIL_STATIC) $(LIBS)

tests/test_split_pathname: $(LIBUTIL_STATIC) tests/test_split_pathname.o
	$(LINKER) $(CPPFLAGS) $(CFLAGS) -o tests/test_split_pathname \
		tests/test_split_pathname.o $(LIBUTIL_STATIC) $(LIBS)

tests/test_split_pathname.o: tests/test_split_pathname.c \
	smvm_malloc.h smvm_util.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -c tests/test_split_pathname.c \
		-o tests/test_split_pathname.o


#########################################################################
## BUILD RULES 
#########################################################################

__complex.o: __complex.c __complex.h config.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -c __complex.c

avltree_intpair.o: avltree_intpair.c avltree_intpair.h smvm_malloc.h \
	smvm_util.h config.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -c avltree_intpair.c

get_options.o: get_options.c get_options.h smvm_malloc.h smvm_util.h config.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -c get_options.c

merge_sort.o: merge_sort.c merge_sort.h smvm_malloc.h smvm_util.h config.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -c merge_sort.c

mt19937ar.o: mt19937ar.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c mt19937ar.c

random_number.o: random_number.c random_number.h smvm_malloc.h smvm_util.h config.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -c random_number.c

smvm_malloc.o: smvm_malloc.c smvm_malloc.h smvm_util.h config.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -c smvm_malloc.c

smvm_util.o: smvm_util.c smvm_util.h random_number.h smvm_malloc.h config.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -c smvm_util.c

sort_joint_arrays.o: sort_joint_arrays.c sort_joint_arrays.h random_number.h \
	smvm_malloc.h smvm_util.h config.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -c sort_joint_arrays.c

timer.o: timer.c timer.h config.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -c timer.c

##########################################################################
## UTILITY STEPS
##########################################################################

# FIXME:  Depends on GNU tar features!!!
backup:
	tar --create --gzip --file=bebop_util.tar.gz *.[ch] Makefile Doxyfile README COPYING BUGS tests/

ctags: 
	ctags *.[ch]

etags:
	etags *.[ch]

clean:
	$(RM) $(LIBUTILOBJ)
	$(RM) $(LIBUTIL_STATIC)
	$(RM) $(LIBUTIL_DYNAMIC)
	$(RM) tests/test_random_number tests/test_random_number.o tests/test_split_pathname tests/test_split_pathname.o
