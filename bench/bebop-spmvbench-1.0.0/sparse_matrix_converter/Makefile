# Makefile
# Author: Mark Hoemmen
# Since: 26 May 2005
# Date: 10 Jun 2006
#
# For a sparse matrix utility that converts between various file formats 
# for sparse matrices, and also between various internal storage formats
# for sparse matrices.  Harwell-Boeing, MatrixMarket and Matlab file
# formats are supported.
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
# Linux-specific
#CFLAGS = -fpic -O3 -Wall -std=c99 
# Darwin-specific
CFLAGS = -fPIC -O3 -Wall -std=c99 
##
## Define USING_VALGRIND_CLIENT_REQUESTS and give an include directory 
## which has a valgrind/ subdirectory, if you want to use some features 
## of Valgrind (see http://valgrind.org) to help debug this code.  You
## probably need to use gcc, and Valgrind currently is only (well-) 
## supported on x86 platforms.
##
#CPPFLAGS = -I../bebop_util -I. -DUSING_VALGRIND_CLIENT_REQUESTS=1 -I/home/eecs/mhoemmen/pkg/include
CPPFLAGS = -I../bebop_util -I. 
##
## Version of CPPFLAGS if we are not using Valgrind.
##
#CPPFLAGS = -I../bebop_util -I. 
LIBS = -L../bebop_util -lbebop_util -lm
LINKER = gcc
RANLIB = echo
RM = rm -f

#
# Object files for the sparse matrix converter library.
#
LIBSMCOBJ = bcoo_matrix.o bcsr_matrix.o coo_matrix.o coord_elem.o \
	    csc_matrix.o csr_matrix.o iohb.o jad_matrix.o mmio.o \
	    read_mm.o sparse_matrix.o sparse_matrix_ops.o sparse_vector.o \
            csr_matmatmult.o interface.o
EXEOBJ = main.o  

LIBSMC_STATIC = libsparse_matrix_converter.a
# Linux-specific
#LIBSMC_DYNAMIC = libsparse_matrix_converter.so
# Darwin-specific
LIBSMC_DYNAMIC = libsparse_matrix_converter.dylib
EXEOUT = sparse_matrix_converter

#########################################################################
## OUTPUT RULES
#########################################################################

all: $(EXEOUT)

dynamic: $(LIBSMC_DYNAMIC)

$(EXEOUT): $(LIBSMC_STATIC) $(EXEOBJ)
	$(LINKER) $(CPPFLAGS) $(CFLAGS) -o $(EXEOUT) $(EXEOBJ) \
		$(LIBSMC_STATIC) $(LIBS)

$(LIBSMC_STATIC): $(LIBSMCOBJ) 
	$(ARCHIVER) $(ARFLAGS) $(LIBSMC_STATIC) $(LIBSMCOBJ) 
	$(RANLIB) $(LIBSMC_STATIC)

$(LIBSMC_DYNAMIC): $(LIBSMCOBJ)
	$(CC) -dynamiclib $(CFLAGS) $(LIBSMCOBJ) $(LIBS) -o $(LIBSMC_DYNAMIC)

smc_prompt: $(LIBSMC_STATIC) prompt.tab.o prompt.lexer.o
	$(LINKER) $(CPPFLAGS) $(CFLAGS) -o smc_prompt prompt.tab.o \
		prompt.lexer.o $(LIBSMC_STATIC) $(LIBS) -lfl

prompt.tab.o: prompt.tab.c sparse_matrix.h sparse_matrix_ops.h 
	$(CC) $(CPPFLAGS) $(CFLAGS) -c prompt.tab.c

prompt.lexer.o: prompt.lexer.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c prompt.lexer.c

prompt.lexer.c: prompt.lex
	flex prompt.lex
	mv lex.yy.c prompt.lexer.c

prompt.tab.c: prompt.y
	bison -d prompt.y

interface: interface_wrap.c 
	$(CC) $(CPPFLAGS) $(CFLAGS) -c interface.c interface_wrap.c -I/usr/apps/python/default/include/python2.3
	$(CC) -shared $(LIBS) interface.o interface_wrap.o $(LIBSMCOUT) -o interface.so

# FIXME: needs more dependencies (on includes)!
interface_wrap.c: interface.i $(LIBSMCOUT)
	swig -python interface.i


#########################################################################
## BUILD RULES 
#########################################################################

main.o: main.c sparse_matrix.h sparse_matrix_ops.h 
	$(CC) $(CPPFLAGS) $(CFLAGS) -c main.c

bcoo_matrix.o: bcoo_matrix.c bcoo_matrix.h bcsr_matrix.h 
	$(CC) $(CPPFLAGS) $(CFLAGS) -c bcoo_matrix.c

bcsr_matrix.o: bcsr_matrix.c bcsr_matrix.h 
	$(CC) $(CPPFLAGS) $(CFLAGS) -c bcsr_matrix.c

coo_matrix.o: coo_matrix.c coo_matrix.h csc_matrix.h read_mm.h  
	$(CC) $(CPPFLAGS) $(CFLAGS) -c coo_matrix.c

coord_elem.o: coord_elem.c coord_elem.h 
	$(CC) $(CPPFLAGS) $(CFLAGS) -c coord_elem.c

csc_matrix.o: csc_matrix.c csc_matrix.h coo_matrix.h iohb.h sparse_vector.h 
	$(CC) $(CPPFLAGS) $(CFLAGS) -c csc_matrix.c

csr_matrix.o: csr_matrix.c csr_matrix.h coo_matrix.h csc_matrix.h \
	iohb.h read_mm.h sparse_vector.h csr_matmatmult.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -c csr_matrix.c

csr_matmatmult.o: csr_matmatmult.c csr_matmatmult.h 
	$(CC) $(CPPFLAGS) $(CFLAGS) -c csr_matmatmult.c

interface.o: interface.c interface.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -c interface.c

iohb.o: iohb.c iohb.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -c iohb.c

jad_matrix.o: jad_matrix.c jad_matrix.h csr_matrix.h 
	$(CC) $(CPPFLAGS) $(CFLAGS) -c jad_matrix.c

mmio.o: mmio.c mmio.h 
	$(CC) $(CPPFLAGS) $(CFLAGS) -c mmio.c

read_mm.o: read_mm.c read_mm.h coo_matrix.h coord_elem.h csc_matrix.h \
	csr_matrix.h mmio.h 
	$(CC) $(CPPFLAGS) $(CFLAGS) -c read_mm.c

sparse_matrix.o: sparse_matrix.c sparse_matrix.h bcoo_matrix.h bcsr_matrix.h \
	coo_matrix.h csc_matrix.h csr_matrix.h jad_matrix.h 
	$(CC) $(CPPFLAGS) $(CFLAGS) -c sparse_matrix.c

sparse_matrix_ops.o: sparse_matrix_ops.c sparse_matrix_ops.h bcoo_matrix.h \
	bcsr_matrix.h coo_matrix.h csc_matrix.h csr_matrix.h jad_matrix.h \
	read_mm.h sparse_matrix.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -c sparse_matrix_ops.c

sparse_vector.o: sparse_vector.c sparse_vector.h 
	$(CC) $(CPPFLAGS) $(CFLAGS) -c sparse_vector.c

##########################################################################
## UTILITY STEPS
##########################################################################

# FIXME:  Depends on GNU tar features!!!
backup:
	tar --create --gzip --file=sparse_matrix_converter.tar.gz *.[ch] *.lex *.y Makefile TODO BUGS README Doxyfile COPYING 

ctags: 
	ctags *.[ch]

etags:
	etags *.[ch]

clean:
	$(RM) $(LIBSMCOBJ) 
	$(RM) $(LIBSMCOUT) 
	$(RM) $(EXEOBJ)
	$(RM) $(EXEOUT)
	$(RM) prompt.tab.c prompt.tab.h prompt.tab.o
	$(RM) prompt.lexer.c prompt.lexer.o
	$(RM) smc_prompt
	$(RM) interface_wrap.*
	$(RM) interface.so
