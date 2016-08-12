# Makefile
# Author: mfh
# Since: 14 June 2005
# 
# Makefile for matrix generator library, and tester executable.
##########################################################################

# Platform-specific configuration options.
ARCHIVER = ar
ARFLAGS = r
CC ?= gcc
CFLAGS = -O3 -Wall -std=c99 
#CFLAGS = -pg -g -O2 -Wall -std=c99 
#CFLAGS = -g -Wall -std=c99 
#CFLAGS = -g -O2 -Wall -std=c99 
CPPFLAGS = -I../sparse_matrix_converter -I../bebop_util -I.
LINKER = gcc
LIBS = -lm -L../sparse_matrix_converter -lsparse_matrix_converter -L../bebop_util -lbebop_util
RANLIB = echo
RM = rm -f

THELIB = lib_matrix_generators.a
TESTEREXE = matrix_generator_tester
LIBOBJ = create_rand.o stencil.o
TESTEROBJ = main.o dict.o


#########################################################################
## OUTPUT RULES
#########################################################################

all: $(THELIB) $(TESTEREXE)

$(THELIB): $(LIBOBJ)
	$(ARCHIVER) $(ARFLAGS) $(THELIB) $(LIBOBJ) 
	$(RANLIB) $(THELIB)

#	$(LINKER) $(CPPFLAGS) $(CFLAGS) -fPIC -shared -o $(THELIB) $(LIBOBJ)

$(TESTEREXE): $(LIBOBJ) $(TESTEROBJ)
	$(LINKER) $(CPPFLAGS) $(CFLAGS) -o $(TESTEREXE) $(TESTEROBJ) $(THELIB) $(LIBS)


#	$(LINKER) $(CPPFLAGS) $(CFLAGS) -rdynamic -o $(TESTEREXE) $(TESTEROBJ) -ldl


#########################################################################
## BUILD RULES 
#########################################################################

### Tester executable object files ###

main.o: main.c create_rand.h stencil.h random_number.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

dynamic_library.o: dynamic_library.c dynamic_library.h 
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

dict.o: dict.c dict.h
	$(CC) $(CPPFLAGS) -DNDEBUG=1 $(CFLAGS) -o $@ -c $<

### Library object files ###

create_rand.o: create_rand.c create_rand.h 
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

stencil.o: stencil.c stencil.h 
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

##########################################################################
## UTILITY STEPS
##########################################################################

# FIXME:  Depends on GNU tar features!!!
backup:
	tar --create --gzip --file=matrix_generator.tar.gz *.[ch] Makefile Doxyfile BUGS COPYING

ctags: 
	ctags *.[ch]

etags:
	etags *.[ch]

clean:
	$(RM) $(LIBOBJ) 
	$(RM) $(TESTEROBJ) 
	$(RM) $(THELIB) 
	$(RM) $(TESTEREXE)

