PROG       = hpcc_spmv_benchmark
INSTALLDIR = ./bin

CC         = gcc
LIBDIR     = 
INCLUDE    = 
SRC        = benchmark_subset.c hpcc_spmv_benchmark.c main.c
# NOTE:  the library loading order is backwards from what you expect!  
LIBS       = -lm -l_bsr_spmv -l_matrix_generators -lsparse_matrix_converter -lbebop_util 
CFLAGS     = -O3 -Wall -std=c99 -I../bebop_util -I../sparse_matrix_converter -I../matrix_generator -I../spmvbench
#CFLAGS     = -O2 -Wall -std=c99 -I../bebop_util -I../sparse_matrix_converter -I../matrix_generator -I../spmvbench
#CFLAGS     = -g -Wall -std=c99 -I../bebop_util -I../sparse_matrix_converter -I../matrix_generator -I../spmvbench

# Setup the loader.
LD         = /usr/bin/ld
LDFLAGS	   = -L../bebop_util -L../sparse_matrix_converter -L../matrix_generator -L../spmvbench

# Object Modules
OBJ = $(SRC:.c=.o)

# Make all Programs
all: $(PROG) 

$(PROG): $(OBJ)
	$(CC) $(LDFLAGS) $(CFLAGS) -o $(PROG) $(OBJ) $(LIBS)

# These guys can be activated by:
# cmd-line> make clean
# cmd-line> make install
clean:
	-rm	-f	$(PROG)	*.o 
install:
	mkdirhier $(INSTALLDIR)
	-cp	$(PROG) 	$(INSTALLDIR)/.
#
# DO NOT DELETE THIS LINE -- make depend depends on it.
