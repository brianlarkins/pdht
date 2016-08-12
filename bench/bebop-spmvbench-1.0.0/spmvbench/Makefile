# Makefile
# Author: mfh
# Last modified: 2004 Nov 28
#
# For register-blocked CSR sparse matrix-vector multiplication component
# of HPCC benchmark.  This particular Makefile was adapted from the HPCC
# STREAM Makefile.
##########################################################################

arch = UNKNOWN

CPPFLAGS = -DSPMV_STANDALONE -I../matrix_generator -I../bebop_util -I../sparse_matrix_converter

include ./config/Make.$(arch)

#include ../hpl/Make.$(arch)
#CPPFLAGS = -I. -I../hpl/include -I../include $(HPL_DEFS)

##
## HPCC uses the nonstandard "CCFLAGS" for the flags to pass to the C compiler.
## It should be CFLAGS instead.
##
#CFLAGS = -DDEBUG $(CCFLAGS) 
CFLAGS = $(CCFLAGS) 
LIBS = -L../matrix_generator -l_matrix_generators -L../bebop_util -lbebop_util \
	-L../sparse_matrix_converter -lsparse_matrix_converter -lm
THELIB = lib_bsr_spmv.a

BENCHMARK_OBJ = avltree.o block_smvm_code.o fill_with_random_block.o \
	smvm_benchmark.o smvm_get_output_file.o smvm_timing_results.o \
	smvm_timing_run.o smvm_verify_result.o smvm_benchmark_results.o

INVOKER_OBJ = hpcc.o smvm_memory_hierarchy_traits.o smvm_invoker.o \
	smvm_hpcc_interface.o block_size_parser.o

TIMING_OBJ = timing.o
TIMING_HEADER = timing.h
CONFIG_H      = config.h

LIBOBJ = $(INVOKER_OBJ) $(BENCHMARK_OBJ) $(TIMING_OBJ)


#########################################################################
# OUTPUT RULES
#########################################################################

all: $(THELIB) lowlevel_tester hpcc_interface_tester

$(THELIB): $(LIBOBJ)
	$(ARCHIVER) $(ARFLAGS) $(THELIB) $(LIBOBJ)
	$(RANLIB) $(THELIB)

lowlevel_tester: $(THELIB) benchmark_tester.o
	$(LINKER) $(CPPFLAGS) $(MPinc) $(CFLAGS) -o lowlevel_tester \
		benchmark_tester.o $(THELIB) $(MPlib) $(LIBS)

hpcc_interface_tester: $(THELIB) smvm_hpcc_interface_tester.o
	$(LINKER) $(CPPFLAGS) $(MPinc) $(CFLAGS) -lm -o hpcc_interface_tester \
		hpcc.o smvm_hpcc_interface_tester.o $(THELIB) $(LIBS)


	
#########################################################################
# BUILD RULES 
#########################################################################

### Invoker files ###

block_size_parser.o: block_size_parser.c block_size_parser.h $(CONFIG_H)
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

hpcc.o: hpcc.c hpcc.h smvm_benchmark_results.h $(CONFIG_H)
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

smvm_hpcc_interface.o: smvm_hpcc_interface.c smvm_hpcc_interface.h \
		hpcc.h smvm_benchmark_results.h smvm_invoker.h $(CONFIG_H)
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

smvm_hpcc_interface_tester.o: smvm_hpcc_interface_tester.c \
		smvm_hpcc_interface.h hpcc.h $(CONFIG_H)
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

smvm_invoker.o: smvm_invoker.c smvm_invoker.h avltree.h block_size_parser.h \
		smvm_memory_hierarchy_traits.h smvm_timing_results.h \
		smvm_benchmark.h $(CONFIG_H)
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

smvm_memory_hierarchy_traits.o: smvm_memory_hierarchy_traits.c \
		smvm_memory_hierarchy_traits.h $(CONFIG_H)
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<



####################
### Timing files ###
####################

timing.o: timing.c timing.h

#######################
### Benchmark files ###
#######################

avltree.o: avltree.c avltree.h $(CONFIG_H)
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

benchmark_tester.o: benchmark_tester.c smvm_benchmark.h $(CONFIG_H)
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

block_smvm_code.o: block_smvm_code.c block_smvm_code.h $(CONFIG_H)
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

fill_with_random_block.o: fill_with_random_block.c fill_with_random_block.h \
		$(CONFIG_H)
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

smvm_benchmark.o: smvm_benchmark.c smvm_benchmark.h $(TIMING_HEADER) \
		block_smvm_code.h $(CONFIG_H)
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

smvm_benchmark_results.o: smvm_benchmark_results.c smvm_benchmark_results.h \
		$(CONFIG_H)
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

smvm_get_output_file.o: smvm_get_output_file.c smvm_get_output_file.h \
		$(CONFIG_H)
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

smvm_timing_results.o: smvm_timing_results.c smvm_timing_results.h \
		$(CONFIG_H)
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

smvm_timing_run.o: smvm_timing_run.c smvm_timing_run.h smvm_get_output_file.h \
		smvm_timing_results.h $(TIMING_HEADER) $(CONFIG_H)
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

smvm_verify_result.o: smvm_verify_result.c smvm_verify_result.h $(CONFIG_H)
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

##########################################################################
# UTILITY STEPS
##########################################################################

etags:
	etags *.[ch]

clean:
	$(RM) $(LIBOBJ) $(THELIB) *~ lowlevel_tester hpcc_interface_tester 
	$(RM) benchmark_tester.o smvm_hpcc_interface_tester.o hpcc.o 
	$(RM) get_options.o
