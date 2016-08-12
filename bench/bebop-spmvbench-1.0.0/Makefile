# Makefile
#
# Overall makefile for standalone SpMV benchmark
##########################################################################

arch = UNKNOWN
include makes/Make.$(arch)

all: libbebop libconverter libgenerator libspmvbench hpccbench

libbebop:
	$(CD) bebop_util && make -f Make.standalone arch=$(arch)

libconverter:
	$(CD) sparse_matrix_converter && make -f Make.standalone arch=$(arch)

libgenerator:
	$(CD) matrix_generator && make -f Make.standalone arch=$(arch)

libspmvbench:
	$(CD) spmvbench && make -f Make.standalone arch=$(arch)

hpccbench:
	$(CD) hpcc_spmv_benchmark && make -f Make.standalone arch=$(arch)

clean:
	$(RM) bebop_util/*.o bebop_util/*.a
	$(RM) bebop_util/tests/*.o
	$(RM) sparse_matrix_converter/*.o sparse_matrix_converter/*.a
	$(RM) matrix_generator/*.o matrix_generator/*.a
	$(RM) spmvbench/*.o spmvbench/*.a
	$(RM) hpcc_spmv_benchmark/*.o hpcc_spmv_benchmark/hpcc_spmv_benchmark
