TOPDIR=..
include $(TOPDIR)/Makefile.include
BUILD_AS_SINGLE_PROGRAM=0
ifeq ($(BUILD_AS_SINGLE_PROGRAM),0)
NPROCS = 1
else
NPROCS = 4
endif

# List all object files here (except the one for the main program)
WRK_OBJS = CInt.o CDouble.o CMatrix.o CreateMatrix.o MatrixAdd.o MatrixMultiply.o CreateMatrixFromIds.o GetMatrixElement.o

# List all header files here
HEADER_FILES = CDouble.h CInt.h CMatrix.h CreateMatrixFromIds.h CreateMatrix.h GetMatrixElement.h MatrixAdd.h MatrixElementValues.h MatrixMultiply.h

test_matrix: test_matrix_manager cht_worker
	

test_matrix_manager: test_matrix.o $(WRK_OBJS) $(CHTDIR)/libcht.a
	$(CC) $(CFLAGS) $(CHTINCL) -o $@ $^ 
cht_worker: $(WRK_OBJS) $(CHTDIR)/libcht.a
	$(CC) $(CFLAGS) $(CHTINCL) -o $@ $^

%.o: %.cc $(HEADER_FILES)
	$(CC) $(CFLAGS) $(CHTINCL) -c $< -o $@

check: test_matrix
	#srun -N $(NPROCS) ./test_matrix_manager 800 3 > test_matrix_output.txt && grep -q 'test_matrix finished OK' test_matrix_output.txt || exit 1; 
	mpirun -np $(NPROCS) ./test_matrix_manager 800 3 > test_matrix_output.txt && grep -q 'test_matrix finished OK' test_matrix_output.txt || exit 1; 
	rm test_matrix_output.txt output.txt;
	@stars="***************************************************" ; \
	printf "\n$$stars\n test_matrix finished OK! \n$$stars\n\n"

clean:
	rm -f *.o test_matrix_manager cht_worker
