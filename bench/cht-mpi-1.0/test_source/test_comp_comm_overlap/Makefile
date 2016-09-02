TOPDIR=..
include $(TOPDIR)/Makefile.include
BUILD_AS_SINGLE_PROGRAM=0
ifeq ($(BUILD_AS_SINGLE_PROGRAM),0)
NPROCS = 1
else
NPROCS = 3
endif

# TEST OVERLAP
test_overlap: test_overlap_manager cht_worker
	
test_overlap_manager: test_overlap.cc test_overlap_task_implementations.cc ChunkObjDummyData.h ChunkObjDummyData.cc ChunkObjInputParams.h ChunkObjInputParams.cc $(CHTDIR)/libcht.a
	$(CC) $(CFLAGS) $(CHTINCL) -o $@ $^ 
cht_worker: test_overlap_task_implementations.cc ChunkObjDummyData.h ChunkObjDummyData.cc ChunkObjInputParams.h ChunkObjInputParams.cc $(CHTDIR)/libcht.a
	$(CC)  $(CFLAGS) $(CHTINCL) -o $@ $^ 
	
check: test_overlap
	mpirun -np $(NPROCS) ./test_overlap_manager > test_overlap_output.txt && grep -q 'test_comp_comm_overlap OK' test_overlap_output.txt || exit 1; 
	rm test_overlap_output.txt output.txt;
	@stars="***************************************************" ; \
	printf "\n$$stars\n test_overlap finished OK! \n$$stars\n\n"
	
clean:
	rm -f test_overlap_manager cht_worker
