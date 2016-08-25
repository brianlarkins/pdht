TOPDIR=..
include $(TOPDIR)/Makefile.include
BUILD_AS_SINGLE_PROGRAM=0
ifeq ($(BUILD_AS_SINGLE_PROGRAM),0)
NPROCS = 1
else
NPROCS = 6
endif

# TEST SERVICE
test_service: test_service_manager cht_worker
	
test_service_manager: test_service.cc test_service_chunk_obj_implementations.cc test_tasksched_task_implementations_ABC.cc $(CHTDIR)/libcht.a
	$(CC) $(CFLAGS) $(CHTINCL) -o $@ $^
cht_worker: test_service_chunk_obj_implementations.cc test_tasksched_task_implementations_ABC.cc $(CHTDIR)/libcht.a
	$(CC) $(CFLAGS) $(CHTINCL) -o $@ $^

check: test_service
	mpirun -np $(NPROCS) ./test_service_manager > test_service_output.txt && grep -q 'test_service: DONE' test_service_output.txt || exit 1; 
	rm test_service_output.txt output.txt;
	@stars="***************************************************" ; \
	printf "\n$$stars\n test_service finished OK! \n$$stars\n\n"

clean:
	rm -f test_service_manager cht_worker
