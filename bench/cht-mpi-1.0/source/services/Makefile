BUILD_AS_SINGLE_PROGRAM=0
CC=mpiCC
CC=mpic++

CFLAGS= -O3

INCLUDES=-I..

MAKEDEPEND= $(CC) -M -DBUILD_AS_SINGLE_PROGRAM=$(BUILD_AS_SINGLE_PROGRAM) $(CFLAGS) $(INCLUDES) -o $*.d $<

VPATH = ./output_service ./output_test_service/ ./chunk_service ./chunk_obj_service ./chunk_test_service ./task_scheduler_service ./service_manager 

SRCS = services_utils.cc \
	determine_process_roles.cc \
	service_worker_main.cc \
	Service_worker.cc \
	OutputService_worker.cc \
	OutputTestService_worker.cc \
	ChunkObjService_worker.cc \
	ChunkService_worker.cc \
	ChunkTestService_worker.cc \
	TaskSchedulerService_worker.cc \
	service_worker_func.cc \
	services_registration_worker.cc \
	Task.cc \
	ChunkObject.cc \
	spawn_workers_func.cc \
	ParentInterface.cc \
	chunks_and_tasks_parent.cc \
	ServiceManager.cc \
	Service_parent.cc \
	OutputService_parent.cc \
	OutputTestService_parent.cc \
	ChunkObjService_parent.cc \
	ChunkService_parent.cc \
	ChunkTestService_parent.cc \
	TaskSchedulerService_parent.cc \
	services_registration_parent.cc \
	chunk_type_macros.cc	\
	task_type_macros.cc

OBJS = $(SRCS:.cc=.o) 

all: libs

libs: libservice.a 

libservice.a: $(OBJS) 
	ar rcs libservice.a $(OBJS) 

%.o: %.cc
	@$(MAKEDEPEND)
	$(CC) -DBUILD_AS_SINGLE_PROGRAM=$(BUILD_AS_SINGLE_PROGRAM) $(CFLAGS) $(INCLUDES) -c $< -o $@

-include $(SRCS:.cc=.d)

clean:
	rm -f *.a *.o *.d

check:
	
