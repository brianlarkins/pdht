#
# makefile overlord for portals distributed hash table.
# author: d. brian larkins
# created: 3/20/16
#

PDHT_TOP=../..



include $(PDHT_TOP)/pdht.mk



diffMPI: CC = $(MPICC)
diffMPI: CFLAGS = $(CFLAGSMPI) -DMPI
diffSHMEM: CC = $(OSHCC)
diffSHMEM: CFLAGS = $(CFLAGSOSHMEM)

OBJS =  alloc.o          \
				init.o           \
				math.o           \
				tensor.o         \
				tree.o           \
        # line eater
.PHONY: all
all: diff

diff: pdhtlibs $(OBJS) diff3d.o 
	$(CC) $(CFLAGS) -lrt -o diff diff3d.o $(OBJS) $(PDHT_LIBS)

diffSHMEM: pdhtoshmemlibs $(OBJS) diff3d.o 
	$(OSHCC) $(CFLAGSOSHMEM) -o diffSHMEM diff3d.o $(OBJS) $(PDHT_OSHMEMLIBS)

diffMPI: pdhtmpilibs $(OBJS) diff3d.o
	$(MPICC) $(CFLAGSMPI) -DMPI=1 -o diffMPI diff3d.o $(OBJS) $(PDHT_MPILIBS)
