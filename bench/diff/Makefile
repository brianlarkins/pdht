#
# makefile overlord for portals distributed hash table.
# author: d. brian larkins
# created: 3/20/16
#

PDHT_TOP=../..



include $(PDHT_TOP)/pdht.mk



diffMPI: CC = $(MPICC)
diffMPI: CFLAGS = $(CFLAGSMPI)

OBJS =  alloc.o          \
				init.o           \
				math.o           \
				tensor.o         \
				tree.o           \
        # line eater
.PHONY: all
all: diff

diff: pdhtlibs $(OBJS) diff3d.o 
	$(CC) $(CFLAGS) -o diff diff3d.o $(OBJS) $(PDHT_LIBS)

diffMPI: pdhtmpilibs $(OBJS) diff3d.o
	$(MPICC) $(CFLAGSMPI) -o diffMPI diff3d.o $(OBJS) $(PDHT_MPILIBS)
