#
# makefile overlord for portals distributed hash table.
# author: d. brian larkins
# created: 3/20/16
#

PDHT_TOP=../..

include $(PDHT_TOP)/pdht.mk

OBJS =  init.o           \
				math.o           \
				tensor.o         \
				tree.o           \
        # line eater

.PHONY: all
all: diff

diff: pdhtlibs $(OBJS) diff3d.o 
	$(CC) $(CFLAGS) -o diff diff3d.o $(OBJS) $(PDHT_LIBS)
