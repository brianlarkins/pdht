#
# makefile overlord for portals distributed hash table.
# author: d. brian larkins
# created: 3/20/16
#

PDHT_TOP=..

include $(PDHT_TOP)/pdht.mk

.PHONY: all
all: simple


barrier: pdhtlibs barrier.c
	$(CC) $(CFLAGS) -o barrier barrier.c $(PDHT_LIBS)

simple: pdhtlibs simple.c
	$(CC) $(CFLAGS) -o simple simple.c $(PDHT_LIBS)
