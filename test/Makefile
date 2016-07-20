#
# makefile overlord for portals distributed hash table.
# author: d. brian larkins
# created: 3/20/16
#

PDHT_TOP=..

include $(PDHT_TOP)/pdht.mk

.PHONY: all
all: trigmult


barrier: pdhtlibs barrier.c
	$(CC) $(CFLAGS) -o barrier barrier.c $(PDHT_LIBS)

simple: pdhtlibs simple.c
	$(CC) $(CFLAGS) -o simple simple.c $(PDHT_LIBS)

notfound: pdhtlibs notfound.c
	$(CC) $(CFLAGS) -o notfound notfound.c $(PDHT_LIBS)

collision: pdhtlibs collision.c
	$(CC) $(CFLAGS) -o collision collision.c $(PDHT_LIBS)

trig: pdhtlibs trig.c
	$(CC) $(CFLAGS) -o trig trig.c $(PDHT_LIBS)

trigmult: pdhtlibs trigmult.c
	$(CC) $(CFLAGS) -o trigmult trigmult.c $(PDHT_LIBS)
