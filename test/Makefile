#
# makefile overlord for portals distributed hash table.
# author: d. brian larkins
# created: 3/20/16
#

PDHT_TOP=..

include $(PDHT_TOP)/pdht.mk

.PHONY: all
all: scaling


barrier: pdhtlibs barrier.c
	$(CC) $(CFLAGS) -o barrier barrier.c $(PDHT_LIBS)

collision: pdhtlibs collision.c
	$(CC) $(CFLAGS) -o collision collision.c $(PDHT_LIBS)

csort: pdhtlibs csort.c
	$(CC) $(CFLAGS) -o csort csort.c $(PDHT_LIBS)

latency: pdhtlibs latency.c
	$(CC) $(CFLAGS) -o latency latency.c $(PDHT_LIBS)

matchlength: pdhtlibs matchlength.c
	$(CC) $(CFLAGS) -o matchlength matchlength.c $(PDHT_LIBS)

notfound: pdhtlibs notfound.c
	$(CC) $(CFLAGS) -o notfound notfound.c $(PDHT_LIBS)

optimes: pdhtlibs optimes.c
	$(CC) $(CFLAGS) -o optimes optimes.c $(PDHT_LIBS)

oshbench: pdhtlibs oshbench.c
	$(CC) $(CFLAGS) -o oshbench oshbench.c $(PDHT_LIBS)

scaling: pdhtlibs scaling.c
	$(CC) $(CFLAGS) -o scaling scaling.c $(PDHT_LIBS)

simple: pdhtlibs simple.c
	$(CC) $(CFLAGS) -o simple simple.c $(PDHT_LIBS)

trig: pdhtlibs trig.c
	$(CC) $(CFLAGS) -o trig trig.c $(PDHT_LIBS)

trigmult: pdhtlibs trigmult.c
	$(CC) $(CFLAGS) -o trigmult trigmult.c $(PDHT_LIBS)
