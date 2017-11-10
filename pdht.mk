#
# makefile overlord for portals distributed hash table.
# author: d. brian larkins
# created: 3/20/16
#

PORTALS_INCLUDEDIR = ../opt/include
PORTALS_LIBDIR     = ../opt/lib
#PORTALS_INCLUDEDIR = /opt/hpctools/include
#PORTALS_LIBDIR     = /opt/hpctools/lib

#CC = clang
CC = gcc
MPICC = mpicc
#GCFLAGS = --std=c99 -g -O3 -D_POSIX_C_SOURCE=199309L -msse4.2 # development
#GCFLAGS = -std=c99 -g -D_POSIX_C_SOURCE=199309L
#GCFLAGS = -g -Wall
#GCFLAGS = -g -Wno-pointer-to-int-cast -Wno-int-to-pointer-cast
#GCFLAGS = --std=c99 -pg -O3 -D_POSIX_C_SOURCE=199309L -msse4.2 # profiling
GCFLAGS = --std=c99 -O3 -D_POSIX_C_SOURCE=199309L -msse4.2 # performance
#GCFLAGS = -O3
CFLAGS = $(GCFLAGS) -I. -I$(PDHT_TOP)/include -I$(PDHT_TOP)/$(PORTALS_INCLUDEDIR)
CFLAGSMPI = $(GCFLAGS) -I. -I$(PDHT_TOP)/includempi -I$(PDHT_TOP)/$(PORTALS_INCLUDEDIR)

#LDFLAGS=-L$(PORTALS_LIBDIR)
MATH_LIB            = -lm
RT_LIB              = -lrt
PMI_LIB             = -lpmi
PORTALS_LIB         = -lportals
PTHREAD_LIB         = -lpthread
PDHT_INSTALL_LIBDIR = $(PDHT_TOP)/lib
PDHT_LIBPDHT        = $(PDHT_INSTALL_LIBDIR)/libpdht.a
PDHT_LIBMPIPDHT     = $(PDHT_INSTALL_LIBDIR)/libmpipdht.a

PDHT_LIBDIRS = $(PDHT_TOP)/libpdht
PDHT_LIBMPIDIRS = $(PDHT_TOP)/libmpipdht

PDHT_LIBS = -L$(PDHT_TOP)/$(PORTALS_LIBDIR) -Wl,-rpath=$(PDHT_TOP)/$(PORTALS_LIBDIR) $(PDHT_LIBPDHT) $(PTHREAD_LIB) $(PMI_LIB) $(PORTALS_LIB) $(MATH_LIB)
PDHT_MPILIBS = $(PDHT_LIBMPIPDHT) $(MATH_LIB)

.PHONY: all

default: all


checkflags:
ifndef PORTALS_INCLUDEDIR
	@echo You must define PORTALS_INCLUDEDIR with the path to your Portals4 headers
	@false
endif
ifndef PORTALS_LIBDIR
	@echo You must define PORTALS_LIBDIR with the path to your Portals4 libraries
	@false
endif
ifndef PDHT_TOP
	@echo You must define PDHT_TOP with the path of the top-level directory
	@false
endif

.PHONY: pdhtlibs $(PDHT_LIBDIRS)
pdhtlibs: pdhtheaders $(PDHT_LIBDIRS)

.PHONY: pdhtlibs $(PDHT_LIBMPIDIRS)
pdhtmpilibs: pdhtheaders $(PDHT_LIBMPIDIRS)

.PHONY: checkflags pdhtheaders
pdhtheaders: 
	for dir in $(PDHT_LIBDIRS); do \
    $(MAKE) -C $$dir headers; \
  done; \
	for dir in $(PDHT_LIBMPIDIRS); do \
		$(MAKE) -C $$dir headers; \
	done

$(PDHT_LIBDIRS):
	$(MAKE) -C $@ GCFLAGS="$(GCFLAGS)" CC="$(CC)" PORTALS_INCLUDEDIR="$(PORTALS_INCLUDEDIR)" PORTALS_LIBDIR="$(PORTALS_LIBDIR)"

$(PDHT_LIBMPIDIRS):
	$(MAKE) -C $@ GCFLAGS="$(GCFLAGS)" CC="$(MPICC)" PORTALS_INCLUDEDIR="$(PORTALS_INCLUDEDIR)" PORTALS_LIBDIR="$(PORTALS_LIBDIR)"

clean:
	for dir in $(PDHT_LIBDIRS); do \
    $(MAKE) -C $$dir clean; \
  done;\
	for dir in $(PDHT_LIBMPIDIRS); do \
		$(MAKE) -C $$dir clean; \
	done;\
	rm -f *~ *.o gmon.out $(LOCAL_EXECS)
