#
# makefile overlord for portals distributed hash table.
# author: d. brian larkins
# created: 3/20/16
#

PORTALS_INCLUDEDIR = /opt/hpctools/include
PORTALS_LIBDIR     = /opt/hpctools/lib

CC = clang
GCFLAGS = --std=c99 -g
#GCFLAGS = -g -Wall
#GCFLAGS = -g -Wno-pointer-to-int-cast -Wno-int-to-pointer-cast
#GCFLAGS = -pg -O3
#GCFLAGS = -O3
CFLAGS = $(GCFLAGS) -I$(PDHT_TOP)/include -I$(PORTALS_INCLUDE)

#LDFLAGS=-L$(PORTALS_LIBDIR)
MATH_LIB            = -lm
PDHT_INSTALL_LIBDIR = $(PDHT_TOP)/lib
PDHT_LIBPDHT        = $(PDHT_TOP)/libpdht.a

PDHT_LIBS = $(PDHT_LIBPDHT)

.PHONY: all
checkflags:
ifndef PORTALS_INCLUDEDIR
  @echo You must define PORTALS_INCLUDEDIR with the path to your Portals4 headers
  @false
endif
ifndef PORTALS_LIBIDR
  @echo You must define PORTALS_LIBDIR with the path to your Portals4 libraries
  @false
endif
ifndef PDHT_TOP
  @echo You must define PDHT_TOP with the path of the top-level directory
  @false
endif

.PHONY: pdhtlibs $(PDHT_LIBDIRS)
pdhtlibs: pdhtheaders $(PDHT_LIBDIRS)

.PHONY: checkflags pdhtheaders
pdhtheaders: 
	for dir in $(PDHT_LIBDIRS); do \
    $(MAKE) -C $$dir headers; \
  done

$(PDHT_LIBDIRS):
	$(MAKE) -C $@ GCFLAGS="$(GCFLAGS)" CC="$(CC)"

clean:
	for dir in $(PDHT_LIBDIRS); do \
    $(MAKE) -C $$dir clean; \
  done
	rm -f *~ *.o gmon.out $(LOCAL_EXECS)
