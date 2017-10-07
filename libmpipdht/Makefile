
#pdht w/ mpi make file
#its shitty cause jack wrote it




CFLAGS=$(GCFLAGS) -I. -I../include -I../$(PORTALS_INCLUDEDIR)

RANLIB = /usr/bin/ranlib

HDRS = pdht.h

OBJS = city.o	\
       commsynch.o \
       hash.o \
       init.o \
       putget.o \
       util.o \
       # line eater?

.PHONY: all
all: libmpipdht.a

libmpipdht.a : headers $(OBJS)
	mkdir -p ../lib
	$(AR) rcs ../lib/libmpipdht.a $(OBJS)
	$(RANLIB) ../lib/libmpipdht.a

$(OBJS): $(HDRS)

.PHONY: headers
headers:
	mkdir -p ../includempi
	cp city.h citycrc.h uthash.h pdht.h ../includempi/.

.PHONY: clean
clean:
	rm -f *~ *.o ../include/pdht.h ../lib/libmpipdht.h




