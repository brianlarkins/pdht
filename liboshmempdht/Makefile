#pdht w/ oshmem w/ AM

CFLAGS=$(GCFLAGS) -I. -I../include
RANLIB = /usr/bin/ranlib

HDRS = pdht.h

OBJS = pdht.o \
			 #line eater

.PHONY: all
all: liboshmempdht.a

liboshmempdht.a : headers $(OBJS)
	mkdir -p ../lib
	$(AR) rcs ../lib/liboshmempdht.a $(OBJS)
	$(RANLIB) ../lib/liboshmempdht.a

$(OBJS): $(HDRS)

.PHONY: headers
headers:
	mkdir -p ../includeoshmem
	cp city.h uthash.h pdht.h ../includeoshmem/.

.PHONY: clean
clean:
	rm -f *~ *.o ../include/pdht.h ../lib/liboshmempdht.h
