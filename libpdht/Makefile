#
# makefile overlord for portals distributed hash table.
# author: d. brian larkins
# created: 3/20/16
#

CFLAGS=$(GCFLAGS) -I../include -I$(PORTALS_INCLUDEDIR)

#LDFLAGS=-lportals -lpmi

RANLIB = /usr/bin/ranlib


OBJS =  checker.o \
        pmi.o     \
        util.o    \
        # line eater

.PHONY: all
all: libpdht.a

libpdht.a: headers $(OBJS)
	$(AR) rcs ../lib/libpdht.a $(OBJS)
	$(RANLIB) ../lib/libpdht.a

.PHONY: headers
headers:
	cp *.h ../include/.

.PHONY: clean
clean:
	rm -f *~ *.o ../include/pdht.h ../lib/libpdht.a
