#
# makefile overlord for portals distributed hash table.
# author: d. brian larkins
# created: 3/20/16
#

PORTALS_INCLUDEDIR=/opt/hpctools/include

CFLAGS=$(GCFLAGS) -I. -I../include -I$(PORTALS_INCLUDEDIR)

#LDFLAGS=-lportals -lpmi

RANLIB = /usr/bin/ranlib


OBJS =  assoc.o      \
        commsynch.o  \
        hash.o       \
        init.o       \
        nbputget.o   \
        putget.o     \
        pmi.o        \
        util.o       \
        # line eater

.PHONY: all
all: libpdht.a

libpdht.a: headers $(OBJS)
	$(AR) rcs ../lib/libpdht.a $(OBJS)
	$(RANLIB) ../lib/libpdht.a

.PHONY: headers
headers:
	cp pdht_inline.h pdht.h ../include/.

.PHONY: clean
clean:
	rm -f *~ *.o ../include/pdht.h ../lib/libpdht.a
