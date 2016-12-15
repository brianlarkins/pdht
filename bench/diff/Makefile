# 
# makefile for transparent memory adaptor
# 
# id: $Id: Makefile 614 2007-03-06 21:28:26Z brian $
# 

### Berkeley UPC
#CC = upcc
#GCFLAGS = -network=vapi -D__BERKELEY_UPC__
#GCFLAGS = -pthreads -network=smp -D__BERKELEY_UPC__
#GCFLAGS = -network=mpi -D__BERKELEY_UPC__

#GCFLAGS = -Wall -Wp,"-g" -pthreads
#GCFLAGS = -g -T=2 -pthreads

### Intrepid UPC
CC = mpicc
#GCFLAGS = -g -Wall
#GCFLAGS = -03
#GCFLAGS = -pg -03
GCFLAGS = -O0 -std=c99 -g -Wall -Wno-pointer-to-int-cast -Wno-int-to-pointer-cast

CFLAGS = $(GCFLAGS) -I../include

LIBDIRS = ../libgcl  ../libgt ../libtc
INCDIRS = $(LIBDIRS)

CLEANLIBS = $(LIBDIRS)

OBJS =  init.o           \
	math.o           \
	tensor.o         \
	tree.o	         \
        # line eater

MATH_LIB    = -lm
LIBDIR      = ../lib
LIBGCL      = $(LIBDIR)/libgcl.a
LIBTC       = $(LIBDIR)/libtc.a
LIBTREE     = $(LIBDIR)/libgt.a
LIBSYNCH    = $(LIBDIR)/libsynch.a
LIBARMCI    = $(LIBDIR)/libarmci.a
LIBS        = $(LIBTREE) $(LIBGCL) $(LIBTC) $(LIBSYNCH) $(LIBARMCI) $(MATH_LIB)


.PHONY: all
all: libs diff3d

.PHONY: libs $(LIBDIRS)
libs: headers $(LIBDIRS)

.PHONY: headers
headers: 
	cp *.h ../include/.
	for dir in $(LIBDIRS); do \
	  $(MAKE) -C $$dir headers; \
	done


$(LIBDIRS):
	$(MAKE) -C $@ GCFLAGS="$(GCFLAGS)" CC="$(CC)"

diff3d: $(OBJS) diff3d.o
	$(CC) $(CFLAGS) -o diff diff3d.o $(OBJS) $(LIBS)

seqdiff: $(OBJS) seqdiff.c
	$(CC) $(CFLAGS) -o seqdiff seqdiff.c $(OBJS)

test-chunk: libs test-chunk.o
	$(CC) $(CFLAGS) -o test-chunk test-chunk.o $(LIBS)

test-task: libs test-task.o
	$(CC) $(CFLAGS) -o test-task test-task.o $(LIBS)

ring_test: ring_test.c ring.c collection.c timer.c
	$(CC) $(CFLAGS) -o ring_test ring_test.c ring.c collection.c timer.c

ring_n_test: ring_n_test.c ring.c collection.c timer.c
	$(CC) $(CFLAGS) -o ring_n_test ring_n_test.c ring.c collection.c timer.c

clean: 
	for dir in $(LIBDIRS); do \
	  $(MAKE) -C $$dir clean; \
	done
	rm -f *~ *.o gmon.out diff task test-task ring_test ring_n_test test-chunk

