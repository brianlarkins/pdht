CC=mpiCC

CFLAGS= -O3

INCLUDES=-I..

MAKEDEPEND= $(CC) -M $(CFLAGS) $(INCLUDES) -o $*.d $<

VPATH = ./threads

SRCS = cht_threads.cc Cond.cc Manager.cc Mutex.cc Thread.cc ThreadGroup.cc

OBJS= $(SRCS:.cc=.o)

all: libs

libs: libutils.a

libutils.a: $(OBJS)
	ar rcs libutils.a $(OBJS)

%.o: %.cc
	@$(MAKEDEPEND)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

-include $(SRCS:.cc=.d)

check:
	
clean:
	rm -f *.a *.o *.d

