CC=mpiCC
CC=mpic++

CFLAGS= -O3

MAKEDEPEND= $(CC) -M $(CFLAGS) $(INCLUDES) -o $*.d $<

SRCS = MPIWrapperMultiple.cc MPIWrapperSerialized.cc MPIWrapperFunneled.cc 

OBJS= $(SRCS:.cc=.o) 

all: libs test_mpiwrapper

libs: libmpiwrapper.a

libmpiwrapper.a: $(OBJS)
	ar rcs libmpiwrapper.a $(OBJS)

test_mpiwrapper: test_mpiwrapper.cc libmpiwrapper.a
	$(CC) $(CFLAGS) -o $@ $^  

%.o: %.cc
	@$(MAKEDEPEND)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

-include $(SRCS:.cc=.d)

clean:
	rm -f *.a *.o *.d test_mpiwrapper
check:
	
