CC=gcc
CFLAGS=--std=c99 -g -I/opt/hpctools/include
LDFLAGS=-L/opt/hpctools/lib

OBJS =  checker.o pmi.o util.o
      
checker: $(OBJS)
	$(CC) $(LDFLAGS) -o checker $(OBJS) -lportals -lpmi

checker.o: checker.c
	$(CC) -c $(CFLAGS) checker.c
