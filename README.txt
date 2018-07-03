One day this will tell a nice story.

For now, we put all of the DHT library code in libpdht and any testing drivers
in test.

To run openSHMEM version, you need the shmem-am of snyjm-18 fork of 
osss-ucx. To build that you need the ucx-shmem-dev branch of snyjm-18
fork of UCX. See those two README 's for further instructions ons building and
running those

To build:

First need to build Portals 4:
  $ cd cfg
  $ run configure for the type of system you are on
  $ cd ../portals4
  $ make -j 20 
  $ make check # should mostly pass (i had 10 faults out of 785 tests)
  $ make install

Then need to build/run a PDHT program:
  $ cd test
  $ make clean
  $ make scaling # or whatever

To run on senna (2 Nodes, 4 processor cores):
  $ srun -N 2 -n 4 ./scaling 
