#include <pdht.h>

extern pdht_context_t *c;

void pdht_barrier(void){
  MPI_Barrier(c->barrier_comm);
}

void pdht_fence(pdht_t *dht){
  MPI_Barrier(c->barrier_comm);
}
