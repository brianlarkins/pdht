#include <pdht.h>


void pdht_barrier(void){
  MPI_Barrier(MPI_COMM_WORLD);

}

void pdht_fence(pdht_t *dht){
  MPI_Barrier(MPI_COMM_WORLD);

}
