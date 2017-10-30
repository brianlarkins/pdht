#include <pdht.h>

#define PDHT_COUNTER_HOLDER 0

extern pdht_context_t *c;

void pdht_barrier(void){
  MPI_Barrier(c->barrier_comm);
}

void pdht_fence(pdht_t *dht){
  MPI_Barrier(c->barrier_comm);
}



int pdht_counter_init(pdht_t *ht, int initval){

  assert(ht->counter_count < PDHT_MAX_COUNTERS);

  if(c->rank == PDHT_COUNTER_HOLDER){
    ht->counters[ht->counter_count] = initval;
  }
  return ht->counter_count++;
}

void pdht_counter_reset(pdht_t *ht, int counter){
  
  char sendBuf[sizeof(message_t) + sizeof(int)];
  
  message_t *msg_reset = (message_t *)sendBuf;
  msg_reset->type = pdhtCounterReset;
  msg_reset->rank = c->rank;
  
  int i;
  for(i = 0; i < c->size; i++){
    if(c->hts[i] == ht) break;
  }
  

  msg_reset->ht_index = i;
  msg_reset->mbits = 0;
  memcpy(&msg_reset->key, &counter, sizeof(int));
  
  MPI_Ssend(msg_reset, sizeof(sendBuf), MPI_CHAR, PDHT_COUNTER_HOLDER, PDHT_TAG_COMMAND, MPI_COMM_WORLD);
  

}

uint64_t pdht_counter_inc(pdht_t *ht, int counter, uint64_t val){
  //preparing message
  char sendBuf[sizeof(message_t) + sizeof(int)];
  MPI_Status status;

  message_t *inc_message;
  
  inc_message = (message_t *)sendBuf;
  inc_message->type = pdhtCounterInc;
  inc_message->rank = c->rank;
  

  int i;
  for(i = 0; i < c->size; i++){
    if(c->hts[i] == ht) break;
  }
  
  inc_message->ht_index = i;
  inc_message->mbits = val; //yea i get this isn't the best but it'll work without wasting that space
  
  memcpy(inc_message->key, &counter, sizeof(int));

  MPI_Send(inc_message, sizeof(sendBuf), MPI_CHAR, PDHT_COUNTER_HOLDER, PDHT_TAG_COMMAND, MPI_COMM_WORLD);
  
  int counter_val;

  MPI_Recv(&counter, sizeof(int), MPI_INT, PDHT_COUNTER_HOLDER, PDHT_COUNTER_REPLY, MPI_COMM_WORLD, &status);
  
  return counter;
}


