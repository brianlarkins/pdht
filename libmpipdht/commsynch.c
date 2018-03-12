#include <pdht.h>
  
#define PDHT_COUNTER_HOLDER 1

extern pdht_context_t *c;

void pdht_barrier(void){
  MPI_Barrier(c->split_comm);
}



void pdht_fence(pdht_t *dht){
  MPI_Barrier(c->split_comm);
}



int pdht_counter_init(pdht_t *ht, uint64_t initval){
  char sendBuf[sizeof(message_t) + sizeof(uint64_t)];
  message_t *msg_init = (message_t *)sendBuf;
  assert(ht->counter_count < PDHT_MAX_COUNTERS);

  msg_init->type = pdhtCounterInit;
  memcpy(&msg_init->key, &initval, sizeof(uint64_t));
  msg_init->rank = c->rank;
  int i;
  for(i = 0; i < c->dhtcount; i++){
    if(c->hts[i] == ht) break;
  }
  
  msg_init->ht_index = i;
  if(c->rank == PDHT_COUNTER_HOLDER){
    ht->counters[ht->counter_count] = initval;
  }
  MPI_Ssend(msg_init, sizeof(message_t) + sizeof(uint64_t), MPI_CHAR, PDHT_COUNTER_HOLDER, PDHT_TAG_COMMAND, MPI_COMM_WORLD);
  
  return ht->counter_count++;
}



void pdht_counter_reset(pdht_t *ht, int counter){
  char sendBuf[sizeof(message_t) + sizeof(int)];
  
  if(c->rank == 0){
    message_t *msg_reset = (message_t *)sendBuf;
    msg_reset->type = pdhtCounterReset;
    msg_reset->rank = c->rank;
    
    int i;
    for(i = 0; i < c->dhtcount; i++){
      if(c->hts[i] == ht) break;
    }
    

    msg_reset->ht_index = i;
    msg_reset->mbits = 0;
    memcpy(&msg_reset->key, &counter, sizeof(int));
    
    MPI_Ssend(msg_reset, sizeof(sendBuf), MPI_CHAR, PDHT_COUNTER_HOLDER, PDHT_TAG_COMMAND, MPI_COMM_WORLD);
  }
  pdht_barrier();
}



uint64_t pdht_counter_inc(pdht_t *ht, int counter, uint64_t val){
  //preparing message
  char sendBuf[sizeof(message_t) + sizeof(int)];
  long unsigned counter_val;
  MPI_Status status;

  message_t *inc_message;
  
  inc_message = (message_t *)sendBuf;
  inc_message->type = pdhtCounterInc;
  inc_message->rank = c->rank;
  
  int i;
  for(i = 0; i < c->dhtcount; i++){
    if(c->hts[i] == ht) break;
  }
  inc_message->ht_index = i;
  inc_message->mbits = val; //yea i get this isn't the best but it'll work without wasting that space
  memcpy(inc_message->key, &counter, sizeof(int));
  
  MPI_Send(inc_message, sizeof(sendBuf), MPI_CHAR, PDHT_COUNTER_HOLDER, PDHT_TAG_COMMAND, MPI_COMM_WORLD);
  
  MPI_Recv(&counter_val, 1, MPI_UNSIGNED_LONG, PDHT_COUNTER_HOLDER, PDHT_COUNTER_REPLY, MPI_COMM_WORLD, &status);

  return counter_val;
}


