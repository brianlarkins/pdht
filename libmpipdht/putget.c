#include <pdht.h>

extern pdht_context_t *c;

/**
 * pdht_put - puts or overwrites an entry in the global hash table
 *   @param key - hash table key
 *   @param ksize - size of key
 *   @param value - value for table entry
 *   @returns status of operation
 */
pdht_status_t pdht_get(pdht_t *dht, void *key, void *value){
  ptl_process_t rank; //rank for hash
  uint64_t mbits; //mbits for hash
  uint32_t ptindex;
  char buf[sizeof(message_t) + PDHT_MAXKEYSIZE + dht->elemsize]; 
  char rbuf[sizeof(reply_t) + dht->elemsize]; 
  message_t *msg;
  reply_t *reply;
  ht_t *instance;
  MPI_Status status;
  int ret;
  int target_rank;

  dht->hashfn(dht,key,&mbits,&ptindex,&rank); //hashing

  
  target_rank = rank.rank + (1 - rank.rank % 2);

  msg = (message_t *)buf;
  msg->type = pdhtGet;

  for(int i = 0;i < c->dhtcount;i++){
    if(dht == c->hts[i]){
      msg->ht_index = i;
    }
  }
  msg->rank = c->rank;
  msg->mbits = mbits;
  // go ask remote for the element
  MPI_Send(msg, sizeof(message_t), MPI_CHAR, target_rank, PDHT_TAG_COMMAND, MPI_COMM_WORLD); 

  //MPI_Send(&mbits,sizeof(mbits),MPI_UNSIGNED_LONG_LONG,rank.rank,2,MPI_COMM_WORLD);//matchbits of the thing i want

  ret = MPI_Recv(rbuf, sizeof(rbuf), MPI_CHAR, target_rank, PDHT_TAG_REPLY,
                 MPI_COMM_WORLD,&status);

  assert(ret == MPI_SUCCESS);

  reply = (reply_t *)rbuf;
  if (reply->status == 0) 
    return PdhtStatusNotFound;

  if (memcmp(&reply->key,key,dht->keysize) != 0)
    return PdhtStatusCollision;

  memcpy(value,&reply->value,dht->elemsize); 

  return PdhtStatusOK;
}



/**
 * pdht_put - adds an entry to the global hash table
 *   @param key - hash table key
 *   @param ksize - size of key
 *   @param value - value for table entry
 *   @returns status of operation
 */
pdht_status_t pdht_put(pdht_t *dht, void *key, void *value){
  ptl_process_t rank;
  uint64_t mbits;
  uint32_t ptindex;
  ht_t *instance;
  message_t *msg;
  MPI_Status status;
  int flag;
  int target_rank;
  char sbuf[sizeof(message_t) + PDHT_MAXKEYSIZE + dht->elemsize];

  dht->hashfn(dht,key,&mbits,&ptindex,&rank);


  // prepare command message
  msg = (message_t *)sbuf;
  msg->type = pdhtPut;
  msg->rank = c->rank;
  msg->mbits= mbits;

  target_rank = rank.rank + (1 - rank.rank % 2);
  for(int i = 0;i < c->dhtcount;i++){
    if(dht == c->hts[i]){
      msg->ht_index = i;
    }
  }
  memcpy(msg->key, key, dht->keysize);
  memcpy(msg->key + PDHT_MAXKEYSIZE, value, dht->elemsize);
  // send command message to target
  MPI_Send(sbuf, sizeof(sbuf), MPI_CHAR, target_rank,
           PDHT_TAG_COMMAND ,MPI_COMM_WORLD);
  MPI_Recv(c->reply_buf, 1, MPI_INT, target_rank, PDHT_TAG_REPLY, MPI_COMM_WORLD, &status);
  return PdhtStatusOK;
}


pdht_status_t pdht_persistent_get(pdht_t *dht, void *key, void *value){
  while(pdht_get(dht, key, value) != PdhtStatusOK);
  return PdhtStatusOK;
}


pdht_status_t pdht_update(pdht_t *dht, void *key, void *value){
  return pdht_put(dht,key,value);
}
