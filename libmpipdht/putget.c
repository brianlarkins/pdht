#include <pdht.h>

extern pdht_context_t *c;

void pdht_get(pdht_t *dht, void *key, void *value){
  ptl_process_t rank; //rank for hash
  uint64_t mbits; //mbits for hash
  uint32_t ptindex;
  dht->hashfn(dht,key,&mbits,&ptindex,&rank); //hashing
  char buf[PDHT_MAXKEYSIZE + dht->elemsize]; //creating buffer to recieve 
  ht_t *instance;
  char *iter;
  //if i have the memory just get it
  if (rank.rank == c->rank){
    

    HASH_FIND_INT(dht->ht,&mbits,instance);
    

    
    if (instance == NULL){
      dht->fuckups++;
      //printf("key : %d \n",*(int*)key);
      //printf("not found\n");
      return;
    }
    
    iter = (char*)instance->value;
    
    memcpy(buf,iter,PDHT_MAXKEYSIZE + dht->elemsize);

    //make sure it wasnt a collision
    if (memcmp(buf,key,dht->keysize) != 0){
      //printf("local found collision :(\n");
      dht->fuckups++;
      return;
    }
    //copy over info if no collision

    memcpy(value,buf+PDHT_MAXKEYSIZE,dht->elemsize);
  }
  //if i dont, then go get it from who does
  else{
    
    msg_type type = pdhtGet;
    MPI_Status status;
    
    message_t msg;
    for(int i = 0;i < c->dhtcount;i++){
      if(dht == c->hts[i]){
        msg.ht_index = i;

      }

    }


    msg.type = type;
    msg.rank = c->rank;

    MPI_Send(&msg,sizeof(message_t),c->msgType,rank.rank,1,MPI_COMM_WORLD);//what i want to do, and where to send it
    
    MPI_Send(&mbits,sizeof(mbits),MPI_UNSIGNED_LONG_LONG,rank.rank,1,MPI_COMM_WORLD);//matchbits of the thing i want
    
    char recvbuf[PDHT_MAXKEYSIZE + dht->elemsize + sizeof(int)];
    if(MPI_Recv(recvbuf,sizeof(recvbuf),MPI_CHAR,rank.rank,0,MPI_COMM_WORLD,&status) != MPI_SUCCESS){
      printf("Mpi fucked up \n");
    }
    if(*(int*)(recvbuf + PDHT_MAXKEYSIZE + dht->elemsize) == 0){
      return;

    }
    if(memcmp(recvbuf,key,dht->keysize) != 0){

      dht->fuckups++;
      //printf("remote found collision :(\n\n\n\n\n");
      return;

    }
    memcpy(value,recvbuf + PDHT_MAXKEYSIZE,dht->elemsize); 
  }

  return;
}



void pdht_put(pdht_t *dht, void *key, void *value){
  
  ptl_process_t rank;
  uint64_t mbits;
  uint32_t ptindex;
  dht->hashfn(dht,key,&mbits,&ptindex,&rank);
  if (rank.rank == c->rank){
    
    ht_t *instance;
    HASH_FIND_INT(dht->ht,&mbits,instance);
    if (instance == NULL){
      instance = (ht_t*)malloc(sizeof(ht_t));
      instance->key = mbits;
      HASH_ADD_INT(dht->ht,key,instance);  
      instance->value = malloc(PDHT_MAXKEYSIZE + dht->elemsize);
    }
    

    memcpy(instance->value,key,PDHT_MAXKEYSIZE);
    memcpy((char*)instance->value + PDHT_MAXKEYSIZE, value, dht->elemsize);
  }
  else{
    msg_type type = pdhtPut;
    MPI_Status status;
    message_t msg;
    msg.type = type;
    msg.rank = c->rank;
    
    for(int i = 0;i < c->dhtcount;i++){
      if(dht == c->hts[i]){
        msg.ht_index = i;

      }

    }

    MPI_Send(&msg,sizeof(msg),c->msgType,rank.rank,1,MPI_COMM_WORLD);//intial message saying what i want to do
    char sendbuf[sizeof(mbits) + PDHT_MAXKEYSIZE + dht->elemsize];
    memcpy(sendbuf,&mbits,sizeof(mbits));
    memcpy(sendbuf + sizeof(mbits),key,PDHT_MAXKEYSIZE);
    memcpy(sendbuf + sizeof(mbits) + PDHT_MAXKEYSIZE,value,dht->elemsize);
    MPI_Send(sendbuf,sizeof(sendbuf),MPI_CHAR,rank.rank,1,MPI_COMM_WORLD); //second message with information on how to do the put
    
    int flag;
    //reciveing confirmation
    MPI_Recv(&flag,sizeof(int),MPI_INT,rank.rank,0,MPI_COMM_WORLD,&status);
  }

}
