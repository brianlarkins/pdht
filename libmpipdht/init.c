#include <pdht.h>

void *pdht_comm(void *arg);

pdht_context_t *c = NULL;

void pdht_init(){
  c = (pdht_context_t *)malloc(sizeof(pdht_context_t));
  c->thread_active = 1;
  int result;
  MPI_Init_thread(NULL,NULL,MPI_THREAD_MULTIPLE,&result);
  int my_rank;
  int size;
  MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
   
  c->size = size;
  c->rank = my_rank;

  //making message datatype
  const int nitems = 2;
  int blocklengths[2] = {1,1};
  MPI_Datatype types[3] = {MPI_INT,MPI_INT,MPI_INT};
  MPI_Datatype mpi_msg_type;
  MPI_Aint offsets[3];

  offsets[0] = offsetof(message_t,type);
  offsets[1] = offsetof(message_t,rank);
  offsets[2] = offsetof(message_t,ht_index);

  MPI_Type_create_struct(nitems,blocklengths, offsets, types, &mpi_msg_type);
  MPI_Type_commit(&mpi_msg_type);

  c->msgType = mpi_msg_type;
}


pdht_t *pdht_create(int keysize, int elemsize,pdht_mode_t mode){
  pdht_t *dht;
  ht_t *ht = NULL;
  dht = (pdht_t *)malloc(sizeof(pdht_t));
  memset(dht,0,sizeof(pdht_t));
  dht->ht = ht;
  dht->elemsize = elemsize;
  dht->hashfn = pdht_hash;
  dht->keysize = keysize;




  if (!c){
    pdht_init();
  }
  c->hts[c->dhtcount] = dht;
  c->dhtcount++;
  


  if (c->dhtcount == 1){
    pthread_create(&c->tid,NULL,pdht_comm,NULL);
  }
  return dht;


}






void *pdht_comm(void *arg){
  
  pdht_t *dht = c->hts[0];
  message_t msg;
  MPI_Status status;
  int requester;
  while(c->thread_active){
    MPI_Recv(&msg,sizeof(msg),c->msgType,MPI_ANY_SOURCE,1,MPI_COMM_WORLD,&status);
    pdht_t *dht = c->hts[msg.ht_index];
    
    requester = msg.rank;
    char *iter;
    if (msg.type == pdhtStop) break;

    if(msg.type == pdhtGet){

      uint64_t recvBuf;
      MPI_Recv(&recvBuf,sizeof(recvBuf),MPI_UNSIGNED_LONG_LONG,requester,1,MPI_COMM_WORLD,&status);
      
      
      uint64_t mbits = recvBuf;
      
      
      ht_t *instance;
      
      
      
      char buf[PDHT_MAXKEYSIZE + dht->elemsize + sizeof(int)];
      HASH_FIND_INT(dht->ht,&mbits,instance);
      if (instance == NULL){
        int failure = 0;
        //printf("not found\n");
        memcpy(buf + PDHT_MAXKEYSIZE + dht->elemsize,&failure,sizeof(int));

      }
      else{
        iter = (char *)instance->value;
        int failure = 1;
        memcpy(buf,instance->value,PDHT_MAXKEYSIZE + dht->elemsize);
        memcpy(buf + PDHT_MAXKEYSIZE + dht->elemsize,&failure,sizeof(int));
      }

      MPI_Send(buf,sizeof(buf),MPI_CHAR,requester,0,MPI_COMM_WORLD);


    }
    else if(msg.type == pdhtPut){
      uint64_t mbits;
      char recvBuf[sizeof(mbits) + PDHT_MAXKEYSIZE + dht->elemsize];
      MPI_Recv(recvBuf,sizeof(recvBuf),MPI_CHAR,requester,1,MPI_COMM_WORLD,&status);  

      memcpy(&mbits,recvBuf,sizeof(mbits));
      
      ht_t *instance;
      HASH_FIND_INT(dht->ht,&mbits,instance);
      if (instance == NULL){
        instance = (ht_t*)malloc(sizeof(ht_t));
        instance->key = mbits;
        HASH_ADD_INT(dht->ht,key,instance);  
        instance->value = malloc(PDHT_MAXKEYSIZE + dht->elemsize);
      }
      

      memcpy(instance->value,recvBuf+sizeof(mbits),PDHT_MAXKEYSIZE + dht->elemsize);
      int flag = 1;
      MPI_Send(&flag,sizeof(int),MPI_INT,requester,0,MPI_COMM_WORLD);

    }
    


  }
}

void pdht_fini(){
  
  pthread_kill(c->tid,0);
  
  MPI_Finalize();
  free(c);
  


}

void free_entries(pdht_t *dht){
  ht_t *cur, *tmp;
  HASH_ITER(hh,dht->ht,cur,tmp){
    HASH_DEL(dht->ht,cur);
    free(cur);

  }

  

}


void pdht_free(pdht_t *dht){
  



  int dht_index;
  int i;
  for(i = 0;i < c->dhtcount;i++){
    if(c->hts[i] == dht){
      break;
    }
  }
  while(i < c->dhtcount){
    c->hts[i] = c->hts[i+1];
    i++;
  }


  free_entries(dht);
  free(dht);
  c->dhtcount--;
  if(c->dhtcount == 0){
    pdht_fini();



  }
}








