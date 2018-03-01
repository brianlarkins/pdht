#include <pdht.h>
#include <execinfo.h>
#include <fcntl.h>
#include <signal.h>

#define MAIN_COLOR 0
#define SERVER_COLOR 1

/* local protos */
void *pdht_comm(void *arg);
void free_entries(pdht_t *dht);


// global pdht context
pdht_context_t *c = NULL;

/**
 * pdht_init - initializes PDHT system
 */
void pdht_init() {
  int result;
  int my_rank;
  int size;

  // init MPI, need MPI_THREAD_MULTIPLE 
#ifdef THREAD_MULTIPLE
  MPI_Init_thread(NULL,NULL,MPI_THREAD_MULTIPLE,&result);
  assert(result == MPI_THREAD_MULTIPLE);
#else
  MPI_Init(NULL, NULL);
#endif
  MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  char pname[MPI_MAX_PROCESSOR_NAME];
  int nlen;
  MPI_Get_processor_name(pname, &nlen);
  printf("%d: %s\n", my_rank, pname);
  fflush(stdout);
  // setup global context
  c = (pdht_context_t *)malloc(sizeof(pdht_context_t));
  memset(c, 0, sizeof(pdht_context_t));
  c->thread_active = 1;
  c->dhtcount = 0;
  c->size = size;
  c->rank = my_rank;
  c->maxbufsize = 0;
  c->pid = getpid();
  printf("c->rank : %d pid : %d \n", c->rank, getpid());
#if 0
  // define  message datatype for MPI
  const int nitems = 3;
  int blocklengths[4] = {1,1,1,1};
  MPI_Datatype types[4] = {MPI_INT,MPI_INT,MPI_INT,MPI_UNSIGNED_LONG};
  MPI_Datatype mpi_msg_type;
  MPI_Aint offsets[4];

  offsets[0] = offsetof(message_t,type);
  offsets[1] = offsetof(message_t,rank);
  offsets[2] = offsetof(message_t,ht_index);
  offsets[3] = offsetof(message_t,mbits);
  MPI_Type_create_struct(nitems,blocklengths, offsets, types, &mpi_msg_type);
  MPI_Type_commit(&mpi_msg_type);

  c->msgType = mpi_msg_type;
#endif 
}



/**
 * pdht_create -- allocates a new dht
 * @returns the newly minted dht
 */
pdht_t *pdht_create(int keysize, int elemsize, pdht_mode_t mode) {
  pdht_t *dht;

  ht_t *ht = NULL;
  int htbuflen;

  pthread_mutex_t *lock = malloc(sizeof(pthread_mutex_t));
  
  pthread_mutex_init(lock,NULL);

  dht = (pdht_t *)calloc(1,sizeof(pdht_t));
  dht->ht = ht;
  dht->elemsize = elemsize;
  dht->hashfn = pdht_hash;
  dht->keysize = keysize;
  dht->ptl.nptes = 1;
  dht->uthash_lock = lock;

  if (!c){
    pdht_init();
  }

#ifndef THREAD_MULTIPLE
  if(c->dhtcount == 1){
    //send request to my comm process to make me a new ht
    int target_rank = c->rank + 1;
    int buflen = sizeof(message_t) + 2 * (sizeof(int));
    char buf[buflen];

    message_t *msg;
    
    msg = (message_t *)buf;
    
    msg->type = pdhtCreateHt;
    msg->rank = c->rank;
    msg->ht_index = 0;
    
    memcpy(msg->key, &keysize, sizeof(int));
    memcpy(msg->key + sizeof(int), &elemsize, sizeof(int));

    MPI_Ssend(buf, buflen, MPI_CHAR, target_rank, PDHT_TAG_COMMAND, MPI_COMM_WORLD);

  }
#endif

  c->hts[c->dhtcount] = dht;
  c->dhtcount++;
  htbuflen = sizeof(message_t) + PDHT_MAXKEYSIZE + elemsize;
  c->maxbufsize = c->maxbufsize >= htbuflen ?  c->maxbufsize : htbuflen;
#ifndef THREAD_MULTIPLE
  if (c->rank % 2){
    MPI_Comm s_comm;

    MPI_Comm_split(MPI_COMM_WORLD, SERVER_COLOR, c->rank, &s_comm);
    c->split_comm = s_comm;
    pdht_comm(NULL);
  }
#endif

  if (c->dhtcount == 1) {
    
#ifdef THREAD_MULTIPLE
    pthread_create(&c->comm_tid,NULL,pdht_comm,NULL);
#endif
    MPI_Comm_split(MPI_COMM_WORLD, MAIN_COLOR,c->rank, &(c->split_comm));
  }
  return dht;
}



/**
 * pdht_comm -- target-side communication thread for MPI PDHT
 * @param arg - unused
 */
void *pdht_comm(void *arg) {
  pdht_t *dht = NULL;
  //char msgbuf[c->maxbufsize];
  MPI_Status status;
  int requester;
  char *iter;
  ht_t *instance;
  uint64_t mbits;
  int flag = 1;
  char *buf = NULL, *tbuf = NULL;
  message_t *msg = NULL;
  reply_t *reply = NULL;
  int buflen;
  int need;
  int offset;
  int index;
  int *counter_index;
  uint64_t increment;
  unsigned long counter_value;
  uint64_t initval;
  
  int target;
  int *keysize;
  int *elemsize;
  int htbuflen;
  int last;
  int put_reply;

  MPI_Request *requests = calloc(sizeof(MPI_Request), c->size / 2);
  char *bufs = calloc(c->maxbufsize, c->size / 2);
  
  char *msgbuf;// = calloc(c->maxbufsize, 1);

  for(int i = 0; i < c->size / 2; i++){
    offset = i * c->maxbufsize;
    target = i * 2;
    MPI_Irecv(bufs + offset, c->maxbufsize, MPI_CHAR, target, PDHT_TAG_COMMAND, MPI_COMM_WORLD, &requests[i]);
  }
  while(c->thread_active) {
    MPI_Waitany(c->size / 2, requests, &index, &status);
    offset = index * c->maxbufsize;
    last = index;
    msgbuf = bufs + offset;
    msg = (message_t *)msgbuf; // cast to access message fields
    dht = c->hts[msg->ht_index];
    requester = msg->rank;
    switch (msg->type) {
      case pdhtStop:
        // game over, go home
#ifdef THREAD_MULTIPLE
        goto done;
        break;
#else
        for(int i = 0; i < c->size / 2; i++){
          if(i != index){
            MPI_Cancel(&requests[i]);
          }
        }
        MPI_Finalize();
        goto done;
#endif

      case pdhtGet:
        // get request, search for entry and send reply to requestor
        // make sure MPI send buffer is big enough
        need = sizeof(reply_t) + dht->elemsize;
        if ((!buf) || (buflen < need)) {
          tbuf = realloc(buf, need);
          if (!tbuf) {
            printf("%d: realloc failure. game over.\n", c->rank); fflush(stdout);
            MPI_Abort(MPI_COMM_WORLD, -1);
          }
          buf = tbuf;

        }
        buflen = need;

        reply = (reply_t *)buf; // cast so we can set header values
#ifdef THREAD_MULTIPLE
        pthread_mutex_lock(dht->uthash_lock);
#endif  
       
        HASH_FIND_INT(dht->ht,&(msg->mbits),instance);

#ifdef THREAD_MULTIPLE
        pthread_mutex_unlock(dht->uthash_lock);
#endif
        if (instance) {
      // found entry
          reply->status = 1; 
          memcpy(&reply->key, instance->value, PDHT_MAXKEYSIZE + dht->elemsize);
        } else {
          // not found
          reply->status = 0;
        }
        MPI_Send(buf,buflen,MPI_CHAR,requester,PDHT_TAG_REPLY,MPI_COMM_WORLD);
        break;


      case pdhtPut:
        // put request, check for existence and add/overwrite as needed
        

#ifdef THREAD_MULTIPLE
        pthread_mutex_lock(dht->uthash_lock);
#endif         
        HASH_FIND_INT(dht->ht,&msg->mbits,instance);
#ifdef THREAD_MULTIPLE       
        pthread_mutex_unlock(dht->uthash_lock);
#endif
        if (!instance) {
          // new entry -- create new HT entry
          instance = (ht_t*)calloc(1,sizeof(ht_t));
          instance->key = msg->mbits;
          instance->value = malloc(PDHT_MAXKEYSIZE + dht->elemsize);

#ifdef THREAD_MULTIPLE
          pthread_mutex_lock(dht->uthash_lock); 
#endif

          HASH_ADD_INT(dht->ht,key,instance);  

#ifdef THREAD_MULTIPLE
          pthread_mutex_unlock(dht->uthash_lock);
#endif   

        }
        // update HT entry with PUT data

        memcpy(instance->value,msg->key,PDHT_MAXKEYSIZE + dht->elemsize);
        MPI_Send(&put_reply, 1, MPI_INT, requester, PDHT_TAG_REPLY, MPI_COMM_WORLD);
        // send ack to requestor
        //MPI_Send(&flag,sizeof(int),MPI_INT,requester,PDHT_TAG_ACK,MPI_COMM_WORLD);
        break;


      case pdhtCounterReset:
        //reset counter
        counter_index = (int *)(msg->key);
        dht->counters[*counter_index] = 0;
        break;

      case pdhtCounterInc:
        increment = msg->mbits;
        counter_index = (int *)(msg->key);
        counter_value = dht->counters[*counter_index];
        MPI_Send(&counter_value, 1, MPI_UNSIGNED_LONG, msg->rank, PDHT_COUNTER_REPLY, MPI_COMM_WORLD);
        dht->counters[*counter_index] += increment;
        break;


      case pdhtCounterInit:
        if(msg->rank == 0){
          initval = *(uint64_t *)(msg->key);
          dht->counters[dht->counter_count] = initval;
          dht->counter_count++;
        }
        break;
      case pdhtCreateHt:

        
        keysize  = (int *)(msgbuf + sizeof(message_t));
        elemsize = (int *)(msgbuf + sizeof(message_t) + sizeof(int));


        dht = (pdht_t *)calloc(1,sizeof(pdht_t));


        ht_t *ht = NULL;

        dht->ht = ht;
        dht->elemsize = *elemsize;
        dht->hashfn = pdht_hash;
        dht->keysize = *keysize;
        dht->ptl.nptes = 1;

        c->hts[c->dhtcount] = dht;
        c->dhtcount++;


        htbuflen = sizeof(message_t) + PDHT_MAXKEYSIZE + *elemsize;
        c->maxbufsize = c->maxbufsize >= htbuflen ?  c->maxbufsize : htbuflen;
        bufs = realloc(bufs, c->size / 2 * c->maxbufsize);
        //msgbuf = realloc(msgbuf, c->maxbufsize);
        break;
    }
    MPI_Irecv(bufs + offset, c->maxbufsize, MPI_CHAR, last * 2, PDHT_TAG_COMMAND, MPI_COMM_WORLD, &requests[last]);
  }
done:
  if (buf) free(buf);


#ifndef THREAD_MULTIPLE
//  MPI_Finalize();
  exit(0);
#endif

}



/**
 * pdht_fini - clean up everything
 */
void pdht_fini() {
  message_t msg;
  int target_rank;
  // force comm thread out of service loop
  c->thread_active = 0;
  // send comm thread game over message
  msg.type = pdhtStop;
  msg.ht_index = 0;
  msg.rank = c->rank;
#ifdef THREAD_MULTIPLE
  target_rank = c->rank;
#else
  target_rank = c->rank + 1;
#endif
  MPI_Ssend(&msg, sizeof(message_t), MPI_CHAR, target_rank, 1, MPI_COMM_WORLD);
  MPI_Finalize();
  free(c);
}



/*
 * free_entries - clean up target-side HT for a PDHT
 * @param dht - ht to clean out
 */
void free_entries(pdht_t *dht){
  ht_t *cur, *tmp;
  HASH_ITER(hh,dht->ht,cur,tmp){
    HASH_DEL(dht->ht,cur);
    free(cur);
  }
}


/**
 * pdht_free - release all resources with a PDHT
 * @param dht - ht to free
 */
void pdht_free(pdht_t *dht) {
  int dht_index;
  int i;

  // bookkeep on ht list
  for(i = 0;i < c->dhtcount;i++){
    if(c->hts[i] == dht){
      break;
    }
  }
  while(i < c->dhtcount){
    c->hts[i] = c->hts[i+1];
    i++;
  }

  // clean out target side entries
  free_entries(dht);
  free(dht);
  c->dhtcount--;
  if(c->dhtcount == 0){
    pdht_fini();
  }
}


/**
 *  pdht_tune - non-functioning stub
 */
void pdht_tune(unsigned opts, pdht_config_t *config) {
  ;
}
void pdht_bthandler(int sig) {
  void *a[100];
  size_t size;

  size = backtrace(a, 100);
  fprintf(stderr, "Error: signal: %d:\n", sig);
  backtrace_symbols_fd(a,size, STDERR_FILENO);
  exit(1);
}
