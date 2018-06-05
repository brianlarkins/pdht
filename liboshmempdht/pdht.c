#include "pdht.h"
#include "city.c"

pdht_context_t *c = NULL;

void pdht_hash(pdht_t *dht, void *key, ptl_match_bits_t *mbits, uint32_t *ptindex, ptl_process_t *rank){
  *mbits = CityHash64((char *)key, dht->keysize);
  (*rank).rank = *mbits % c->size;
}

void func_put(void *item, void *args, int index){
  am_args_t *entry = args;
  ht_t *instance;
  pdht_t *ht = c->hts[entry->ht_index];
  HASH_FIND(hh, ht->ht, &entry->key, sizeof(uint64_t), instance);
  if(!instance){
      instance = calloc(1, sizeof(ht_t));
      instance->key = entry->key;
      instance->value = malloc(sizeof(uint64_t) + ht->elemsize);
      HASH_ADD(hh, ht->ht, key, sizeof(uint64_t), instance);
  }
  memcpy(instance->value, args + c->data_offset, sizeof(ht->elemsize));
}

void func_get(void *item, void *args, int index){
  get_args_t get_args = *(get_args_t *)args;
  pdht_t *ht = c->hts[get_args.ht_index];
  ht_t *instance;
  //TODO Collision checking... maybe
  HASH_FIND(hh, ht->ht, &(get_args.key), sizeof(uint64_t), instance);
  memcpy(item, instance->value, ht->elemsize);
}

void pdht_init(){
  shmem_init();
  shmem_init_am();
  c = (pdht_context_t *)calloc(sizeof(pdht_context_t), 1);

  c->dhtcount = 0;
  c->rank = shmem_my_pe();
  c->size = shmem_n_pes();

  c->data_offset = offsetof(am_args_t, val); 
  c->value_offset = offsetof(ht_t, value);

  c->get_handle = shmem_insert_cb(SHMEM_AM_GET, func_get);
  c->put_handle = shmem_insert_cb(SHMEM_AM_PUT, func_put);

}

//setup and tear down
pdht_t *pdht_create(size_t keysize, size_t elemsize, pdht_mode_t mode){
  pdht_t *dht;
  ht_t *ht = NULL;
  if(!c){
    pdht_init();
  }

  dht = (pdht_t *)calloc(sizeof(pdht_t), 1);
  dht->ht = ht;
  dht->elemsize = elemsize;
  dht->keysize = keysize;

  pdht_sethash(dht, pdht_hash);
  /*XXX this should be changed to elemsize * 2 at some point because we will only have two threads calling ucp_progress */
  /* maybe ? b/c initiator needs to know */
  /* Change this to support returning PdhtStatus (elemsize + sizeof(int))*/
  dht->ht_shmem_space = shmem_malloc(elemsize * c->size); 

  c->hts[c->dhtcount] = dht;
  c->dhtcount++;
  
  shmem_barrier_all();

  return dht;

}

void pdht_fini(){
  shmem_finalize();
  free(c);
}

void pdht_free(pdht_t *dht){
  int i;

  for(i = 0; i < c->dhtcount; i++){
    if(c->hts[i] = dht){
      break;
    }
  }
  while(i < c->dhtcount){
    c->hts[i] = c->hts[i+1];
    i++;
  }
  
  shmem_free(dht->ht_shmem_space);
  ht_t *cur, *tmp;
  HASH_ITER(hh, dht->ht, cur, tmp){
    HASH_DEL(dht->ht, cur);
    free(cur);
  }
  free(dht);
  c->dhtcount--;

  if(c->dhtcount == 0){
    pdht_fini();
  }
}

//synchonization
void pdht_barrier(void){
  shmem_barrier_all();
}
void pdht_fence(pdht_t *dht){
  shmem_fence();
  shmem_fence_am();
}
pdht_status_t   pdht_reduce(void *in, void *out, pdht_reduceop_t op, pdht_datatype_t type, int elems);
pdht_status_t   pdht_allreduce(void *in, void *out, pdht_reduceop_t op, pdht_datatype_t type, int elems);

/*TODO optimize for local puts */
// communication operations
pdht_status_t pdht_put(pdht_t *dht, void *key, void *value){
  am_args_t *args;
  void *target_addr;
  ptl_process_t target;
  uint32_t ptindex; //gross portals hold over, we can probably get rid of this
  char msg[sizeof(am_args_t) + dht->elemsize];
  
  /* can get rid of all of these */
  int ht_index;
  ptl_match_bits_t mbits;
  
  for(ht_index = 0; ht_index < c->dhtcount; ht_index++){
    if(dht = c->hts[ht_index])
      break;
  }

  dht->hashfn(dht, key, &mbits, &ptindex, &target);
  
  args = (am_args_t *)msg;
  args->ht_index = ht_index;
  args->key = mbits;
  
  memcpy(args->val, value, dht->elemsize); 
  target_addr = ((char *)dht->ht_shmem_space) + (c->rank * dht->elemsize);
  shmem_put_am(target_addr, 1, dht->elemsize, target.rank, c->put_handle, args, sizeof(am_args_t) + dht->elemsize);
  return PdhtStatusOK;
}

pdht_status_t pdht_update(pdht_t *dht, void *key, void *value){
  return pdht_put(dht, key, value);
}

/*TODO optimize for local lookups */
pdht_status_t pdht_get(pdht_t *dht, void *key, void *value){
  get_args_t args;
  int ht_index;
  ptl_process_t target;
  ptl_match_bits_t mbits;
  void *target_addr;
  uint32_t ptindex; //another dumb portals thing

  for(ht_index = 0; ht_index < c->dhtcount; ht_index++){
    if(dht = c->hts[ht_index])
      break;
  }

  dht->hashfn(dht, key, &mbits, &ptindex, &target);

  args.ht_index = ht_index;
  args.key = mbits;

  target_addr = ((char *)(dht->ht_shmem_space)) + (c->rank * dht->elemsize);
 
  shmem_get_am(value, target_addr, 1, dht->elemsize, target.rank, c->get_handle, &args, sizeof(get_args_t));
  //TODO We have to send the actuall key and value back. This is not a big deal like Put because its not
  //an AM thing.
  return PdhtStatusOK;
}
pdht_status_t   pdht_persistent_get(pdht_t *dht, void *key, void *value){
  return pdht_get(dht, key, value);
}

//Hash function operations
void pdht_sethash(pdht_t *dht, pdht_hashfunc hfun){
  dht->hashfn = hfun;
}

//utility stuff
void            pdht_print_stats(pdht_t *dht){;}
double          pdht_average_time(pdht_t *dht, pdht_timer_t timer){;}
void            pdht_tune(unsigned opts, pdht_config_t *config){;}
