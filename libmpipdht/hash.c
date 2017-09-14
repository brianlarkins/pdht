#include <pdht.h>

extern pdht_context_t *c;

void pdht_hash(pdht_t *dht, void *key, ptl_match_bits_t *mbits, uint32_t *ptindex, ptl_process_t *rank){
  *mbits = CityHash64((char  *)key,dht->keysize);
  
  (*rank).rank = *mbits % c->size;
}



void pdht_sethash(pdht_t *dht, pdht_hashfunc hfun){

  dht->hashfn = hfun;


}
