/********************************************************/
/*                                                      */
/*  hash.c - PDHT hash function operations              */
/*                                                      */
/*  author: d. brian larkins                            */
/*  created: 3/20/16                                    */
/*                                                      */
/********************************************************/

#include <pdht_impl.h>

/**
 * @file
 * 
 * portals distributed hash function utilities
 */

/** 
 * pdht_hash() - associative accumulate operation into a hashed object
 *  @param dht hash table structure
 *  @param key key of entry to hash
 *  @returns match bits for portals request
 */
void pdht_hash(pdht_t *dht, void *key, ptl_match_bits_t *mbits, uint32_t *ptindex, ptl_process_t *rank) {
   // THIS IS ALL TOTAL BULLSHIT
   *mbits = *(ptl_match_bits_t *)key;
   //printf("mbits: %lu %lu\n", *mbits, *(unsigned long *)key);
   *ptindex = *mbits % dht->nptes;
   (*rank).rank  = 1; // always assume that 0 is where the hash table entries are
}



/**
 * pdht_sethash() - use an alternate hash function over the default
 *   @param dht hash table structure
 *   @param hfun pointer to hash function
 */
void pdht_sethash(pdht_t *dht, pdht_hashfunc hfun) {
  dht->hashfn = hfun;
}
