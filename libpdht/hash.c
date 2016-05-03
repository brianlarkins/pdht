/********************************************************/
/*                                                      */
/*  hash.c - PDHT hash function operations              */
/*                                                      */
/*  author: d. brian larkins                            */
/*  created: 3/20/16                                    */
/*                                                      */
/********************************************************/

#include <pdht.h>

/**
 * @file
 * 
 * portals distributed hash function utilities
 */

static long _fakebits=0;

/** 
 * pdht_hash() - associative accumulate operation into a hashed object
 *  @param dht hash table structure
 *  @param key key of entry to hash
 *  @returns match bits for portals request
 */
void pdht_hash(pdht_t *dht, void *key, ptl_match_bits_t *mbits, ptl_process_t *rank) {
   // THIS IS ALL TOTAL BULLSHIT
   *mbits = _fakebits++;
   (*rank).rank  = 0; // always assume that 0 is where the hash table entries are
}
