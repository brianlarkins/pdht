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
ptl_match_bits pdht_hash(pdht_t *dht, void *key) {
   return _fakebits++;
}
