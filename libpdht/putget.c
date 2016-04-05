/********************************************************/
/*                                                      */
/*  putget.c - PDHT put / get operations                */
/*                                                      */
/*  author: d. brian larkins                            */
/*  created: 3/20/16                                    */
/*                                                      */
/********************************************************/

#include <pdht.h>

/**
 * @file
 * 
 * portals distributed hash table put/get ops
 */

// XXX i think we need value size too
//  - maybe make keys a specific struct with void* + size to reduce clutter

/**
 * pdht_put - puts or overwrites an entry in the global hash table
 *   @param key - hash table key
 *   @param ksize - size of key
 *   @param value - value for table entry
 *   @returns status of operation
 */
pdht_status_t pdht_put(pdht_t *dht, void *key, int ksize, void *value) {

  // hash key -> rank + match bits
	
  return PdhtStatusOK;
}



/**
 * pdht_get - gets an entry from the global hash table
 *   @param key - hash table key
 *   @param ksize - size of key
 *   @param value - value for table entry
 *   @returns status of operation
 */
pdht_status_t pdht_get(pdht_t *dht, void *key, int ksize, void **value) {
  return PdhtStatusOK;
}

/**
 * pdht_insert - manually inserts a hash table entry into the global hash table
 *  @param bits - Portals match bits for the table entry
 *  @param value - value for table entry
 *  @returns status of operation
 */
 pdht_status_t pdht_insert(pdht *dht, ptl_match_bits_t bits, void *value) {

 }