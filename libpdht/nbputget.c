/********************************************************/
/*                                                      */
/*  nbputget.c - PDHT asynch non-blocking operations    */
/*                                                      */
/*  author: d. brian larkins                            */
/*  created: 3/20/16                                    */
/*                                                      */
/********************************************************/

#include <pdht.h>

/**
 * @file
 * 
 * portals distributed hash table non-blocking ops
 */


/**
 * pdht_nbput - asynchronously puts or overwrites an entry in the global hash table
 *   @param key - hash table key
 *   @param value - value for table entry
 *   @returns handle for completion operations
 */
pdht_handle_t pdht_nbput(pdht_t *dht, void *key, void *value) {
  // worry about spending too much time finding a free NB handle here

  // hash key to get dest rank + match bits
  // find new entry in pending[rank] array
  // fill out entry
  // PtlPut()
  // return index
  return PDHT_NULL_HANDLE;
}



/**
 * pdht_nbget - asynchronously gets an entry from the global hash table
 *   @param key - hash table key
 *   @param value - value of the global entry
 *   @returns handle for completion operations
 */
pdht_handle_t pdht_nbget(pdht_t *dht, void *key, void **value) {
  // worry about spending too much time finding a free NB handle here

  // hash key to get dest rank + match bits
  // find new entry in pending[rank] array
  // fill out entry
  // PtlPut()
  // return index
  return PDHT_NULL_HANDLE;
}
