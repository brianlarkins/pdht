/********************************************************/
/*                                                      */
/*  putget.c - PDHT put / get operations                */
/*                                                      */
/*  author: d. brian larkins                            */
/*  created: 3/20/16                                    */
/*                                                      */
/********************************************************/

#include <pdht_impl.h>

/**
 * @file
 * 
 * portals distributed hash table put/get ops
 */

//  - maybe make keys a specific struct with void* + size to reduce clutter

/**
 * pdht_put - puts or overwrites an entry in the global hash table
 *   @param key - hash table key
 *   @param ksize - size of key
 *   @param value - value for table entry
 *   @returns status of operation
 */
pdht_status_t pdht_put(pdht_t *dht, void *key, void *value) {

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
pdht_status_t pdht_get(pdht_t *dht, void *key, void **value) {
  int ret;
  ptl_match_bits_t mbits = pdht_hash(dht, key);

  ret = PtlGet(dht->ptl.md, 
  return PdhtStatusOK;
}


/**
 * pdht_insert - manually inserts a hash table entry into the global hash table
 *  @param bits - Portals match bits for the table entry
 *  @param value - value for table entry
 *  @returns status of operation
 */
 pdht_status_t pdht_insert(pdht_t *dht, ptl_match_bits_t bits, void *value) {
  _pdht_ht_entry_t *hte;
  int ret;
  
  // find our next spot -- pointer math
  hte = (_pdht_ht_entry_t *)((dht->nextfree * dht->entrysize) + (char *)dht->ht);
 
  // initialize value
  memcpy(hte->data, value, dht->elemsize);
  
  dht->nextfree++;
   
  // create counter for our entry
  ret = PtlCTAlloc(dht->ptl.lni, &hte->ct);
  if (ret != PTL_OK) {
     pdht_dprintf("pdht_insert: counter allocation failed\n");
     goto error;
  } 

  // update our ME template with new match bits
  dht->ptl.me.ct_handle = hte->ct;
  dht->ptl.me.match_bits = bits;

  ret = PtlMEAppend(dht->ptl.lni, dht->ptl.ptindex, &dht->ptl.me, PTL_PRIORITY_LIST, NULL, &hte->me);
  if (ret != PTL_OK) {
     pdht_dprintf("pdht_insert: match-list insertion failed\n");
     goto error;
  } 
  return PdhtStatusOK;

error:
  return PdhtStatusError;
 }
