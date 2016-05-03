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
  ptl_match_bits_t mbits; 
  ptl_process_t rank;
  ptl_size_t loffset = (ptl_size_t)(value); // this might be totally wrong
  ptl_ct_event_t ctevent;
  ptl_event_t fault;
  int ret;
  unsigned which;

  // 1. hash key -> rank + match bits + element
  pdht_hash(dht, key, &mbits, &rank);

  // 2. put hash entry on target
  ret = PtlPut(dht->ptl.lmd, loffset, dht->elemsize, PTL_CT_ACK_REQ, rank, dht->ptl.putindex, 
               mbits, 0, value, 0);
  if (ret != PTL_OK) {
     pdht_dprintf("pdht_get: PtlPut() failed\n");
     goto error;
  }

  // increment our counter for local put/gets from our MD
  dht->ptl.lcount++;

  // 3. need to check for fail event or success count (200ms timeout)
  ret = PtlCTPoll(&dht->ptl.lmdct, &dht->ptl.lcount, 1, 200, &ctevent, &which);
  if (ret == PTL_OK) {
    return PdhtStatusOK;
  } else if (ret == PTL_CT_NONE_REACHED) {
    // timed out, repeatdly check for fault or success
    while (1) {
      // check failure EQ for flow control
      ret = PtlEQPoll(&dht->ptl.lmdeq, 1, 100, &fault, &which);
      if (ret == PTL_OK) {
        // flow control on remote NI
        if (fault.ni_fail_type == PTL_NI_PT_DISABLED) {
          pdht_dprintf("pdht_put: flow control on remote rank: %d\n", rank);
          goto error;
        }

      } else if (ret != PTL_EQ_EMPTY) {
        pdht_dprintf("pdht_put: PtlEQPoll() error: %s\n", pdht_ptl_error(ret));
      }

      // check for success again
      ret = PtlCTPoll(&dht->ptl.lmdct, &dht->ptl.lcount, 1, 200, &ctevent, &which);
      if (ret == PTL_OK) {
        return PdhtStatusOK;
      }
    }
    return PdhtStatusOK; // never gets here due to infinite loop
  } else {
     pdht_dprintf("pdht_put: PtlCTPoll() error\n");
     goto error;
  }

error:
  return PdhtStatusError;
}



/**
 * pdht_get - gets an entry from the global hash table
 *   @param key - hash table key
 *   @param ksize - size of key
 *   @param value - value for table entry
 *   @returns status of operation
 */
pdht_status_t pdht_get(pdht_t *dht, void *key, void **value) {
  ptl_match_bits_t mbits; 
  int roffset = sizeof(_pdht_ht_entry_t);
  ptl_ct_event_t ctevent;
  ptl_process_t rank;
  ptl_size_t loffset = (ptl_size_t)(*value); // this might be totally wrong
  int ret;

#if 0
  pdht_hash(dht, key, &mbits, &rank);

  // assumes: that loffset = address of *value
  ret = PtlGet(dht->ptl.lmd, loffset, dht->elemsize, rank, dht->ptl.getindex, mbits, roffset, NULL);
  if (ret != PTL_OK) {
     pdht_dprintf("pdht_get: PtlGet() failed\n");
     goto error;
  } 

  dht->ptl.lcount++;

  ret = PtlCTWait(dht->ptl.lct, dht->ptl.lcount, &ctevent);
  if (ret != PTL_OK) {
     pdht_dprintf("pdht_get: PtlCTWait() failed\n");
     goto error;
  }
  if (ctevent.success != dht->ptl.lcount)
    dht->ptl.lcount = ctevent.success;

#endif
  return PdhtStatusOK;

error:
  return PdhtStatusError;
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
  
#if 0
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
  
  // XXX - need to upate ME with actual start, length fields to correctly locate hte

  ret = PtlMEAppend(dht->ptl.lni, dht->ptl.getindex, &dht->ptl.me, PTL_PRIORITY_LIST, NULL, &hte->me);
  if (ret != PTL_OK) {
     pdht_dprintf("pdht_insert: match-list insertion failed\n");
     goto error;
  } 
#endif 
  return PdhtStatusOK;

error:
  return PdhtStatusError;
 }
