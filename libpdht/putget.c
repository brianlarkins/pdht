/********************************************************/
/*                                                      */
/*  putget.c - PDHT put / get operations                */
/*                                                      */
/*  author: d. brian larkins                            */
/*  created: 3/20/16                                    */
/*                                                      */
/********************************************************/

#include <alloca.h>
#include <pdht_impl.h>

/**
 * @file
 * 
 * portals distributed hash table put/get ops
 */


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
  ptl_size_t loffset;
  ptl_size_t lsize;
  ptl_ct_event_t ctevent;
  ptl_event_t fault;
  int ret;
  unsigned which;
  ptl_me_t *mep;
  char *valp, *kval;
  pdht_status_t rval = PdhtStatusOK;

  PDHT_START_TIMER(dht, ptimer);
  PDHT_START_TIMER(dht, t3);
  PDHT_START_TIMER(dht, t4);
  dht->stats.puts++;

  // 1. hash key -> rank + match bits + element
  dht->hashfn(dht, key, &mbits, &rank);


  // 1.5 figure out what we need to send to far end
  //   - send just HT element usually, need to also send
  //     match bits for triggered-only updates
  switch (dht->pmode) {
    case PdhtPendingTrig:
      // XXX - THIS CODE IS HIGHLY DEPENDENDENT ON PORTALS ME DATA DEFINITION - XXX
      // to workaround, could issue two puts, one for match_bits, one for data
      //  trigger on +2 event counts
      mep = alloca(sizeof(ptl_me_t) + PDHT_MAXKEYSIZE + dht->elemsize); // no need to free later.
      valp = (char *)mep + sizeof(ptl_me_t);
      memcpy(valp, key, PDHT_MAXKEYSIZE); // ugh. copying. copy key into element.
      memcpy(valp + PDHT_MAXKEYSIZE, value, dht->elemsize); // double ugh. -- pointer math
      mep->match_bits = mbits;
      mep->ignore_bits = 0;
      mep->min_free = 0;

      loffset = (ptl_size_t)&mep->match_bits; // only send from match_bits field and beyond
      lsize = (sizeof(ptl_me_t) - offsetof(ptl_me_t, match_bits)) + PDHT_MAXKEYSIZE + dht->elemsize; 
      break;

    case PdhtPendingPoll:
    default:
      kval = alloca(PDHT_MAXKEYSIZE + dht->elemsize); // no need to free later.
      memcpy(kval,key,PDHT_MAXKEYSIZE);
      memcpy(kval + PDHT_MAXKEYSIZE,value,dht->elemsize);
      loffset = (ptl_size_t)kval;
      lsize = PDHT_MAXKEYSIZE + dht->elemsize;
      //pdht_dprintf("put: key: %lu\n", *(u_int64_t *)kval);
      //pdht_dprintf("put: value: %f\n", *(double *)((char *)kval + PDHT_MAXKEYSIZE));
      break;
  }

#ifdef PDHT_DEBUG_TRACE
  pdht_dprintf("put: key: %lu val: %lu onto pending queue of %d\n", *(unsigned long *)key, *(unsigned long *)value, rank);
#endif  
  PtlCTGet(dht->ptl.lmdct, &dht->ptl.curcounts);
  //pdht_dprintf("pdht_put: pre: success: %lu fail: %lu\n", dht->ptl.curcounts.success, dht->ptl.curcounts.failure);


  // 2. put hash entry on target
  ret = PtlPut(dht->ptl.lmd, loffset, lsize, PTL_CT_ACK_REQ, rank, dht->ptl.putindex, 
      mbits, 0, value, 0);
  if (ret != PTL_OK) {
    pdht_dprintf("pdht_put: PtlPut() failed\n");
    goto error;
  }

  PDHT_STOP_TIMER(dht, t3);

  ret = PtlCTWait(dht->ptl.lmdct, dht->ptl.curcounts.success+1, &ctevent);
  if (ret != PTL_OK) {
    pdht_dprintf("pdht_put: PtlCTWait() failed\n");
    goto error;
  }

  PDHT_STOP_TIMER(dht, t4);

  //pdht_dprintf("pdht_put: post: success: %lu fail: %lu\n", ctevent.success, ctevent.failure);

  if (ctevent.failure > dht->ptl.curcounts.failure) {
    ret = PtlEQWait(dht->ptl.lmdeq, &fault);
    if (ret == PTL_OK) {
      if (fault.ni_fail_type == PTL_NI_PT_DISABLED) {
        pdht_dprintf("pdht_put: flow control on remote rank: %d\n", rank);
        goto error;
      } else {
        pdht_dprintf("pdht_get: found fail event: %s\n", pdht_event_to_string(fault.type));   
        pdht_dump_event(&fault);
      }
    } else {
      pdht_dprintf("pdht_put: PtlEQWait() error: %s\n", pdht_ptl_error(ret));
    }
  }

#if 0
  // 3. need to check for fail event or success count (200ms timeout)
  ret = PtlCTPoll(&dht->ptl.lmdct, &dht->ptl.curcounts, 1, 200, &ctevent, &which);
  pdht_dprintf("pdht_put: success: %lu failure: %lu lcount: %lu\n", ctevent.success, ctevent.failure, dht->ptl.lcount);
  if (ret == PTL_OK) {
    pdht_dprintf("returning\n");
    goto done;

  } else if (ret == PTL_CT_NONE_REACHED) {
    pdht_dprintf("going again\n");
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
        goto done;
      }
    }
    goto done;

  } else {
    pdht_dprintf("pdht_put: PtlCTPoll() error\n");
    goto error;
  }
#endif

done:
  PDHT_STOP_TIMER(dht, ptimer);
  return rval;

error:
  PDHT_STOP_TIMER(dht, ptimer);
  return PdhtStatusError;
}



/**
 * pdht_get - gets an entry from the global hash table
 *   @param key - hash table key
 *   @param ksize - size of key
 *   @param value - value for table entry
 *   @returns status of operation
 */
pdht_status_t pdht_get(pdht_t *dht, void *key, void *value) {
  ptl_match_bits_t mbits; 
  unsigned long roffset = 0;
  ptl_ct_event_t ctevent;
  ptl_process_t rank;
  char buf[PDHT_MAXKEYSIZE + dht->elemsize];
  ptl_event_t ev;
  int ret;
  pdht_status_t rval = PdhtStatusOK;

  PDHT_START_TIMER(dht, gtimer);
  PDHT_START_TIMER(dht, t1);
  PDHT_START_TIMER(dht, t2);
  dht->stats.gets++;

  dht->hashfn(dht, key, &mbits, &rank);

  PtlCTGet(dht->ptl.lmdct, &dht->ptl.curcounts);
  //pdht_dprintf("pdht_get: pre: success: %lu fail: %lu\n", dht->ptl.curcounts.success, dht->ptl.curcounts.failure);

#ifdef PDHT_DEBUG_TRACE
  pdht_dprintf("pdht_get: key: %lu from active queue of %d with match: %lu\n", *(unsigned long *)key, rank, mbits);
#endif  


  ret = PtlGet(dht->ptl.lmd, (ptl_size_t)buf, PDHT_MAXKEYSIZE + dht->elemsize, rank, dht->ptl.getindex, mbits, roffset, NULL);
  if (ret != PTL_OK) {
    pdht_dprintf("pdht_get: PtlGet() failed\n");
    goto error;
  }
  PDHT_STOP_TIMER(dht, t1);

  // check for completion or failure
  ret = PtlCTWait(dht->ptl.lmdct, dht->ptl.curcounts.success+1, &ctevent);
  if (ret != PTL_OK) {
    pdht_dprintf("pdht_get: PtlCTWait() failed\n");
    goto error;
  }

  PDHT_STOP_TIMER(dht, t2);
  //pdht_dprintf("pdht_get: event counter: success: %lu failure: %lu\n", ctevent.success, ctevent.failure);


  if (ctevent.failure > dht->ptl.curcounts.failure) {
    ret = PtlEQWait(dht->ptl.lmdeq, &ev);
    if (ret == PTL_OK) {
      if (ev.type == PTL_EVENT_REPLY) {
#ifdef PDHT_DEBUG_TRACE
        pdht_dprintf("pdht_get: key: %lu not found\n", *(unsigned long *)key);
#endif  
        ctevent.success = 0;
        ctevent.failure = -1;
        PtlCTInc(dht->ptl.lmdct, ctevent);
        dht->stats.notfound++;
        rval = PdhtStatusNotFound;
        goto done;
      } else {
        pdht_dprintf("pdht_get: found fail event: %s\n", pdht_event_to_string(ev.type));   
        pdht_dump_event(&ev);
      }
    } else {
      pdht_dprintf("pdht_get: PtlEQWait() failed\n");
      goto error;
    }
  }

  //pdht_dprintf("collision checking: %lu == %lu (size: %d)\n", *(u_int64_t *)buf, *(u_int64_t *)key, dht->keysize);

  // fetched entry has key + value concatenated, validate key
  if (memcmp(buf, key, dht->keysize) != 0) {
    // keys don't match, this must be a collision
    dht->stats.collisions++;
    //pdht_dprintf("get: found collision between: %lu and %lu\n", *(u_int64_t *)key, *(u_int64_t *)buf);
    rval = PdhtStatusCollision;
    goto done;
  }

  // looks good, copy value to application buffer
  // skipping over the embedded key data (for collision detection)
  memcpy(value, buf + dht->keysize, dht->elemsize); // pointer math

done:
  // get of non-existent entry should hit fail counter + PTL_EVENT_REPLY event
  // in PTL_EVENT_REPLY event, we should get ni_fail_type
  // ni_fail_type should be: PTL_NI_DROPPED

  PDHT_STOP_TIMER(dht, gtimer);
  return rval;

error:
  PDHT_STOP_TIMER(dht, gtimer);
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
