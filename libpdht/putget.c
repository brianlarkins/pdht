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


// local-only discriminator for add/update/put operations
typedef enum { PdhtPTQPending, PdhtPTQActive } pdht_ptq_t;
static inline pdht_status_t pdht_do_put(pdht_t *dht, void *key, void *value, pdht_ptq_t which);
static void pdht_keystr(void *key, char* str);
static void pdht_dump_entry(pdht_t *dht, void *exp, void *act);

/**
 * pdht_put - puts or overwrites an entry in the global hash table
 *   @param key - hash table key
 *   @param ksize - size of key
 *   @param value - value for table entry
 *   @returns status of operation
 */
static inline pdht_status_t pdht_do_put(pdht_t *dht, void *key, void *value, pdht_ptq_t which) {
  ptl_match_bits_t mbits; 
  ptl_process_t rank;
  uint32_t ptindex;
  ptl_size_t loffset;
  ptl_size_t lsize;
  ptl_ct_event_t ctevent, current, reset;
  ptl_event_t fault;
  ptl_pt_index_t ptl_pt_index;
  int toobusy = 1, ret = 0, again = 0;
  ptl_me_t *mep;
  char *valp, *kval;
  pdht_status_t rval = PdhtStatusOK;
  struct timespec ts;

  PDHT_START_TIMER(dht, ptimer);

  dht->stats.puts++;
  // 1. hash key -> rank + match bits + element
  dht->hashfn(dht, key, &mbits, &ptindex, &rank);

  dht->stats.rankputs[rank.rank]++;

  // 1.5 figure out what we need to send to far end
  //   - send just HT element usually, need to also send
  //     match bits for triggered-only updates
  switch (dht->pmode) {
    case PdhtPendingTrig:

      // need to setup different ME start/lengths for pending/active puts
      if (which == PdhtPTQPending) {

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

      } else {
        // if putting to active queue, don't worry about match bits madness
        kval = alloca(PDHT_MAXKEYSIZE + dht->elemsize); // no need to free later.
        memcpy(kval,key,PDHT_MAXKEYSIZE);
        memcpy(kval + PDHT_MAXKEYSIZE,value,dht->elemsize);
        loffset = (ptl_size_t)kval;
        lsize = PDHT_MAXKEYSIZE + dht->elemsize;
      }
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

  //#define PDHT_DEBUG_TRACE
#ifdef PDHT_DEBUG_TRACE
  pdht_dprintf("put: key: %lu val: %lu onto pending queue of %d\n", *(unsigned long *)key, *(unsigned long *)value, rank);
#endif  

  // find out which ME queue we're off too...
  if (which == PdhtPTQPending) {
    ptl_pt_index = dht->ptl.putindex[ptindex];  // put/add to pending
    dht->stats.pendputs++;
  } else {
    ptl_pt_index = dht->ptl.getindex[ptindex];  // update to active
  }

  PtlCTGet(dht->ptl.lmdct, &dht->ptl.curcounts);
  current = dht->ptl.curcounts;
  //pdht_dprintf("pdht_put: pre: success: %lu fail: %lu\n", dht->ptl.curcounts.success, dht->ptl.curcounts.failure);

  // may have to re-attempt put if we run into flow control, repeat until successful
  do { 

    toobusy = 0; // default is to only repeat once

    // put hash entry on target
    ret = PtlPut(dht->ptl.lmd, loffset, lsize, PTL_ACK_REQ, rank, ptl_pt_index,
        mbits, 0, value, 0);

    if (ret != PTL_OK) {
      pdht_dprintf("pdht_put: PtlPut(key: %lu, rank: %d, ptindex: %d) failed: %s\n",
          *(long *)key, rank.rank, dht->ptl.putindex[ptindex], pdht_ptl_error(ret));
      goto error;
    }



    // wait for local completion
    ret = PtlCTWait(dht->ptl.lmdct, current.success+1, &ctevent);
    if (ret != PTL_OK) {
      pdht_dprintf("pdht_put: PtlCTWait() failed\n");
      goto error;
    }


    //pdht_dprintf("pdht_put: post: success: %lu fail: %lu\n", ctevent.success, ctevent.failure);

    // check for errors
    if (ctevent.failure > current.failure) {
      ret = PtlEQWait(dht->ptl.lmdeq, &fault);

      if (ret == PTL_OK) {

        if (fault.ni_fail_type == PTL_NI_PT_DISABLED) {
          // flow control event generated only on initial drop
          pdht_dprintf("pdht_put: flow control on remote rank: %d : %d\n", rank, dht->stats.puts);
          PtlCTGet(dht->ptl.lmdct, &current);
          ts.tv_sec = 0;
          ts.tv_nsec = 10000000; // 10ms
          nanosleep(&ts, NULL);
          again = 1;
          toobusy = 1; // reset loop sentinel
          reset.success = 0;
          reset.failure = -1;
          PtlCTInc(dht->ptl.lmdct, reset); // reset failure count
        } else { // if (fault.ni_fail_type != PTL_NI_OK) {
          pdht_dprintf("pdht_put: found fail event: %s\n", pdht_event_to_string(fault.type));
          pdht_dump_event(&fault);
        }

        // vim auto indenting is dumb as shit.
        } else {
          pdht_dprintf("pdht_put: PtlEQWait() error: %s\n", pdht_ptl_error(ret));
        }
      }

    if (again && (ctevent.success == dht->ptl.curcounts.success)) {
        pdht_dprintf("pdht_put: (again) flow control on remote rank: %d : %d\n", rank, dht->stats.puts);
        PtlCTGet(dht->ptl.lmdct, &current);
        //pdht_dprintf("pdht_put: post (again): success: %lu fail: %lu\n", current.success, current.failure);
        nanosleep(&ts, NULL);
        toobusy = 1; // reset loop sentinel
      }

    } while (toobusy);

done:
    PDHT_STOP_TIMER(dht, ptimer);
    return rval;

error:
    PDHT_STOP_TIMER(dht, ptimer);
    return PdhtStatusError;
  }



  /**
   * pdht_add - adds an entry in the global hash table
   *   @param key - hash table key
   *   @param ksize - size of key
   *   @param value - value for table entry
   *   @returns status of operation
   */
  pdht_status_t pdht_add(pdht_t *dht, void *key, void *value) {
    return pdht_do_put(dht,key,value, PdhtPTQPending);
  }



  /**
   * pdht_put - adds an entry to the global hash table
   *   @param key - hash table key
   *   @param ksize - size of key
   *   @param value - value for table entry
   *   @returns status of operation
   */
  pdht_status_t pdht_put(pdht_t *dht, void *key, void *value) {
    return pdht_do_put(dht,key,value,PdhtPTQPending);
  }



  /**
   * pdht_update - overwrites an entry in the global hash table
   *   @param key - hash table key
   *   @param ksize - size of key
   *   @param value - value for table entry
   *   @returns status of operation
   */
  pdht_status_t pdht_update(pdht_t *dht, void *key, void *value) {
    return pdht_do_put(dht,key,value,PdhtPTQActive);
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
  uint32_t ptindex;
  ptl_ct_event_t ctevent;
  ptl_process_t rank;
  char buf[PDHT_MAXKEYSIZE + dht->elemsize];
  ptl_event_t ev;
  int ret;
  pdht_status_t rval = PdhtStatusOK;

  PDHT_START_TIMER(dht, gtimer);


  dht->stats.gets++;

  dht->hashfn(dht, key, &mbits, &ptindex, &rank);

  dht->stats.ptcounts[ptindex]++;
  
  PtlCTGet(dht->ptl.lmdct, &dht->ptl.curcounts);
  //pdht_dprintf("pdht_get: pre: success: %lu fail: %lu\n", dht->ptl.curcounts.success, dht->ptl.curcounts.failure);

#ifdef PDHT_DEBUG_TRACE
  pdht_dprintf("pdht_get: key: %lu from active queue of %d with match: %lu\n", *(unsigned long *)key, rank, mbits);
#endif  

  ret = PtlGet(dht->ptl.lmd, (ptl_size_t)buf, PDHT_MAXKEYSIZE + dht->elemsize, rank, dht->ptl.getindex[ptindex], mbits, roffset, NULL);

  if (ret != PTL_OK) {
    pdht_dprintf("pdht_get: PtlGet(key: %lu, rank: %d, ptindex: %d/%d) failed: (%s) : %d\n", *(long *)key, rank.rank, ptindex, dht->ptl.putindex[ptindex], pdht_ptl_error(ret), dht->stats.gets);
    goto error;
  }
  
  // check for completion or failure
  ret = PtlCTWait(dht->ptl.lmdct, dht->ptl.curcounts.success+1, &ctevent);
  if (ret != PTL_OK) {
    pdht_dprintf("pdht_get: PtlCTWait() failed\n");
    goto error;
  }

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
      } else if (ev.ni_fail_type != PTL_NI_OK) {
        pdht_dprintf("pdht_get: found fail event: %s\n", pdht_event_to_string(ev.type));   
        pdht_dump_event(&ev);
      } 
    } else {
      pdht_dprintf("pdht_get: PtlEQWait() failed\n");
      goto error;
    }
  }
  // fetched entry has key + value concatenated, validate key
  if (memcmp(buf, key, dht->keysize) != 0) {
    // keys don't match, this must be a collision
    dht->stats.collisions++;
    pdht_dprintf("pdht_get: found collision. %d\n");
    pdht_dump_entry(dht, key, buf);
    rval = PdhtStatusCollision;
    PDHT_STOP_TIMER(dht,t4);
    goto done;
  }
  // looks good, copy value to application buffer
  // skipping over the embedded key data (for collision detection)
  memcpy(value, buf + PDHT_MAXKEYSIZE, dht->elemsize); // pointer math
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
pdht_status_t pdht_insert(pdht_t *dht, ptl_match_bits_t bits, uint32_t ptindex, void *key, void *value) {
  char *index;
  _pdht_ht_entry_t *hte;
  ptl_me_t me;
  int ret;

  // find our next spot 
  index = (char *)dht->ht;
  index += (dht->nextfree * dht->entrysize); // pointer math
  hte = (_pdht_ht_entry_t *)index;
  assert(dht->nextfree == pdht_find_bucket(dht, hte));

  // setup ME to append to active list
  me.start         = &hte->key;
  me.length        = PDHT_MAXKEYSIZE + dht->elemsize; // storing HT key _and_ HT entry in each elem.
  me.ct_handle     = PTL_CT_NONE;
  me.uid           = PTL_UID_ANY;
  // disable auto-unlink events, we just check for PUT completion
  me.options       = PTL_ME_OP_GET | PTL_ME_IS_ACCESSIBLE 
    | PTL_ME_EVENT_UNLINK_DISABLE | PTL_ME_EVENT_LINK_DISABLE;
  me.match_id.rank = PTL_RANK_ANY;
  me.match_bits    = bits;
  me.ignore_bits   = 0;

  memcpy(&hte->key, key, PDHT_MAXKEYSIZE); // fucking shoot me.
  memcpy(&hte->data, value, dht->elemsize);

  //pdht_dprintf("inserting val: %lu on rank %d ptindex %d [%d] matchbits %lu\n", 
  //     *(unsigned long *)value, c->rank, ptindex, dht->ptl.getindex[ptindex], bits);
  ret = PtlMEAppend(dht->ptl.lni, dht->ptl.getindex[ptindex], &me, PTL_PRIORITY_LIST, hte, &hte->ame);
  if (ret != PTL_OK) {
    pdht_dprintf("pdht_insert: ME append failed (active) : %s\n", pdht_ptl_error(ret));
    exit(1);
  }

  dht->nextfree++;

  return PdhtStatusOK;

error:
  return PdhtStatusError;
}


static void pdht_dump_entry(pdht_t *dht, void *exp, void *act) {
  char *cp;
  printf("expected: ");
  cp = exp;
  for (int i=0; i<dht->keysize; i++)  printf("%2hhx ", cp[i]);
  printf("\nactual:   ");
  cp = act;
  for (int i=0; i<dht->keysize; i++)  printf("%2hhx ", cp[i]);
  printf("\n");
}

  static void pdht_keystr(void *key, char* str) {
    long *k = (long *)key;
    // sprintf(str, "<%ld,%ld,%ld> @ %ld", k[0], k[1], k[2], k[3]); // MADNESS
    sprintf(str, "%ld", k[0]); // tests
  }
