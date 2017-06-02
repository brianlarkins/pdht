/*******************************************************/
/*                                                      */
/*  trig.c - PDHT triggered put setup                   */
/*                                                      */
/*  author: d. brian larkins                            */
/*  created: 5/20/16                                    */
/*                                                      */
/********************************************************/

#include <pdht_impl.h>

/**
 * @file
 * 
 * portals distributed hash table triggered ops code
 */


void *pdht_trig_progress(void *arg);

/**
 * pdht_trig_init -- initializes triggered operations for pending puts
 * @param dht - hash table data structure
 */
void pdht_trig_init(pdht_t *dht) {
  int ret;
  _pdht_ht_trigentry_t *hte;
  char *iter; // used for pointer math
  ptl_event_t ev;
  ptl_ct_event_t inc;
  unsigned hdrsize;

  // only xfer latter half of ME entry in put to pending ME
  hdrsize = sizeof(ptl_me_t) - offsetof(ptl_me_t, match_bits); 

  // setup triggered CT ops to increment success field by one
  inc.success = 1;
  inc.failure = 0;

  for (int ptindex = 0; ptindex < dht->ptl.nptes; ptindex++) {

    // allocate event queue for pending puts
    ret = PtlEQAlloc(dht->ptl.lni, dht->pendq_size, &dht->ptl.eq[ptindex]);
    if (ret != PTL_OK) {
      pdht_dprintf("pdht_trig_init: PtlEQAlloc failure\n");
      exit(1);
    }

    // allocate PTE for pending put
    ret = PtlPTAlloc(dht->ptl.lni, PTL_PT_ONLY_USE_ONCE | PTL_PT_FLOWCTRL,
        dht->ptl.eq[ptindex], dht->ptl.putindex_base+ptindex, &dht->ptl.putindex[ptindex]);
    if (ret != PTL_OK) {
      pdht_dprintf("pdht_trig_init: PtlPTAlloc failure\n");
      exit(1);
    }

    // iterator = ht[PTE * QSIZE] (i.e. PENDINGQ_SIZE per PTE)
    iter = (char *)dht->ht + ((dht->pendq_size * ptindex) * dht->entrysize);
    //pdht_dprintf("trig init append: ptindex: %d %d ht[%d] userp: %p\n", ptindex, dht->ptl.putindex[ptindex], pdht_find_bucket(dht, iter), iter);

    // append one-time match entres to the put PTE to catch incoming puts
    for (int i=0; i < dht->pendq_size; i++) {
      hte = (_pdht_ht_trigentry_t *)iter;

      // allocate per-pending elem trigger event counter
      ret = PtlCTAlloc(dht->ptl.lni, &hte->tct);
      if (ret != PTL_OK) {
        pdht_dprintf("pdht_trig_init: PtlCTAlloc failure: %s (%d)\n", pdht_ptl_error(ret), i);
        exit(1);
      }

      // set pending ME params / options
      hte->me.start         = &hte->me.match_bits; // each entry has a unique memory buffer
      hte->me.length        = hdrsize + PDHT_MAXKEYSIZE + dht->elemsize;
      hte->me.uid           = PTL_UID_ANY;
      hte->me.options       = PTL_ME_OP_PUT 
                            | PTL_ME_USE_ONCE 
                            | PTL_ME_EVENT_CT_COMM 
                            | PTL_ME_IS_ACCESSIBLE | PTL_ME_EVENT_UNLINK_DISABLE 
                            | PTL_ME_EVENT_LINK_DISABLE;
      hte->me.match_id.rank = PTL_RANK_ANY;
      hte->me.match_bits    = __PDHT_PENDING_MATCH;
      hte->me.ignore_bits   = 0xffffffffffffffff; // ignore it all
      hte->me.ct_handle     = hte->tct;

      // append ME to the pending ME list
      ret = PtlMEAppend(dht->ptl.lni, dht->ptl.putindex[ptindex], &hte->me, PTL_PRIORITY_LIST, hte, &hte->pme);
      if (ret != PTL_OK) {
        pdht_dprintf("pdht_trig_init: PtlMEAppend error (%d:%d)\n", i,hte->me.length);
        exit(1);
      }

      // fix up ME entry data for future triggered append
      hte->me.start         = &hte->key;
      hte->me.length        = PDHT_MAXKEYSIZE + dht->elemsize;
      hte->me.options       = PTL_ME_OP_GET 
                            | PTL_ME_IS_ACCESSIBLE 
                            | PTL_ME_EVENT_COMM_DISABLE
                            | PTL_ME_EVENT_UNLINK_DISABLE;
      hte->me.ignore_bits   = 0;

      // once match bits have been copied, append to active match list
      ret = PtlTriggeredMEAppend(dht->ptl.lni, dht->ptl.getindex[ptindex], &hte->me, PTL_PRIORITY_LIST,
          hte, &hte->ame, hte->tct, 1);
      if (ret != PTL_OK) {
        pdht_dprintf("pdht_trig_init: PtlTriggeredMEAppend error (iteration %d)\n",i );
        exit(1);
      }

      iter += dht->entrysize; // pointer math, danger.
    }
  }

  // nextfree points to the first empty hash entry that doesn't have a pending trigger setup
  dht->nextfree = dht->ptl.nptes * dht->pendq_size; // free = DEFAULT_TABLE_SIZE - PENDINGQ_SIZE

  if (c->dhtcount == 1)
    pthread_create(&c->progress_tid, NULL, pdht_trig_progress, NULL);
}



/**
 * pdht_trig_fini - cleans up triggered structures
 * @param dht - hash table data structure
 */
void pdht_trig_fini(pdht_t *dht) {
  struct timespec ts;
  _pdht_ht_trigentry_t *hte;
  char *iter;
  int ret;

  // disable any new messages arriving on the portals table put entry
  for (int ptindex=0; ptindex < dht->ptl.nptes; ptindex++) 
    PtlPTDisable(dht->ptl.lni, dht->ptl.putindex[ptindex]);

  if (c->dhtcount == 1) 
    pthread_cancel(c->progress_tid);

  // remove all match entries from the table
  iter = (char *)dht->ht;

  // XXX - TODO probably want to cancel all triggered ops on all pending entries

  for (int i=0; i<dht->maxentries; i++) {
    hte = (_pdht_ht_trigentry_t *)iter;

    // pending/put ME entries
    if (!PtlHandleIsEqual(hte->pme, PTL_INVALID_HANDLE)) {
      ret = PtlMEUnlink(hte->pme);
      while (ret == PTL_IN_USE) {
        ts.tv_sec = 0;
        ts.tv_nsec = 20000;  // 20ms
        nanosleep(&ts, NULL);
        ret = PtlMEUnlink(hte->pme);
      }
    }

    // active/get ME entries
    if (!PtlHandleIsEqual(hte->ame, PTL_INVALID_HANDLE)) {
      ret = PtlMEUnlink(hte->ame);
      while (ret == PTL_IN_USE) {
        ts.tv_sec = 0;
        ts.tv_nsec = 20000;  // 20ms
        nanosleep(&ts, NULL);
        ret = PtlMEUnlink(hte->ame);
      }
    }
    iter += dht->entrysize; // pointer math, danger.
  }

  // free our table entry
  for (int ptindex=0; ptindex < dht->ptl.nptes; ptindex++)  {
    PtlPTFree(dht->ptl.lni, dht->ptl.getindex[ptindex]); 
    // kill off event queue after PT free'd
    if (!PtlHandleIsEqual(dht->ptl.eq[ptindex], PTL_EQ_NONE)) {
      PtlEQFree(dht->ptl.eq[ptindex]);
    }
  }

  // release all storage for ht objects
  free(dht->ht);
}



/**
 * pdht_trig_progress - progress thread for triggered mode
 */ 
void *pdht_trig_progress(void *arg) {
  pdht_t *dht;
  _pdht_ht_trigentry_t *hte;
  char *index; // used for pointer math
  ptl_me_t me;
  unsigned hdrsize;
  int ret;

  // XXX TODO
  //   - this may vomit due to limits max CTs is pending q size
  //     not local table size

  while (1) {

    /* iterate over all active tables */
    for (int cur=0; cur < c->dhtcount; cur++) {
      dht = c->hts[cur]; 

      // count up link events from appending things to the active queue
      pthread_mutex_lock(&dht->completion_mutex);
      pdht_finalize_puts(dht);
      pthread_mutex_unlock(&dht->completion_mutex);

      index = (char *)dht->ht;
      index += (dht->nextfree * dht->entrysize);

      // only xfer latter half of ME entry in put to pending ME
      hdrsize = sizeof(ptl_me_t) - offsetof(ptl_me_t, match_bits); 

      // for each active PTE in this table,
      for (int ptindex=0; ptindex < dht->ptl.nptes; ptindex++) {

        // check to see if we've exhausted pending ME entries
        if (dht->stats.tappends[ptindex] >= dht->pendq_size) {
          
          // if so, refill the pending queue
          for (int i=0; i < dht->pendq_size; i++) {
            hte = (_pdht_ht_trigentry_t *)index;

            // allocate per-pending elem trigger event counter
            ret = PtlCTAlloc(dht->ptl.lni, &hte->tct);
            if (ret != PTL_OK) {
              pdht_dprintf("pdht_trig_init: PtlCTAlloc failure: %s (%d)\n", pdht_ptl_error(ret), i);
              exit(1);
            }

            // set pending ME params / options
            hte->me.start         = &hte->me.match_bits; // each entry has a unique memory buffer
            hte->me.length        = hdrsize + PDHT_MAXKEYSIZE + dht->elemsize;
            hte->me.uid           = PTL_UID_ANY;
            hte->me.options       = PTL_ME_OP_PUT 
                                  | PTL_ME_USE_ONCE 
                                  | PTL_ME_EVENT_CT_COMM 
                                  | PTL_ME_IS_ACCESSIBLE | PTL_ME_EVENT_UNLINK_DISABLE 
                                  | PTL_ME_EVENT_LINK_DISABLE;
            hte->me.match_id.rank = PTL_RANK_ANY;
            hte->me.match_bits    = __PDHT_PENDING_MATCH;
            hte->me.ignore_bits   = 0xffffffffffffffff; // ignore it all
            hte->me.ct_handle     = hte->tct;

            // append ME to the pending ME list
            ret = PtlMEAppend(dht->ptl.lni, dht->ptl.putindex[ptindex], &hte->me, PTL_PRIORITY_LIST, hte, &hte->pme);
            if (ret != PTL_OK) {
              pdht_dprintf("pdht_trig_init: PtlMEAppend error (%d:%d)\n", i,hte->me.length);
              exit(1);
            }

            // fix up ME entry data for future triggered append
            hte->me.start         = &hte->key;
            hte->me.length        = PDHT_MAXKEYSIZE + dht->elemsize;
            hte->me.options       = PTL_ME_OP_GET 
                                  | PTL_ME_IS_ACCESSIBLE 
                                  | PTL_ME_EVENT_COMM_DISABLE
                                  | PTL_ME_EVENT_UNLINK_DISABLE;
            hte->me.ignore_bits   = 0;

            // once match bits have been copied, append to active match list
            ret = PtlTriggeredMEAppend(dht->ptl.lni, dht->ptl.getindex[ptindex], 
                                       &hte->me, PTL_PRIORITY_LIST,
                                       hte, &hte->ame, hte->tct, 1);
            if (ret != PTL_OK) {
              pdht_dprintf("pdht_trig_init: PtlTriggeredMEAppend error (iteration %d)\n",i );
              exit(1);
            }

            index += dht->entrysize; // pointer math, danger.
            dht->nextfree++;
          } // refill

          dht->stats.tappends[ptindex] = 0; // reset the number of consumed pending entries

        } // if exhausted
      } // PTE loop
    } // HT loop
  } // forever loop
  return NULL;
}




void print_count(pdht_t *dht, char *msg) {
  if (c->rank == 1) {
    _pdht_ht_trigentry_t *hte = (_pdht_ht_trigentry_t *)dht->ht;
    ptl_ct_event_t ce;
    PtlCTGet(hte->tct, &ce);
    printf("%s\n", msg);
    printf("  counter: %"PRIu64" %"PRIu64"\n", ce.success, ce.failure);
    printf("  start: %p\n", hte->me.start);
    printf("  length: %"PRIu64" [elemsize: %d]\n", hte->me.length, dht->elemsize);
    printf("  match: %"PRIx64" %"PRIx64"\n", hte->me.match_bits, hte->me.ignore_bits);
    printf("  min_free: %"PRIu64"\n", hte->me.min_free);
  }
}
