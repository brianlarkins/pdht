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


/**
 * pdht_trig_init -- initializes triggered operations for pending puts
 * @param dht - hash table data structure
 */
void pdht_trig_init(pdht_t *dht) {
  int ret;
  _pdht_ht_trigentry_t *hte;
  char *iter; // used for pointer math
  ptl_event_t ev;
  unsigned hdrsize;

  // only xfer latter half of ME entry in put to pending ME
  hdrsize = sizeof(ptl_me_t) - offsetof(ptl_me_t, match_bits); 

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

    iter = (char *)dht->ht + ((dht->pendq_size * ptindex) * dht->entrysize);

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
      hte->me.options       = PTL_ME_OP_PUT | PTL_ME_USE_ONCE | PTL_ME_EVENT_CT_COMM 
                            | PTL_ME_IS_ACCESSIBLE | PTL_ME_EVENT_UNLINK_DISABLE | PTL_ME_EVENT_LINK_DISABLE;
      hte->me.match_id.rank = PTL_RANK_ANY;
      hte->me.match_bits    = __PDHT_PENDING_MATCH;
      hte->me.ignore_bits   = 0xffffffffffffffff; // ignore it all
      hte->me.ct_handle     = hte->tct;

      // append ME to the pending ME list
      ret = PtlMEAppend(dht->ptl.lni, dht->ptl.putindex_base+ptindex, &hte->me, PTL_PRIORITY_LIST, hte, &hte->pme);
      if (ret != PTL_OK) {
        pdht_dprintf("pdht_trig_init: PtlMEAppend error\n");
        exit(1);
      }

      // fix up ME entry data for future triggered append
      hte->me.start         = &hte->key;
      hte->me.length        = PDHT_MAXKEYSIZE + dht->elemsize;
      hte->me.options       = PTL_ME_OP_GET | PTL_ME_IS_ACCESSIBLE | PTL_ME_EVENT_UNLINK_DISABLE;
      hte->me.ignore_bits   = 0;

      // clean out the LINK events from the event queue
      //PtlEQWait(dht->ptl.eq[ptindex], &ev); // could set PTL_ME_EVENT_LINK_DISABLE on ME

      // once match bits have been copied, append to active match list
      ret = PtlTriggeredMEAppend(dht->ptl.lni, __PDHT_ACTIVE_INDEX+ptindex, &hte->me, PTL_PRIORITY_LIST,
          hte, &hte->ame, hte->tct, 1);
    }
    iter += dht->entrysize; // pointer math, danger.
  }

  dht->nextfree = dht->ptl.nptes * dht->pendq_size; // free = DEFAULT_TABLE_SIZE - PENDINGQ_SIZE
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

  // kill off event queue
  for (int ptindex=0; ptindex < dht->ptl.nptes; ptindex++)  {
    if (!PtlHandleIsEqual(dht->ptl.eq[ptindex], PTL_EQ_NONE)) {
      PtlEQFree(dht->ptl.eq[ptindex]);
    }
  }

  // remove all match entries from the table
  iter = (char *)dht->ht;

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
  for (int ptindex=0; ptindex < dht->ptl.nptes; ptindex++) 
    PtlPTFree(dht->ptl.lni, dht->ptl.getindex[ptindex]);

  // release all storage for ht objects
  free(dht->ht);
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
