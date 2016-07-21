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
  int tablesize, ret;
  _pdht_ht_trigentry_t *hte;
  char *iter; // used for pointer math
  ptl_event_t ev;
  unsigned hdrsize;


  // allocate event queue for pending puts
  ret = PtlEQAlloc(dht->ptl.lni, PDHT_PENDINGQ_SIZE, &dht->ptl.eq);
  if (ret != PTL_OK) {
    pdht_dprintf("pdht_trig_init: PtlEQAlloc failure\n");
    exit(1);
  }

  // allocate PTE for pending put
  ret = PtlPTAlloc(dht->ptl.lni, PTL_PT_ONLY_USE_ONCE | PTL_PT_FLOWCTRL,
      dht->ptl.eq, __PDHT_PENDING_INDEX, &dht->ptl.putindex);
  if (ret != PTL_OK) {
    pdht_dprintf("pdht_trig_init: PtlPTAlloc failure\n");
    exit(1);
  }

  // allocate array for hash table data
  dht->entrysize = (sizeof(_pdht_ht_trigentry_t)) + dht->elemsize;

  dht->ht = calloc(PDHT_DEFAULT_TABLE_SIZE, dht->entrysize);
  if (!dht->ht) {
    pdht_dprintf("pdht_trig_init: calloc %lu bytes (entrysize: %lu) yields error: %s\n", PDHT_DEFAULT_TABLE_SIZE * dht->entrysize, dht->entrysize, strerror(errno));
    exit(1);
  }

  // use a byte pointer as iterator over variable-sized element array
  iter = (char *)dht->ht;

  // only xfer latter half of ME entry in put to pending ME
  hdrsize = sizeof(ptl_me_t) - offsetof(ptl_me_t, match_bits); 

  // append one-time match entres to the put PTE to catch incoming puts
  for (int i=0; i < PDHT_DEFAULT_TABLE_SIZE; i++) {
    hte = (_pdht_ht_trigentry_t *)iter;

    /* these statements initialize _all_ DHT entries */

    // setup match list entry template
    hte->me.length      = dht->elemsize + hdrsize; // need to include latter half of ptl_me_t in xfer
    hte->me.uid         = PTL_UID_ANY;
    // disable AUTO unlink events, we just check for PUT completion
    hte->me.options     = PTL_ME_OP_PUT | PTL_ME_USE_ONCE | PTL_ME_EVENT_CT_COMM
      | PTL_ME_IS_ACCESSIBLE | PTL_ME_EVENT_UNLINK_DISABLE;
    hte->me.match_id.rank = PTL_RANK_ANY;
    hte->me.match_bits  = __PDHT_PENDING_MATCH; // this is ignored, each one of these is a wildcard
    hte->me.ignore_bits = 0xffffffffffffffff; // ignore it all

    hte->pme = PTL_INVALID_HANDLE; // initialize pending put ME as invalid
    hte->ame = PTL_INVALID_HANDLE; // initialize active/get ME as invalid

    /* these statements initialize _pending_ DHT buckets */

    if (i<PDHT_PENDINGQ_SIZE) {

      // allocate per-pending elem trigger event counter
      ret = PtlCTAlloc(dht->ptl.lni, &hte->tct);
      if (ret != PTL_OK) {
        pdht_dprintf("pdht_trig_init: PtlCTAlloc failure\n");
        exit(1);
      }


      hte->me.start  = &hte->me.match_bits; // each entry has a unique memory buffer
      //hte->me.length = dht->elemsize; // done with ME/match_bits + data stuff
      hte->me.length = (sizeof(ptl_me_t) - offsetof(ptl_me_t, match_bits)) + PDHT_MAXKEYSIZE + dht->elemsize;
      hte->me.ct_handle   = hte->tct;

      //if ((c->rank == 1) && (i==0)) {
      //   print_count(dht, "pending:");
      //}

      // append ME to the pending ME list
      ret = PtlMEAppend(dht->ptl.lni, __PDHT_PENDING_INDEX, &hte->me, PTL_PRIORITY_LIST, hte, &hte->pme);
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
      PtlEQWait(dht->ptl.eq, &ev);

      //if ((c->rank == 1) && (i==0)) {
      // print_count(dht, "trigger:");
      //}

      // once match bits have been copied, append to active match list
      ret = PtlTriggeredMEAppend(dht->ptl.lni, __PDHT_ACTIVE_INDEX, &hte->me, PTL_PRIORITY_LIST,
          hte, &hte->ame, hte->tct, 1);
    }
    iter += dht->entrysize; // pointer math, danger.
  }

  dht->nextfree = PDHT_PENDINGQ_SIZE; // free = DEFAULT_TABLE_SIZE - PENDINGQ_SIZE
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
  PtlPTDisable(dht->ptl.lni, dht->ptl.putindex);

  // kill off event queue
  if (!PtlHandleIsEqual(dht->ptl.eq, PTL_EQ_NONE)) {
    PtlEQFree(dht->ptl.eq);
  }

  // remove all match entries from the table
  iter = (char *)dht->ht;

  for (int i=0; i<PDHT_DEFAULT_TABLE_SIZE; i++) {
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
  PtlPTFree(dht->ptl.lni, dht->ptl.getindex);

  // release all storage for ht objects
  free(dht->ht);
}


void print_count(pdht_t *dht, char *msg) {
  if (c->rank == 1) {
    _pdht_ht_trigentry_t *hte = (_pdht_ht_trigentry_t *)dht->ht;
    ptl_ct_event_t ce;
    PtlCTGet(hte->tct, &ce);
    printf("%s\n", msg);
    printf("  counter: %lu %lu\n", ce.success, ce.failure);
    printf("  start: %p\n", hte->me.start);
    printf("  length: %lu [elemsize: %d]\n", hte->me.length, dht->elemsize);
    printf("  match: %lx %lx\n", hte->me.match_bits, hte->me.ignore_bits);
    printf("  min_free: %lu\n", hte->me.min_free);
  }
}
