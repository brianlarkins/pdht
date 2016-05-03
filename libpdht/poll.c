/*******************************************************/
/*                                                      */
/*  nbputget.c - PDHT asynch non-blocking operations    */
/*                                                      */
/*  author: d. brian larkins                            */
/*  created: 3/20/16                                    */
/*                                                      */
/********************************************************/

#include <pdht_impl.h>

/**
 * @file
 * 
 * portals distributed hash table polling task code
 */


/**
 * pdht_polling_init -- initializes polling thread
 * @param dht - hash table data structure
 */
void pdht_polling_init(pdht_t *dht) {
  int tablesize, ret;
  _pdht_ht_entry_t *hte;
  char *iter; // used for pointer math
  ptl_me_t me;
  ptl_event_t ev;


  // allocate event queue for pending puts
  ret = PtlEQAlloc(dht->ptl.lni, PDHT_PENDINGQ_SIZE, &dht->ptl.eq);
  if (ret != PTL_OK) {
    pdht_dprintf("pdht_polling_init: PtlEQAlloc failure\n");
    exit(1);
  }

  // allocate PTE for pending put
  ret = PtlPTAlloc(dht->ptl.lni, PTL_PT_ONLY_USE_ONCE | PTL_PT_FLOWCTRL,
                   dht->ptl.eq, __PDHT_PUT_INDEX, &dht->ptl.putindex);
  if (ret != PTL_OK) {
    pdht_dprintf("pdht_polling_init: PtlPTAlloc failure\n");
    exit(1);
  }

  // allocate array for hash table data
  dht->entrysize = (sizeof(_pdht_ht_entry_t)) + dht->elemsize;

  dht->ht = calloc(PDHT_DEFAULT_TABLE_SIZE, dht->entrysize);
  if (!dht->ht) {
    pdht_dprintf("pdht_polling_init: calloc error: %s\n", strerror(errno));
    exit(1);
  }

  // use a byte pointer as iterator over variable-sized element array
  iter = (char *)dht->ht;

  // default match-list entry values
  me.length      = dht->elemsize; // storing ME handle _and_ HT entry in each elem.
  me.ct_handle   = PTL_CT_NONE;
  me.uid         = PTL_UID_ANY;
  // disable AUTO unlink events, we just check for PUT completion
  me.options     = PTL_ME_OP_PUT | PTL_ME_USE_ONCE 
                 | PTL_ME_IS_ACCESSIBLE | PTL_ME_EVENT_UNLINK_DISABLE;
  me.match_id.rank = PTL_RANK_ANY;
  me.match_bits  = __PDHT_PUT_MATCH; // this is ignored, each one of these is a wildcard
  me.ignore_bits = 0xffffffffffffffff; // ignore it all

  dht->nextfree = PDHT_PENDINGQ_SIZE; // free = DEFAULT_TABLE_SIZE - PENDINGQ_SIZE

  // append one-time match entres to the put PTE to catch incoming puts
  for (int i=0; i < PDHT_DEFAULT_TABLE_SIZE; i++) {
    hte = (_pdht_ht_entry_t *)iter;
    hte->pme = PTL_INVALID_HANDLE; // initialize pending put ME as invalid
    hte->gme = PTL_INVALID_HANDLE; // initialize active/get ME as invalid

    if (i<PDHT_PENDINGQ_SIZE) {
      me.start  = &hte->data; // each entry has a unique memory buffer
      ret = PtlMEAppend(dht->ptl.lni, __PDHT_PUT_INDEX, &me, PTL_PRIORITY_LIST, hte, &hte->pme);
      if (ret != PTL_OK) {
        pdht_dprintf("pdht_polling_init: PtlMEAppend error\n");
        exit(1);
      } 
      // clean out the LINK events from the event queue
      PtlEQWait(dht->ptl.eq, &ev);
    }

    iter += dht->entrysize; // pointer math, danger.

  }

}



/**
 * pdht_polling_fini - cleans up polling structures
 * @param dht - hash table data structure
 */
void pdht_polling_fini(pdht_t *dht) {
  struct timespec ts;
  _pdht_ht_entry_t *hte;
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
    hte = (_pdht_ht_entry_t *)iter;

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
    if (!PtlHandleIsEqual(hte->gme, PTL_INVALID_HANDLE)) {
      ret = PtlMEUnlink(hte->gme);
      while (ret == PTL_IN_USE) {
        ts.tv_sec = 0;
        ts.tv_nsec = 20000;  // 20ms
        nanosleep(&ts, NULL);
        ret = PtlMEUnlink(hte->gme);
      }
    }
    iter += dht->entrysize; // pointer math, danger.
  }

  // free our table entry
  PtlPTFree(dht->ptl.lni, dht->ptl.getindex);

  // release all storage for ht objects
  free(dht->ht);
}



/**
 * pdht_poll - runs periodic polling events
 * @param dht distributed hash table data structure
 */
void pdht_poll(pdht_t *dht) {
  _pdht_ht_entry_t *hte;
  ptl_event_t ev;
  ptl_me_t me;
  char *index;
  int init, ret;

  init = 0;

  // need to run through the event queue for the putindex
  while ((ret = PtlEQGet(dht->ptl.eq,&ev)) == PTL_OK) {

    // only initialize ME data if a) have pending events and b) we haven't already
    if (init == 0) {
      // default match-list entry values
      me.length        = dht->elemsize; // storing ME handle _and_ HT entry in each elem.
      me.ct_handle     = PTL_CT_NONE;
      me.uid           = PTL_UID_ANY;
      // disable auto-unlink events, we just check for PUT completion
      me.options       = PTL_ME_OP_GET | PTL_ME_IS_ACCESSIBLE | PTL_ME_EVENT_UNLINK_DISABLE;
      me.match_id.rank = PTL_RANK_ANY;

      init = 1; // remember this day, my son.
    }


    // found something to do, is it something we care about?
    if (ev.type == PTL_EVENT_PUT) {
      hte = (_pdht_ht_entry_t *)ev.user_ptr;

      // if get ME is inactive, then this is a new put()
      if (PtlHandleIsEqual(hte->gme, PTL_INVALID_HANDLE)) {

        me.start = &hte->data; // hte points to entire entry
        me.options       = PTL_ME_OP_GET | PTL_ME_IS_ACCESSIBLE | PTL_ME_EVENT_UNLINK_DISABLE;
        me.match_bits = ev.match_bits; // copy over match_bits from put
        me.ignore_bits   = 0;
    
        // append this pending put to the active PTE match list
        ret = PtlMEAppend(dht->ptl.lni, __PDHT_GET_INDEX, &me, PTL_PRIORITY_LIST, hte, &hte->gme);
        if (ret != PTL_OK) {
          pdht_dprintf("pdht_poll: ME append failed (get)\n");
          exit(1);
        }

        // setup ME entry to replace the one that was just consumed
        index =  (char *)dht->ht;
        index += (dht->nextfree * dht->entrysize); // pointer math, caution

        hte = (_pdht_ht_entry_t *)index; // index into HT entry table
        me.start       = &hte->data;
        me.options     = PTL_ME_OP_PUT | PTL_ME_USE_ONCE 
                       | PTL_ME_IS_ACCESSIBLE | PTL_ME_EVENT_UNLINK_DISABLE;
        me.match_bits  = __PDHT_PUT_MATCH; // this is ignored, each one of these is a wildcard
        me.ignore_bits = 0xffffffffffffffff; // ignore it all

        dht->nextfree++; // update next free space in local HT table
        assert(dht->nextfree < PDHT_DEFAULT_TABLE_SIZE);

        // add replacement entry to put/pending ME
        ret = PtlMEAppend(dht->ptl.lni, __PDHT_PUT_INDEX, &me, PTL_PRIORITY_LIST, hte, &hte->pme);
        if (ret != PTL_OK) {
           pdht_dprintf("pdht_poll: PtlMEAppend error (put)\n");
        }
        // clean out the LINK events from the event queue
        PtlEQWait(dht->ptl.eq, &ev);

      // otherwise this is an overwrite of an existing value
      } else {
         // XXX - what to do about overwrites?
         ; // do nothing -- overwrites are potential race conditions
      }

    } else {
      pdht_dprintf("pdht_poll: got event for %s\n", pdht_event_to_string(ev.type));
    }
  }

  // check for other problems
  if (ret != PTL_EQ_EMPTY) {
    pdht_dprintf("pdht_poll: event queue issue: %s\n", pdht_ptl_error(ret));
  }
}


