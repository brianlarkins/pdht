/*******************************************************/
/*                                                      */
/*  nbputget.c - PDHT asynch non-blocking operations    */
/*                                                      */
/*  author: d. brian larkins                            */
/*  created: 3/20/16                                    */
/*                                                      */
/********************************************************/

#include <pdht_impl.h>

static char *pdht_event_to_string(ptl_event_kind_t evtype);

/**
 * @file
 * 
 * portals distributed hash table polling task code
 */



void pdht_polling_init(pdht_t *dht) {
  int tablesize, ret;
  _pdht_ht_entry_t *hte;
  char *iter; // used for pointer math
  ptl_me_t me;


  // allocate event queue for pending puts
  ret = PtlEQAlloc(dht->ptl.lni, PDHT_EVENTQ_SIZE, &dht->ptl.eq);
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

  // append one-time match entres to the put PTE to catch incoming puts
  for (int i=0; i < PDHT_DEFAULT_TABLE_SIZE; i++) {
    hte = (_pdht_ht_entry_t *)iter;
    me.start  = &hte->data; // each entry has a unique memory buffer
    ret = PtlMEAppend(dht->ptl.lni, __PDHT_PUT_INDEX, &me, PTL_PRIORITY_LIST, hte, &hte->me);
    iter += dht->entrysize; // pointer math, danger.
  }
}


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

    if (!PtlHandleIsEqual(hte->me, PTL_INVALID_HANDLE)) {
      ret = PtlMEUnlink(hte->me);
      while (ret == PTL_IN_USE) {
        ts.tv_sec = 0;
        ts.tv_nsec = 20000;  // 20ms
        nanosleep(&ts, NULL);
        ret = PtlMEUnlink(hte->me);
      }
    }
    iter += dht->entrysize; // pointer math, danger.
  }

  // free our table entry
  PtlPTFree(dht->ptl.lni, dht->ptl.getindex);

  // release all storage for ht objects
  free(dht->ht);
}



void pdht_poll(pdht_t *dht) {
  _pdht_ht_entry_t *hte;
  ptl_event_t ev;
  ptl_me_t me;
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
      // disable AUTO unlink events, we just check for PUT completion
      me.options       = PTL_ME_OP_PUT | PTL_ME_IS_ACCESSIBLE | PTL_ME_EVENT_UNLINK_DISABLE;
      me.match_id.rank = PTL_RANK_ANY;
      me.ignore_bits   = 0;

      init = 1; // remember this day, my son.
    }

    // found something to do, is it something we care about?
    if (ev.type == PTL_EVENT_PUT) {
      hte = (_pdht_ht_entry_t *)ev.user_ptr;

      me.start = &hte->data; // hte points to entire entry
      me.match_bits = ev.match_bits; // copy over match_bits from put
    
      ret = PtlMEAppend(dht->ptl.lni, __PDHT_GET_INDEX, &me, PTL_PRIORITY_LIST, hte, &hte->me);
      if (ret != PTL_OK) {
        pdht_dprintf("pdht_poll: ME append failed\n");
        exit(1);
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


static char *pdht_event_to_string(ptl_event_kind_t evtype) {
  char *ret = "(unmatched)";
  switch (evtype) {
    case PTL_EVENT_GET:
      ret = "PTL_EVENT_GET";
      break;
    case PTL_EVENT_GET_OVERFLOW:
      ret = "PTL_EVENT_GET_OVERFLOW";
      break;
    case PTL_EVENT_PUT:
      ret = "PTL_EVENT_PUT";
      break;
    case PTL_EVENT_PUT_OVERFLOW:
      ret = "PTL_EVENT_PUT_OVERFLOW";
      break;
    case PTL_EVENT_ATOMIC:
      ret = "PTL_EVENT_ATOMIC";
      break;
    case PTL_EVENT_ATOMIC_OVERFLOW:
      ret = "PTL_EVENT_ATOMIC_OVERFLOW";
      break;
    case PTL_EVENT_FETCH_ATOMIC:
      ret = "PTL_EVENT_FETCH_ATOMIC";
      break;
    case PTL_EVENT_FETCH_ATOMIC_OVERFLOW:
      ret = "PTL_EVENT_FETCH_ATOMIC_OVERFLOW";
      break;
    case PTL_EVENT_REPLY:
      ret = "PTL_EVENT_REPLY";
      break;
    case PTL_EVENT_SEND:
      ret = "PTL_EVENT_SEND";
      break;
    case PTL_EVENT_ACK:
      ret = "PTL_EVENT_ACK";
      break;
    case PTL_EVENT_PT_DISABLED:
      ret = "PTL_EVENT_PT_DISABLED";
      break;
    case PTL_EVENT_LINK:
      ret = "PTL_EVENT_LINK";
      break;
    case PTL_EVENT_AUTO_UNLINK:
      ret = "PTL_EVENT_AUTO_UNLINK";
      break;
    case PTL_EVENT_AUTO_FREE:
      ret = "PTL_EVENT_AUTO_FREE";
      break;
    case PTL_EVENT_SEARCH:
      ret = "PTL_EVENT_SEARCH";
      break;
  }
  return ret;
}
