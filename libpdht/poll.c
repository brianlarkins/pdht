/********************************************************/
/*                                                      */
/*  poll.c - PDHT polling thread operations             */
/*                                                      */
/*  author: d. brian larkins                            */
/*  created: 3/20/16                                    */
/*                                                      */
/********************************************************/

#include <pthread.h>
#include <time.h>

#include <pdht_impl.h>

/**
 * @file
 * 
 * portals distributed hash table polling task code
 */

static pthread_t _pdht_poll_tid;

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


  // allocate array for hash table data
  dht->entrysize = (sizeof(_pdht_ht_entry_t)) + dht->elemsize;

  eprintf("polling init: hash table entry size: %lu (%d + %d)\n", dht->entrysize, sizeof(_pdht_ht_entry_t), dht->elemsize);

  dht->ht = calloc(PDHT_DEFAULT_TABLE_SIZE, dht->entrysize);
  if (!dht->ht) {
    pdht_dprintf("pdht_polling_init: calloc error: %s\n", strerror(errno));
    exit(1);
  }

  // use a byte pointer as iterator over variable-sized element array
  iter = (char *)dht->ht;


  // initialize entire hash table array
  for (int i=0; i < PDHT_DEFAULT_TABLE_SIZE; i++) {
    hte = (_pdht_ht_entry_t *)iter;
    hte->pme = PTL_INVALID_HANDLE; // initialize pending put ME as invalid
    hte->ame = PTL_INVALID_HANDLE; // initialize active ME as invalid
    iter += dht->entrysize; // pointer math, danger.
  }

  // default match-list entry values
  me.length      = PDHT_MAXKEYSIZE + dht->elemsize; // storing key _and_ value for each entry
  me.ct_handle   = PTL_CT_NONE;
  me.uid         = PTL_UID_ANY;
  // disable AUTO unlink events, we just check for PUT completion
  me.options     = PTL_ME_OP_PUT | PTL_ME_USE_ONCE 
    | PTL_ME_IS_ACCESSIBLE | PTL_ME_EVENT_UNLINK_DISABLE | PTL_ME_EVENT_LINK_DISABLE;
  me.match_id.rank = PTL_RANK_ANY;
  me.match_bits  = __PDHT_PENDING_MATCH; // this is ignored, each one of these is a wildcard
  me.ignore_bits = 0xffffffffffffffff; // ignore it all


  // deal with multiple PTEs per hash table to handle match list length
  for (int ptindex = 0; ptindex < dht->nptes; ptindex++) {

    // allocate event queue for pending puts
    ret = PtlEQAlloc(dht->ptl.lni, PDHT_PENDINGQ_SIZE, &dht->ptl.eq[ptindex]);
    if (ret != PTL_OK) {
      pdht_dprintf("pdht_polling_init: PtlEQAlloc failure\n");
      exit(1);
    }

    // allocate PTE for pending put
    ret = PtlPTAlloc(dht->ptl.lni, PTL_PT_ONLY_USE_ONCE | PTL_PT_FLOWCTRL,
        dht->ptl.eq[ptindex], __PDHT_PENDING_INDEX+ptindex, &dht->ptl.putindex[ptindex]);
    if (ret != PTL_OK) {
      pdht_dprintf("pdht_polling_init: PtlPTAlloc failure [%d] : %s\n", ptindex, pdht_ptl_error(ret));
      exit(1);
    }
    //pdht_dprintf("%d: %d %d\n", ptindex, __PDHT_PENDING_INDEX+ptindex, dht->ptl.putindex[ptindex]);

    // iterator = ht[PTE * QSIZE] (i.e. PENDINGQ_SIZE per PTE)
    iter = (char *)dht->ht + ((PDHT_PENDINGQ_SIZE * ptindex) * dht->entrysize); 

    // append one-time match entres to the put PTE to catch incoming puts
    for (int i=0; i < PDHT_PENDINGQ_SIZE; i++) {
      hte = (_pdht_ht_entry_t *)iter;
      me.start  = &hte->key; // each entry has a unique memory buffer (starts with key)

      //pdht_dprintf("init append: %d %d userp: %p\n", i, pdht_find_bucket(dht, hte), hte);
      ret = PtlMEAppend(dht->ptl.lni, dht->ptl.putindex[ptindex], &me, PTL_PRIORITY_LIST, hte, &hte->pme);
      if (ret != PTL_OK) {
        pdht_dprintf("pdht_polling_init: [%d] PTE: %d PtlMEAppend error: %s\n", i, dht->ptl.putindex[ptindex], pdht_ptl_error(ret));
        exit(1);
      } 
      // clean out the LINK events from the event queue
      // PtlEQWait(dht->ptl.eq, &ev);

      iter += dht->entrysize; // pointer math, danger.
    }
  }

  dht->nextfree = dht->nptes * PDHT_PENDINGQ_SIZE; // free = DEFAULT_TABLE_SIZE - PENDINGQ_SIZE

  pthread_create(&_pdht_poll_tid, NULL, pdht_poll, dht);
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
  for (int ptindex=0; ptindex < dht->nptes; ptindex++)
    PtlPTDisable(dht->ptl.lni, dht->ptl.putindex[ptindex]);

  // kill polling thread
  pthread_cancel(_pdht_poll_tid);

  // kill off event queues
  for (int ptindex=0; ptindex < dht->nptes; ptindex++) {
    if (!PtlHandleIsEqual(dht->ptl.eq[ptindex], PTL_EQ_NONE)) {
      PtlEQFree(dht->ptl.eq[ptindex]);
    }
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

  // free our table entries
  for (int ptindex=0; ptindex < dht->nptes; ptindex++)
    PtlPTFree(dht->ptl.lni, dht->ptl.putindex[ptindex]);

  // release all storage for ht objects
  free(dht->ht);
}



/**
 * pdht_poll - runs periodic polling events
 * @param dht distributed hash table data structure
 */
void *pdht_poll(void *arg) {
  pdht_t *dht = (pdht_t *)arg;
  _pdht_ht_entry_t *hte;
  ptl_event_t ev;
  ptl_me_t me;
  char *index;
  int ret;
  unsigned int ptindex;


  // NOTE: we have one thread per DHT, we could re-write to handle _all_ 
  //   polling activity in one thread, need to keep a list of "active" hts
  //   and switch over to PtlEQPoll()

  // default match-list entry values
  me.length        = PDHT_MAXKEYSIZE + dht->elemsize; // storing HT key _and_ HT entry in each elem.
  me.ct_handle     = PTL_CT_NONE;
  me.uid           = PTL_UID_ANY;
  // disable auto-unlink events, we just check for PUT completion
  me.options       = PTL_ME_OP_GET | PTL_ME_IS_ACCESSIBLE | PTL_ME_EVENT_UNLINK_DISABLE;
  me.match_id.rank = PTL_RANK_ANY;

  eprintf("Polling thread is active\n");


  // need to run through the event queue for the putindex
  while (1) {
    //if ((ret = PtlEQWait(dht->ptl.eq,&ev)) == PTL_OK)  {

    if ((ret = PtlEQPoll(dht->ptl.eq,dht->nptes, PTL_TIME_FOREVER, &ev, &ptindex)) == PTL_OK)  {
      PDHT_START_TIMER(dht,t5);

      // found something to do, is it something we care about?
      if (ev.type == PTL_EVENT_PUT) {
        hte = (_pdht_ht_entry_t *)ev.user_ptr;
        //eprintf("+");
#ifdef PDHT_DEBUG_TRACE
        pdht_dprintf("poll: key: %lu is pending queue bound for %d, new pending is: %d\n", *(unsigned long *)hte->key, pdht_find_bucket(dht, hte), dht->nextfree);
#endif  

        // if get ME is inactive, then this is a new put()
        if (PtlHandleIsEqual(hte->ame, PTL_INVALID_HANDLE)) {

          //eprintf("processing key: %lu\n", (u_int64_t)be64toh(hte->key));
          me.start         = &hte->key; // hte points to entire entry, start at key+val
          me.options       = PTL_ME_OP_GET | PTL_ME_IS_ACCESSIBLE 
            | PTL_ME_EVENT_UNLINK_DISABLE | PTL_ME_EVENT_LINK_DISABLE;
          me.match_bits = ev.match_bits; // copy over match_bits from put
          me.ignore_bits   = 0;

          // append this pending put to the active PTE match list
          ret = PtlMEAppend(dht->ptl.lni, dht->ptl.getindex[ptindex], &me, PTL_PRIORITY_LIST, hte, &hte->ame);
          if (ret != PTL_OK) {
            pdht_dprintf("pdht_poll: ME append failed (active): %s\n", pdht_ptl_error(ret));
            exit(1);
          }
        } 

        // setup ME entry to replace the one that was just consumed
        index =  (char *)dht->ht;
        index += (dht->nextfree * dht->entrysize); // pointer math, caution

        hte = (_pdht_ht_entry_t *)index; // index into HT entry table
        me.start       = &hte->key; // put stores key+val into ht
        me.options     = PTL_ME_OP_PUT | PTL_ME_USE_ONCE 
          | PTL_ME_IS_ACCESSIBLE | PTL_ME_EVENT_UNLINK_DISABLE 
          | PTL_ME_EVENT_LINK_DISABLE;
        me.match_bits  = __PDHT_PENDING_MATCH; // this is ignored, each one of these is a wildcard
        me.ignore_bits = 0xffffffffffffffff; // ignore it all


        assert(dht->nextfree <= PDHT_DEFAULT_TABLE_SIZE);

        //pdht_dprintf("append: %d %d userp: %p\n", dht->nextfree, pdht_find_bucket(dht, hte), hte);

        // add replacement entry to put/pending ME
        ret = PtlMEAppend(dht->ptl.lni, dht->ptl.putindex[ptindex], &me, PTL_PRIORITY_LIST, hte, &hte->pme);
        if (ret != PTL_OK) {
          pdht_dprintf("pdht_poll: PtlMEAppend error (pending): %s\n", pdht_ptl_error(ret));
        }

        dht->nextfree++; // update next free space in local HT table

      } else {
        pdht_dprintf("pdht_poll: got event for %s\n", pdht_event_to_string(ev.type));
        pdht_dump_event(&ev);
      }
    } else {
      // check for other problems
      pdht_dprintf("pdht_poll: event queue issue: %s\n", pdht_ptl_error(ret));
    }
    PDHT_STOP_TIMER(dht,t5);
  }
  return NULL;
  }
