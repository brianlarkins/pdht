/********************************************************/
/*                                                      */
/*  atomics.c - PDHT atomic and counter operations      */
/*                                                      */
/*  author: d. brian larkins                            */
/*  created: 2/2/17                                     */
/*                                                      */
/********************************************************/

#include <pdht_impl.h>

struct _pdht_atomic_data_s {
  int64_t old;
  int64_t new;
  int64_t compare;
};
typedef struct _pdht_atomic_data_s _pdht_atomic_data_t;

/*
 * pdht_atomic_init - initialize data structures for atomic operations
 * @param ht - a PDHT hash table
 */
int pdht_atomic_init(pdht_t *ht) {
  ptl_md_t md;
  int ret;
  _pdht_atomic_data_t *atomicptr = NULL;

  // get CT ready for local counter MD events
  ret = PtlCTAlloc(ht->ptl.lni, &ht->ptl.atomic_ct);
  if (ret != PTL_OK) {
    pdht_dprintf("pdht_atomic_init: unable to create CT for atomics. -- %s\n", pdht_ptl_error(ret));
    return -1;
  }

  ret = posix_memalign((void **)&atomicptr, sizeof(int64_t), sizeof(_pdht_atomic_data_t));
  if (ret != 0) {
    pdht_dprintf("unable to get aligned memory block for atomic ops - %s\n", strerror(errno));
    return -1;
  }
  ht->ptl.atomic_scratch = atomicptr; // save this to use for the actual atomic op later
  memset(atomicptr, 0, sizeof(_pdht_atomic_data_t));

  md.start     = ht->ptl.atomic_scratch;
  md.length    = sizeof(_pdht_atomic_data_t);
  md.options   = PTL_MD_EVENT_CT_REPLY;
  md.eq_handle = PTL_EQ_NONE;
  md.ct_handle = ht->ptl.atomic_ct;

  ret = PtlMDBind(ht->ptl.lni, &md, &ht->ptl.atomic_md);
  if (ret != PTL_OK) {
    pdht_dprintf("pdht_atomic_init: unable to create MD for atomics");
    return -1;
  }
  return 0;
}



/*
 * pdht_atomic_free - cleanup after atomics 
 * @param ht - a PDHT hash table
 */
void pdht_atomic_free(pdht_t *ht) {
  PtlCTFree(ht->ptl.atomic_ct);    // free counter
  PtlMDRelease(ht->ptl.atomic_md); // free memory descriptor
  free(ht->ptl.atomic_scratch);    // free scratch space
}



/* 
 * pdht_atomic_cswap - atomically compare and swap an int64 inside a HT entry
 * @note this is hard-coded to int64_t, but could take a type parameter and 
 *        support all portals atomic types (and cswap variants)
 * @param ht - a PDHT hash table
 * @param key - the key of the value to atomically update
 * @param offset - offset inside the hash table entry to modify
 * @param old - copyout of the remote value _prior_ to the swap
 * @param new - new value to swap into HT entry
 */
pdht_status_t pdht_atomic_cswap(pdht_t *ht, void *key, size_t offset, int64_t *old, int64_t new) {
  ptl_match_bits_t mbits;
  uint32_t ptindex;
  ptl_pt_index_t ptl_ptindex;
  ptl_ct_event_t ctevent, ct2;
  ptl_process_t rank;
  _pdht_atomic_data_t *as;
  ptl_size_t oldoff, newoff;
  int retries = 5;
  int ret;

  // find out where we need to perform operation
  ht->hashfn(ht, key, &mbits, &ptindex, &rank);
  ptl_ptindex = ht->ptl.getindex[ptindex];

  // setup scratch space and get current counter value
  as = (_pdht_atomic_data_t *)ht->ptl.atomic_scratch; 
  as->old = 0;
  as->new = new;
  as->compare = *old;
  oldoff = offsetof(_pdht_atomic_data_t, old);
  newoff = offsetof(_pdht_atomic_data_t, new);

  do {
    PtlCTGet(ht->ptl.atomic_ct, &ctevent);
    //eprintf("pdht_atomic_cswap: pre: success: %lu fail: %lu rank: %lu\n", ctevent.success, ctevent.failure, rank.rank);

    // perform atomic cswap
    ret = PtlSwap(ht->ptl.atomic_md, oldoff, ht->ptl.atomic_md, newoff,
        sizeof(int64_t), rank, ptl_ptindex, mbits, offset + PDHT_MAXKEYSIZE,
        NULL, 0, &as->compare, PTL_CSWAP, PTL_INT64_T);
    if (ret != PTL_OK) {
      pdht_dprintf("pdht_atomic_cswap: PtlSwap failed\n");
      return PdhtStatusError;
    }

#define RELIABLE_TARGETS
#ifdef RELIABLE_TARGETS
    // wait for completion
    ret = PtlCTWait(ht->ptl.atomic_ct, ctevent.success+1, &ct2);
    if (ret != PTL_OK) {
      pdht_dprintf("pdht_atomic_cswap: PtlCTWait failed\n");
      return PdhtStatusError;
    }
#else
    ptl_size_t splusone = ctevent.success+1;
    int which;
    ret = PtlCTPoll(&ht->ptl.atomic_ct, &splusone, 1, 4, &ct2, &which);
    if (ret == PTL_CT_NONE_REACHED) {
      pdht_dprintf("pdht_atomic_cswap: timed out waiting for reply\n");
      return PdhtStatusError;
    } else if (ret != PTL_OK) {
      pdht_dprintf("pdht_atomic_cswap: PtlCTPoll failed\n");
      return PdhtStatusError;
    }
#endif

    if (ct2.failure > ctevent.failure) {
      ct2.success = 0;
      ct2.failure = -1;
      PtlCTInc(ht->ptl.atomic_ct, ct2);
      retries--;
    } else {
      *old = as->old;
      retries = -1;
      return PdhtStatusOK;
    }
  } while (retries > 0);

  //printf("oldoff: %d newoff: %d offset: %d %12"PRIx64" %d, %d %ld\n", oldoff, newoff, offset, mbits, ptl_ptindex, rank, as->old);
  //pdht_dprintf("pdht_atomic_cswap: failure : %d : %d\n", c->rank, rank.rank);

  return PdhtStatusError;
}

/*
 * pdht_counter_init - initializes a new atomic counter for a hash table
 * @param ht - a PDHT hash table
 * @param initval - initial counter value
 * @returns index of the new counter
 */
int pdht_counter_init(pdht_t *ht, int initval) {
  int cindex, ret;
  ptl_me_t me;
  ptl_md_t md;

  // XXX these structures are leaked and not cleaned up on hash table removal

  cindex = ht->countercount++;

  if (c->rank == 0) {
    // create new counter md in array of counter objects in pdht_t

    // create MD for the target side counter array
    ht->counters[cindex] = initval;

    me.start         = &ht->counters[cindex];
    me.length        = sizeof(ht->counters[cindex]);
    me.ct_handle     = PTL_CT_NONE; // no counter on target side
    me.uid           = PTL_UID_ANY;
    me.options       = PTL_ME_OP_PUT | PTL_ME_OP_GET;  // | PTL_ME_EVENT_CT_COMM;
    me.match_id.rank = PTL_RANK_ANY;
    me.match_bits    = cindex;
    me.ignore_bits   = 0;

    ret = PtlMEAppend(ht->ptl.lni, __PDHT_COUNTER_INDEX, &me, PTL_PRIORITY_LIST, NULL, &ht->ptl.centries[cindex]);
    if (ret != PTL_OK) {
      pdht_dprintf("pdht_counter_init: unable to append ME for counter.\n");
      return -1;
    }
  }

  // get CT ready for local counter MD events
  ret = PtlCTAlloc(ht->ptl.lni, &ht->ptl.countcts[cindex]);
  if (ret != PTL_OK) {
    pdht_dprintf("pdht_counter_init: unable to create CT for counter (%d). -- %s\n", cindex, pdht_ptl_error(ret));
    return -1;
  }

  ht->lcounts[cindex] = 0;
  md.start     = &ht->lcounts[cindex];
  md.length    = sizeof(ht->lcounts[cindex]);
  md.options   = PTL_MD_EVENT_CT_REPLY;
  md.eq_handle = PTL_EQ_NONE;
  md.ct_handle = ht->ptl.countcts[cindex];

  ret = PtlMDBind(ht->ptl.lni, &md, &ht->ptl.countmds[cindex]);
  if (ret != PTL_OK) {
    pdht_dprintf("pdht_counter_init: unable to create MD for counter (%d).\n", cindex);
    return -1;
  }

  return cindex;
}



/**
 * pdht_counter_reset - collectively reset an HT atomic counter
 * @param ht - a hash table
 * @param counter - which counter to reset
 */
void pdht_counter_reset(pdht_t *ht, int counter) {
  ptl_process_t r0 = { .rank = 0 };
  ptl_ct_event_t ctevent;
  int ret;

  PtlCTGet(ht->ptl.countcts[counter], &ctevent);

  ht->lcounts[counter] = 0; // set our value to zero

  // swap with rank 0's value.
  if (c->rank == 0) {

  // use PtlSwap if we want one-sided reset/set (add parameter to initialize value)
#if 0
  if (c->rank == 1)
    ret = PtlSwap(ht->ptl.countmds[counter], 0, 
                  ht->ptl.countmds[counter], 0, 
                  sizeof(uint64_t), r0, __PDHT_COUNTER_INDEX,
                  counter, 0, NULL, 0, PTL_SWAP, PTL_UINT64_T);
    if (ret != PTL_OK) {
      pdht_dprintf("pdht_counter_reset: swap error\n");
      return;
    }

    printf("rank 0 waiting\n");
    ret = PtlCTWait(ht->ptl.countcts[counter], ctevent.success+1, &ctevent);
    if (ret != PTL_OK) {
      pdht_dprintf("pdht_counter_reset: PtlCTWait failed\n");
      return;
    }
#endif
    ht->counters[counter] = 0;
  } 
  pdht_barrier();
  // not handling atomic failure (ctevent.failure)
}



/**
 * pdht_counter_inc - increments a counter by a given value
 * @param ht a hash table
 * @param counter index of the counter to modify
 * @param val amount to increment counter by
 * @returns existing counter value
 */
uint64_t pdht_counter_inc(pdht_t *ht, int counter, uint64_t val) {
  ptl_process_t r0 = { .rank = 0 };
  ptl_ct_event_t ctevent;
  int ret;


  // get the current counter values
  PtlCTGet(ht->ptl.countcts[counter], &ctevent);

  // set target value to the parameter
  ht->lcounts[counter] = val;

  // fetch and add to counter on rank zero
  ret = PtlFetchAtomic(ht->ptl.countmds[counter], 0, ht->ptl.countmds[counter], 0,
      sizeof(ht->lcounts[counter]), r0, __PDHT_COUNTER_INDEX,
      counter, 0, NULL, 0, PTL_SUM, PTL_UINT64_T);
  if (ret != PTL_OK) {
    pdht_dprintf("pdht_counter_inc: fetch add error\n");
    return -1;
  }

  ret = PtlCTWait(ht->ptl.countcts[counter], ctevent.success+1, &ctevent);
  if (ret != PTL_OK) {
    pdht_dprintf("pdht_counter_inc: PtlCTWait failed\n");
    return -1;
  }
  // not handling atomic failure (ctevent.failure)
  return ht->lcounts[counter];
}
