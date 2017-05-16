/********************************************************/
/*                                                      */
/*  atomics.c - PDHT atomic counter operations          */
/*                                                      */
/*  author: d. brian larkins                            */
/*  created: 2/2/17                                     */
/*                                                      */
/********************************************************/

#include <pdht_impl.h>

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
 * pdht_counter_reset - reset an HT atomic counter
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
  ret = PtlAtomic(ht->ptl.countmds[counter], 0, sizeof(uint64_t),
		  PTL_CT_ACK_REQ, r0, __PDHT_COUNTER_INDEX, counter, 0,
		  NULL, 0, PTL_SWAP, PTL_UINT64_T);
  if (ret != PTL_OK) {
    pdht_dprintf("pdht_counter_reset: swap error\n");
    return;
  }

  ret = PtlCTWait(ht->ptl.countcts[counter], ctevent.success+1, &ctevent);
  if (ret != PTL_OK) {
    pdht_dprintf("pdht_counter_reset: PtlCTWait failed\n");
    return;
  }

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
