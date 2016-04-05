/********************************************************/
/*                                                      */
/*  commsynch.c - PDHT communication completion ops     */
/*                                                      */
/*  author: d. brian larkins                            */
/*  created: 3/20/16                                    */
/*                                                      */
/********************************************************/

#include <pdht_impl.h>

#define __PDHT_BARRIER_INDEX 23

/**
 * @file
 * 
 * portals distributed hash table synchronization ops
 */


/**
 * pdht_fence - ensures completion of put/get operations to destination rank
 * @param rank rank of process to ensure completion
 */
void pdht_fence(int rank) {
}



/**
 * pdht_allfence - ensures completion of all pending put/get operations
 */
void pdht_allfence(void) {
}



/**
 * pdht_test - checks status of an asynchronous put/get operation
 * @param h handle of pending operation
 * @returns status of operation
 */
pdht_status_t pdht_test(pdht_handle_t h) {
  return PdhtStatusOK;
}



/**
 * pdht_wait - blocks process until an asynchronous operation completes
 * @param h handle of pending operation
 * @returns status of operation
 */
pdht_status_t pdht_wait(pdht_handle_t h) {
  return PdhtStatusOK;
}



/**
 * pdht_waitrank - blocks process until all asynchronous operations complete wrt one process rank
 * @param h handle of pending operation
 * @returns status of operation
 */
pdht_status_t pdht_waitrank(int rank) {
  return PdhtStatusOK;
}



/**
 * pdht_waitall - blocks process until all asynchronous operations complete 
 * @returns status of operation
 */
pdht_status_t pdht_waitall(void) {
  return PdhtStatusOK;
}


/**
 * pdht_barrier_init() - initializes the barrier universe
 */
void pdht_barrier_init(pdht_context_t *c) {
  ptl_md_t md;
  ptl_handle_md_t md_handle;
  ptl_pt_index_t index;
  ptl_le_t le;
  ptl_handle_le_t le_handle;
  int ret;


  // initialize our barrier count to one (for ourself)
  c->ptl.barrier_count = 1;

  // create/bind send side MD
  md.start     = NULL;
  md.length    = 0;
  md.options   = PTL_MD_UNORDERED;
  md.eq_handle = PTL_EQ_NONE;
  md.ct_handle = PTL_CT_NONE;

  ret = PtlMDBind(c->ptl.lni, &md, &c->ptl.barrier_md);
  if (ret != PTL_OK) {
    pdht_dprintf("barrier_init: MD Bind failed\n");
    exit(1);
  }

  // create portals table entry
  ret = PtlPTAlloc(c->ptl.lni, 0, PTL_EQ_NONE, __PDHT_BARRIER_INDEX, &index);
  if ((ret != PTL_OK) || (index != __PDHT_BARRIER_INDEX)) {
    pdht_dprintf("barrier_init: PTAlloc failed\n");
    exit(1);
  } 

  // create  LE with counter for barrier message send/receive
  c->ptl.barrier_ct = PTL_INVALID_HANDLE;

  ret = PtlCTAlloc(c->ptl.lni, &c->ptl.barrier_ct);
  if (ret != PTL_OK)  {
    pdht_dprintf("barrier_init: CTAlloc failed\n");
    exit(1);
  } 

  le.start     = NULL;
  le.length    = 0;
  le.uid       = PTL_UID_ANY;
  le.options   = PTL_LE_OP_PUT | PTL_LE_ACK_DISABLE | PTL_LE_EVENT_CT_COMM | PTL_ME_EVENT_LINK_DISABLE;
  le.ct_handle = c->ptl.barrier_ct;
  ret = PtlLEAppend(c->ptl.lni, __PDHT_BARRIER_INDEX, &le, PTL_PRIORITY_LIST, NULL, &le_handle);
  if (ret != PTL_OK)  {
    pdht_dprintf("barrier_init: PtlLEAppend failed\n");
    exit(1);
  } 
} 



/**
 * pdht_barrier - blocks process until all processes reach the barrier
 */ 
void pdht_barrier(void) {
  ptl_process_t  p, l, r;
  ptl_size_t     test;
  ptl_ct_event_t cval;
  int ret; 

  // use binary tree of processes
  p.rank = ((c->rank + 1) >> 1) - 1;
  l.rank = ((c->rank + 1) << 1) - 1;
  r.rank = l.rank + 1;

  //pdht_dprintf("p: %lu l: %lu r: %lu\n",p.rank, l.rank, r.rank);

  // wait for children to enter barrier
  if (l.rank < c->size) {
    test = c->ptl.barrier_count++;
    if (r.rank < c->size) 
      test = c->ptl.barrier_count++;
    //pdht_dprintf("waiting for %d messages from l: %lu r: %lu\n", test, l.rank, r.rank);
    ret = PtlCTWait(c->ptl.barrier_ct, test, &cval);
    if (ret != PTL_OK) {
      //pdht_dprintf("barrier: CTWait failed (children)\n");
      exit(1);
    }
    //pdht_dprintf("children have entered barrier\n");
  }

  // children tell parents that they have entered
  if (c->rank > 0)  {
    //pdht_dprintf("notifying %lu\n", p.rank);
    ret = PtlPut(c->ptl.barrier_md, 0, 0, PTL_NO_ACK_REQ, p, __PDHT_BARRIER_INDEX, 0, 0, NULL, 0);
    test = c->ptl.barrier_count++;

    ret = PtlCTWait(c->ptl.barrier_ct, test, &cval);
    if (ret != PTL_OK) {
      pdht_dprintf("barrier: CTWait failed (parent)\n");
      exit(1);
    }
  }

  // wake up waiting children
  if (l.rank < c->size) {
    ret = PtlPut(c->ptl.barrier_md, 0, 0, PTL_NO_ACK_REQ, l, __PDHT_BARRIER_INDEX, 0, 0, NULL, 0);
    if (ret != PTL_OK) {
      pdht_dprintf("barrier: Put failed (wake left)\n");
      exit(1);
    }

    if (r.rank < c->size) {
      ret = PtlPut(c->ptl.barrier_md, 0, 0, PTL_NO_ACK_REQ, r, __PDHT_BARRIER_INDEX, 0, 0, NULL, 0);
      if (ret != PTL_OK) {
        pdht_dprintf("barrier: Put failed (wake right)\n");
        exit(1);
      }
    }
  }
} 
