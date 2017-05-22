/********************************************************/
/*                                                      */
/*  commsynch.c - PDHT communication completion ops     */
/*                                                      */
/*  author: d. brian larkins                            */
/*  created: 3/20/16                                    */
/*                                                      */
/********************************************************/

#include <pdht_impl.h>

/**
 * @file
 * 
 * portals distributed hash table synchronization ops
 */

static void reduce_zip(void *dest, void *src, pdht_datatype_t ty, int s);
static int pdht_collective_size(pdht_datatype_t type);

/**
 * pdht_collective_init() - initializes the collective universe
 */
void pdht_collective_init(pdht_context_t *c) {
  ptl_md_t md;
  ptl_handle_md_t md_handle;
  ptl_pt_index_t index;
  ptl_me_t me;
  ptl_handle_me_t me_handle;
  int ret;

  // initialize our collective count to one (for ourself)
  c->ptl.collective_count = 1; // keep general collective op counts separate from barriers
  c->ptl.barrier_count = 1;

  // allocate some scratch place for reductions (assumes longs/doubles are 64-bits each)
  c->ptl.collective_lscratch = calloc(2*PDHT_MAX_REDUCE_ELEMS, sizeof(uint64_t));
  // set start to middle of scratch space (danger, pointer math)
  c->ptl.collective_rscratch = ((char *)(c->ptl.collective_lscratch)) + (PDHT_MAX_REDUCE_ELEMS * sizeof(uint64_t));

  // create/bind send side MD
  md.start     = NULL;
  md.length    = PTL_SIZE_MAX;
  md.options   = PTL_MD_UNORDERED;
  md.eq_handle = PTL_EQ_NONE;
  md.ct_handle = PTL_CT_NONE;

  ret = PtlMDBind(c->ptl.lni, &md, &c->ptl.collective_md);
  if (ret != PTL_OK) {
    pdht_dprintf("collective_init: MD Bind failed: %s\n", pdht_ptl_error(ret));
    exit(1);
  }

  // create portals table entry
  ret = PtlPTAlloc(c->ptl.lni, 0, PTL_EQ_NONE, __PDHT_COLLECTIVE_INDEX, &index);
  if ((ret != PTL_OK) || (index != __PDHT_COLLECTIVE_INDEX)) {
    pdht_dprintf("collective_init: PTAlloc failed\n");
    exit(1);
  } 

  // create ME with counter for collective message send/receive
  c->ptl.collective_ct = PTL_INVALID_HANDLE;
  c->ptl.barrier_ct    = PTL_INVALID_HANDLE;

  // get counter to attach to match entry
  ret = PtlCTAlloc(c->ptl.lni, &c->ptl.collective_ct);
  if (ret != PTL_OK)  {
    pdht_dprintf("collective_init: CTAlloc failed\n");
    exit(1);
  } 

  // get another just for barrier ops (to avoid interleaving issues)
  ret = PtlCTAlloc(c->ptl.lni, &c->ptl.barrier_ct);
  if (ret != PTL_OK)  {
    pdht_dprintf("collective_init: CTAlloc failed\n");
    exit(1);
  } 

  // allocate matchlist entry for barriers
  me.start      = NULL;
  me.length     = 0;
  me.ct_handle  = c->ptl.barrier_ct;
  me.uid        = PTL_UID_ANY;
  me.options    = PTL_ME_OP_PUT | PTL_ME_ACK_DISABLE | PTL_ME_EVENT_CT_COMM | PTL_ME_EVENT_LINK_DISABLE;
  me.match_id.rank = PTL_RANK_ANY;
  me.match_bits  = __PDHT_BARRIER_MATCH;
  me.ignore_bits = 0;
  ret = PtlMEAppend(c->ptl.lni, __PDHT_COLLECTIVE_INDEX, &me, PTL_PRIORITY_LIST, NULL, &c->ptl.barrier_me);
  if (ret != PTL_OK)  {
    pdht_dprintf("collective_init: PtlMEAppend failed\n");
    exit(1);
  } 

  // allocate left matchlist entry for reductions
  me.start      = c->ptl.collective_lscratch;
  me.length     = sizeof(uint64_t)*PDHT_MAX_REDUCE_ELEMS;
  me.ct_handle  = c->ptl.collective_ct;
  me.match_bits  = __PDHT_REDUCE_LMATCH;
  ret = PtlMEAppend(c->ptl.lni, __PDHT_COLLECTIVE_INDEX, &me, PTL_PRIORITY_LIST, NULL, &c->ptl.reduce_lme);
  if (ret != PTL_OK)  {
    pdht_dprintf("collective_init: PtlMEAppend failed\n");
    exit(1);
  } 

  // allocate right matchlist entry for reductions
  me.start      = c->ptl.collective_rscratch;
  me.match_bits  = __PDHT_REDUCE_RMATCH;
  ret = PtlMEAppend(c->ptl.lni, __PDHT_COLLECTIVE_INDEX, &me, PTL_PRIORITY_LIST, NULL, &c->ptl.reduce_rme);
  if (ret != PTL_OK)  {
    pdht_dprintf("collective_init: PtlMEAppend failed\n");
    exit(1);
  } 
} 



/**
 * pdht_collective_fini - cleans up collective op infrastructure
 */
void pdht_collective_fini(void) {
  if (c->ptl.collective_lscratch) 
    free(c->ptl.collective_lscratch);

  PtlMEUnlink(c->ptl.barrier_me);
  PtlMEUnlink(c->ptl.reduce_lme);
  PtlMEUnlink(c->ptl.reduce_rme);
  PtlCTFree(c->ptl.collective_ct);
  PtlCTFree(c->ptl.barrier_ct);
  PtlPTFree(c->ptl.lni, __PDHT_COLLECTIVE_INDEX);
  PtlMDRelease(c->ptl.collective_md);
}



/**
 * pdht_barrier - blocks process until all processes reach the barrier
 */ 
void pdht_barrier(void) {
  ptl_process_t  p, l, r;
  ptl_size_t     test;
  ptl_ct_event_t cval, cval2;
  int ret; 

  // use binary tree of processes
  p.rank = ((c->rank + 1) >> 1) - 1;
  l.rank = ((c->rank + 1) << 1) - 1;
  r.rank = l.rank + 1;

  //pdht_dprintf("**** p: %lu l: %lu r: %lu\n", p.rank, l.rank, r.rank);
  //PtlCTGet(c->ptl.barrier_ct, &cval2);
  //pdht_dprintf("before: %lu\n", cval2.success);

  // wait for children to enter barrier
  if (l.rank < c->size) {
    test = c->ptl.barrier_count++;
    if (r.rank < c->size) 
      test = c->ptl.barrier_count++;
    pdht_lprintf(PDHT_DEBUG_VERBOSE, "waiting for %d messages from l: %lu r: %lu\n", test, l.rank, r.rank);
    ret = PtlCTWait(c->ptl.barrier_ct, test, &cval);
    if (ret != PTL_OK) {
      pdht_dprintf("barrier: CTWait failed (children)\n");
      exit(1);
    }
    pdht_lprintf(PDHT_DEBUG_VERBOSE, "children have entered barrier\n");
  }

  // children tell parents that they have entered
  if (c->rank > 0)  {
    pdht_lprintf(PDHT_DEBUG_VERBOSE, "notifying %lu\n", p.rank);
    ret = PtlPut(c->ptl.collective_md, 0, 0, PTL_NO_ACK_REQ, p, __PDHT_COLLECTIVE_INDEX, 
        __PDHT_BARRIER_MATCH, 0, NULL, 0);
    test = c->ptl.barrier_count++;

    ret = PtlCTWait(c->ptl.barrier_ct, test, &cval);
    if (ret != PTL_OK) {
      pdht_dprintf("barrier: CTWait failed (parent)\n");
      exit(1);
    }
  }

  // wake up waiting children
  if (l.rank < c->size) {
    pdht_lprintf(PDHT_DEBUG_VERBOSE, "notifying %lu\n", l.rank);
    ret = PtlPut(c->ptl.collective_md, 0, 0, PTL_NO_ACK_REQ, l, __PDHT_COLLECTIVE_INDEX, 
        __PDHT_BARRIER_MATCH, 0, NULL, 0);
    if (ret != PTL_OK) {
      pdht_dprintf("barrier: Put failed (wake left)\n");
      exit(1);
    }
  }

  if (r.rank < c->size) {
    pdht_lprintf(PDHT_DEBUG_VERBOSE, "notifying %lu\n", r.rank);
    ret = PtlPut(c->ptl.collective_md, 0, 0, PTL_NO_ACK_REQ, r, __PDHT_COLLECTIVE_INDEX, 
        __PDHT_BARRIER_MATCH, 0, NULL, 0);
    if (ret != PTL_OK) {
      pdht_dprintf("barrier: Put failed (wake right)\n");
      exit(1);
    }
  }
}



/**
 * pdht_reduce - blocks process until all processes reach the barrier
 * @param in buffer containing this processes source data
 * @param out buffer containing the final reduction value (only valid on root process)
 * @param op operation to perform
 * @param type type data type of elements
 * @param elems number of elements
 * @returns status of operation
 */ 
pdht_status_t pdht_reduce(void *in, void *out, pdht_reduceop_t op, pdht_datatype_t type, int elems) {
  ptl_process_t  p, l, r;
  ptl_size_t     test;
  ptl_ct_event_t cval, cval2;
  ptl_match_bits_t destme;
  int kids = 0;
  int ret; 
  int tysize = 0;

  tysize = pdht_collective_size(type);

  if (op != PdhtReduceOpSum) {
    pdht_dprintf("pdht_reduce: unsupported reduction operation\n");
    exit(1);
  }

  // figure out if we are the left or right child of our parent
  destme = ((c->rank % 2) != 0) ? __PDHT_REDUCE_LMATCH : __PDHT_REDUCE_RMATCH;

  // use binary tree of processes
  l.rank = ((c->rank + 1) << 1) - 1;
  r.rank = l.rank + 1;
  p.rank = ((c->rank + 1) >> 1) - 1;

  // parents wait for children to send reduce buffers
  if (l.rank < c->size) {
    kids++;
    test = c->ptl.collective_count++;
    if (r.rank < c->size)  {
      test = c->ptl.collective_count++;
      kids++;
    }
    //pdht_dprintf("waiting for %d messages from l: %lu r: %lu\n", test, l.rank, r.rank);
    ret = PtlCTWait(c->ptl.collective_ct, test, &cval);
    if (ret != PTL_OK) {
      pdht_dprintf("pdht_reduce: CTWait failed (children)\n");
      exit(1);
    }
  }

  if (kids > 0) {
    reduce_zip(in, c->ptl.collective_lscratch, type, elems);
    if (kids == 2) 
      reduce_zip(in, c->ptl.collective_rscratch, type, elems);
  }

  // send our reduced buffer to our parent
  if (c->rank > 0)  {
    ret = PtlPut(c->ptl.collective_md, (ptl_size_t)in, tysize*elems, PTL_NO_ACK_REQ, p, __PDHT_COLLECTIVE_INDEX, 
        destme, 0, NULL, 0);
  } else {
    // we're the root our local data has been zipped with each l/r subtree
    memcpy(out, in, tysize*elems);
  }
} 



/**
 * pdht_broadcast - broadcasts data from rank zero to all other processes
 * @param buf buffer to brodast from (rank 0)/to (other ranks)
 * @param type type of data to broadcast
 * @param elems number of data elements to broadcast
 * @returns status of operation
 */
pdht_status_t pdht_broadcast(void *buf, pdht_datatype_t type, int elems) {
  ptl_process_t  p, l, r;
  ptl_size_t     test;
  ptl_ct_event_t cval, cval2;
  ptl_match_bits_t destme;
  int ret; 
  int tysize = 0;

  // always use left ME for broadcast
  destme = __PDHT_REDUCE_LMATCH; 

  // use binary tree of processes
  l.rank = ((c->rank + 1) << 1) - 1;
  r.rank = l.rank + 1;

  tysize = pdht_collective_size(type);

  // wait for parent to send message
  if (c->rank > 0) {
    test = c->ptl.collective_count++;
    ret = PtlCTWait(c->ptl.collective_ct, test, &cval);
    if (ret != PTL_OK) {
      pdht_dprintf("pdht_broadcast: CTWait failed (parent)\n");
      exit(1);
    }
    // copy recieved data to our buffer
    memcpy(buf, c->ptl.collective_lscratch, tysize*elems);
  }

  // parents send the message to children (if any)
  if (l.rank < c->size) {
    ret = PtlPut(c->ptl.collective_md, (ptl_size_t)buf, tysize*elems, PTL_NO_ACK_REQ, l, 
                 __PDHT_COLLECTIVE_INDEX, destme, 0, NULL, 0);
    if (ret != PTL_OK) {
      pdht_dprintf("pdht_broadcast: Put failed (send left)\n");
      exit(1);
    }
  }

  if (r.rank < c->size) {
    ret = PtlPut(c->ptl.collective_md, (ptl_size_t)buf, tysize*elems, PTL_NO_ACK_REQ, r, 
                 __PDHT_COLLECTIVE_INDEX, destme, 0, NULL, 0);
    if (ret != PTL_OK) {
      pdht_dprintf("pdht_broadcast: Put failed (send left)\n");
      exit(1);
    }
  }
}



/**
 * pdht_allreduce - blocks process until all processes reach the barrier
 * @param in buffer containing this processes source data
 * @param out buffer containing the final reduction value (only valid on root process)
 * @param op operation to perform
 * @param type type data type of elements
 * @param elems number of elements
 * @returns status of operation
 */ 
pdht_status_t pdht_allreduce(void *in, void *out, pdht_reduceop_t op, pdht_datatype_t type, int elems) {
   pdht_status_t ret;
   
   // this is a lazy implementation. i get it.
   ret = pdht_reduce(in, out, op, type, elems);
   if (ret != PdhtStatusOK)
      return ret;
   pdht_barrier(); 
   ret = pdht_broadcast(out, type, elems);
   return ret;
}



/**
 * pdht_fence - ensures completion of put/get operations to destination rank
 * @param rank rank of process to ensure completion
 */
void pdht_fence(int rank) {
  // Q: this is only important for NB ops, put/get are blocking/synch, no?
  // Q: full events or counting? (more vs. less accounting data)
  // Q: what happens to pdht_handles for outstanding nb xfers after fence?
  //   - would be nice to invalidate them all after fencing (done!)
  //   - otherwise, need to mandate explicit call to a wait() op

  // requires: bookeeping on pending put/get/atomics
  //  - use an array of pending to each rank (100 for each process * 1K * size + status <= 1MB)
  //    - pending queues should self-flush when full (with call to this routine)
  //  - put needs to match all outstanding puts with PTL_EVENT_ACKs at Initiator
  //  - get needs to match all outstanding gets with PTL_EVENT_REPLYs at Initiator
  //  - fetchatomic matches PTL_EVENT_REPLY too (distinguishing?) (use atomic_op/type in event)

  // for i in pending[rank]:
  //   sum count or bytes transferred
  // call update bookkeeping from EQ
  //   - once EQ entry is retrieved, it's gone from the EQ
  //   - do all processing on EQ somewhere else, then use most
  //     recent metadata to perform fence

  // check global accounting state for outstanding communication counts/sizes
  // if pending > completed
  //   call handler which calls PtlEQWait and handles a single event (may not be one we're looking for)
}



/**
 * pdht_allfence - ensures completion of all pending put/get operations
 *    
 */
void pdht_allfence(void) {
  // same as above, but handle all ranks
  // only deal with EQ once for all ranks
}



/**
 * pdht_test - checks status of an asynchronous put/get operation
 * @param h handle of pending operation
 * @returns status of operation
 */
pdht_status_t pdht_test(pdht_handle_t h) {
  // handle is index into pending communication array
  // simply return status of entry
  return PdhtStatusOK;
}



/**
 * pdht_wait - blocks process until an asynchronous operation completes
 * @param h handle of pending operation
 * @returns status of operation
 */
pdht_status_t pdht_wait(pdht_handle_t h) {
  // handle is index into pending communication array
  // if status is complete, then return OK
  // if status is pending, then deal with EQ shit until our prince has arrived
  // anybody who wants our EQ handling to be re-entrant should be killed.
  return PdhtStatusOK;
}



/**
 * pdht_waitrank - blocks process until all asynchronous operations complete wrt one process rank
 * @param h handle of pending operation
 * @returns status of operation
 */
pdht_status_t pdht_waitrank(int rank) {
  // same as above, but for all outstanding entries in pending array for a given rank
  return PdhtStatusOK;
}



/**
 * pdht_waitall - blocks process until all asynchronous operations complete 
 * @returns status of operation
 */
pdht_status_t pdht_waitall(void) {
  // same as above, but for all outstanding entries in pending array for all ranks
  return PdhtStatusOK;
}



/**
 * reduce_zip - helper function to accumulate (sum-only) reduced data
 * @param dest  accumulation buffer
 * @param src source data
 * @param ty data type
 * @param s number of elements
 */
static void reduce_zip(void *dest, void *src, pdht_datatype_t ty, int s) {
  int    *is, *id;
  long   *ls, *ld;
  double *ds, *dd;
  char   *cs, *cd;

  switch(ty) {
    case IntType:
      is = src;
      id = dest;
      for (int i=0;i<s;i++) id[i] += is[i];
      break;
    case LongType:
      ls = src;
      ld = dest;
      for (int i=0;i<s;i++) ld[i] += ls[i]; 
      break;
    case DoubleType:
      ds = src;
      dd = dest;
      for (int i=0;i<s;i++) dd[i] += ds[i]; 
      break;
    case CharType:
    case BoolType:
      cs = src;
      cd = dest;
      for (int i=0;i<s;i++) cd[i] += cs[i]; 
      break;
    default:
      pdht_dprintf("pdht_reduce: unsupported reduction datatype\n");
      exit(1);
  }
}



/**
 * pdht_collective_size - helper function to return data type size for collective ops
 * @param type PDHT datatype
 * @returns size of datatype
 */
static int pdht_collective_size(pdht_datatype_t type) {
  int tysize;
  switch(type) {
    case IntType:
      tysize = sizeof(int);
      break;
    case LongType:
      tysize = sizeof(long);
      break;
    case DoubleType:
      tysize = sizeof(double);
      break;
    case CharType:
    case BoolType:
      tysize = sizeof(char);
      break;
    default:
      pdht_dprintf("pdht_reduce: unsupported reduction datatype\n");
      exit(1);
  }
  return tysize;
}

