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

static int barrier_landing = 0;

static void reduce_zip(void *dest, void *src, pdht_reduceop_t op, pdht_datatype_t ty, int s);
static void reduce_sum(void *dest, void *src, pdht_datatype_t ty, int s);
static void reduce_min(void *dest, void *src, pdht_datatype_t ty, int s);
static void reduce_max(void *dest, void *src, pdht_datatype_t ty, int s);
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
  c->ptl.barrier_count    = 0;
  c->ptl.bcast_count      = 0;
  c->ptl.reduce_count     = 0;

  // allocate some scratch place for reductions (assumes longs/doubles are 64-bits each)
  c->ptl.collective_lscratch = calloc(3*PDHT_MAX_REDUCE_ELEMS, sizeof(uint64_t));
  // set start to middle of scratch space (danger, pointer math)
  c->ptl.collective_rscratch = ((char *)(c->ptl.collective_lscratch)) + (PDHT_MAX_REDUCE_ELEMS * sizeof(uint64_t));
  c->ptl.collective_scratch = ((char *)(c->ptl.collective_rscratch)) + (PDHT_MAX_REDUCE_ELEMS * sizeof(uint64_t));

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
  c->ptl.barrier_ct    = PTL_INVALID_HANDLE;
  c->ptl.bcast_ct      = PTL_INVALID_HANDLE;
  c->ptl.reduce_ct     = PTL_INVALID_HANDLE;

  // get counter to attach to match entry
  ret = PtlCTAlloc(c->ptl.lni, &c->ptl.bcast_ct);
  if (ret != PTL_OK)  {
    pdht_dprintf("collective_init: CTAlloc failed\n");
    exit(1);
  } 

  // get another for reductions
  ret = PtlCTAlloc(c->ptl.lni, &c->ptl.reduce_ct);
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
  me.start      = &barrier_landing;
  me.length     = sizeof(barrier_landing);
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

  // allocate matchlist entry for broadcasts
  me.start      = c->ptl.collective_lscratch; // re-use common scratch space
  me.length     = sizeof(uint64_t)*PDHT_MAX_REDUCE_ELEMS;
  me.ct_handle  = c->ptl.bcast_ct;
  me.match_bits  = __PDHT_BCAST_MATCH;
  ret = PtlMEAppend(c->ptl.lni, __PDHT_COLLECTIVE_INDEX, &me, PTL_PRIORITY_LIST, NULL, &c->ptl.bcast_me);
  if (ret != PTL_OK)  {
    pdht_dprintf("collective_init: PtlMEAppend failed\n");
    exit(1);
  } 

  // allocate left matchlist entry for reductions
  me.start      = c->ptl.collective_lscratch;
  me.length     = sizeof(uint64_t)*PDHT_MAX_REDUCE_ELEMS;
  me.ct_handle  = c->ptl.reduce_ct;
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
  PtlCTFree(c->ptl.bcast_ct);
  PtlCTFree(c->ptl.reduce_ct);
  PtlCTFree(c->ptl.barrier_ct);
  PtlPTFree(c->ptl.lni, __PDHT_COLLECTIVE_INDEX);
  PtlMDRelease(c->ptl.collective_md);
}



/**
 * pdht_barrier - blocks process until all processes reach the barrier
 */ 
void pdht_barrier(void) {
  ptl_process_t  p, l, r;
  int nchildren = 0;
  int nparent   = (c->rank != 0) ? 1 : 0; // root has no parent
  ptl_size_t     count_base;
  ptl_ct_event_t cval;
  int ret; 
  int x = 1;

  // use binary tree of processes
  p.rank = ((c->rank + 1) >> 1) - 1;
  l.rank = ((c->rank + 1) << 1) - 1;
  r.rank = l.rank + 1;

  //PtlCTGet(c->ptl.barrier_ct, &cval2);
  //pdht_dprintf("**** p: %lu l: %lu r: %lu barrier_count: %d counter: %lu\n", p.rank, l.rank, r.rank, c->ptl.barrier_count, cval2.success);
  

  if (l.rank < c->size)
    nchildren++;
  if (r.rank < c->size)
    nchildren++;

  // barrier_count tracks barrier events, the number of received messages is related to
  // our tree connectivity
  count_base = c->ptl.barrier_count * (nparent + nchildren);

  // if we have children, wait for them to send us a barrier entry message
  if (nchildren > 0) {
    //pdht_dprintf("waiting for %d messages from l: %lu r: %lu \n", nchildren, l.rank, r.rank);
    ret = PtlCTWait(c->ptl.barrier_ct, count_base + nchildren, &cval);
    if (ret != PTL_OK) {
      pdht_dprintf("pdht_barrier: CTWait failed (children)\n");
      exit(1);
    }

    if (cval.failure > 0) {
      pdht_dprintf("pdht_barrier: found failure event waiting on child messages\n");
      exit(1);
    }
  }


  // if we are a child, send a message to our parent (our subtree has all arrived at barrier)
  if (c->rank > 0)  {

    //pdht_dprintf("notifying %lu : ct: %d count: %d\n", p.rank, cval2.success, c->ptl.barrier_count);
    ret = PtlPut(c->ptl.collective_md, (ptl_size_t)&x, sizeof(x), PTL_NO_ACK_REQ, p, __PDHT_COLLECTIVE_INDEX,
        __PDHT_BARRIER_MATCH, 0, NULL, 0);
    if (ret != PTL_OK) {
      pdht_dprintf("barrier: Put failed (put parent): %s\n", pdht_ptl_error(ret));
      exit(1);
    }

    // wait for downward broadcast from parent
    ret = PtlCTWait(c->ptl.barrier_ct, count_base + (nchildren + nparent), &cval);
    if (ret != PTL_OK) {
      pdht_dprintf("barrier: CTWait failed (parent)\n");
      exit(1);
    }
    if (cval.failure > 0) {
      pdht_dprintf("pdht_barrier: found failure event waiting on child messages\n");
      exit(1);
    }
  }

  // wake up waiting children
  if (nchildren > 0) {
    //pdht_dprintf("waking left: %lu\n", l.rank);
    ret = PtlPut(c->ptl.collective_md, (ptl_size_t)&x, sizeof(x), PTL_NO_ACK_REQ, l, __PDHT_COLLECTIVE_INDEX,
        __PDHT_BARRIER_MATCH, 0, NULL, 0);
    if (ret != PTL_OK) {
      pdht_dprintf("barrier: Put failed (wake left): %s\n", pdht_ptl_error(ret));
      exit(1);
    }
  }

  if (nchildren > 1) {
    //pdht_dprintf("waking right: %lu\n", l.rank);
    ret = PtlPut(c->ptl.collective_md, (ptl_size_t)&x, sizeof(x), PTL_NO_ACK_REQ, r, __PDHT_COLLECTIVE_INDEX,
        __PDHT_BARRIER_MATCH, 0, NULL, 0);
    if (ret != PTL_OK) {
      pdht_dprintf("barrier: Put failed (wake right): %s\n", pdht_ptl_error(ret));
      exit(1);
    }
  }

  // keep track of how many barriers have occurred
  c->ptl.barrier_count++;
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
  int nchildren = 0;
  int nparent   = (c->rank != 0) ? 1 : 0; // root has no parent
  ptl_size_t     count_base;
  ptl_ct_event_t cval;
  ptl_match_bits_t destme;
  int ret = 0, tysize = 0, x = 0;

  switch (op) {
  case PdhtReduceOpSum:
  case PdhtReduceOpMin:
  case PdhtReduceOpMax:
    break;
  default:
    pdht_dprintf("pdht_reduce: unsupported reduction operation\n");
    exit(1);
  }

  tysize = pdht_collective_size(type);

  // figure out if we are the left or right child of our parent
  destme = ((c->rank % 2) != 0) ? __PDHT_REDUCE_LMATCH : __PDHT_REDUCE_RMATCH;

  // use binary tree of processes
  p.rank = ((c->rank + 1) >> 1) - 1;
  l.rank = ((c->rank + 1) << 1) - 1;
  r.rank = l.rank + 1;

  if (l.rank < c->size)
    nchildren++;
  if (r.rank < c->size)
    nchildren++;

  // reduce_count tracks reduce events, the number of received messages is related to
  // our tree connectivity
  count_base = c->ptl.reduce_count * (nparent + nchildren);

  // if we have children, wait for them to send us reduce message
  if (nchildren > 0) {
    //pdht_dprintf("pdht_reduce: waiting for %d messages from l: %lu r: %lu \n", nchildren, l.rank, r.rank);
    ret = PtlCTWait(c->ptl.reduce_ct, count_base + nchildren, &cval);
    if (ret != PTL_OK) {
      pdht_dprintf("pdht_reduce: CTWait failed (children)\n");
      exit(1);
    } else if (cval.failure > 0) {
      pdht_dprintf("pdht_reduce: found failure event waiting on child messages\n");
      exit(1);
    }
  }

	memcpy(c->ptl.collective_scratch, in, (tysize*elems));
  if (nchildren > 0) {
    reduce_zip(in, c->ptl.collective_lscratch, op, type, elems);
    if (nchildren == 2) 
      reduce_zip(in, c->ptl.collective_rscratch, op, type, elems);
  }

  // send our reduced buffer to our parent
  if (c->rank > 0)  {
    //pdht_dprintf("notifying %lu : ct: %d count: %d\n", p.rank, cval2.success, c->ptl.barrier_count);
    ret = PtlPut(c->ptl.collective_md, (ptl_size_t)c->ptl.collective_scratch, tysize*elems, PTL_NO_ACK_REQ, p, __PDHT_COLLECTIVE_INDEX,
        destme, 0, NULL, 0);
    if (ret != PTL_OK) {
      pdht_dprintf("pdht_reduce: Put failed (put parent): %s\n", pdht_ptl_error(ret));
      exit(1);
    }

    // wait for downward broadcast from parent
    ret = PtlCTWait(c->ptl.reduce_ct, count_base + (nchildren + nparent), &cval);
    if (ret != PTL_OK) {
      pdht_dprintf("barrier: CTWait failed (parent)\n");
      exit(1);
    }
    if (cval.failure > 0) {
      pdht_dprintf("pdht_barrier: found failure event waiting on child messages\n");
      exit(1);
    }
  } else {
    // we're the root our local data has been zipped with each l/r subtree
    memcpy(out, c->ptl.collective_scratch, tysize*elems);
  }

  // wake up waiting children
  if (nchildren > 0) {
    //pdht_dprintf("waking left: %lu\n", l.rank);
    ret = PtlPut(c->ptl.collective_md, (ptl_size_t)&x, sizeof(x), PTL_NO_ACK_REQ, l, __PDHT_COLLECTIVE_INDEX,
        __PDHT_REDUCE_LMATCH, 0, NULL, 0);
    if (ret != PTL_OK) {
      pdht_dprintf("barrier: Put failed (wake left): %s\n", pdht_ptl_error(ret));
      exit(1);
    }
  }

  if (nchildren > 1) {
    //pdht_dprintf("waking right: %lu\n", l.rank);
    ret = PtlPut(c->ptl.collective_md, (ptl_size_t)&x, sizeof(x), PTL_NO_ACK_REQ, r, __PDHT_COLLECTIVE_INDEX,
        __PDHT_REDUCE_LMATCH, 0, NULL, 0);
    if (ret != PTL_OK) {
      pdht_dprintf("barrier: Put failed (wake right): %s\n", pdht_ptl_error(ret));
      exit(1);
    }
  }

  // keep track of how many reductions have occurred
  c->ptl.reduce_count++;

  return PdhtStatusOK;
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
  ptl_size_t     bcast_base;
  ptl_ct_event_t cval;
  int nparent   = (c->rank != 0) ? 1 : 0; // root has no parent
  int tysize = 0, nchildren = 0, ret, x;

  // use binary tree of processes
  p.rank = ((c->rank + 1) >> 1) -1;
  l.rank = ((c->rank + 1) << 1) - 1;
  r.rank = l.rank + 1;

  tysize = pdht_collective_size(type);

  if (l.rank < c->size)
    nchildren++;
  if (r.rank < c->size)
    nchildren++;

  // bcast_count tracks broadcast events, compute the right number of
  //    expected messages for our counter threshold
  bcast_base = c->ptl.bcast_count * (nparent + nchildren);

  // wait for parent to send message
  if (c->rank > 0) {
    ret = PtlCTWait(c->ptl.bcast_ct, bcast_base + nparent, &cval);
    if (ret != PTL_OK) {
      pdht_dprintf("pdht_broadcast: CTWait failed (parent)\n");
      exit(1);
    }
    // copy recieved data to our buffer
    memcpy(buf, c->ptl.collective_lscratch, tysize*elems);
  }

  // parents send the message to children (if any)
  if (nchildren > 0) {
    ret = PtlPut(c->ptl.collective_md, (ptl_size_t)buf, tysize*elems, PTL_NO_ACK_REQ, l,
                 __PDHT_COLLECTIVE_INDEX, __PDHT_BCAST_MATCH, 0, NULL, 0);
    if (ret != PTL_OK) {
      pdht_dprintf("pdht_broadcast: Put failed (send left)\n");
      exit(1);
    }
  }

  if (nchildren > 1) {
    ret = PtlPut(c->ptl.collective_md, (ptl_size_t)buf, tysize*elems, PTL_NO_ACK_REQ, r,
                 __PDHT_COLLECTIVE_INDEX, __PDHT_BCAST_MATCH, 0, NULL, 0);
    if (ret != PTL_OK) {
      pdht_dprintf("pdht_broadcast: Put failed (send right)\n");
      exit(1);
    }
  }

  // now we need to wait for confirmation from our children
  if (nchildren != 0)  {
    ret = PtlCTWait(c->ptl.bcast_ct, bcast_base + (nparent + nchildren), &cval);
    if (ret != PTL_OK) {
      pdht_dprintf("pdht_broadcast: CTWait failed (parent)\n");
      exit(1);
    }
  }

  // notify our parent that we've received the broadcast
  if (c->rank > 0) {
    x = 0;
    ret = PtlPut(c->ptl.collective_md, (ptl_size_t)&x, sizeof(x), PTL_NO_ACK_REQ, p, __PDHT_COLLECTIVE_INDEX,
        __PDHT_BCAST_MATCH, 0, NULL, 0);
    if (ret != PTL_OK) {
      pdht_dprintf("pdht_broadcast: Put failed (parent confirmation)\n");


      exit(1);
    }
  }

  c->ptl.bcast_count++;
  return PdhtStatusOK;
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
 * pdht_fence - ensures completion of put/get operations
 * @param dht hash table
 */
void pdht_fence(pdht_t *dht) {
  int ret;
  int sbuf[2], rbuf[2];


  do {
    // count up link events from appending things to the active queue
    pthread_mutex_lock(&dht->completion_mutex);
    pdht_finalize_puts(dht);
    pthread_mutex_unlock(&dht->completion_mutex);

    sbuf[0] = dht->stats.pendputs;
    sbuf[1] = dht->stats.appends;
    pdht_allreduce(sbuf, rbuf, PdhtReduceOpSum, IntType, 2);
    
    //pdht_eprintf(PDHT_DEBUG_NONE, "expected: %d actual: %d\n", rbuf[0], rbuf[1]);
  } while (rbuf[0] > rbuf[1]);

  // reset all the pending counters
  dht->stats.pendputs = 0;
  dht->stats.appends = 0;
}



/**
 * pdht_finalize_puts - tally completion of put operations for fence
 * @param dht a hash table
 * @returns status of operation
 */
pdht_status_t pdht_finalize_puts(pdht_t *dht) {
  int ret = PTL_OK;
  ptl_event_t ev;
  unsigned int ptindex;
  ptl_time_t timeout;

  // NOTE: this function may be called by the polling/trigger progress threads 
  // _or_ during the fence operation _or_ in get when something  is local. 
  // it should be protected by a mutex in the caller

  timeout = 10; // only block the polling progress thread or fence operation for up to 10 ms
     
  while ((ret = PtlEQPoll(dht->ptl.aeq,dht->ptl.nptes, timeout, &ev, &ptindex)) == PTL_OK) {
    //pdht_dump_event(&ev);
    if (ev.type == PTL_EVENT_LINK) {
      dht->stats.appends++;
      dht->stats.tappends[ptindex]++;
    }
    else if (ev.type == PTL_EVENT_SEARCH){
      pthread_mutex_lock(&dht->local_gets_flag_mutex);
      dht->local_get_flag = 1;
      if(ev.ni_fail_type == PTL_NI_NO_MATCH){
        dht->local_get_flag = -1;
      }
      else{
        *(ptl_event_t **)ev.user_ptr = ev.start;
      }
      pthread_mutex_unlock(&dht->local_gets_flag_mutex);
      return ret;
    }
  }
  if ((ret != PTL_EQ_EMPTY) && (ret != PTL_INTERRUPTED)) {
    printf("pdht_finalize_puts: PtlEQPoll error %s\n", pdht_ptl_error(ret));
  }
  return ret;
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
static void reduce_zip(void *dest, void *src, pdht_reduceop_t op, pdht_datatype_t ty, int s) {

  switch(op) {
  case PdhtReduceOpSum:
    reduce_sum(dest, src, ty, s);
    break;
  case PdhtReduceOpMin:
    reduce_min(dest, src, ty, s);
    break;
  case PdhtReduceOpMax:
    reduce_max(dest, src, ty, s);
    break;
  }
}



/**
 * reduce_sum - helper function to accumulate (sum-only) reduced data
 * @param dest  accumulation buffer
 * @param src source data
 * @param ty data type
 * @param s number of elements
 */
static void reduce_sum(void *dest, void *src, pdht_datatype_t ty, int s) {
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
 * reduce_min - helper function to accumulate reduced data
 * @param dest  accumulation buffer
 * @param src source data
 * @param ty data type
 * @param s number of elements
 */
static void reduce_min(void *dest, void *src, pdht_datatype_t ty, int s) {
  int    *is, *id;
  long   *ls, *ld;
  double *ds, *dd;
  char   *cs, *cd;

  switch(ty) {
    case IntType:
      is = src;
      id = dest;
      for (int i=0;i<s;i++) id[i] = id[i] < is[i] ? id[i] : is[i];
      break;
    case LongType:
      ls = src;
      ld = dest;
      for (int i=0;i<s;i++) ld[i] = ld[i] < ls[i] ? ld[i] : ls[i];
      break;
    case DoubleType:
      ds = src;
      dd = dest;
      for (int i=0;i<s;i++) dd[i] = dd[i] < ds[i] ? dd[i] : ds[i];
      break;
    case CharType:
    case BoolType:
      cs = src;
      cd = dest;
      for (int i=0;i<s;i++) cd[i] = cd[i] < cs[i] ? cd[i] : cs[i];
      break;
    default:
      pdht_dprintf("pdht_reduce: unsupported reduction datatype\n");
      exit(1);
  }
}



/**
 * reduce_max - helper function to accumulate reduced data
 * @param dest  accumulation buffer
 * @param src source data
 * @param ty data type
 * @param s number of elements
 */
static void reduce_max(void *dest, void *src, pdht_datatype_t ty, int s) {
  int    *is, *id;
  long   *ls, *ld;
  double *ds, *dd;
  char   *cs, *cd;

  switch(ty) {
    case IntType:
      is = src;
      id = dest;
      for (int i=0;i<s;i++) id[i] = id[i] > is[i] ? id[i] : is[i];
      break;
    case LongType:
      ls = src;
      ld = dest;
      for (int i=0;i<s;i++) ld[i] = ld[i] > ls[i] ? ld[i] : ls[i];
      break;
    case DoubleType:
      ds = src;
      dd = dest;
      for (int i=0;i<s;i++) dd[i] = dd[i] > ds[i] ? dd[i] : ds[i];
      break;
    case CharType:
    case BoolType:
      cs = src;
      cd = dest;
      for (int i=0;i<s;i++) cd[i] = cd[i] > cs[i] ? cd[i] : cs[i];
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

