/********************************************************/
/*                                                      */
/*  init.c - portals distributed hash table init        */
/*                                                      */
/*  author: d. brian larkins                            */
/*  created: 3/20/16                                    */
/*                                                      */
/********************************************************/

#define _XOPEN_SOURCE 700
#include <time.h>
#include <pdht_impl.h>

pdht_context_t *c = NULL; // a global variable. for shame.

static void print_fucking_mapping(void);

/**
 * @file
 * 
 * portals distributed hash table startup/shutdown
 */


/**
 * pdht_create -- allocates a new dht
 * @returns the newly minted dht
 */
pdht_t *pdht_create(int keysize, int elemsize, pdht_mode_t mode) {
  pdht_t *dht;
  ptl_md_t md;
  int ret;

  if (!c) {
    pdht_init();
  }

  dht = (pdht_t *)malloc(sizeof(pdht_t));
  memset(dht, 0, sizeof(pdht_t));

  c->dhtcount++; // register ourselves globally on this process

  dht->keysize = keysize;
  dht->elemsize = elemsize;
  dht->mode = mode;

  // portals info
  dht->ptl.lni = c->ptl.lni;

  // setup polling data structures for pending puts
  pdht_polling_init(dht);

  // allocate event counter
  ret = PtlCTAlloc(dht->ptl.lni, &dht->ptl.lmdct);
  if (ret != PTL_OK) {
    pdht_dprintf("pdht_create: PtlCTAlloc failure\n");
    exit(1);
  }

  // allocate event queue
  ret = PtlEQAlloc(dht->ptl.lni, PDHT_EVENTQ_SIZE, &dht->ptl.lmdeq);
  if (ret != PTL_OK) {
    pdht_dprintf("pdht_create: PtlEQAlloc failure\n");
    exit(1);
  }

  // create memory descriptor (MD) to allow for remote access of our memory
  md.start  = NULL;
  md.length = PTL_SIZE_MAX; 
  md.options = PTL_MD_EVENT_SUCCESS_DISABLE | PTL_MD_EVENT_CT_ACK | PTL_MD_EVENT_CT_SEND | PTL_MD_EVENT_CT_REPLY;
  md.eq_handle = dht->ptl.lmdeq;
  md.ct_handle = dht->ptl.lmdct;

  // bind MD to our matching interface
  ret = PtlMDBind(c->ptl.lni, &md, &dht->ptl.lmd);
  if (ret != PTL_OK) {
    pdht_dprintf("pdht_create: PtlMDBind failure\n");
    exit(1);
  }

  // create PTE for matching gets, will be populated by pending put poller
  ret = PtlPTAlloc(dht->ptl.lni, 0, PTL_EQ_NONE, __PDHT_GET_INDEX, &dht->ptl.getindex);
  if (ret != PTL_OK) {
    pdht_dprintf("pdht_create: PtlPTAlloc failure\n");
    exit(1);
  }

  return dht;
}



/**
 * pdht_free -- frees a new dht
 * @param dht - the dht to free
 */
void pdht_free(pdht_t *dht) {
  assert(dht);
  struct timespec ts;
  int ret;

  c->dhtcount--;

  // disable incoming gets
  PtlPTDisable(dht->ptl.lni, dht->ptl.getindex);

  // cleans up from pending put MEs 
  // -- also removes all MEs from both put/get PTEs
  pdht_polling_fini(dht);

  // free our table entry
  PtlPTFree(dht->ptl.lni, dht->ptl.getindex);

  // free our memory descriptor
  PtlMDRelease(dht->ptl.lmd);

  // release event queue and counter
  PtlCTFree(dht->ptl.lmdct);
  PtlEQFree(dht->ptl.lmdeq);

  // clean up everything if we're last out the door
  if (c->dhtcount <= 0) {
    pdht_fini();
  }

  free(dht);
}



/**
 * pdht_clear - resets a DHT
 * @param dht - the dht to clear
 */
void pdht_clear(pdht_t *dht) {
  // memset the chunk(s) of memory reserved for a DHT
}



/**
 * pdht_init - initializes PDHT system
 */
void pdht_init(void) {
  ptl_ni_limits_t ni_req_limits;
  ptl_process_t me;
  int ret;

  c = (pdht_context_t *)malloc(sizeof(pdht_context_t));
  memset(c,0,sizeof(pdht_context_t));

  ret = PtlInit();
  if (ret != PTL_OK) {
    pdht_dprintf("portals initialization error\n");
    exit(-1);
  }

  ret = PtlNIInit(PTL_IFACE_DEFAULT,
      PTL_NI_NO_MATCHING | PTL_NI_PHYSICAL,
      PTL_PID_ANY, NULL, NULL, &(c->ptl.phy));
  if (ret != PTL_OK) {
    pdht_dprintf("Portals physical NI initialization problem. (return=%d)\n", ret);
    exit(-1);
  }

  init_pmi();

  eprintf("Initializing Portals 4\n");
  eprintf("Initializing Network Interface\n");

  // request portals NI limits
  ni_req_limits.max_entries = 1024;
  ni_req_limits.max_unexpected_headers = 1024;
  ni_req_limits.max_mds = 1024;
  ni_req_limits.max_eqs = 1024;
  ni_req_limits.max_cts = 1024;
  ni_req_limits.max_pt_index = 64;
  ni_req_limits.max_iovecs = 1024;
  ni_req_limits.max_list_size = 1024;
  ni_req_limits.max_triggered_ops = 1024;
  ni_req_limits.max_msg_size = LONG_MAX;
  ni_req_limits.max_atomic_size = 512;
  ni_req_limits.max_fetch_atomic_size = 512;
  ni_req_limits.max_waw_ordered_size = 512;
  ni_req_limits.max_war_ordered_size = 512;
  ni_req_limits.max_volatile_size = 512;
#ifdef PTL_TOTAL_DATA_ORDERING
  eprintf(" (configured for total data ordering on short get/puts)\n");
  ni_req_limits.features = PTL_TOTAL_DATA_ORDERING;
#else
  ni_req_limits.features = 0;
#endif


  // we are not using non-matching right now
#if USE_NON_MATCHING
  // create non-matching logical NI
  ret = PtlNIInit(PTL_IFACE_DEFAULT,
      PTL_NI_NO_MATCHING | PTL_NI_LOGICAL,
      PTL_PID_ANY,
      &ni_req_limits,
      &(c->ptl.ni_limits),
      &(c->ptl.lni));

#else
  // create matching logical NI
  ret = PtlNIInit(PTL_IFACE_DEFAULT,
      PTL_NI_MATCHING | PTL_NI_LOGICAL,
      PTL_PID_ANY,
      &ni_req_limits,
      &(c->ptl.ni_limits),
      &(c->ptl.lni));
#endif // matching or non-matching?

  if (ret != PTL_OK) {
    eprintf("Portals logical NI initialization error\n");
    goto error;
  }

  //print_fucking_mapping();

  ret = PtlSetMap(c->ptl.lni, c->size, c->ptl.mapping);
  if (ret != PTL_OK) {
    eprintf("Portals physical/logical mapping failed : %s.\n", pdht_ptl_error(ret));
    goto error;
  }
  /*
   * have to call one more PMI/MPI barrier to guarantee we have our own counters 
   * ready for our internal barrier operation... 
   */
  pdht_barrier_init(c);
  init_only_barrier(); // safe to use pdht_barrier() after this

  return;

error:
  if (c->ptl.mapping)
    free(c->ptl.mapping);
  exit(-1);
}



/**
 * pdht_fini - initializes PDHT system
 */
void pdht_fini(void) {

  // free up barrier initialization stuff (PT Entry, MD)


  // barrier structures
  PtlMEUnlink(c->ptl.barrier_me);
  PtlCTFree(c->ptl.barrier_ct);
  PtlPTFree(c->ptl.lni, __PDHT_BARRIER_INDEX);
  PtlMDRelease(c->ptl.barrier_md);

  PtlNIFini(c->ptl.lni);
  PtlNIFini(c->ptl.phy);
  if (c->ptl.mapping)
    free(c->ptl.mapping);
  PtlFini();

  free(c);
  c = NULL;
}


static void print_fucking_mapping() {
  int i;
  for (i=0;i<c->size;i++) {
    printf("process %d: map[%u] nid: %u pid: %u\n", c->rank, i, c->ptl.mapping[i].phys.nid, c->ptl.mapping[i].phys.pid);
    fflush(stdout);
  }
}
