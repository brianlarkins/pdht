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

pdht_context_t *c = NULL;       // a global variable. for shame.
pdht_config_t   *__pdht_config = NULL; // another one. i'm over it. (used only during init)

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
  _pdht_ht_entry_t *hte;
  pdht_config_t cfg;
  char *iter; // used for pointer math
  int ret;

  if (!__pdht_config) {
     cfg.nptes      = PDHT_DEFAULT_NUM_PTES;
     cfg.pendmode   = PDHT_DEFAULT_PMODE;
     cfg.maxentries = PDHT_DEFAULT_TABLE_SIZE;
     cfg.pendq_size = PDHT_PENDINGQ_SIZE;
  } else {
    memcpy(&cfg, __pdht_config, sizeof(pdht_config_t));
  }

  if (!c) {
    pdht_init(&cfg);
  }

  dht = (pdht_t *)malloc(sizeof(pdht_t));
  memset(dht, 0, sizeof(pdht_t));

  c->dhtcount++; // register ourselves globally on this process

  if (keysize > PDHT_MAXKEYSIZE) {
    pdht_dprintf("pdht_create: keysize greater than PDHT_MAXKEYSIZE: %d > %d\n", keysize, PDHT_MAXKEYSIZE);
    pdht_dprintf("\t (update value in pdht_impl.h and recompile)\n");
  }

  dht->keysize = keysize;
  dht->elemsize = elemsize;
  dht->maxentries = cfg.maxentries;
  dht->pendq_size = cfg.pendq_size;
  dht->mode = mode;
  dht->pmode = cfg.pendmode;

  dht->hashfn = pdht_hash;

  // portals info
  dht->ptl.nptes = cfg.nptes;
  assert(dht->ptl.nptes < PDHT_MAX_PTES);
  dht->ptl.putindex_base = __PDHT_ACTIVE_INDEX + dht->ptl.nptes;
  dht->ptl.lni           = c->ptl.lni;

  // allocate array for hash table data
  if (dht->pmode == PdhtPendingPoll) 
    dht->entrysize = (sizeof(_pdht_ht_entry_t)) + dht->elemsize;
  else if (dht->pmode == PdhtPendingTrig)
    dht->entrysize = (sizeof(_pdht_ht_trigentry_t)) + dht->elemsize;
  else
    pdht_eprintf(PDHT_DEBUG_NONE, "pdht_create: illegal polling mode\n");

  // print runtime settings
  pdht_eprintf(PDHT_DEBUG_WARN, "pdht_create: hash table entry size: %lu (%d + %d)\n", 
         dht->entrysize, sizeof(_pdht_ht_entry_t), dht->elemsize);
  pdht_eprintf(PDHT_DEBUG_WARN, "\tcontext: %lu bytes ht: %lu bytes table: %lu\n", 
        sizeof(pdht_context_t), sizeof(pdht_t), dht->maxentries * dht->entrysize);
  pdht_eprintf(PDHT_DEBUG_WARN, "\tmax table size: %d pending q size: %d\n", dht->maxentries, dht->pendq_size);
  pdht_eprintf(PDHT_DEBUG_WARN, "\tPT Entries: %d initial pending entries: %d\n", dht->ptl.nptes, dht->ptl.nptes*dht->pendq_size);
  if (dht->pmode == PdhtPendingPoll) {
    pdht_eprintf(PDHT_DEBUG_WARN, "\tpending PTE mode: polling\n");
  } else {
    pdht_eprintf(PDHT_DEBUG_WARN, "\tpending PTE mode: triggered\n");
  }

  dht->ht = calloc(dht->maxentries, dht->entrysize);
  if (!dht->ht) {
    pdht_dprintf("pdht_create: calloc error: %s\n", strerror(errno));
    exit(1);
  }

  // use a byte pointer as iterator over variable-sized element array
  iter = (char *)dht->ht;

  // initialize entire hash table array
  for (int i=0; i < dht->maxentries; i++) {
    hte = (_pdht_ht_entry_t *)iter;
    hte->pme = PTL_INVALID_HANDLE; // initialize pending put ME as invalid
    hte->ame = PTL_INVALID_HANDLE; // initialize active ME as invalid
    iter += dht->entrysize; // pointer math, danger.
  }

  // setup data structures for pending puts
  if (dht->pmode == PdhtPendingPoll) {
    pdht_polling_init(dht);
  } else {
    pdht_trig_init(dht);
  }

  // allocate event counter
  ret = PtlCTAlloc(dht->ptl.lni, &dht->ptl.lmdct);
  if (ret != PTL_OK) {
    pdht_dprintf("pdht_create: PtlCTAlloc failure\n");
    exit(1);
  }

  // allocate event queue
  ret = PtlEQAlloc(dht->ptl.lni, dht->pendq_size, &dht->ptl.lmdeq);
  if (ret != PTL_OK) {
    pdht_dprintf("pdht_create: PtlEQAlloc failure\n");
    exit(1);
  }

  // create memory descriptor (MD) to allow for remote access of our memory
  md.start  = NULL;
  md.length = PTL_SIZE_MAX; 
  md.options = PTL_MD_EVENT_SUCCESS_DISABLE | PTL_MD_EVENT_CT_ACK | PTL_MD_EVENT_CT_REPLY;
  md.eq_handle = dht->ptl.lmdeq;
  md.ct_handle = dht->ptl.lmdct;

  // bind MD to our matching interface
  ret = PtlMDBind(c->ptl.lni, &md, &dht->ptl.lmd);
  if (ret != PTL_OK) {
    pdht_dprintf("pdht_create: PtlMDBind failure\n");
    exit(1);
  }

  for (int ptindex=0; ptindex < dht->ptl.nptes; ptindex++) {
    // create PTE for matching gets, will be populated by pending put poller
    ret = PtlPTAlloc(dht->ptl.lni, dht->ptl.ptalloc_opts, PTL_EQ_NONE, __PDHT_ACTIVE_INDEX+ptindex, &dht->ptl.getindex[ptindex]);
    if (ret != PTL_OK) {
      pdht_dprintf("pdht_create: PtlPTAlloc failure [%d]\n", ptindex);
      exit(1);
    }
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
  for (int ptindex=0; ptindex < dht->ptl.nptes; ptindex++) 
    PtlPTDisable(dht->ptl.lni, dht->ptl.getindex[ptindex]);

  // cleans up from pending put MEs 
  // -- also removes all MEs from both put/get PTEs
  pdht_polling_fini(dht);

  // free our table entries
  for (int ptindex=0; ptindex < dht->ptl.nptes; ptindex++) 
    PtlPTFree(dht->ptl.lni, dht->ptl.getindex[ptindex]);

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
 ** pdht_tune - sets tunable parameters for PDHT
 * @param opts - bit flags marking modified parameters
 * @param config - config tunable structure
 */
void pdht_tune(unsigned opts, pdht_config_t *config) {
  if (!config)
     return;

  if (!__pdht_config) {
     __pdht_config = (pdht_config_t *)calloc(1, sizeof(pdht_config_t));
     __pdht_config->nptes        = PDHT_DEFAULT_NUM_PTES;
     __pdht_config->pendmode     = PDHT_DEFAULT_PMODE;
     __pdht_config->maxentries   = PDHT_DEFAULT_TABLE_SIZE;
     __pdht_config->pendq_size   = PDHT_PENDINGQ_SIZE;
     __pdht_config->ptalloc_opts = PDHT_PTALLOC_OPTIONS;
  }
  if (opts & PDHT_TUNE_NPTES) 
    __pdht_config->nptes        = config->nptes;
  if (opts & PDHT_TUNE_PMODE) 
    __pdht_config->pendmode     = config->pendmode;
  if (opts & PDHT_TUNE_ENTRY) 
    __pdht_config->maxentries   = config->maxentries;
  if (opts & PDHT_TUNE_PENDQ)
    __pdht_config->pendq_size   = config->pendq_size;
  if (opts & PDHT_TUNE_PTOPT)
    __pdht_config->ptalloc_opts = config->ptalloc_opts;

  // copy back tunables, so app can see
  memcpy(config,__pdht_config, sizeof(pdht_config_t));
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
void pdht_init(pdht_config_t *cfg) {
  ptl_ni_limits_t ni_req_limits;
  ptl_process_t me;
  int ret;

  // turn off output buffering for everyone's sanity
  setbuf(stdout, NULL);

  c = (pdht_context_t *)malloc(sizeof(pdht_context_t));
  memset(c,0,sizeof(pdht_context_t));

  ret = PtlInit();
  if (ret != PTL_OK) {
    pdht_dprintf("portals initialization error\n");
    exit(-1);
  }
#if 0
  ret = PtlNIInit(PTL_IFACE_DEFAULT,
      PTL_NI_NO_MATCHING | PTL_NI_PHYSICAL,
      PTL_PID_ANY, NULL, NULL, &(c->ptl.phy));
  if (ret != PTL_OK) {
    pdht_dprintf("Portals physical NI initialization problem. (return=%d)\n", ret);
    exit(-1);
  }
#endif

  // request portals NI limits
  ni_req_limits.max_entries = cfg->maxentries;
  ni_req_limits.max_unexpected_headers = 1024;
  ni_req_limits.max_mds = 1024;
  ni_req_limits.max_eqs = (cfg->nptes)+2;
  ni_req_limits.max_cts = (cfg->nptes*cfg->pendq_size)+2;
  //ni_req_limits.max_eqs = PDHT_DEFAULT_TABLE_SIZE;
  //ni_req_limits.max_cts = PDHT_DEFAULT_TABLE_SIZE;
  //ni_req_limits.max_pt_index = 64;
  ni_req_limits.max_pt_index = 2*cfg->nptes + 1;
  ni_req_limits.max_iovecs = 1024;
  ni_req_limits.max_list_size = cfg->maxentries;
  ni_req_limits.max_triggered_ops = (cfg->nptes*cfg->pendq_size)+100;
  ni_req_limits.max_msg_size = LONG_MAX;
  ni_req_limits.max_atomic_size = 512;
  ni_req_limits.max_fetch_atomic_size = 512;
  ni_req_limits.max_waw_ordered_size = 512;
  ni_req_limits.max_war_ordered_size = 512;
  ni_req_limits.max_volatile_size = 512;
#undef PDHT_WANT_DATA_ORDERING
  //#ifdef PTL_TOTAL_DATA_ORDERING
#ifdef PDHT_WANT_DATA_ORDERING
  pdht_eprintf(PDHT_DEBUG_WARN, " (configured for total data ordering on short get/puts)\n");
  ni_req_limits.features = PTL_TOTAL_DATA_ORDERING;
#else
  ni_req_limits.features = 0;
#endif

  // create matching logical NI
  ret = PtlNIInit(PTL_IFACE_DEFAULT,
      PTL_NI_MATCHING | PTL_NI_LOGICAL,
      PTL_PID_ANY,
      &ni_req_limits,
      &(c->ptl.ni_limits),
      &(c->ptl.lni));

  init_pmi();

  pdht_eprintf(PDHT_DEBUG_WARN, "Initializing Portals 4\n");
  pdht_eprintf(PDHT_DEBUG_WARN, "Initializing Network Interface\n");

  if (ret != PTL_OK) {
    pdht_eprintf(PDHT_DEBUG_NONE, "Portals logical NI initialization error\n");
    goto error;
  }


  c->dbglvl = PDHT_DEBUG_WARN;

  pdht_eprintf(PDHT_DEBUG_WARN, "\tmax_entries: %d\n", c->ptl.ni_limits.max_entries);
  pdht_eprintf(PDHT_DEBUG_WARN, "\tmax_unexpected_headers: %d\n", c->ptl.ni_limits.max_unexpected_headers);
  pdht_eprintf(PDHT_DEBUG_WARN, "\tmax_mds: %d\n", c->ptl.ni_limits.max_mds);
  pdht_eprintf(PDHT_DEBUG_WARN, "\tmax_eqs: %d\n", c->ptl.ni_limits.max_eqs);
  pdht_eprintf(PDHT_DEBUG_WARN, "\tmax_cts: %d\n", c->ptl.ni_limits.max_cts);
  pdht_eprintf(PDHT_DEBUG_WARN, "\tmax_pt_index: %d\n", c->ptl.ni_limits.max_pt_index);
  pdht_eprintf(PDHT_DEBUG_WARN, "\tmax_iovecs: %d\n", c->ptl.ni_limits.max_iovecs);
  pdht_eprintf(PDHT_DEBUG_WARN, "\tmax_list_size: %d\n", c->ptl.ni_limits.max_list_size);
  pdht_eprintf(PDHT_DEBUG_WARN, "\tmax_triggered_ops: %d\n", c->ptl.ni_limits.max_triggered_ops);
  pdht_eprintf(PDHT_DEBUG_WARN, "\tmax_msg_size: %d\n", c->ptl.ni_limits.max_msg_size);
  pdht_eprintf(PDHT_DEBUG_WARN, "\tmax_atomic_size: %d\n", c->ptl.ni_limits.max_atomic_size);
  pdht_eprintf(PDHT_DEBUG_WARN, "\tmax_fetch_atomic_size: %d\n", c->ptl.ni_limits.max_fetch_atomic_size);
  pdht_eprintf(PDHT_DEBUG_WARN, "\tmax_waw_ordered_size: %d\n", c->ptl.ni_limits.max_waw_ordered_size);
  pdht_eprintf(PDHT_DEBUG_WARN, "\tmax_war_ordered_size: %d\n", c->ptl.ni_limits.max_war_ordered_size);
  pdht_eprintf(PDHT_DEBUG_WARN, "\tmax_volatile_size: %d\n", c->ptl.ni_limits.max_volatile_size);

  //print_fucking_mapping();

  ret = PtlSetMap(c->ptl.lni, c->size, c->ptl.mapping);
  if (ret != PTL_OK) {
    pdht_eprintf(PDHT_DEBUG_NONE, "Portals physical/logical mapping failed : %s.\n", pdht_ptl_error(ret));
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
  //PtlNIFini(c->ptl.phy);
  if (c->ptl.mapping)
    free(c->ptl.mapping);
  PtlFini();

  free(c);
  c = NULL;
  if (__pdht_config)
    free(__pdht_config);
}


static void print_fucking_mapping() {
  int i;
  for (i=0;i<c->size;i++) {
    printf("process %d: map[%u] nid: %u pid: %u\n", c->rank, i, c->ptl.mapping[i].phys.nid, c->ptl.mapping[i].phys.pid);
    fflush(stdout);
  }
}
