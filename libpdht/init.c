/********************************************************/
/*                                                      */
/*  init.c - portals distributed hash table init        */
/*                                                      */
/*  author: d. brian larkins                            */
/*  created: 3/20/16                                    */
/*                                                      */
/********************************************************/

#include <pdht.h>

pdht_context_t *c = NULL; // a global variable. for shame.

/**
 * @file
 * 
 * portals distributed hash table startup/shutdown
 */


/**
 * pdht_create -- allocates a new dht
 * @returns the newly minted dht
 */
pdht_t *pdht_create(void) {
  pdht_t *dht;

  if (!c) {
    pdht_init();
  }

  dht = (pdht_t *)malloc(sizeof(pdht_t));
  memset(dht, 0, sizeof(pdht_t));

  dht->ctx = c;
  dht->ctx->dhtcount++; // register ourselves globally on this process

  return dht;
}



/**
 * pdht_free -- frees a new dht
 * @param dht - the dht to free
 */
void dht_free(pdht_t *dht) {
  assert(dht);

  dht->ctx->dhtcount--;
  if (dht->ctx->dhtcount <= 0) {
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


  eprintf("Initializing Portals 4\n");

  ret = PtlInit();
  if (ret != PTL_OK) {
    pdht_dprintf("portals initialization error\n");
    exit(-1);
  }

  ret = PtlNIInit(PTL_IFACE_DEFAULT,
      PTL_NI_NO_MATCHING | PTL_NI_PHYSICAL,
      PTL_PID_ANY, NULL, NULL, &c->ptl.phy);
  if (ret != PTL_OK) {
    pdht_dprintf("Portals physical NI initialization problem. (return=%d)\n", ret);
    exit(-1);
  }

  init_pmi(c);

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


  // create logical NI
  ret = PtlNIInit(PTL_IFACE_DEFAULT,
      PTL_NI_NO_MATCHING | PTL_NI_LOGICAL,
      PTL_PID_ANY,
      &ni_req_limits,
      &c->ptl.ni_limits,
      &c->ptl.lni);
  if (ret != PTL_OK) {
    eprintf("Portals logical NI initialization error\n");
    goto error;
  }

  ret = PtlSetMap(c->ptl.lni, c->size, c->ptl.mapping);
  if (ret != PTL_OK) {
    eprintf("Portals physical/logical mapping failed.\n");
    goto error;
  }

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
  free(c);
  c = NULL;
}



