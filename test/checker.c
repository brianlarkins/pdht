/****************************************************/
/*                                                  */
/*                                                  */
/*                                                  */
/*                                                  */
/*                                                  */
/****************************************************/
#define _XOPEN_SOURCE
#include "cfg.h"

#define SHARED_BUF_SIZE 1000000

// global variables
cfg_t c;

// local protos
static void init(void);
int main(int argc, char **argv);


void init(void) {
  ptl_ni_limits_t ni_req_limits;
  ptl_process_t me;
  int ret;

  printf("Initializing Portals 4\n");
  ret = PtlInit();
  if (ret != PTL_OK) {
    printf("portals initialization error\n");
    exit(-1);
  }

  ret = PtlNIInit(PTL_IFACE_DEFAULT,
      PTL_NI_NO_MATCHING | PTL_NI_PHYSICAL,
      PTL_PID_ANY, NULL, NULL, &c.phy);
  if (ret != PTL_OK) {
    printf("Portals physical NI initialization problem. (return=%d)\n", ret);
    exit(-1);
  }

  init_pmi();

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
      &c.ni_limits,
      &c.lni);
  if (ret != PTL_OK) {
    eprintf("Portals logical NI initialization error\n");
    goto error;
  }

  ret = PtlSetMap(c.lni, c.size, c.mapping);
  if (ret != PTL_OK) {
    eprintf("Portals physical/logical mapping failed.\n");
    goto error;
  }

  return;

error:
  if (c.mapping)
    free(c.mapping);
  exit(-1);
}


void fini(void) {

  PtlNIFini(c.lni);
  PtlNIFini(c.phy);
  if (c.mapping) 
    free(c.mapping);
  PtlFini();
}

/* main */
int main(int argc, char **argv) {
  ptl_md_t pmd;
  ptl_le_t ple;
  ptl_handle_md_t pmdh;
  ptl_handle_md_t pleh;
  ptl_pt_index_t want, pti;
  ptl_process_t them;
  int ret;
  char *shared = NULL;

  putenv("PTL_DISABLE_MEM_REG_CACHE=1");

  init();

  //ret = PtlSetMap(nihdl, 

  printf("%d saying hello out of %d\n", c.rank, c.size);

  want = 0;
  ret = PtlPTAlloc(c.lni, 0, PTL_EQ_NONE, want, &pti);
  if (ret != PTL_OK)
    printf("PtlPTAlloc error: %d\n", ret);
  // unmatched I/O (traditional PGAS)

  // process 0
  // - do portal bullshit to expose memory via portals
  //    - create memory descirptor
  // - write to known memory location in pgas shared space
  // does local need to use PtlGet() to access memory in shared window?
  if (c.rank == 0) {
    them.rank = 1;

    // initialize shared memory buffer
    shared = (char *)malloc(SHARED_BUF_SIZE * sizeof(char));
    for (int i=0; i<SHARED_BUF_SIZE; i++) 
      shared[i] = (15+i)%256;

    // setup non-matching list element
    ple.start     = shared;
    ple.length    = SHARED_BUF_SIZE * sizeof(char);
    ple.uid       = PTL_UID_ANY;
    ple.options   = PTL_LE_OP_PUT | PTL_LE_OP_GET | PTL_LE_EVENT_CT_COMM;

    // allocate counter for tracking completion of put/get events
    ret = PtlCTAlloc(c.lni, &ple.ct_handle);
    if (ret != PTL_OK)
      printf("PtlCTAlloc error: %d\n", ret);

    // append LE to priority list
    ret = PtlLEAppend(c.lni, pti, &ple, PTL_PRIORITY_LIST, NULL, &pleh);
    if (ret != PTL_OK)
      printf("PtlLEAppend error: %d\n", ret);

    printf("%d: sharedbuf[0] contains: %d\n", c.rank, shared[0]);
    PMI_Barrier(); // wait for everybody to get ready

    fflush(stdout);


    PMI_Barrier(); // wait for all done
    printf("%d: sharedbuf[0] contains: %d (post-update)\n", c.rank, shared[0]);

    ret = PtlLEUnlink(pleh);
    if (ret != PTL_OK)
      printf("PtlLEUnlink error: %d\n", ret);
    ret = PtlPTFree(c.lni, pti);
    if (ret != PTL_OK)
      printf("PtlPTFree error: %d\n", ret);

  } else if (c.rank == 1) {
    // process 1
    // - get memory from p0 - print
    // - put memory to p0
    // - get memory from p0 - print
    // matched I/O (use to select from multiple windows?)

    int roffset  = 0;
    int loffset  = 0;
    ptl_hdr_data_t hdr_data = 0xdeadbeeffeedbeef;
    ptl_ct_event_t ctc;

    them.rank = 0;

    // create local MD for data that we get
    char *lbuf  = (char *)malloc(SHARED_BUF_SIZE * sizeof(char));
    // initialize local memory buffer
    for (int i=0; i<SHARED_BUF_SIZE; i++) 
      lbuf[i] = 127; 

    pmd.start   = lbuf;
    pmd.length  = SHARED_BUF_SIZE * sizeof(char);
    pmd.options = PTL_MD_EVENT_CT_REPLY | PTL_MD_EVENT_CT_ACK; // sets counter/event behaviors upon read/writes to this md
    pmd.eq_handle = PTL_EQ_NONE; // pointer to event queue handle

    // allocate counter for tracking completion of put/get events
    ret = PtlCTAlloc(c.lni, &pmd.ct_handle);
    if (ret != PTL_OK)
      printf("PtlCTAlloc error: %d\n", ret);

    // bind local MD
    ret = PtlMDBind(c.lni, &pmd, &pmdh);
    if (ret != PTL_OK)
      printf("PtlMDBind error: %d\n", ret);

    printf("%d: localbuf[0] initially contains: %d\n", c.rank, lbuf[0]);

    PMI_Barrier(); // wait for setup to complete
    fflush(stdout);

    // get memory from p0
    ret = PtlGet(pmdh, loffset, sizeof(char), them, pti, 0, roffset, NULL);
    if (ret != PTL_OK)
      printf("PtlGet error: %d\n", ret);

    ret = PtlCTWait(pmd.ct_handle, 1, &ctc);

    printf("%d: localbuf[0] contains: %d (from get)\n", c.rank, lbuf[0]);
    lbuf[0] -= 1;
    printf("%d: localbuf[0] contains: %d (post local modification)\n", c.rank, lbuf[0]);

    // put memory to p0
    ret = PtlPut(pmdh, 0, sizeof(char), PTL_CT_ACK_REQ, them, pti, 0, roffset, NULL, hdr_data);
    if (ret != PTL_OK)
      printf("PtlPut error: %d\n", ret);
    ret = PtlCTWait(pmd.ct_handle, 2, &ctc);

    // get memory from p0
    ret = PtlGet(pmdh, loffset, sizeof(char), them, pti, 0, roffset, NULL);
    if (ret != PTL_OK)
      printf("PtlGet error: %d\n", ret);
    ret = PtlCTWait(pmd.ct_handle, 3, &ctc);

    PMI_Barrier();  // wait for all done

    printf("%d: localbuf[0] contains: %d (from PtlGet)\n", c.rank, lbuf[0]);

    PtlCTGet(pmd.ct_handle, &ctc);
    printf("%d: successful portals operations: %d\n", c.rank, (int)ctc.success);

    PtlMDRelease(pmdh);
  }



  fini();

  return 0;

error:
  fini();
  return -1;
}
