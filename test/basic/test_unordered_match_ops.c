#include <portals4.h>
#include <../support.h>

#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <sched.h>

#include "testing.h"

#if INTERFACE == 1
# define ENTRY_T  ptl_me_t
# define HANDLE_T ptl_handle_me_t
# define NI_TYPE  PTL_NI_MATCHING
# define OPTIONS  (PTL_ME_OP_GET | PTL_ME_EVENT_CT_COMM)
# define APPEND   PtlMEAppend
# define UNLINK   PtlMEUnlink
#else
# define ENTRY_T  ptl_le_t
# define HANDLE_T ptl_handle_le_t
# define NI_TYPE  PTL_NI_NO_MATCHING
# define OPTIONS  (PTL_LE_OP_GET | PTL_LE_EVENT_CT_COMM)
# define APPEND   PtlLEAppend
# define UNLINK   PtlLEUnlink
#endif /* if INTERFACE == 1 */

int main(int   argc,
         char *argv[])
{
    ptl_handle_ni_t ni_logical;
    ptl_process_t   myself;
    ptl_pt_index_t  logical_pt_index;
    uint64_t        value, readval;
    ENTRY_T         value_e;
    HANDLE_T        value_e_handle;
    ptl_md_t        read_md;
    ptl_handle_md_t read_md_handle;
    int             num_procs;

#if INTERFACE != 1
    assert("This test is only valid for MEs" == 0);
#endif

    CHECK_RETURNVAL(PtlInit());

    CHECK_RETURNVAL(libtest_init());

    num_procs = libtest_get_size();

    CHECK_RETURNVAL(PtlNIInit(PTL_IFACE_DEFAULT, NI_TYPE | PTL_NI_LOGICAL,
                              PTL_PID_ANY, NULL, NULL, &ni_logical));

    CHECK_RETURNVAL(PtlSetMap(ni_logical, num_procs,
                              libtest_get_mapping(ni_logical)));

    CHECK_RETURNVAL(PtlGetId(ni_logical, &myself));

    CHECK_RETURNVAL(PtlPTAlloc(ni_logical, PTL_PT_MATCH_UNORDERED, PTL_EQ_NONE, PTL_PT_ANY,
                               &logical_pt_index));
    assert(logical_pt_index == 0);
    /* Now do the initial setup on ni_logical */
    value          = myself.rank + 0xdeadbeefc0d1f1ed;
    value_e.start  = &value;
    value_e.length = sizeof(uint64_t);
    value_e.uid    = PTL_UID_ANY;
#if INTERFACE == 1
    value_e.match_id.rank = PTL_RANK_ANY;
    value_e.match_bits    = 1;
    value_e.ignore_bits   = 0;
#endif
    value_e.options = OPTIONS;
    CHECK_RETURNVAL(PtlCTAlloc(ni_logical, &value_e.ct_handle));
    CHECK_RETURNVAL(APPEND(ni_logical, 0, &value_e, PTL_PRIORITY_LIST, NULL,
                           &value_e_handle));
    /* Now do a barrier (on ni_physical) to make sure that everyone has their
     * logical interface set up */
    libtest_barrier();

    /* now I can communicate between ranks with ni_logical */

    /* set up the landing pad so that I can read others' values */
    read_md.start     = &readval;
    read_md.length    = sizeof(uint64_t);
    read_md.options   = PTL_MD_EVENT_CT_REPLY;
    read_md.eq_handle = PTL_EQ_NONE;   // i.e. don't queue send events
    CHECK_RETURNVAL(PtlCTAlloc(ni_logical, &read_md.ct_handle));
    CHECK_RETURNVAL(PtlMDBind(ni_logical, &read_md, &read_md_handle));

    /* read rank 0's value */
    {
        ptl_ct_event_t ctc;
        ptl_process_t  r0 = { .rank = 0 };
        CHECK_RETURNVAL(PtlGet(read_md_handle, 0, sizeof(uint64_t), r0,
                               logical_pt_index, 1, 0, NULL));
        CHECK_RETURNVAL(PtlCTWait(read_md.ct_handle, 1, &ctc));
        assert(ctc.failure == 0);
    }
    /*printf("%i readval: %llx\n", (int)myself.rank,
     *     (unsigned long long)readval);*/
    assert(readval == 0xdeadbeefc0d1f1ed);
    if (myself.rank == 0) {
        NO_FAILURES(value_e.ct_handle, num_procs);
    }
    CHECK_RETURNVAL(PtlMDRelease(read_md_handle));
    CHECK_RETURNVAL(PtlCTFree(read_md.ct_handle));
    CHECK_RETURNVAL(UNLINK(value_e_handle));
    CHECK_RETURNVAL(PtlCTFree(value_e.ct_handle));

    /* cleanup */
    CHECK_RETURNVAL(PtlPTFree(ni_logical, logical_pt_index));
    CHECK_RETURNVAL(PtlNIFini(ni_logical));
    CHECK_RETURNVAL(libtest_fini());
    PtlFini();

    return 0;
}

/* vim:set expandtab: */
