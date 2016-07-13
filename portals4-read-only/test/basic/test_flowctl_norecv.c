/*
 * Test flow control where we have no receive space
 */

#include <portals4.h>
#include <support.h>

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
# define OPTIONS  (PTL_ME_OP_PUT | PTL_ME_OP_GET | PTL_ME_EVENT_CT_COMM | PTL_ME_NO_TRUNCATE)
# define APPEND   PtlMEAppend
# define UNLINK   PtlMEUnlink
#else
# define ENTRY_T  ptl_le_t
# define HANDLE_T ptl_handle_le_t
# define NI_TYPE  PTL_NI_NO_MATCHING
# define OPTIONS  (PTL_LE_OP_PUT | PTL_LE_OP_GET | PTL_LE_EVENT_CT_COMM | PTL_LE_EVENT_COMM_DISABLE | PTL_LE_USE_ONCE)
# define APPEND   PtlLEAppend
# define UNLINK   PtlLEUnlink
#endif /* if INTERFACE == 1 */

#define ITERS 16

int main(int   argc,
         char *argv[])
{
    ptl_handle_ni_t ni_handle;
    ptl_process_t   *procs;
    int             rank;
    ptl_pt_index_t  pt_index, signal_pt_index;
    ENTRY_T         value_e;
    HANDLE_T        signal_e_handle;
    HANDLE_T        signal_e2_handle;
    int             num_procs;
    ptl_handle_eq_t eq_handle;
    ptl_handle_ct_t ct_handle;
    ptl_handle_md_t md_handle;

    CHECK_RETURNVAL(PtlInit());

    CHECK_RETURNVAL(libtest_init());

    rank = libtest_get_rank();
    num_procs = libtest_get_size();
    if (num_procs < 2) {
        fprintf(stderr, "test_flowctl_noeq requires at least two processes\n");
        return 77;
    }

    int iters;

    if (num_procs < ITERS)
        iters = ITERS*2+1;
    else
        iters = ITERS;

    CHECK_RETURNVAL(PtlNIInit(PTL_IFACE_DEFAULT, NI_TYPE | PTL_NI_LOGICAL,
                              PTL_PID_ANY, NULL, NULL, &ni_handle));
    procs = libtest_get_mapping(ni_handle);
    CHECK_RETURNVAL(PtlSetMap(ni_handle, num_procs, procs));


    if (0 == rank) {
        
        /* create data PT space */
        CHECK_RETURNVAL(PtlEQAlloc(ni_handle, (num_procs - 1) * iters + 64, &eq_handle));
        CHECK_RETURNVAL(PtlPTAlloc(ni_handle, PTL_PT_FLOWCTRL, eq_handle, 5,
                                   &pt_index));

        /* create signal ME */
        CHECK_RETURNVAL(PtlCTAlloc(ni_handle, &ct_handle));
        CHECK_RETURNVAL(PtlPTAlloc(ni_handle, 1, eq_handle, 6,
                                   &signal_pt_index));
        value_e.start = malloc(iters/2);
#if INTERFACE == 1
        value_e.length = iters/2;
#endif
#if INTERFACE == 0
        value_e.length = 1;
#endif
        value_e.ct_handle = ct_handle;
        value_e.uid = PTL_UID_ANY;
        value_e.options = OPTIONS | PTL_LE_EVENT_CT_COMM;
#if INTERFACE == 1
        value_e.min_free = 1;
        value_e.match_id.rank = PTL_RANK_ANY;
        value_e.match_bits = 0;
        value_e.ignore_bits = 0;
#endif
        CHECK_RETURNVAL(APPEND(ni_handle, 5, &value_e, PTL_OVERFLOW_LIST, NULL, &signal_e_handle));
#if INTERFACE == 0
        for (int i=1; i < iters/2; i++) {
            value_e.start += i;
            CHECK_RETURNVAL(APPEND(ni_handle, 5, &value_e, PTL_OVERFLOW_LIST, NULL, &signal_e_handle));
        }
#endif
    } else {
        ptl_md_t        md;

        /* 16 extra just in case... */
        CHECK_RETURNVAL(PtlEQAlloc(ni_handle, iters * 2 + 16, &eq_handle));

        int send_it = 1;
        md.start = &send_it;
        md.length = 1;
        md.options = 0;
        md.eq_handle = eq_handle;
        md.ct_handle = PTL_CT_NONE;

        CHECK_RETURNVAL(PtlMDBind(ni_handle, &md, &md_handle));
    }

    libtest_barrier();

    if (0 == rank) {
        ptl_ct_event_t  ct;
        ptl_event_t ev;
        int ret, count = 0, saw_flowctl = 0;

        /* wait for signal counts */
        CHECK_RETURNVAL(PtlCTWait(ct_handle, iters/2, &ct));
        if (ct.success + ct.failure != iters/2) {
            abort();
        }
        /* wait for event entries */
#if INTERFACE == 1
        while (count+saw_flowctl < iters) {
#endif
#if INTERFACE == 0
        while (count+saw_flowctl < iters) {
#endif
            ret = PtlEQGet(eq_handle, &ev);
            if (PTL_OK == ret) {
                count++;
            } else if (ret == PTL_EQ_EMPTY) {
                continue;  
            } else if (ret == PTL_EVENT_GET_OVERFLOW){
                count++;
                fprintf(stderr,"found EQ entry #%i in overflow list\n",count);
            } else {
                fprintf(stderr, "0: Unexpected return code from EQGet: %d\n", ret);
                return 1;
            }

            if (ev.type == PTL_EVENT_PT_DISABLED) {
                saw_flowctl++;
                break;
            }
        }

        fprintf(stderr, "0: Saw %d flowctl\n", saw_flowctl);
        if (saw_flowctl == 0) {
            abort();
        }

#if INTERFACE == 1
        value_e.start=malloc(iters*30);
        value_e.length=iters*30;
        CHECK_RETURNVAL(APPEND(ni_handle, 5, &value_e, PTL_PRIORITY_LIST, NULL, &signal_e2_handle));
        ret = PTL_OK;
#endif
#if INTERFACE == 0
        for (int i=1; i <= iters/2; i++) {
            CHECK_RETURNVAL(APPEND(ni_handle, 5, &value_e, PTL_PRIORITY_LIST, NULL, &signal_e2_handle));
        }
#endif
        while (ret != PTL_EQ_EMPTY)
            ret = PtlEQGet(eq_handle, &ev);
    } else {
        ptl_process_t target;
        ptl_event_t ev;
        int ret, count = 0, fails = 0;
        int i;

        target.rank = 0;
        for (i = 0 ; i < iters ; ++i) {
            CHECK_RETURNVAL(PtlPut(md_handle,
                                   0,
                                   1,
                                   PTL_ACK_REQ,
                                   target,
                                   5,
                                   0,
                                   i,
                                   NULL,
                                   0));
        }

        while (count < iters) {
            ret = PtlEQGet(eq_handle, &ev);
            if (PTL_EQ_EMPTY == ret) {
                continue;
            } else if (PTL_OK != ret) {
                fprintf(stderr, "%d: PtlEQGet returned %d\n", rank, ret);
                return 1;
            }

            if (ev.ni_fail_type == PTL_NI_OK) {
                if (ev.type == PTL_EVENT_SEND) {
                    continue;
                } else if (ev.type == PTL_EVENT_ACK) {
                    count++;
                } else {
                    fprintf(stderr, "%d: Unexpected event type %d\n", rank, ev.type);
                }
            } else if (ev.ni_fail_type == PTL_NI_PT_DISABLED) {
                count++;
                fails++;
            } else {
                fprintf(stderr, "%d: Unexpected fail type: %d\n", rank, ev.ni_fail_type);
                return 1;
            }
        }

        fprintf(stderr, "%d: Saw %d of %d ACKs as fails\n", rank, fails, count);

    }
    fflush(stderr);
    libtest_barrier();

    if (0 == rank) {
#if INTERFACE == 1
        CHECK_RETURNVAL(UNLINK(signal_e_handle));
        CHECK_RETURNVAL(UNLINK(signal_e2_handle));
#endif
        CHECK_RETURNVAL(PtlPTFree(ni_handle, signal_pt_index));
        CHECK_RETURNVAL(PtlCTFree(ct_handle));
        CHECK_RETURNVAL(PtlPTFree(ni_handle, pt_index));
        CHECK_RETURNVAL(PtlEQFree(eq_handle));
    } else {
        CHECK_RETURNVAL(PtlMDRelease(md_handle));
        CHECK_RETURNVAL(PtlEQFree(eq_handle));
    }

    CHECK_RETURNVAL(PtlNIFini(ni_handle));
    CHECK_RETURNVAL(libtest_fini());
    PtlFini();

    return 0;
}

/* vim:set expandtab: */
