#include <portals4.h>
#include <support.h>

#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>

#define LOOPS 1000000

#define CHECK_RETURNVAL(x) do {                                                                                                               \
        int ret;                                                                                                                              \
        switch (ret = x) {                                                                                                                    \
            case PTL_IGNORED: case PTL_OK: break;                                                                                             \
            case PTL_FAIL: fprintf(stderr, "=> %s returned PTL_FAIL (line %u)\n", # x, (unsigned int)__LINE__); abort(); break;               \
            case PTL_NO_SPACE: fprintf(stderr, "=> %s returned PTL_NO_SPACE (line %u)\n", # x, (unsigned int)__LINE__); abort(); break;       \
            case PTL_ARG_INVALID: fprintf(stderr, "=> %s returned PTL_ARG_INVALID (line %u)\n", # x, (unsigned int)__LINE__); abort(); break; \
            case PTL_NO_INIT: fprintf(stderr, "=> %s returned PTL_NO_INIT (line %u)\n", # x, (unsigned int)__LINE__); abort(); break;         \
            default: fprintf(stderr, "=> %s returned failcode %i (line %u)\n", # x, ret, (unsigned int)__LINE__); abort(); break;             \
        } } while (0)

#if INTERFACE == 1
# define ENTRY_T  ptl_me_t
# define HANDLE_T ptl_handle_me_t
# define NI_TYPE  PTL_NI_MATCHING
# define OPTIONS  (PTL_ME_OP_PUT)
# define APPEND   PtlMEAppend
# define UNLINK   PtlMEUnlink
#else
# define ENTRY_T  ptl_le_t
# define HANDLE_T ptl_handle_le_t
# define NI_TYPE  PTL_NI_NO_MATCHING
# define OPTIONS  (PTL_LE_OP_PUT)
# define APPEND   PtlLEAppend
# define UNLINK   PtlLEUnlink
#endif /* if INTERFACE == 1 */

int main(int   argc,
         char *argv[])
{
    ptl_handle_ni_t ni_logical;
    ptl_pt_index_t  logical_pt_index;
    ptl_process_t   myself;
    struct timeval  start, stop;
    int             potato = 0;
    ENTRY_T         potato_catcher;
    HANDLE_T        potato_catcher_handle;
    ptl_md_t        potato_launcher;
    ptl_handle_md_t potato_launcher_handle;
    int             num_procs;
    ptl_handle_eq_t pt_eq_handle;

    CHECK_RETURNVAL(PtlInit());

    CHECK_RETURNVAL(libtest_init());

    num_procs = libtest_get_size();

    if (NULL != getenv("MAKELEVEL") && num_procs > 2) {
        return 77;
    }

    CHECK_RETURNVAL(PtlNIInit
                        (PTL_IFACE_DEFAULT, NI_TYPE | PTL_NI_LOGICAL, PTL_PID_ANY,
                        NULL, NULL, &ni_logical));

    CHECK_RETURNVAL(PtlSetMap(ni_logical, num_procs,
                              libtest_get_mapping(ni_logical)));

    CHECK_RETURNVAL(PtlGetId(ni_logical, &myself));
    CHECK_RETURNVAL(PtlEQAlloc(ni_logical, 100, &pt_eq_handle));
    CHECK_RETURNVAL(PtlPTAlloc
                        (ni_logical, 0, pt_eq_handle, PTL_PT_ANY,
                        &logical_pt_index));
    assert(logical_pt_index == 0);
    /* Now do the initial setup on ni_logical */
    potato_catcher.start   = &potato;
    potato_catcher.length  = sizeof(potato);
    potato_catcher.uid     = PTL_UID_ANY;
    potato_catcher.options = OPTIONS;
#if INTERFACE == 1
    potato_catcher.match_id.rank = PTL_RANK_ANY;
    potato_catcher.match_bits    = 1;
    potato_catcher.ignore_bits   = ~potato_catcher.match_bits;
#endif
    potato_catcher.ct_handle = PTL_CT_NONE;
    CHECK_RETURNVAL(APPEND
                        (ni_logical, logical_pt_index, &potato_catcher,
                        PTL_PRIORITY_LIST, NULL, &potato_catcher_handle));
    {
        ptl_event_t event;
        CHECK_RETURNVAL(PtlEQWait(pt_eq_handle, &event));       // wait for link event
        assert(event.type == PTL_EVENT_LINK);
    }
    /* Now do a barrier (on ni_physical) to make sure that everyone has their
     * logical interface set up */
    libtest_barrier();

    /* now I can communicate between ranks with ni_logical */

    /* set up the potato launcher */
    potato_launcher.start     = &potato;
    potato_launcher.length    = sizeof(potato);
    potato_launcher.options   = PTL_MD_EVENT_CT_ACK | PTL_MD_EVENT_CT_SEND;
    potato_launcher.eq_handle = PTL_EQ_NONE;    // i.e. don't queue send events
    CHECK_RETURNVAL(PtlCTAlloc(ni_logical, &potato_launcher.ct_handle));
    CHECK_RETURNVAL(PtlMDBind
                        (ni_logical, &potato_launcher, &potato_launcher_handle));

    /* rank 0 starts the potato going */
    if (myself.rank == 0) {
        ptl_process_t nextrank;
        nextrank.rank  = myself.rank + 1;
        nextrank.rank *= (nextrank.rank <= num_procs - 1);
        gettimeofday(&start, NULL);
        CHECK_RETURNVAL(PtlPut(potato_launcher_handle, 0, potato_launcher.length,
                               (LOOPS == 1) ? PTL_OC_ACK_REQ : PTL_NO_ACK_REQ,
                               nextrank, logical_pt_index, 1, 0, NULL, 1));
    }

    {                                  /* the potato-passing loop */
        size_t         waitfor;
        ptl_ct_event_t ctc;
        ptl_process_t  nextrank;
        nextrank.rank  = myself.rank + 1;
        nextrank.rank *= (nextrank.rank <= num_procs - 1);
        for (waitfor = 1; waitfor <= LOOPS; ++waitfor) {
            do {
                ptl_event_t event;
                CHECK_RETURNVAL(PtlEQWait(pt_eq_handle, &event));       // wait for potato
                if (event.type != PTL_EVENT_PUT) {
                    printf("unexpected event: %i\n", (int)event.type);
                } else {
                    break;
                }
            } while (1);
            /* I have the potato! */
            ++potato;
            if (potato < LOOPS * (num_procs)) { // otherwise, the recipient may have exited
                /* Bomb's away! */
                if (myself.rank == 0) {
                    CHECK_RETURNVAL(PtlPut(potato_launcher_handle, 0,
                                           potato_launcher.length,
                                           (waitfor == (LOOPS - 1)) ? PTL_OC_ACK_REQ : PTL_NO_ACK_REQ,
                                           nextrank, logical_pt_index, 3, 0, NULL, 2));
                } else {
                    CHECK_RETURNVAL(PtlPut(potato_launcher_handle, 0,
                                           potato_launcher.length,
                                           (waitfor == LOOPS) ? PTL_OC_ACK_REQ : PTL_NO_ACK_REQ,
                                           nextrank, logical_pt_index, 3, 0, NULL, 2));
                }
            }
        }
        // make sure that last send completed before exiting
        CHECK_RETURNVAL(PtlCTWait(potato_launcher.ct_handle, LOOPS + 1, &ctc));
        assert(ctc.failure == 0);
        if (myself.rank == 0) {
            printf("Final value of potato = %i\n", potato);
        }
    }
    if (myself.rank == 0) {
        double accumulate = 0.0;
        gettimeofday(&stop, NULL);
        accumulate =
            (stop.tv_sec + stop.tv_usec * 1e-6) - (start.tv_sec +
                                                   start.tv_usec * 1e-6);
        /* calculate the average time waiting */
        printf("Total time: %g secs\n", accumulate);
        accumulate /= LOOPS;
        printf("Average time around the loop: %g microseconds\n",
               accumulate * 1e6);
        accumulate /= num_procs;
        printf("Average catch-to-toss latency: %g microseconds\n",
               accumulate * 1e6);
    }

    /* cleanup */
    CHECK_RETURNVAL(PtlMDRelease(potato_launcher_handle));
    CHECK_RETURNVAL(PtlCTFree(potato_launcher.ct_handle));
    CHECK_RETURNVAL(UNLINK(potato_catcher_handle));

    /* major cleanup */
    CHECK_RETURNVAL(PtlPTFree(ni_logical, logical_pt_index));
    CHECK_RETURNVAL(PtlEQFree(pt_eq_handle));
    CHECK_RETURNVAL(PtlNIFini(ni_logical));
    CHECK_RETURNVAL(libtest_fini());

    PtlFini();

    return 0;
}

/* vim:set expandtab: */
