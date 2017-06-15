/********************************************************/
/*                                                      */
/*  pdht.h - portals distributed hash table             */
/*                                                      */
/*  author: d. brian larkins                            */
/*  created: 3/20/16                                    */
/*                                                      */
/********************************************************/


/**
 * @file
 * 
 * portals distributed hash table specification
 */

#pragma once

#include <limits.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/types.h>

#include <pthread.h>

#ifdef __APPLE__
#include <portals4/pmi.h>
#else
#include <slurm/pmi.h>
#endif // pmi.h
#include <portals4.h>

#define PDHT_MAX_TABLES        20
#define PDHT_MAX_PTES          25
#define PDHT_MAX_COUNTERS      20
#define PDHT_MAX_REDUCE_ELEMS 128

/**********************************************/
/* statistics/performance data                */
/**********************************************/

/**
 *  timer counter/accumulator
 */
struct pdht_timer_s {
  struct timespec total;  //!< total time accumulated
  struct timespec last;   //!< keeps last start time
  struct timespec temp;   //!< for accumulation
};
typedef struct pdht_timer_s pdht_timer_t;

struct pdht_stats_s {
  u_int64_t    puts;
  u_int64_t    pendputs;      // track PtlPuts to pending q for fence
  u_int64_t    appends;       // track complete appends to active q
  u_int64_t    tappends[PDHT_MAX_PTES];      // track complete appends to active q
  u_int64_t    gets;
  u_int64_t    collisions;
  u_int64_t    notfound;
  u_int64_t    ptcounts[PDHT_MAX_PTES];
  pdht_timer_t ptimer; // put timer
  pdht_timer_t gtimer; // get timer
  pdht_timer_t t1; // utility timer 1
  pdht_timer_t t2; // utility timer 2
  pdht_timer_t t3; // utility timer 3
  pdht_timer_t t4; // utility timer 4
  pdht_timer_t t5; // utility timer 5
  pdht_timer_t t6; // utility timer 6
};
typedef struct pdht_stats_s pdht_stats_t;


/**********************************************/
/* sub structures contained in global context */
/**********************************************/
struct pdht_s; // forward ref

// polling queue - ME append list entry
struct pdht_append_s {
  void            *start;
  ptl_size_t       length;
  ptl_match_bits_t bits;
  char             entry[0];
};
typedef struct pdht_append_s pdht_append_t;
   
// per-process queue of pending ME appends
struct pdht_pollqent_s {
  ptl_handle_md_t  md;            // memory descriptor for queue
  ptl_handle_ct_t  ct;            // event counter for queue
  u_int64_t        completed;     // polling thread counts completion 
  pdht_append_t   *appendlist;    // ME entries to process
};
typedef struct pdht_pollqent_s pdht_pollqent_t;

// overall ME append queue
struct pdht_pollq_s {
  u_int32_t          qlen;              // entries per process queue
  u_int32_t          qentrysize;        // max elemsize for an ht entry
  pdht_pollqent_t   *q;                 // queue data
};
typedef struct pdht_pollq_s pdht_pollq_t;


/* portals specific data for global context */
struct pdht_portals_s {
  ptl_handle_ni_t  phy;                 //!< physical NI
  ptl_handle_ni_t  lni;                 //!< logical NI
  ptl_ni_limits_t  ni_limits;           //!< logical NI limits
  ptl_process_t   *mapping;             //!< physical/logical NI mapping
  ptl_handle_md_t  collective_md;       //!< collective MD handle
  ptl_handle_me_t  barrier_me;          //!< barrier ME handle
  ptl_handle_me_t  reduce_lme;          //!< reduce left ME handle
  ptl_handle_me_t  reduce_rme;          //!< reduce right ME handle
  ptl_handle_ct_t  barrier_ct;          //!< barrier CT handle
  ptl_handle_ct_t  collective_ct;       //!< collective CT handle
  ptl_size_t       collective_count;    //!< collective count
  ptl_size_t       barrier_count;       //!< barrier  count
  void            *collective_lscratch; //!< scratch space for collective ops
  void            *collective_rscratch; //!< scratch space for collective ops
  u_int32_t        pt_entries;          //!< number of portals table entries per hash table
};
typedef struct pdht_portals_s pdht_portals_t;

/**********************************************/
/* global context data structure              */
/**********************************************/
struct pdht_context_s {
  int              dhtcount;     //!< DHTs that have been created
  struct pdht_s   *hts[PDHT_MAX_TABLES]; //!< array of active hash tables
  int              rank;         //!< process rank
  int              size;         //!< process count
  int              dbglvl;       //!< debug level for error printing
  pdht_portals_t   ptl;          //!< Portals 4 ADTs
  pthread_t       progress_tid; //!< progress thread id
};
typedef struct pdht_context_s pdht_context_t;

extern pdht_context_t *c;


/*************************************************/
/* sub structures contained in per-DHT structure */
/*************************************************/

/* communication mode */
enum pdht_mode_e {
  PdhtModeStrict,     // blocking,synchronous
  PdhtModeBundled,    // batched puts, fence at end of bundled puts
  PdhtModeAsync       // anything goes
};
typedef enum pdht_mode_e pdht_mode_t;
#define PDHT_DEFAULT_MODE PdhtModeStrict

/* pending mode */
enum pdht_pmode_e {
  PdhtPendingPoll,
  PdhtPendingTrig
};
typedef enum pdht_pmode_e pdht_pmode_t;
#define PDHT_DEFAULT_PMODE PdhtPendingTrig
//#define PDHT_DEFAULT_PMODE PdhtPendingPoll

/* DHT operatation status */
enum pdht_status_e {
  PdhtStatusOK,
  PdhtStatusError,
  PdhtStatusNotFound,
  PdhtStatusCollision
};
typedef enum pdht_status_e pdht_status_t;

#define PDHT_NULL_HANDLE -1
typedef int pdht_handle_t;

typedef void (*pdht_hashfunc)(struct pdht_s *dht, void *key, ptl_match_bits_t *bits, uint32_t *ptindex, ptl_process_t *rank);


/* portals-specific data structures */
struct pdht_htportals_s {
  ptl_handle_ni_t lni;                          //!< portals logical NI
  unsigned        ptalloc_opts;                 //!< options to PtlPTAlloc (for unordered matching)
  unsigned        nptes;                        //!< number of pending / active PTE pairs
  ptl_pt_index_t  putindex_base;                //!< actve_base (1) + dht->nptes
  ptl_pt_index_t  getindex[PDHT_MAX_PTES];      //!< portal table entry index
  ptl_pt_index_t  putindex[PDHT_MAX_PTES];      //!< portal table entry index
  ptl_handle_eq_t eq[PDHT_MAX_PTES];            //!< event queue for put PT entry
  ptl_handle_eq_t aeq[PDHT_MAX_PTES];           //!< event queue for get PT entry (fence/triggered)
  ptl_me_t        me;                           //!< default match entry for ht
  ptl_handle_md_t lmd;                          //!< memory descriptor for any outgoing put/gets
  ptl_handle_eq_t lmdeq;                        //!< event queue for local MD
  ptl_handle_ct_t lmdct;                        //!< counter for local MD
  ptl_handle_me_t centries[PDHT_MAX_COUNTERS];  //!< ME entries for atomic counters (target, rank 0 only)
  ptl_handle_md_t countmds[PDHT_MAX_COUNTERS];  //!< MDs for initiator counter ops (initiator, all ranks)
  ptl_handle_ct_t countcts[PDHT_MAX_COUNTERS];  //!< CTs for initiator counter ops (initiator, all ranks)
  ptl_ct_event_t  curcounts;                    //!< current fail/success counts for local MD state (tracks progress)
  ptl_size_t      lfail;                        //!< number of strict messages received
};
typedef struct pdht_htportals_s pdht_htportals_t;


/**********************************************/
/* main DHT data structure                    */
/**********************************************/
struct pdht_s {
  void             *ht;  
  unsigned          keysize;
  unsigned          elemsize;
  unsigned          entrysize;
  unsigned          maxentries;
  unsigned          pendq_size;
  pdht_hashfunc     hashfn;
  unsigned          nextfree;
  pdht_mode_t       mode;
  pdht_pmode_t      pmode;
  pdht_stats_t      stats;
  uint64_t          counters[PDHT_MAX_COUNTERS]; // rank 0 target (master) counters
  uint64_t          lcounts[PDHT_MAX_COUNTERS];  // initiator side buffers
  int               countercount; // :)
  pthread_mutex_t   completion_mutex;    //!< thread mutex to synch between progress thread and fence
  pdht_htportals_t  ptl;
  pdht_status_t   (*put)(struct pdht_s *dht, void *k, void *v);
  pdht_status_t   (*get)(struct pdht_s *dht, void *k, void **v);
  pdht_handle_t   (*nbput)(struct pdht_s *dht, void *k, void *v);
  pdht_handle_t   (*nbget)(struct pdht_s *dht, void *k, void **v);
};
typedef struct pdht_s pdht_t;

/* application-specified tunable parameters */
#define PDHT_TUNE_NPTES      0x01
#define PDHT_TUNE_PMODE      0x02
#define PDHT_TUNE_ENTRY      0x04
#define PDHT_TUNE_PENDQ      0x08
#define PDHT_TUNE_PTOPT      0x10
#define PDHT_TUNE_ALL        0xffffffff
struct pdht_config_s {
  unsigned      nptes;
  pdht_pmode_t  pendmode;
  unsigned      maxentries;
  unsigned      pendq_size;
  unsigned      ptalloc_opts;
};
typedef struct pdht_config_s pdht_config_t;

/* atomic associative operators */
enum pdht_oper_e {
  AssocOpAdd
};
typedef enum pdht_oper_e pdht_oper_t;

/* atomic operation data types */
enum pdht_datatype_e {
  IntType,
  LongType,
  DoubleType,
  CharType,
  BoolType
};
typedef enum pdht_datatype_e pdht_datatype_t;

enum pdht_reduceop_e {
  PdhtReduceOpSum
};
typedef enum pdht_reduceop_e pdht_reduceop_t;


/* DHT iterators (UNSUPPORTED) */
struct pdht_iter_s {
  // XXX iteration state stuff needs added.
  int   (*hasnext)(struct pdht_iter_s *it);
  void *(*next)(struct pdht_iter_s *it);
};
typedef struct pdht_iter_s pdht_iter_t;

/********************************************************/
/* portals distributed hash table public interface      */
/********************************************************/

// create / destroy single DHT -- init.c
pdht_t              *pdht_create(int keysize, int elemsize, pdht_mode_t mode);
void                 pdht_free(pdht_t *dht);
void                 pdht_tune(unsigned opts, pdht_config_t *config);


// Communication Completion Operations -- commsynch.c
void                 pdht_barrier(void);
void                 pdht_fence(pdht_t *dht);
pdht_status_t        pdht_reduce(void *in, void *out, pdht_reduceop_t op, pdht_datatype_t type, int elems);
pdht_status_t        pdht_allreduce(void *in, void *out, pdht_reduceop_t op, pdht_datatype_t type, int elems);
pdht_status_t        pdht_broadcast(void *buf, pdht_datatype_t type, int elems);

pdht_status_t        pdht_test(pdht_handle_t h);
pdht_status_t        pdht_wait(pdht_handle_t h);
pdht_status_t        pdht_waitrank(int rank);
pdht_status_t        pdht_waitall(void);

// Put / Get Operations -- putget.c
pdht_status_t        pdht_put(pdht_t *dht, void *key, void *value);
pdht_status_t        pdht_add(pdht_t *dht, void *key, void *value);
pdht_status_t        pdht_update(pdht_t *dht, void *key, void *value);
pdht_status_t        pdht_get(pdht_t *dht, void *key, void *value);
pdht_status_t        pdht_insert(pdht_t *dht, ptl_match_bits_t bits, uint32_t ptindex, void * key, void *value);

// Asynchronous Put / Get Operations -- nbputget.c
pdht_handle_t        pdht_nbput(pdht_t *dht, void *key, void *value);
pdht_handle_t        pdht_nbget(pdht_t *dht, void *key, void **value);

// Associative Update Operations -- assoc.c
pdht_status_t        pdht_acc(pdht_t *dht, void *key, pdht_datatype_t type, pdht_oper_t op, void *value);
pdht_handle_t        pdht_nbacc(pdht_t *dht, void *key, pdht_datatype_t type, pdht_oper_t op, void *value);

// Hash Function Operations -- hash.c
void                 pdht_sethash(pdht_t *dht, pdht_hashfunc hfun);

  
// Iteration operations -- iter.c
pdht_status_t        pdht_iterate(pdht_t *dht, pdht_iter_t *it);
pdht_status_t        pdht_iterate_single(pdht_t *dht, pdht_iter_t *it);
int                  pdht_hasnext(pdht_iter_t *it);
void                *pdht_getnext(pdht_iter_t *it);


// atomics / counter support atomics.c
int                  pdht_counter_init(pdht_t *ht, int initval);
uint64_t             pdht_counter_inc(pdht_t *ht, int counter, uint64_t val);
void                 pdht_counter_reset(pdht_t *ht, int counter);

//trig.c - temp
void print_count(pdht_t *dht, char *msg);


//util.c
void   pdht_print_stats(pdht_t *dht);
void   pdht_print_active(pdht_t *dht, void nprinter(void *node));

#include <pdht_inline.h>
