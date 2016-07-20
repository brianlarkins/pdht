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
#include <unistd.h>
#include <sys/types.h>

#include <slurm/pmi.h>
#include <portals4.h>
#include <pdht.h>


/**********************************************/
/* statistics/performance data                */
/**********************************************/
struct pdht_stats_s {
  u_int64_t puts;
  u_int64_t gets;
  u_int64_t collisions;
  u_int64_t notfound;
  double    ptime;
  double    gtime;
};
typedef struct pdht_stats_s pdht_stats_t;


/**********************************************/
/* sub structures contained in global context */
/**********************************************/

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
  ptl_handle_ni_t  phy;           //!< physical NI
  ptl_handle_ni_t  lni;           //!< logical NI
  ptl_ni_limits_t  ni_limits;     //!< logical NI limits
  ptl_process_t   *mapping;       //!< physical/logical NI mapping
  ptl_handle_md_t  barrier_md;    //!< barrier MD handle
  ptl_handle_me_t  barrier_me;    //!< barrier ME handle
  ptl_handle_ct_t  barrier_ct;    //!< barrier CT handle
  ptl_size_t       barrier_count; //!< barrier count
};
typedef struct pdht_portals_s pdht_portals_t;

/**********************************************/
/* global context data structure              */
/**********************************************/
struct pdht_context_s {
  int              dhtcount;     //!< DHTs that have been created
  int              rank;         //!< process rank
  int              size;         //!< process count
  pdht_portals_t   ptl;          //!< Portals 4 ADTs
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
#define PDHT_DEFAULT_PMODE PdhtPendingPoll

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

struct pdht_s;
typedef void (*pdht_hashfunc)(struct pdht_s *dht, void *key, ptl_match_bits_t *bits, ptl_process_t *rank);


/* portals-specific data structures */
struct pdht_htportals_s {
  ptl_handle_ni_t lni;           //!< portals logical NI
  unsigned        getindex;      //!< portal table entry index
  unsigned        putindex;      //!< portal table entry index
  ptl_handle_eq_t eq;            //!< event queue for put PT entry
  ptl_me_t        me;            //!< default match entry for ht
  ptl_handle_md_t lmd;           //!< memory descriptor for any outgoing put/gets
  ptl_handle_eq_t lmdeq;         //!< event queue for local MD
  ptl_handle_ct_t lmdct;         //!< counter for local MD
  ptl_size_t      lcount;        //!< number of strict messages received
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
  pdht_hashfunc     hashfn;
  unsigned          nextfree;
  pdht_mode_t       mode;
  pdht_pmode_t      pmode;
  pdht_stats_t      stats;
  pdht_htportals_t  ptl;
  pdht_status_t   (*put)(struct pdht_s *dht, void *k, void *v);
  pdht_status_t   (*get)(struct pdht_s *dht, void *k, void **v);
  pdht_handle_t   (*nbput)(struct pdht_s *dht, void *k, void *v);
  pdht_handle_t   (*nbget)(struct pdht_s *dht, void *k, void **v);
};
typedef struct pdht_s pdht_t;

enum pdht_oper_e {
  AssocOpAdd
};
typedef enum pdht_oper_e pdht_oper_t;

enum pdht_datatype_e {
  IntType,
  DoubleType,
  CharType,
  BoolType
};
typedef enum pdht_datatype_e pdht_datatype_t;

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


// Communication Completion Operations -- commsynch.c
void                 pdht_fence(int rank);
void                 pdht_allfence(void);
pdht_status_t        pdht_test(pdht_handle_t h);
pdht_status_t        pdht_wait(pdht_handle_t h);
pdht_status_t        pdht_waitrank(int rank);
pdht_status_t        pdht_waitall(void);
void                 pdht_barrier(void);

// Put / Get Operations -- putget.c
pdht_status_t        pdht_put(pdht_t *dht, void *key, void *value);
pdht_status_t        pdht_get(pdht_t *dht, void *key, void *value);

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

//trig.c - temp
void print_count(pdht_t *dht, char *msg);

#include <pdht_inline.h>
