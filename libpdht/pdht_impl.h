/***********************************************************/
/*                                                         */
/*  pdht_imple.h - PDHT private implementation protos/ADTs */
/*                                                         */
/*  author: d. brian larkins                               */
/*  created: 3/20/16                                       */
/*                                                         */
/***********************************************************/

#pragma once
#define _XOPEN_SOURCE 700

#include <assert.h>
#include <errno.h>
#include <stdarg.h>
#include <time.h>
#include <pdht.h>

#define PDHT_DEBUG
//#define PDHT_DEBUG_TRACE

#ifdef PDHT_DEBUG
  #define pdht_dprintf(...) pdht_dbg_printf(__VA_ARGS__)
#else
  #define pdht_dprintf(...) ;
#endif

// default table size should be bigger than 2x pending queue size
#define PDHT_DEFAULT_TABLE_SIZE 300000
#define PDHT_PENDINGQ_SIZE      100000

#define PDHT_MAXKEYSIZE 8

#define __PDHT_PENDING_INDEX 22
#define __PDHT_PENDING_MATCH 0xcafef00d
#define __PDHT_ACTIVE_INDEX 23

#define __PDHT_BARRIER_INDEX 24
#define __PDHT_BARRIER_MATCH 0xdeadbeef

#define PDHT_START_TIMER(HT,TMR) HT->stats.TMR.last   = pdht_get_wtime();
#define PDHT_STOP_TIMER(HT,TMR)  HT->stats.TMR.total += pdht_get_wtime() - HT->stats.TMR.last;
#define PDHT_READ_TIMER(HT,TMR)  HT->stats.TMR.total

/**
 * @file
 * 
 * portals distributed hash table implementations ADTs
 */

// overlay struct for casting
struct _pdht_ht_entry_s {
   ptl_handle_me_t   pme;  // pending ME handle
   ptl_handle_me_t   ame;  // active ME handle
   char              key[PDHT_MAXKEYSIZE]; // fixed data size for key right now.
   char              data[0]; // opaque payload
};
typedef struct _pdht_ht_entry_s _pdht_ht_entry_t;

// overlay struct for casting
struct _pdht_ht_trigentry_s {
   ptl_handle_me_t   pme;  // pending ME handle
   ptl_handle_me_t   ame;  // active ME handle
   ptl_handle_ct_t   tct;  // trigger counter for each entry
   ptl_me_t          me;   // ME buffer for copying match bits over
   char              key[PDHT_MAXKEYSIZE]; // fixed data size for key right now.
   char              data[0]; // opaque payload
};
typedef struct _pdht_ht_trigentry_s _pdht_ht_trigentry_t;


/********************************************************/
/* portals distributed hash table prototypes            */
/********************************************************/

// Initialization / Finalization -- init.c
void                 pdht_init(void);
void                 pdht_fini(void);
void                 pdht_clearall(void);

// commsynch.c
void                 pdht_barrier_init(pdht_context_t *c);

// pmi.c
void init_pmi();
void init_only_barrier(void);

// util.c
int  eprintf(const char *format, ...);
int  pdht_dbg_printf(const char *format, ...);
char *pdht_ptl_error(int error_code);
char *pdht_event_to_string(ptl_event_kind_t evtype);
void pdht_dump_event(ptl_event_t *ev);


// hash.c - PDHT hash function operations
void pdht_hash(pdht_t *dht, void *key, ptl_match_bits_t *bits, ptl_process_t *rank);

// poll.c - PDHT polling tasks
void pdht_polling_init(pdht_t *dht);
void pdht_polling_fini(pdht_t *dht);
void pdht_poll(pdht_t *dht);

// trig.c - PDHT triggered tasks
void pdht_trig_init(pdht_t *dht);
void pdht_trig_fini(pdht_t *dht);
