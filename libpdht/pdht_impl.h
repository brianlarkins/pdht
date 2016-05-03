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

#ifdef PDHT_DEBUG
  #define pdht_dprintf(...) pdht_dbg_printf(__VA_ARGS__)
#else
  #define pdht_dprintf(...) ;
#endif

#define PDHT_DEFAULT_TABLE_SIZE 1000
#define PDHT_EVENTQ_SIZE         500

#define __PDHT_PUT_INDEX 22
#define __PDHT_PUT_MATCH 0xcafef00d
#define __PDHT_GET_INDEX 23

#define __PDHT_BARRIER_INDEX 24
#define __PDHT_BARRIER_MATCH 0xdeadbeef

/**
 * @file
 * 
 * portals distributed hash table implementations ADTs
 */

// overlay struct for casting
struct _pdht_ht_entry_s {
   ptl_handle_me_t   me;
   char              data[0];
};
typedef struct _pdht_ht_entry_s _pdht_ht_entry_t;


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

// hash.c - PDHT hash function operations
void pdht_hash(pdht_t *dht, void *key, ptl_match_bits_t *bits, ptl_process_t *rank);

// poll.c - PDHT polling tasks
void pdht_polling_init(pdht_t *dht);
void pdht_polling_fini(pdht_t *dht);
void pdht_poll(pdht_t *dht);
