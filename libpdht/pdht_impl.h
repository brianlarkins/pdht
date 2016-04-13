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
#include <stdarg.h>
#include <time.h>
#include <pdht.h>

#define PDHT_DEBUG

#ifdef PDHT_DEBUG
  #define pdht_dprintf(...) pdht_dbg_printf(__VA_ARGS__)
#else
  #define pdht_dprintf(...) ;
#endif

#define PDHT_DEFAULT_TABLE_SIZE 100000
#define PDHT_EVENTQ_SIZE        400

#define __PDHT_BARRIER_INDEX 23
#define __PDHT_BARRIER_MATCH 0xdeadbeef

/**
 * @file
 * 
 * portals distributed hash table implementations ADTs
 */


// overlay struct for casting
struct _pdht_ht_entry_s {
   ptl_handle_me_t   me;
   ptl_handle_ct_t   ct;
   char              data[0];
};
typedef struct _pdht_ht_entry_s _pdht_ht_entry_t;

// global (per-process private) data structures
extern pdht_context_t *c;

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
