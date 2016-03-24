/***********************************************************/
/*                                                         */
/*  pdht_imple.h - PDHT private implementation protos/ADTs */
/*                                                         */
/*  author: d. brian larkins                               */
/*  created: 3/20/16                                       */
/*                                                         */
/***********************************************************/

#pragma once
#include <assert.h>
#include <stdarg.h>
#include <pdht.h>


#define PDHT_DEBUG

#ifdef PDHT_DEBUG
  #define pdht_dprintf(...) pdht_dbg_printf(__VA_ARGS__)
#else
  #define pdht_dprintf(...)
#endif

/**
 * @file
 * 
 * portals distributed hash table implementations ADTs
 */


// global (per-process private) data structures
extern pdht_context_t *c;

/********************************************************/
/* portals distributed hash table prototypes            */
/********************************************************/

// Initialization / Finalization -- init.c
void                 pdht_init(void);
void                 pdht_fini(void);
void                 pdht_clearall(void);

// pmi.c
void init_pmi(pdht_context_t *c);

// util.c
int  eprintf(const char *format, ...);
int  pdht_dbg_printf(const char *format, ...);