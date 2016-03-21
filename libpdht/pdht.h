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

#include <slurm/pmi.h>
#include <portals4.h>


/**
 *  global configuration data structure
 */
struct pdht_context_s {
   int              rank;         //!< process rank
   int              size;         //!< process count
   ptl_handle_ni_t  phy;          //!< physical NI
   ptl_handle_ni_t  lni;          //!< logical NI
   ptl_ni_limits_t  ni_limits;    //!< logical NI limits
   ptl_process_t   *mapping;      //!< physical/logical NI mapping
};
typedef struct pdht_context_s pdht_context_t;


struct pdht_s {
   pdht_context_t *ctx;
};
typedef struct pdht_s pdht_t;

enum pdht_status_e {
   PdhtStatusOK,
   PdhtStatusError
};
typedef enum pdht_status_e pdht_status_t;

typedef int pdht_handle_t;

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


// global (per-process private) data structures
extern pdht_context_t c;

/********************************************************/
/* portals distributed hash table prototypes            */
/********************************************************/

// Initialization / Finalization -- init.c
void                 pdht_init(void);
void                 pdht_fini(void);
void                 pdht_clearall(void);


// create / destroy single DHT
pdht_t              *pdht_create(pdht_context_t *ctx);
void                 pdht_free(pdht_t *dht);


// Communication Completion Operations -- commsynch.c
void                 pdht_fence(int rank);
void                 pdht_allfence(void);
pdht_status_t        pdht_test(pdht_handle_t h);
pdht_status_t        pdht_wait(pdht_handle_t h);
pdht_status_t        pdht_waitrank(int rank);
pdht_status_t        pdht_waitall(void);

// Put / Get Operations -- putget.c
pdht_status_t        pdht_put(void *key, int ksize, void *value);
pdht_status_t        pdht_get(void *key, int ksize, void **value);
// XXX probably should make these inlined for type safety
#define pdht_puti(k,v) pdht_put(k,sizeof(int),v)
#define pdht_geti(k,v) pdht_put(k,sizeof(int),v)
#define pdht_putf(k,v) pdht_put(k,sizeof(double),v)
#define pdht_getf(k,v) pdht_put(k,sizeof(double),v)

// Asynchronous Put / Get Operations -- nbputget.c
pdht_handle_t        pdht_nbput(void *key, int ksize, void *value);
pdht_handle_t        pdht_nbget(void *key, int ksize, void **value);
// XXX probably should make these inlined for type safety
#define pdht_nbputi(k,v) pdht_nbput(k,sizeof(int),v)
#define pdht_nbgeti(k,v) pdht_nbput(k,sizeof(int),v)
#define pdht_nbputf(k,v) pdht_nbput(k,sizeof(double),v)
#define pdht_nbgetf(k,v) pdht_nbput(k,sizeof(double),v)

// Associative Update Operations -- assoc.c
pdht_status_t pdht_acc(void *key, int ksize, pdht_datatype_t type, pdht_oper_t op, void *value);
pdht_handle_t pdht_nbacc(void *key, int ksize, pdht_datatype_t type, pdht_oper_t op, void *value);

// pmi.c
void init_pmi(void);

// util.c
int  eprintf(const char *format, ...);
