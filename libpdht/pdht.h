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
#include <pdht.h>


struct pdht_portals_s {
   ptl_handle_ni_t  phy;          //!< physical NI
   ptl_handle_ni_t  lni;          //!< logical NI
   ptl_ni_limits_t  ni_limits;    //!< logical NI limits
   ptl_process_t   *mapping;      //!< physical/logical NI mapping
};
typedef struct pdht_portals_s pdht_portals_t;

/**
 *  global configuration data structure
 */
struct pdht_context_s {
   int              dhtcount;     //!< DHTs that have been created
   int              rank;         //!< process rank
   int              size;         //!< process count
   pdht_portals_t   ptl;          //!< Portals 4 ADTs
};
typedef struct pdht_context_s pdht_context_t;

enum pdht_status_e {
   PdhtStatusOK,
   PdhtStatusError
};
typedef enum pdht_status_e pdht_status_t;

typedef int pdht_handle_t;

struct pdht_s {
   pdht_context_t *ctx;
   pdht_status_t (*put)(struct pdht_s *dht, void *k, int ksize, void *v);
   pdht_status_t (*get)(struct pdht_s *dht, void *k, int ksize, void **v);
   pdht_handle_t (*nbput)(struct pdht_s *dht, void *k, int ksize, void *v);
   pdht_handle_t (*nbget)(struct pdht_s *dht, void *k, int ksize, void **v);
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

struct pdht_iter_t {
   // XXX teration state stuff needs added.
   int   (*hasnext)(pdht_iter_t *it);
   void *(*next)(pdht_iter_t *it);
}

/********************************************************/
/* portals distributed hash table public interface      */
/********************************************************/

// create / destroy single DHT -- init.c
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
pdht_status_t        pdht_put(pdht_t *dht, void *key, int ksize, void *value);
pdht_status_t        pdht_get(pdht_t *dht, void *key, int ksize, void **value);

// Asynchronous Put / Get Operations -- nbputget.c
pdht_handle_t        pdht_nbput(pdht_t *dht, void *key, int ksize, void *value);
pdht_handle_t        pdht_nbget(pdht_t *dht, void *key, int ksize, void **value);

// Associative Update Operations -- assoc.c
pdht_status_t        pdht_acc(pdht_t *dht, void *key, int ksize, pdht_datatype_t type, pdht_oper_t op, void *value);
pdht_handle_t        pdht_nbacc(pdht_t *dht, void *key, int ksize, pdht_datatype_t type, pdht_oper_t op, void *value);

// Iteration operations -- iter.c
pdht_status_t        pdht_iterate(pdht_t *dht, pdht_iter_t *it);
pdht_status_t        pdht_iterate_single(pdht_t *dht, pdht_iter_t *it);
int                  pdht_hasnext(pdht_iter_t *it);
void                *pdht_getnext(pdht_iter_t *it);

#include <pdht_inline.h>
