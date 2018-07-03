#include <limits.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/types.h>
#include <shmem.h>
#include "uthash.h"

#define PDHT_MAX_TABLES 2
#define PDHT_MAX_COUNTERS 20

static inline struct timespec pdht_get_wtime(){
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC_RAW, &ts);
  return ts;
}

#define PDHT_INIT_ATIMER(TMR) do { TMR.total.tv_sec = 0; TMR.total.tv_nsec = 0; } while (0)
#define PDHT_START_ATIMER(TMR) TMR.last   = pdht_get_wtime();
#define PDHT_STOP_ATIMER(TMR) do {\
                                  TMR.temp = pdht_get_wtime();\
                                  TMR.total.tv_sec += (TMR.temp.tv_sec - TMR.last.tv_sec);\
                                  TMR.total.tv_nsec += (TMR.temp.tv_nsec - TMR.last.tv_nsec);\
                                } while (0)
// PDHT_READ_TIMER returns elapsed time in nanoseconds
#define PDHT_READ_ATIMER(TMR)  ((1000000000L * (TMR.total.tv_sec)) + TMR.total.tv_nsec)
#define PDHT_READ_ATIMER_USEC(TMR)  PDHT_READ_ATIMER(TMR)/1000.0
#define PDHT_READ_ATIMER_MSEC(TMR)  PDHT_READ_ATIMER(TMR)/(double)1e6
#define PDHT_READ_ATIMER_SEC(TMR)   PDHT_READ_ATIMER(TMR)/(double)1e9

#define PdhtPendingTrig 0
#define PTL_PT_MATCH_UNORDERED 0
#define PDHT_DEFAULT_RANK_HINT 0 
#define PDHT_TUNE_ALL 0

typedef struct ht_s{
  uint64_t        key;
  void           *value;
  UT_hash_handle  hh;
} ht_t;

/* TODO Actual key not mbits */
typedef struct am_args_s{
  int       ht_index;
  uint64_t  key;
  char      val[0];
} am_args_t;

typedef struct get_args_s{
  int ht_index;
  uint64_t key;
} get_args_t;

typedef struct pdht_timer_s{
  struct timespec total;
  struct timespec last;
  struct timespec temp;
} pdht_timer_t;

/* Because we are porting from Portals */
typedef struct ptl_proces_s{
  int rank;
} ptl_process_t;

typedef uint64_t ptl_match_bits_t;

typedef enum pdht_status_e {
  PdhtStatusOK,
  PdhtStatusError,
  PdhtStatusNotFound,
  PdhtStatusCollision
} pdht_status_t;

typedef enum pdht_reduceop_e{
  PdhtReduceOpSum,
  PdhtReduceOpMin,
  PdhtReduceOpMax
} pdht_reduceop_t;

typedef enum pdht_mode_e{
  PdhtModeStrict,
  PdhtModeBundled,
  PdhtModeAsync
} pdht_mode_t;

typedef enum pdht_datatype_e{
  IntType,
  LongType,
  DoubleType,
  CharType,
  BoolType
} pdht_datatype_t;

typedef enum pdht_local_gets_e{
  PdhtRegular,
  PdhtSearchLocal
} pdht_local_gets_t;

typedef struct pdht_config_s{ 
  int pendmode;
  long unsigned maxentries;
  int quiet;
  pdht_local_gets_t local_gets;
  long unsigned pendq_size;
  long unsigned ptalloc_opts;
  int nptes;
  int rank;
} pdht_config_t;

struct pdht_s;
typedef void (*pdht_hashfunc)(struct pdht_s *dht, void *key, ptl_match_bits_t *bits, uint32_t *ptindex, ptl_process_t *rank);

typedef struct pdht_s{
  int          *counters;
  int           next_counter;
  unsigned      keysize;
  unsigned      elemsize;
  pdht_hashfunc hashfn;
  ht_t         *ht;
  void         *ht_shmem_space;
} pdht_t;

typedef struct pdht_context_s{
  int               dhtcount;
  struct pdht_s    *hts[PDHT_MAX_TABLES];
  int               rank;
  int               size;
  shmem_am_handle_t get_handle;
  shmem_am_handle_t put_handle;
  size_t            data_offset;
  size_t            value_offset;
  void             *pSync;
  void             *pWrk;
} pdht_context_t;

extern pdht_context_t *c;

//setup and tear down
pdht_t *pdht_create(size_t keysize, size_t elemsize, pdht_mode_t mode);
void pdht_free(pdht_t *dht);
void pdht_tune(unsigned opts, pdht_config_t *config);
//synchonization
void            pdht_barrier(void);
void            pdht_fence(pdht_t *dht);
pdht_status_t   pdht_reduce(void *in, void *out, pdht_reduceop_t op, pdht_datatype_t type, int elems);
pdht_status_t   pdht_allreduce(void *in, void *out, pdht_reduceop_t op, pdht_datatype_t type, int elems);

// communication operations
pdht_status_t   pdht_put(pdht_t *dht, void *key, void *value);
pdht_status_t   pdht_update(pdht_t *dht, void *key, void *value);
pdht_status_t   pdht_get(pdht_t *dht, void *key, void *value);
pdht_status_t   pdht_persistent_get(pdht_t *dht, void *key, void *value);

//Hash function operations
void            pdht_sethash(pdht_t *dht, pdht_hashfunc hfun);

//utility stuff
void            pdht_print_stats(pdht_t *dht);
double          pdht_average_time(pdht_t *dht, pdht_timer_t timer);
int eprintf(const char *format, ...);

int pdht_counter_init(pdht_t *dht, int init_val);
int pdht_counter_inc(pdht_t *dht, int counter, int inc_val);
void pdht_counter_reset(pdht_t *dht, int counter);

