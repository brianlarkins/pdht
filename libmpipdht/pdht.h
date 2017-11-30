#include <assert.h>
#include <pthread.h>
#include <signal.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include <portals4.h>
#include "uthash.h"
#include "mpi.h"

/* timer definitions */
struct pdht_timer_s{
  double total;
  double last;
  double temp;

};
typedef struct pdht_timer_s pdht_timer_t;

//#define THREAD_MULTIPLE 0

#define PDHT_PTALLOC_OPTIONS PTL_PT_MATCH_UNORDERED
#define PDHT_MAX_TABLES 20
#define PDHT_MAXKEYSIZE 32
#define PDHT_MAX_COUNTERS 20

#define PDHT_START_ATIMER(TMR) TMR.last   = MPI_Wtime();
#define PDHT_STOP_ATIMER(TMR) do {\
                                  TMR.temp = MPI_Wtime();\
                                  TMR.total += (TMR.temp - TMR.last);\
                                } while (0)
// PDHT_READ_TIMER returns elapsed time in nanoseconds
#define PDHT_READ_ATIMER(TMR)  1000000000L * (TMR.total)
#define PDHT_READ_ATIMER_USEC(TMR)  PDHT_READ_ATIMER(TMR)/1000.0
#define PDHT_READ_ATIMER_MSEC(TMR)  PDHT_READ_ATIMER(TMR)/(double)1e6
#define PDHT_READ_ATIMER_SEC(TMR)   PDHT_READ_ATIMER(TMR)/(double)1e9
#define PDHT_INIT_ATIMER(TMR) do { TMR.total = 0;} while (0)


/* Portals defs for API exposed details (see hash function) */
#if 0
typedef uint64_t ptl_match_bits_t;
struct ptl_process_s { uint64_t rank; };
typedef struct ptl_process_s ptl_process_t;
#endif
typedef struct pdht_ptl_s { int nptes; } pdht_ptl_t;



/* MPI impementation details */

// options for local access optimizations
enum pdht_local_gets_e{
  PdhtRegular,
  PdhtSearchLocal
};
typedef enum pdht_local_gets_e pdht_local_gets_t;


struct pdht_s; // forward ref

// hash function proto
typedef void (*pdht_hashfunc)(struct pdht_s *dht, void *key, ptl_match_bits_t *mbits, uint32_t *ptindex, ptl_process_t *rank);

// message types
typedef enum { pdhtGet, pdhtPut, pdhtStop, pdhtCounterInit, pdhtCounterReset, pdhtCounterInc, pdhtCreateHt } msg_type;

enum pdht_datatype_e {
  IntType,
  LongType,
  DoubleType,
  CharType,
  BoolType
};
typedef enum pdht_datatype_e pdht_datatype_t;

enum pdht_reduceop_e{
  PdhtReduceOpSum,
  PdhtReduceOpMin,
  PdhtReduceOpMax
};
typedef enum pdht_reduceop_e pdht_reduceop_t;


/* global context structure */
struct pdht_context_s {
  int            rank;
  int            size;  
  struct pdht_s *hts[PDHT_MAX_TABLES]; // active HTs
  int            dhtcount;             // # active HTs
  pthread_t      comm_tid;             // tid of comm thread
  int            thread_active;        // comm thread status
  MPI_Datatype   msgType;              // registered MPI datatype for structured comms
  int            maxbufsize;           // maximum MPI receive buffer size for all HTs
  MPI_Comm       split_comm;
  int            pid;
};
typedef struct pdht_context_s pdht_context_t;

extern pdht_context_t *c;

/* MPI private UTHASH HT entry structure */
struct ht_s{
  uint64_t key;
  void *value;
  UT_hash_handle hh;
  
};
typedef struct ht_s ht_t;


/* request message structure for PDHT ops */
struct message_s {
  msg_type         type;
  int              rank;
  int              ht_index;
  ptl_match_bits_t mbits;
  char             key[0];  // optional payload for put requests
};
typedef struct message_s message_t;


/* reply message structure (for casting) */
struct reply_s {
  char             status;  // found / not found bit
  char             key[PDHT_MAXKEYSIZE];
  char             value[0];
};
typedef struct reply_s reply_t;

/* define PDHT-MPI message tags */
#define PDHT_TAG_COMMAND    1
#define PDHT_TAG_REPLY      2
#define PDHT_TAG_ACK        3
#define PDHT_COUNTER_REPLY  4
#define PDHT_CREATE_HT      5 

/* fake tuning structure to not break regular PDHT benches/apps */
//tuning stuff that is not used
#define PDHT_TUNE_NPTES     0x01
#define PDHT_TUNE_PMODE     0x02
#define PDHT_TUNE_ENTRY     0x04
#define PDHT_TUNE_PENDQ     0x08
#define PDHT_TUNE_PTOPT     0x10
#define PDHT_TUNE_QUIET     0x20
#define PDHT_TUNE_GETS      0x40
#define PDHT_TUNE_ALL       0xffffffff
struct pdht_config_s{
  int pendmode;
  long unsigned maxentries;
  int quiet;
  pdht_local_gets_t local_gets;
  long unsigned pendq_size;
  long unsigned ptalloc_opts;
  int nptes;
};
typedef struct pdht_config_s pdht_config_t;


/* MPI implementation for a PDHT table */
struct pdht_s{
  ht_t          *ht;
  pdht_hashfunc  hashfn;
  unsigned       elemsize;
  unsigned       keysize;
  int            fuckups;
  pdht_ptl_t     ptl;
  pthread_mutex_t *uthash_lock;
  uint64_t       counters[PDHT_MAX_COUNTERS];//only used by rank 0
  int            counter_count;
};
typedef struct pdht_s pdht_t;


/* legacy/compatibility definitions */
enum pdht_mode_e{
  PdhtModeStrict,
  PdhtModeBundled,
  PdhtModeAsync
};
typedef enum pdht_mode_e pdht_mode_t;

enum pdht_pmode_e{
  PdhtPendingPoll,
  PdhtPendingTrig
};
typedef enum pdht_pmode_e pdht_pmode_t;

enum pdht_status_e{
  PdhtStatusOK,
  PdhtStatusError,
  PdhtStatusNotFound,
  PdhtStatusCollision
};
typedef enum pdht_status_e pdht_status_t;

//declaring functions the user can use
pdht_t *pdht_create(int keysize, int elemsize, pdht_mode_t mode);
void pdht_free(pdht_t *dht);
void pdht_tune(unsigned opts, pdht_config_t *config);

//putget ops
pdht_status_t pdht_put(pdht_t *dht, void *key, void *value);
pdht_status_t pdht_get(pdht_t *dht, void *key, void *value);
pdht_status_t pdht_update(pdht_t *dht, void *key, void *value);
pdht_status_t pdht_persistent_get(pdht_t *dht, void *key, void *value);

//commsynch
void pdht_barrier(void);
void pdht_fence(pdht_t *dht);
int pdht_counter_init(pdht_t *ht, uint64_t initval);
void pdht_counter_reset(pdht_t *ht, int counter);
uint64_t pdht_counter_inc(pdht_t *ht, int counter, uint64_t val);


//util
void    pdht_print_all(pdht_t *dht);
void    pdht_print_stats(pdht_t *dht);
int     eprintf(const char *format, ...);
void    pdht_print_active(pdht_t *dht, void kprinter(void *key), void vprinter(void *value));
double  pdht_average_time(pdht_t *dht, pdht_timer_t timer);

//hash
void pdht_sethash(pdht_t *dht,pdht_hashfunc hfun);
void pdht_hash(pdht_t *dht, void *key, ptl_match_bits_t *mbits, uint32_t *ptindex, ptl_process_t *rank);
