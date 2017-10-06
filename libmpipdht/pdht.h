#include <pthread.h>

#include <stdlib.h>
#include <stdio.h>
#include "uthash.h"
#include <string.h>
#include <signal.h>
#include <stdarg.h>
#include <time.h>
#include <portals4.h>

#include "mpi.h"



struct pdht_timer_s{
  double total;
  double last;
  double temp;

};
typedef struct pdht_timer_s pdht_timer_t;



#define PDHT_PTALLOC_OPTIONS PTL_PT_MATCH_UNORDERED
#define PDHT_MAX_TABLES 20
#define PDHT_MAXKEYSIZE 32


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

//typedef ptl stuff so things can still work
typedef uint64_t ptl_match_bits_t;




enum pdht_local_gets_e{
  PdhtOptimized,
  PdhtRegular
};
typedef enum pdht_local_gets_e pdht_local_gets_t;


struct pdht_s;

//delcaring what a pdht_hashfunc is
typedef void (*pdht_hashfunc)(struct pdht_s *dht, void *key, ptl_match_bits_t *mbits, uint32_t *ptindex, ptl_process_t *rank);

//what type of message i am getting
typedef enum {pdhtGet, pdhtPut, pdhtStop} msg_type;






//shhhh im trying to be portals
struct fake_portals_s{
  int nptes;//idk?

};
typedef struct fake_portals_s fake_portals_t;


//slimmed down pdht_context
struct pdht_context_s{
  int rank;
  struct pdht_s *hts[PDHT_MAX_TABLES];
  int dhtcount;
  int size; 
  pthread_t tid;
  int thread_active;

  MPI_Datatype msgType;
};
typedef struct pdht_context_s pdht_context_t;


//uthash struct
struct ht_s{
  uint64_t key;
  void *value;
  UT_hash_handle hh;
  
};
typedef struct ht_s ht_t;


//initial message for passing over mpi
struct message_s{
  msg_type type;
  int rank;
  int ht_index;
};
typedef struct message_s message_t;

struct pdht_config_s{
  int nptes;
  int pendmode;
  long unsigned maxentries;
  int quiet;
  pdht_local_gets_t local_gets;
  long unsigned pendq_size;
  long unsigned ptalloc_opts;

};
typedef struct pdht_config_s pdht_config_t;


//slimmed down pdht
struct pdht_s{
  ht_t *ht;
  pdht_hashfunc hashfn;
  unsigned elemsize;
  unsigned keysize;
  fake_portals_t ptl;

  int fuckups;
};
typedef struct pdht_s pdht_t;


//tuning stuff that is not used
#define PDHT_TUNE_NPTES     0x01
#define PDHT_TUNE_PMODE     0x02
#define PDHT_TUNE_ENTRY     0x04
#define PDHT_TUNE_PENDQ     0x08
#define PDHT_TUNE_PTOPT     0x10
#define PDHT_TUNE_QUIET     0x20
#define PDHT_TUNE_GETS      0x40
#define PDHT_TUNE_ALL       0xffffffff



//need this for declaration, never actually used, honestly dont know what it does in normal pdht
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
void pdht_put(pdht_t *dht, void *key, void *value);
void pdht_get(pdht_t *dht, void *key, void *value);
void pdht_update(pdht_t *dht, void *key, void *value);

//commsynch
void pdht_barrier(void);
void pdht_fence(pdht_t *dht);

//util
void pdht_print_all(pdht_t *dht);
void pdht_print_stats(pdht_t *dht);
int eprintf(const char *format, ...);

//hash
void pdht_sethash(pdht_t *dht,pdht_hashfunc hfun);
void pdht_hash(pdht_t *dht, void *key, ptl_match_bits_t *mbits, uint32_t *ptindex, ptl_process_t *rank);
