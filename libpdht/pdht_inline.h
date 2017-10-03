/********************************************************/
/*                                                      */
/*  pdht_inline.h - PDHT inline functions               */
/*                                                      */
/*  author: d. brian larkins                            */
/*  created: 3/20/16                                    */
/*                                                      */
/********************************************************/


/**
 * @file
 * 
 * portals distributed hash table inline functions
 */

#pragma once

/**
 * pdht_get_rank - gets rank of calling process
 *  @returns rank of calling process
 */
static inline long pdht_get_rank(void) {
  return c->rank;
}


/* timing routines */

/**
 *  pdht_get_wtime - get a wall clock time for performance analysis
 */
static inline struct timespec pdht_get_wtime() {
  struct timespec ts;
#if __APPLE__
  clock_gettime(_CLOCK_MONOTONIC, &ts);
#else
  clock_gettime(CLOCK_MONOTONIC_RAW, &ts);
#endif // clock_gettime
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

// PDHT_READ_TIMER returns elapsed time in nanoseconds
#ifdef PDHT_USE_INTERNAL_TIMERS
  #define PDHT_START_TIMER(HT,TMR) PDHT_START_ATIMER(HT->stats.TMR)
  #define PDHT_STOP_TIMER(HT,TMR)  PDHT_STOP_ATIMER(HT->stats.TMR)
  #define PDHT_READ_TIMER(HT,TMR)  PDHT_READ_ATIMER(HT->stats.TMR)

  #define PDHT_READ_TIMER_USEC(HT,TMR)  PDHT_READ_TIMER(HT,TMR)/1000.0
  #define PDHT_READ_TIMER_MSEC(HT,TMR)  PDHT_READ_TIMER(HT,TMR)/(double)1e6
  #define PDHT_READ_TIMER_SEC(HT,TMR)   PDHT_READ_TIMER(HT,TMR)/(double)1e9
#else
  #define PDHT_START_TIMER(HT,TMR) 
  #define PDHT_STOP_TIMER(HT,TMR)  
  #define PDHT_READ_TIMER(HT,TMR) 0

  #define PDHT_READ_TIMER_USEC(HT,TMR) 0
  #define PDHT_READ_TIMER_MSEC(HT,TMR) 0
  #define PDHT_READ_TIMER_SEC(HT,TMR)  0
#endif // PDHT_USE_INTERNAL_TIMERS


#ifdef DEPRECATED

#define PDHT_START_TIMER(HT,TMR) HT->stats.TMR.last   = pdht_get_wtime();
#define PDHT_STOP_TIMER(HT,TMR) do {\
                                  HT->stats.TMR.temp = pdht_get_wtime();\
                                  HT->stats.TMR.total.tv_sec += (HT->stats.TMR.temp.tv_sec - HT->stats.TMR.last.tv_sec);\
                                  HT->stats.TMR.total.tv_nsec += (HT->stats.TMR.temp.tv_nsec - HT->stats.TMR.last.tv_nsec);\
                                } while (0)
// PDHT_READ_TIMER returns elapsed time in nanoseconds
#define PDHT_READ_TIMER(HT,TMR)  ((1000000000L * (HT->stats.TMR.total.tv_sec)) + HT->stats.TMR.total.tv_nsec)


/* blocking shortcuts */

/**
 * pdht_puti - put an integer keyed value into the DHT
 *  @param k integer key
 *  @param v DHT entry to put
 *  @returns status of the operation
 */
static inline pdht_status_t pdht_puti(pdht_t *dht, int k, void *v) {
   return pdht_put(dht, &k, sizeof(int), v);
}



/**
 * pdht_geti - get an integer keyed value from the DHT
 *  @param k integer key
 *  @param v DHT entry to get
 *  @returns status of the operation
 */
static inline pdht_status_t pdht_geti(pdht_t *dht, int k, void **v) {
   return pdht_get(dht, &k, sizeof(int), v);
}



/**
 * pdht_putf - put an double keyed value into the DHT
 *  @param k double key
 *  @param v DHT entry to put
 *  @returns status of the operation
 */
static inline pdht_status_t pdht_putf(pdht_t *dht, double k, void *v) {
   return pdht_put(dht, &k, sizeof(double), v);
}



/**
 * pdht_getf - get an double/float keyed value from the DHT
 *  @param k double key
 *  @param v DHT entry to get
 *  @returns status of the operation
 */
static inline pdht_status_t pdht_getf(pdht_t *dht, double k, void **v) {
   return pdht_get(dht, &k, sizeof(double), v);
}


/* non-blocking shortcuts */

/**
 * pdht_nbputi - put an integer keyed value into the DHT
 *  @param k integer key
 *  @param v DHT entry to put
 *  @returns status of the operation
 */
static inline pdht_status_t pdht_nbputi(pdht_t *dht, int k, void *v) {
   return pdht_put(dht, &k, sizeof(int), v);
}



/**
 * pdht_nbgeti - get an integer keyed value from the DHT
 *  @param k integer key
 *  @param v DHT entry to get
 *  @returns status of the operation
 */
static inline pdht_status_t pdht_nbgeti(pdht_t *dht, int k, void **v) {
   return pdht_get(dht, &k, sizeof(int), v);
}



/**
 * pdht_nbputf - put an double keyed value into the DHT
 *  @param k double key
 *  @param v DHT entry to put
 *  @returns status of the operation
 */
static inline pdht_status_t pdht_nbputf(pdht_t *dht, double k, void *v) {
   return pdht_nbput(dht, &k, sizeof(double), v);
}



/**
 * pdht_nbgetf - get an double/float keyed value from the DHT
 *  @param k double key
 *  @param v DHT entry to get
 *  @returns status of the operation
 */
static inline pdht_status_t pdht_nbgetf(pdht_t *dht, double k, void **v) {
   return pdht_nbget(dht, &k, sizeof(double), v);
}
#endif

static inline u_int64_t pdht_find_bucket(pdht_t *dht, void *p) {
  char *start = dht->ht;
  char *end   = p;
  u_int64_t ret = (end-start) / dht->entrysize;
  if ((start + (ret * dht->entrysize)) != end) {
    printf("%d: misaligned pointer: expected: %p found %p\n", c->rank, start + (ret * dht->entrysize), end); fflush(stdout);
  }
  return ret;
}

