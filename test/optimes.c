#define _XOPEN_SOURCE 600
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>

#include <pdht.h>


#define NITER 1000

extern pdht_context_t *c;
int eprintf(const char *format, ...);

int main(int argc, char **argv);

//#define START_TIMER(TMR) TMR.last = pdht_get_wtime();
//#define STOP_TIMER(TMR) TMR.total += pdht_get_wtime() - TMR.last;
//#define READ_TIMER(TMR) TMR.total

void localhash(pdht_t *dht, void *key, ptl_match_bits_t *mbits, uint32_t *ptindex, ptl_process_t *rank) {
  (*rank).rank = 0;
  *mbits = *(unsigned long *)key;
  //*ptindex = *(unsigned long *)key % dht->nptes;
  *ptindex = 1;
}

void remotehash(pdht_t *dht, void *key, ptl_match_bits_t *mbits, uint32_t *ptindex, ptl_process_t *rank) {
  (*rank).rank = 1;
  *mbits = *(unsigned long *)key;
  //*ptindex = *(unsigned long *)key % dht->nptes;
  *ptindex = 1;
}


int main(int argc, char **argv) {
  pdht_t *ht;
  pdht_status_t ret;
  size_t elemsize = sizeof(unsigned long);
  pdht_timer_t ltimer;
  unsigned long key = 10; // whatever, just increasing monotonically
  void *val = NULL;
  int opt, maxiters = NITER;
  pdht_timer_t lnotfound, rnotfound;
  pdht_timer_t lget, rget;
  pdht_timer_t lput, rput;
  pdht_timer_t hget, tget;
  pdht_timer_t total;
 
  while ((opt = getopt(argc, argv, "i:s:")) != -1) {
    switch (opt) {
    case 'i':
       maxiters = atoi(optarg);
       break;
    case 's':
       elemsize = atoi(optarg);
       break;
    } 
  }

  val = malloc(elemsize);
  memset(val,0,elemsize);
  

  // create hash table
  ht = pdht_create(sizeof(unsigned long), elemsize, PdhtModeStrict);


  if (c->size != 2) {
    if (c->rank == 0) {
      printf("requires two (and only two) processes to run\n");
    }
    goto done;
  }

  PDHT_START_ATIMER(total);
  sleep(5);
  PDHT_STOP_ATIMER(total);
  printf("%12.5f sec %12.5f ms %12.5f us %12.5f ns\n",
     PDHT_READ_ATIMER_SEC(total),
     PDHT_READ_ATIMER_MSEC(total),
     PDHT_READ_ATIMER_USEC(total),
     (double)PDHT_READ_ATIMER(total));

  PDHT_START_ATIMER(total);

  pdht_sethash(ht, localhash);
    
  // nothing should be in hash, so everything should be not found
  PDHT_START_ATIMER(lnotfound);
  for (int iter=0; iter < maxiters; iter++) {
    if (c->rank == 0) {
      pdht_get(ht, &key, val);
      key++;
    } 
  }
  PDHT_STOP_ATIMER(lnotfound);
 
  pdht_barrier();
  pdht_sethash(ht, remotehash);

  // repeat for remote checks
  PDHT_START_ATIMER(rnotfound);
  for (int iter=0; iter < maxiters; iter++) {
    if (c->rank == 0) {
      pdht_get(ht, &key, val);
      key++;
    } 
  }
  PDHT_STOP_ATIMER(rnotfound);
 
  // get times for local puts / gets
  pdht_barrier();
  key = 10;
  pdht_sethash(ht, localhash);

  PDHT_START_ATIMER(lput);
  for (int iter=0; iter < maxiters; iter++) {
    if (c->rank == 0) {
      pdht_put(ht, &key, val);
      key++;
    } 
  }
  PDHT_STOP_ATIMER(lput);

  pdht_barrier();
  key = 10;
  PDHT_START_ATIMER(lget);
  for (int iter=0; iter < maxiters; iter++) {
    if (c->rank == 0) {
      pdht_get(ht, &key, val);
      key++;
    } 
  }
  PDHT_STOP_ATIMER(lget);

  // times for remote puts 
  pdht_barrier();
  pdht_sethash(ht, remotehash);
  key = 100 + maxiters;

  PDHT_START_ATIMER(rput);
  for (int iter=0; iter < maxiters; iter++) {
    if (c->rank == 0) {
      pdht_put(ht, &key, val);
      key++;
    } 
  }
  PDHT_STOP_ATIMER(rput);

  pdht_barrier();
  key = 100 + maxiters;

  // remote gets
  PDHT_START_ATIMER(rget);
  for (int iter=0; iter < maxiters; iter++) {
    if (c->rank == 0) {
      pdht_get(ht, &key, val);
      key++;
    } 
  }
  PDHT_STOP_ATIMER(rget);

  pdht_barrier();

  // head of matchlist gets
  unsigned long headkey = 100+maxiters;
  // head of ME list gets
  PDHT_START_ATIMER(hget);
  for (int iter=0; iter < maxiters; iter++) {
    if (c->rank == 0) {
      pdht_get(ht, &headkey, val);
    } 
  }
  PDHT_STOP_ATIMER(hget);

  pdht_barrier();

  // tail of matchlist gets
  key--; // key was autoincremented at end of gets/loop
  PDHT_START_ATIMER(tget);
  for (int iter=0; iter < maxiters; iter++) {
    if (c->rank == 0) {
      pdht_get(ht, &key, val);
    } 
  }
  PDHT_STOP_ATIMER(tget);

  pdht_barrier();
  PDHT_STOP_ATIMER(total);

  pdht_print_stats(ht);
  eprintf("\n\nelemsize %lu iterations: %d elapsed time: %12.7f s\n", elemsize, maxiters, PDHT_READ_ATIMER(total)/(double)1e9);

  eprintf("  totl: local  get: %12.7f s  put: %12.7f s   -- notfound: %12.7f s\n", 
     PDHT_READ_ATIMER_SEC(lget),
     PDHT_READ_ATIMER_SEC(lput),
     PDHT_READ_ATIMER_SEC(lnotfound));
  eprintf("  totl: remote get: %12.7f s  put: %12.7f s   -- notfound: %12.7f s \n", 
     PDHT_READ_ATIMER_SEC(rget),
     PDHT_READ_ATIMER_SEC(rput),
     PDHT_READ_ATIMER_SEC(rnotfound));

  eprintf("\n");
  eprintf("  unit: local  get: %12.7f us put: %12.7f us  -- notfound: %12.7f us\n", 
     (PDHT_READ_ATIMER_USEC(lget)/(double)maxiters), 
     (PDHT_READ_ATIMER_USEC(lput)/(double)maxiters), 
     (PDHT_READ_ATIMER_USEC(lnotfound)/(double)maxiters));
  eprintf("  unit: remote get: %12.7f us put: %12.7f us  -- notfound: %12.7f us \n", 
     (PDHT_READ_ATIMER_USEC(rget)/(double)maxiters), 
     (PDHT_READ_ATIMER_USEC(rput)/(double)maxiters), 
     (PDHT_READ_ATIMER_USEC(rnotfound)/(double)maxiters));
  eprintf("  unit: matchlist head: %12.7f us matchlist tail: %12.7f us\n", 
     (PDHT_READ_ATIMER_USEC(hget)/(double)maxiters), 
     (PDHT_READ_ATIMER_USEC(tget)/(double)maxiters));
done:
  pdht_free(ht);
  free(val);
}
