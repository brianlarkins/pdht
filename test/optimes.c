#define _XOPEN_SOURCE 600
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>

#include <pdht.h>


#define NITER 1000

extern pdht_context_t *c;
int eprintf(const char *format, ...);

int main(int argc, char **argv);

#define START_TIMER(TMR) TMR.last = pdht_get_wtime();
#define STOP_TIMER(TMR) TMR.total += pdht_get_wtime() - TMR.last;
#define READ_TIMER(TMR) TMR.total

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

  START_TIMER(total);

  pdht_sethash(ht, localhash);
    
  // nothing should be in hash, so everything should be not found
  START_TIMER(lnotfound);
  for (int iter=0; iter < maxiters; iter++) {
    if (c->rank == 0) {
      pdht_get(ht, &key, val);
      key++;
    } 
  }
  STOP_TIMER(lnotfound);
 
  pdht_barrier();
  pdht_sethash(ht, remotehash);

  // repeat for remote checks
  START_TIMER(rnotfound);
  for (int iter=0; iter < maxiters; iter++) {
    if (c->rank == 0) {
      pdht_get(ht, &key, val);
      key++;
    } 
  }
  STOP_TIMER(rnotfound);
 
  // get times for local puts / gets
  pdht_barrier();
  key = 10;
  pdht_sethash(ht, localhash);

  START_TIMER(lput);
  for (int iter=0; iter < maxiters; iter++) {
    if (c->rank == 0) {
      pdht_put(ht, &key, val);
      key++;
    } 
  }
  STOP_TIMER(lput);

  pdht_barrier();
  key = 10;
  START_TIMER(lget);
  for (int iter=0; iter < maxiters; iter++) {
    if (c->rank == 0) {
      pdht_get(ht, &key, val);
      key++;
    } 
  }
  STOP_TIMER(lget);

  // times for remote puts 
  pdht_barrier();
  pdht_sethash(ht, remotehash);
  key = 100 + maxiters;

  START_TIMER(rput);
  for (int iter=0; iter < maxiters; iter++) {
    if (c->rank == 0) {
      pdht_put(ht, &key, val);
      key++;
    } 
  }
  STOP_TIMER(rput);

  pdht_barrier();
  key = 100 + maxiters;

  // remote gets
  START_TIMER(rget);
  for (int iter=0; iter < maxiters; iter++) {
    if (c->rank == 0) {
      pdht_get(ht, &key, val);
      key++;
    } 
  }
  STOP_TIMER(rget);

  pdht_barrier();

  // head of matchlist gets
  unsigned long headkey = 100+maxiters;
  // head of ME list gets
  START_TIMER(hget);
  for (int iter=0; iter < maxiters; iter++) {
    if (c->rank == 0) {
      pdht_get(ht, &headkey, val);
    } 
  }
  STOP_TIMER(hget);

  pdht_barrier();

  // tail of matchlist gets
  key--; // key was autoincremented at end of gets/loop
  START_TIMER(tget);
  for (int iter=0; iter < maxiters; iter++) {
    if (c->rank == 0) {
      pdht_get(ht, &key, val);
    } 
  }
  STOP_TIMER(tget);

  pdht_barrier();
  STOP_TIMER(total);

  pdht_print_stats(ht);
  eprintf("\n\nelemsize %lu iterations: %d elapsed time: %12.7f s\n", elemsize, maxiters, READ_TIMER(total));

  eprintf("  totl: local  get: %12.7f s  put: %12.7f s   -- notfound: %12.7f s\n", 
     READ_TIMER(lget),
     READ_TIMER(lput),
     READ_TIMER(lnotfound));
  eprintf("  totl: remote get: %12.7f s  put: %12.7f s   -- notfound: %12.7f s \n", 
     READ_TIMER(rget),
     READ_TIMER(rput),
     READ_TIMER(rnotfound));

  eprintf("\n");
  eprintf("  unit: local  get: %12.7f ms put: %12.7f ms  -- notfound: %12.7f ms\n", 
     (READ_TIMER(lget)/(double)maxiters)*1000, 
     (READ_TIMER(lput)/(double)maxiters)*1000, 
     (READ_TIMER(lnotfound)/(double)maxiters)*1000);
  eprintf("  unit: remote get: %12.7f ms put: %12.7f ms  -- notfound: %12.7f ms \n", 
     (READ_TIMER(rget)/(double)maxiters)*1000, 
     (READ_TIMER(rput)/(double)maxiters)*1000, 
     (READ_TIMER(rnotfound)/(double)maxiters)*1000);
  eprintf("  unit: matchlist head: %12.7f ms matchlist tail: %12.7f ms\n", 
     (READ_TIMER(hget)/(double)maxiters)*1000, 
     (READ_TIMER(tget)/(double)maxiters)*1000);
done:
  pdht_free(ht);
  free(val);
}
