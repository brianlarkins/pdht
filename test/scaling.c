#define _XOPEN_SOURCE 600
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>

#include <pdht.h>

/*
 * TO RUN:
 *   - need to make sure that PDHT_DEFAULT_TABLE_SIZE is at least 100K * NUM_PTES
 *   - set PDHT_DEFAULT_NUM_PTES to 2
 *   - set PDHT_DEFAULT_TABLE_SIZE 250000
 */

#define NITER 1000

int maxentries = NITER;

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

void ahash(pdht_t *dht, void *key, ptl_match_bits_t *mbits, uint32_t *ptindex, ptl_process_t *rank) {
  *mbits = *(unsigned long *)key;
  *ptindex = *(unsigned long *)key % dht->nptes;
  (*rank).rank = *(unsigned long *)key % c->size;
}


int main(int argc, char **argv) {
  pdht_t *ht;
  pdht_status_t ret;
  size_t elemsize = sizeof(unsigned long);
  pdht_timer_t ltimer;
  unsigned long key = 0; // whatever, just increasing monotonically
  unsigned long val = 0;
  int opt;
  int iters = 1;
  pdht_timer_t gtimer,total;

  setenv("PTL_DISABLE_MEM_REG_CACHE","1",1);
  while ((opt = getopt(argc, argv, "i:s:v:")) != -1) {
    switch (opt) {
      case 'i':
        iters = atoi(optarg);
      case 's':
        elemsize = atoi(optarg);
        break;
      case 'v':
        maxentries = atoi(optarg);
        break;
    } 
  }


  // create hash table
  ht = pdht_create(sizeof(unsigned long), elemsize, PdhtModeStrict);

  eprintf("starting run with %d processes, each with %d entries (total reads == %d)\n", c->size, maxentries,iters*maxentries);
  pdht_barrier();

  START_TIMER(total);

  pdht_sethash(ht, ahash);

  // each process puts maxentries elements into distributed hash
  key = c->rank;
  for (int iter=0; iter < maxentries; iter++) {
    val = key + 10;
    pdht_insert(ht, key, key % ht->nptes, &key, &val);
    key += c->size;
    if ((iter % 1000) == 0) eprintf(".");
  }
  sleep(5);
  pdht_barrier();
  eprintf(".\n");

  eprintf("starting fetches\n");

  // now we time getting maxentries
  for (int it=0; it<iters; it++) {
    key = (c->rank != 0) ? c->rank-1 : c->size - 1;
    START_TIMER(gtimer);
    for (int iter=0; iter < maxentries; iter++) {
      pdht_get(ht, &key, &val);
      key += c->size;
    }
    STOP_TIMER(gtimer);
    pdht_barrier();
  }

  printf("%d: %12.7f ms\n", c->rank,  (READ_TIMER(gtimer)*1000));

  pdht_barrier();

  STOP_TIMER(total);
  eprintf("total elapsed time: %12.7f\n", READ_TIMER(total));

  pdht_print_stats(ht);

done:
  pdht_free(ht);
}
