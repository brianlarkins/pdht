#define _XOPEN_SOURCE 600
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>

#include <pdht.h>

/*
 * TO RUN:
 */

#define NITER 1000

// entries *per-node*
int maxentries = 1;

extern pdht_context_t *c;
int eprintf(const char *format, ...);

int main(int argc, char **argv);

void ahash(pdht_t *dht, void *key, ptl_match_bits_t *mbits, uint32_t *ptindex, ptl_process_t *rank) {
  *mbits = *(unsigned long *)key;
  *ptindex = *(unsigned long *)key % dht->nptes;
  (*rank).rank = *(unsigned long *)key / maxentries;
}


int main(int argc, char **argv) {
  pdht_t *ht;
  pdht_status_t ret;
  size_t elemsize = sizeof(unsigned long);
  pdht_timer_t ltimer;
  unsigned long key = 0; // whatever, just increasing monotonically
  unsigned long val = 0;
  int opt;
  pdht_timer_t gtimer,total;

  //setenv("PTL_DISABLE_MEM_REG_CACHE","1",1);
  while ((opt = getopt(argc, argv, "v:s:")) != -1) {
    switch (opt) {
      case 's':
        elemsize = atoi(optarg);
        break;
      case 'v':
        maxentries = atoi(optarg);
        break;
    } 
  }

  // create hash table
  ht = pdht_create(sizeof(unsigned long), sizeof(unsigned long), PdhtModeStrict);

  eprintf("starting run with %d processes, each with %d entries\n", c->size, maxentries);
  pdht_barrier();

  PDHT_START_ATIMER(total);

  pdht_sethash(ht, ahash);
  
  key = c->rank;
  val = 10 + key;
  pdht_insert(ht, key, key % ht->nptes, &key, &val);

  pdht_barrier();

  // now we time getting maxentries

  for (int rank=0; rank<c->size; rank++) {
    pdht_get(ht, &rank, &val);
    printf("%d: key: %lu val: %lu\n", c->rank, key, val); fflush(stdout);
  }

  //printf("%d: %12.7f ms\n", c->rank,  (READ_TIMER(gtimer)*1000));

  pdht_barrier();

  PDHT_STOP_ATIMER(total);
  eprintf("total elapsed time: %12.7f\n", PDHT_READ_ATIMER_SEC(total));

  pdht_print_stats(ht);

done:
  pdht_free(ht);
}
