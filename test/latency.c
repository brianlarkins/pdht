#include <sys/time.h>
#include <sys/resource.h>

#include <pdht.h>

#define NITER 10000

// polling
// 20000 = unit latency of 14.183 ms (283s run)
// 10000 = unit latency of 14.173 ms (141s run)
// 1000  = unit latency of  9.84 ms (9.8s run)
// 100   = unit latency of 11.6 ms (1.2s run)
// 10    = unit latency of 15.8 ms (.16s run)
/*
pdht statistics for rank 0
    puts:              20000     gets:          2175542
    collisions:            0     notfound:      2155542
    puttime:       96.8309 sec    gets:       186.7387 sec
    t1:             1.1660 sec    t2:         186.2566 sec
    t3:             0.0151 sec    t4:          96.8288 sec
    t5:             0.0000 sec    t6:         283.6596 sec
*/

// triggered
// 20000 = unit latency of 5.87 ms (117.4s)
// 10000 = unit latency of 5.77 ms (57.7s)
// 1000  = unit latency of 4.77 ms (4.7s)
// 100   = unit latency of 4.49 ms (.45s)
// 10    = unit latency of 5.05 ms (.05s)

/*
pdht statistics for rank 0
    puts:              20000     gets:            20000
    collisions:            0     notfound:            0
    puttime:      113.9631 sec    gets:         3.4240 sec
    t1:             0.0132 sec    t2:           3.4215 sec
    t3:             0.0165 sec    t4:         113.9610 sec
    t5:             0.0000 sec    t6:         117.3898 sec
*/

#define ELEMSIZE  sizeof(unsigned long)

extern pdht_context_t *c;
int eprintf(const char *format, ...);

int main(int argc, char **argv);


void fakehash(pdht_t *dht, void *key, ptl_match_bits_t *mbits, uint32_t *ptindex, ptl_process_t *rank) {
  (*rank).rank = 1;
  *mbits = *(unsigned long *)key;
  *ptindex = *(unsigned long *)key % dht->nptes;
}


int main(int argc, char **argv) {
  pdht_t *ht;
  pdht_status_t ret;
  size_t elemsize = ELEMSIZE;
  pdht_timer_t ltimer;
  unsigned long key = 10; // whatever, just increasing monotonically
  unsigned long val = 0;

  // create hash table
  ht = pdht_create(sizeof(unsigned long), elemsize, PdhtModeStrict);


  if (c->size != 2) {
    if (c->rank == 0) {
      printf("requires two (and only two) processes to run\n");
    }
    goto done;
  }

  pdht_sethash(ht, fakehash);

  for (int iter=0; iter < NITER; iter++) {

    if (c->rank == 0) {

      PDHT_START_TIMER(ht,t6);
      ret = pdht_put(ht, &key, &key);
 
      while (pdht_get(ht, &key, &val) != PdhtStatusOK)
        ;

      PDHT_STOP_TIMER(ht,t6);
      key++;

    } 

    pdht_barrier();
  }


  pdht_barrier();
  pdht_print_stats(ht);
  eprintf("\n\ntotal latency between put and get:\t %12.7f s\n", PDHT_READ_TIMER_SEC(ht,t6));
  eprintf("unit latency between put and get:\t%12.7f ms\n", (PDHT_READ_TIMER_MSEC(ht,t6)/(double)NITER));

done:
  pdht_free(ht);
}
