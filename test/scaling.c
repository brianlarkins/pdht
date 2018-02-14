#define _XOPEN_SOURCE 600
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>

#include <pdht.h>

#define NITER 1000

extern pdht_context_t *c;
int eprintf(const char *format, ...);

int main(int argc, char **argv);

/*
void localhash(pdht_t *dht, void *key, ptl_match_bits_t *mbits, uint32_t *ptindex, ptl_process_t *rank) {
  (*rank).rank = 0;
  *mbits = *(unsigned long *)key;
  ptindex = *(unsigned long *)key % dht->ptl.nptes;
  *ptindex = 0;
}

void remotehash(pdht_t *dht, void *key, ptl_match_bits_t *mbits, uint32_t *ptindex, ptl_process_t *rank) {
  (*rank).rank = 1;
  *mbits = *(unsigned long *)key;
  ptindex = *(unsigned long *)key % dht->ptl.nptes;
  *ptindex = 0;
}
*/


/*
 * bench specific hash function
 *   - key is expected to be a long integer
 *   - key % P == target rank
 */
void ahash(pdht_t *dht, void *key, ptl_match_bits_t *mbits, uint32_t *ptindex, ptl_process_t *rank) {
  *mbits = *(unsigned long *)key;
  *ptindex = *(unsigned long *)key % dht->ptl.nptes;
  (*rank).rank = *(unsigned long *)key % c->size;
}


int main(int argc, char **argv) {
  pdht_t *ht;
  pdht_status_t ret;
  size_t elemsize = sizeof(unsigned long);
  unsigned long key = 0; // whatever, just increasing monotonically
  void *val = NULL;
  int opt;
  int iters = 1, numentries = NITER;
  pdht_timer_t gtimer, ptimer, total;
  double local[3], avg[3], min[3], max[3];
  int using_mpi = 0;

  pdht_config_t cfg;
  cfg.nptes        = 1;
  cfg.pendmode     = PdhtPendingTrig;
  cfg.maxentries   =  102000; // 250000;
  cfg.pendq_size   =  51000; // 100000;
  cfg.ptalloc_opts = 0;
  cfg.local_gets   = 0;
  cfg.rank         = PDHT_DEFAULT_RANK_HINT;

  while ((opt = getopt(argc, argv, "hi:ls:pquv:U")) != -1) {
    switch (opt) {
      case 'd':
        cfg.quiet = 1;
        break;
      case 'h':
        printf("usage: scaling -htu -i <iters> -s <elemsize> -v <numentries>\n");
        printf("\t-d be quieter\n");
        printf("\t-h this message\n");
        printf("\t-i # iterations\n");
        printf("\t-s element size\n");
        printf("\t-p use polling ops (instead of triggered)\n");
        printf("\t-u use unordered matching (hashing) in Portals\n");
        printf("\t-v number of entries each processor inserts/reads\n");
        printf("\t-l optimizes local gets\n");
        break;
      case 'i':
        iters = atoi(optarg);
        break;
      case 'l':
        cfg.local_gets = PdhtSearchLocal;
        break;
      case 'p':
        cfg.pendmode = PdhtPendingPoll;
        break;
      case 'q':
        cfg.quiet = 1;
        break;
      case 's':
        elemsize = atoi(optarg);
        break;
      case 'u':
        cfg.ptalloc_opts = PTL_PT_MATCH_UNORDERED;
        break;
      case 'v':
        numentries = atoi(optarg);
        break;
      case 'U':
        setenv("PTL_IGNORE_UMMUNOTIFY", "1",1);
        break;
    } 
  }

  // setenv("PTL_DISABLE_MEM_REG_CACHE","1",1);
  //setenv("PTL_LOG_LEVEL","3",1);
  setenv("PTL_DEBUG","1",1);
  setenv("PTL_PROGRESS_NOSLEEP","1",1);

  val = malloc(elemsize); // alloc our token value object

  // create hash table
  pdht_tune(PDHT_TUNE_ALL, &cfg);
  
  ht = pdht_create(sizeof(unsigned long), elemsize, PdhtModeStrict);
  
  eprintf("starting run with %d processes, each with %d entries (total reads == %d)\n", 
      c->size, numentries,iters*numentries);
  pdht_barrier();

  PDHT_START_ATIMER(total);

  //pdht_sethash(ht, ahash);

  // each process puts numentries elements into distributed hash
  key = c->rank;
  memset(&ptimer, 0, sizeof(ptimer));
  PDHT_START_ATIMER(ptimer);
  for (int iter=0; iter < numentries; iter++) {
    //val = key + 10;
    pdht_put(ht, &key, val);
    key += c->size;
    if ((iter % 1000) == 0) eprintf(".");
  }
  PDHT_STOP_ATIMER(ptimer);

  pdht_fence(ht);
  eprintf(".\n");

  eprintf("starting fetches\n");

  memset(&gtimer, 0, sizeof(gtimer));
  // now we time getting numentries
  for (int it=0; it<iters; it++) {
    
    // set starting key to be our left neighbor
    key = (c->rank != 0) ? c->rank-1 : c->size - 1;

    PDHT_START_ATIMER(gtimer);
    for (int iter=0; iter < numentries; iter++) {
      pdht_get(ht, &key, val);
      key += c->size;
    }
    PDHT_STOP_ATIMER(gtimer);
    pdht_barrier();
  }
  PDHT_STOP_ATIMER(total);

  local[0] = PDHT_READ_ATIMER_MSEC(ptimer);
  local[1] = PDHT_READ_ATIMER_MSEC(gtimer);
  local[2] = PDHT_READ_ATIMER_MSEC(total);

  printf("rank %d: put: %12.7f get: %12.7f total: %12.7f\n", c->rank, local[0], local[1], local[2]);
  pdht_allreduce(local, &avg, PdhtReduceOpSum, DoubleType, 3);
  pdht_allreduce(local, &min, PdhtReduceOpMin, DoubleType, 3);
  pdht_allreduce(local, &max, PdhtReduceOpMax, DoubleType, 3);

  // fix up averages
  for (int i=0; i<3; i++)  {
#ifndef MPI
    avg[i] = avg[i] / (double)c->size;
#else
    int mpisize = c->size / 2;
    avg[i] = avg[i] / (double)(mpisize);
#endif
  }

  eprintf("put times         : %12.7f / %12.7f / %12.7f  ms (avg/min/max)\n", 
      avg[0], min[0], max[0]);

  eprintf("get times         : %12.7f / %12.7f / %12.7f  ms (avg/min/max)\n", 
      avg[1], min[1], max[1]);

  eprintf("total elapsed time: %12.7f / %12.7f / %12.7f  ms (avg/min/max)\n", 
      avg[2], min[2], max[2]);

  eprintf("put throughput    : %12.3f MBps\n", (c->size*elemsize*numentries)/(avg[0]*10e6));
  eprintf("get throughput    : %12.3f MBps %12.3f reads/sec\n", 
	(c->size*elemsize*numentries*iters)/(avg[1]*10e6),
	(c->size*numentries*iters)/(avg[1]));
  pdht_print_stats(ht);

done:
  pdht_free(ht);
}
