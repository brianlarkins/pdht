
#define _XOPEN_SOURCE 600

#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>

#include <pdht.h>


#define NITER 10000

extern pdht_context_t *c;
int eprintf(const char *format, ...);

int main(int argc, char **argv);


void localhash(pdht_t *dht, void *key, ptl_match_bits_t *mbits, uint32_t *ptindex, ptl_process_t *rank) {
  (*rank).rank = 0;
  *mbits = *(unsigned long *)key;
  *ptindex = *(unsigned long *)key % dht->ptl.nptes;
  //*ptindex = 1;
}

void remotehash(pdht_t *dht, void *key, ptl_match_bits_t *mbits, uint32_t *ptindex, ptl_process_t *rank) {
#if MPI
  (*rank).rank = 2;
#else
  (*rank).rank = 1;
#endif
  
  *mbits = *(unsigned long *)key;
  *ptindex = *(unsigned long *)key % dht->ptl.nptes;
  //*ptindex = 1;
}


int main(int argc, char **argv) {
  pdht_t *ht;
  pdht_status_t ret;
  size_t elemsize = sizeof(unsigned long);
  pdht_timer_t ltimer;
  unsigned long key = 1; // whatever, just increasing monotonically
  void *val = NULL;
  int opt, maxentries, maxiters;
  int rawmode = 0;
  int mlistentry[] = { 1, 10, 100, 1000, 10000, 50000 };
  pdht_timer_t rgets[6];
  pdht_timer_t lgets[6];
  pdht_timer_t lnotfound, rnotfound;
  pdht_timer_t lput, rput;
  pdht_timer_t lupdate, rupdate;
  pdht_timer_t total;
  pdht_config_t cfg;

  maxentries = NITER;
  maxiters   = NITER;
  cfg.nptes        = 1;
  cfg.pendmode     = PdhtPendingTrig;
  cfg.maxentries   = maxentries < 100000 ? 101000 : 2*maxentries;
  cfg.pendq_size   = maxentries < 100000 ? 51000 : maxentries+1;
  cfg.ptalloc_opts = 0;
  cfg.quiet        = 0;
  cfg.local_gets   = PdhtRegular;
  while ((opt = getopt(argc, argv, "dhi:n:s:puUl")) != -1) {
    switch (opt) {
      case 'd':
        rawmode = 1;
        cfg.quiet = 1;
        break;
      case 'h':
        fprintf(stderr,"usage:\t%s [-dhuU] [-i iters] [-s size]\n", argv[0]);
        fprintf(stderr,"\t\t-d - print data in raw form for GNUPlot\n");
        fprintf(stderr,"\t\t-h - help\n");
        fprintf(stderr,"\t\t-i iters - # of iterations to run\n");
        fprintf(stderr,"\t\t-s size - size in bytes of PDHT object\n");
        fprintf(stderr,"\t\t-p - use polling thread for inserts\n");
        fprintf(stderr,"\t\t-u - use unordered matching\n");
        fprintf(stderr,"\t\t-U - ignore ummunotify\n");
        exit(0);
        break;
      case 'i':
        maxiters = atoi(optarg);
        break;
      case 'n':
        maxentries = atoi(optarg);
        cfg.maxentries   = maxentries < 100000 ? 101000 : 2*maxentries;
        cfg.pendq_size   = maxentries < 100000 ? 51000 : maxentries+1;
        break;
      case 's':
        elemsize = atoi(optarg);
        break;
      case 'p':
        cfg.pendmode = PdhtPendingPoll;
        break;
      case 'u':
        cfg.ptalloc_opts = PTL_PT_MATCH_UNORDERED;
        
        break;
      case 'U':
        setenv("PTL_IGNORE_UMMUNOTIFY", "1",1);
        break;
      case 'l':
        cfg.local_gets = PdhtOptimized;
        break;
    } 
  }

  val = malloc(elemsize);
  memset(val,0,elemsize);


  // create hash table
  pdht_tune(PDHT_TUNE_ALL, &cfg);
  ht = pdht_create(sizeof(unsigned long), elemsize, PdhtModeStrict);
  pdht_barrier();
/*
  if (c->size != 2) {
    if (c->rank == 1) {
      printf("requires two (and only two) processes to run\n");
    }
    goto done;
  }
*/

  pdht_sethash(ht, localhash);

  PDHT_INIT_ATIMER(total);
  PDHT_INIT_ATIMER(lnotfound);
  PDHT_INIT_ATIMER(rnotfound);
  PDHT_INIT_ATIMER(lput);
  PDHT_INIT_ATIMER(rput);
  PDHT_INIT_ATIMER(lupdate);
  PDHT_INIT_ATIMER(rupdate);

  PDHT_START_ATIMER(total);

  // TIMING: latency for local not found elements
  // nothing should be in hash, so everything should be not found
  key = 1; // local uses odds
  if (c->rank == 0) {
    PDHT_START_ATIMER(lnotfound);
    for (int iter=0; iter < maxentries; iter++) {
      pdht_get(ht, &key, val);
      key+=2;
    } 
    PDHT_STOP_ATIMER(lnotfound);
  }

  // TIMING: latency for remote not found elements
  pdht_barrier();
  key = 2; // remote uses evens
  pdht_sethash(ht, remotehash);
  // repeat for remote checks
  if (c->rank == 0) {
    PDHT_START_ATIMER(rnotfound);
    for (int iter=0; iter < maxentries; iter++) {
      pdht_get(ht, &key, val);
      key+=2;
    } 
    PDHT_STOP_ATIMER(rnotfound);
  }


  // timing put / update / get operations

  // TIMING: latency for local put operations
  pdht_barrier();
  key = 1; // local uses odds
  pdht_sethash(ht, localhash);
  if (c->rank == 0) {
    PDHT_START_ATIMER(lput);
    for (int iter=0; iter < maxentries; iter++) {
      pdht_put(ht, &key, val);
      key+=2;
    } 
    PDHT_STOP_ATIMER(lput);
    if (!rawmode)
      printf("last local: %ld\n",key-2);
  }

  // TIMING: latency for local update operations
  // NOTE: only checking for updating first matchlist entry
  pdht_barrier();
  key = 1;
  memset(val,1,elemsize);
  if (c->rank == 0) {
    PDHT_START_ATIMER(lupdate);
    for (int iter=0; iter < maxiters; iter++) {
      pdht_update(ht, &key, val);
    } 
    PDHT_STOP_ATIMER(lupdate);
  }

  // TIMING: latency for remote put operations
  pdht_barrier();
  pdht_sethash(ht, remotehash);
  key = 2; // remote uses evens
  memset(val,0,elemsize);
  if (c->rank == 0) {
    PDHT_START_ATIMER(rput);
    for (int iter=0; iter < maxentries; iter++) {
      pdht_put(ht, &key, val);
      key+=2;
    }
    PDHT_STOP_ATIMER(rput);
    if (!rawmode)
      printf("last remote: %ld\n",key-2);
  } 
  // TIMING: latency for remote update operations
  // NOTE: only checking for updating first matchlist entry
  pdht_barrier();
  key = 2;
  memset(val,1,elemsize);
  if (c->rank == 0) {
    PDHT_START_ATIMER(rupdate);
    for (int iter=0; iter < maxiters; iter++) {
      pdht_update(ht, &key, val);
    } 
    PDHT_STOP_ATIMER(rupdate);
  }


  // TIMING: gets for local and remote elements at sample points in matchlist
  pdht_barrier();

  // for 1, 10, 100, 1000, 10000, 100000
  for (int e=0; e<6; e++) {
    PDHT_INIT_ATIMER(rgets[e]);
    PDHT_INIT_ATIMER(lgets[e]);

    // time local gets
    pdht_sethash(ht, localhash);
    key = 2 * mlistentry[e] - 1; // e.g. 100th entry is (2*100-1) == 199
    if (c->rank == 0) {
      PDHT_START_ATIMER(lgets[e]);

      for (int iter=0; iter<maxiters; iter++) {
        pdht_get(ht,&key,val);
      }
      PDHT_STOP_ATIMER(lgets[e]);
    }

    pdht_barrier();

    // time remote gets
    pdht_sethash(ht, remotehash);
    key = 2 * mlistentry[e]; // e.g. 100th entry is (2*100) = 200
    if (c->rank == 0) {
      PDHT_START_ATIMER(rgets[e]);
      for (int iter=0; iter<maxiters; iter++) {
        pdht_get(ht,&key,val);
      }
      PDHT_STOP_ATIMER(rgets[e]);
    }
  }

#if 0
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

#endif 

  pdht_barrier();
  PDHT_STOP_ATIMER(total);

  //pdht_print_stats(ht);

  if (!rawmode) {
    // element size, # of elements, iterations
    eprintf("\n\nelemsize %lu elements: %d iterations: %d elapsed time: %12.7f s\n", 
        elemsize, maxentries, maxiters, PDHT_READ_ATIMER_SEC(total));

    // not found
    eprintf("  unit: local notfound: %12.7f us remote notfound %12.7f us\n",
        (PDHT_READ_ATIMER_USEC(lnotfound)/(double)maxentries),
        (PDHT_READ_ATIMER_USEC(rnotfound)/(double)maxentries));
    eprintf("  unit: local put:      %12.7f us remote put:     %12.7f us\n",
        (PDHT_READ_ATIMER_USEC(lput)/(double)maxentries), 
        (PDHT_READ_ATIMER_USEC(rput)/(double)maxentries));
    eprintf("  unit: local update:   %12.7f us remote update:  %12.7f us\n",
        (PDHT_READ_ATIMER_USEC(lupdate)/(double)maxiters), 
        (PDHT_READ_ATIMER_USEC(rupdate)/(double)maxiters));
    for (int e=0; e<6; e++) {
      eprintf("  [%7d]: local get: %12.7f us remote get:     %12.7f us\n", mlistentry[e],
        (PDHT_READ_ATIMER_USEC(lgets[e])/(double)maxiters), 
        (PDHT_READ_ATIMER_USEC(rgets[e])/(double)maxiters));
    }
  } else {
    // iterations elemsize lput rput lupdate rupdate lhead rhead ltail rtail
#if 0
    eprintf("  %7d %7d %9.5f %9.5f ", maxentries, elemsize,
        (PDHT_READ_ATIMER_USEC(lput)/(double)maxentries), 
        (PDHT_READ_ATIMER_USEC(rput)/(double)maxentries));
    eprintf(" %9.5f %9.5f ", 
        (PDHT_READ_ATIMER_USEC(lupdate)/(double)maxiters), 
        (PDHT_READ_ATIMER_USEC(rupdate)/(double)maxiters));
#endif
    // { 1, 10, 100, 1000, 10000, 50000 }
    for (int e=0; e<6; e++) {
      eprintf(" %7d %9.5f %9.5f\n", mlistentry[e],
        (PDHT_READ_ATIMER_USEC(lgets[e])/(double)maxiters), 
        (PDHT_READ_ATIMER_USEC(rgets[e])/(double)maxiters));
    }
  }
done:
  pdht_free(ht);
  free(val);
}
