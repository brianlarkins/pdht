#define _XOPEN_SOURCE 600
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>

#include <pdht.h>

#define offsetof(type, member)  __builtin_offsetof (type, member)

#define ASIZE 100

void pdht_poll(pdht_t *dht);

extern pdht_context_t *c;

int main(int argc, char **argv);

struct data_s {
  double dval;
  int64_t used;
};
typedef struct data_s data_t;

void mkprinter(void *key) {
  unsigned long k = *(unsigned long *)key;
  printf("%lu ", k);
}

void mvprinter(void *val) {
  data_t *dp = (data_t *)val;
  printf("dval: %f used: %ld", dp->dval, dp->used);
}

int main(int argc, char **argv) {
  pdht_t *ht;
  pdht_config_t cfg;
  size_t off;
  int64_t old;

  cfg.nptes = 1;
  cfg.pendmode = PdhtPendingTrig;
  cfg.maxentries = 50000;
  cfg.pendq_size = 5000;
  cfg.ptalloc_opts = PTL_PT_MATCH_UNORDERED;
  cfg.local_gets = PdhtSearchLocal;
  cfg.rank = PDHT_DEFAULT_RANK_HINT;


  //  setenv("PTL_IGNORE_UMMUNOTIFY", "1",1);

  // create hash table
  pdht_tune(PDHT_TUNE_ALL, &cfg);
  ht = pdht_create(sizeof(unsigned long), ASIZE * sizeof(data_t), PdhtModeStrict);

  data_t val;

  if (c->rank == 0) {
    for (unsigned long key=0; key < ASIZE; key++) {
      val.dval = key + .1;
      val.used = 0;
      pdht_put(ht, &key, &val);
    }
  }

  pdht_fence(ht);
  pdht_barrier();
  pdht_print_active(ht, mkprinter, mvprinter);

  off = offsetof(data_t, used);
  old = 0;
  for (unsigned long key=0; key < ASIZE; key++) {
    pdht_atomic_cswap(ht, &key, off, &old, 1);
    if (old != 1) {
        printf("%d: swapped %lu from %ld to 1\n", c->rank, key, old);
    } else {
        printf("%d: somebody already swapped %ld\n", c->rank, key);
  }
}

  pdht_fence(ht);
  pdht_barrier();
  pdht_print_active(ht, mkprinter, mvprinter);

  pdht_free(ht);
}
