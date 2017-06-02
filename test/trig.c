#include <sys/time.h>
#include <sys/resource.h>

#include <pdht.h>

#define ASIZE 10

void pdht_poll(pdht_t *dht);

void remotehash(pdht_t *dht, void *key, ptl_match_bits_t *mbits, uint32_t *ptindex, ptl_process_t *rank) {
  (*rank).rank = 1;
  *mbits = *(unsigned long *)key;
  *ptindex = *(unsigned long *)key % dht->ptl.nptes;
  *ptindex = 0;
}

extern pdht_context_t *c;

int main(int argc, char **argv);

int main(int argc, char **argv) {
  pdht_t *ht;
  unsigned long key = 10;
  double pbuf[ASIZE], gbuf[ASIZE];

  pdht_config_t cfg;
  cfg.pendmode = PdhtPendingTrig;
  pdht_tune(PDHT_TUNE_PMODE, &cfg);

  // create hash table
  ht = pdht_create(sizeof(unsigned long), ASIZE * sizeof(double), PdhtModeStrict);

  pdht_sethash(ht, remotehash);

  if (c->size != 2) {
    if (c->rank == 0) { 
      printf("requires two (and only two) processes to run\n");
    }
    goto done;
  }


  if (c->rank == 0) {
    for (int i=0; i<ASIZE; i++) {
      pbuf[i] = i * 1.1;
    }

    printf("%d: putting object\n", c->rank);
    pdht_put(ht, &key, pbuf);
    pdht_barrier();

    printf("%d: getting object\n", c->rank);
    //printf("%d: gbuf: %p &gbuf: %p\n", c->rank, gbuf, &gbuf);
    pdht_get(ht, &key, gbuf);

    for (int i=0; i<ASIZE; i++) {
      printf("gbuf[%d] = %5.1f\n", i, gbuf[i]);
    }

  } else {
    //print_count(ht, "before barrier:");
    pdht_barrier();
    //print_count(ht, "after barrier:");
  }

  pdht_fence(ht);
  pdht_barrier();

done:
  pdht_free(ht);
}
