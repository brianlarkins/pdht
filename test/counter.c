#define _XOPEN_SOURCE 600
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>

#include <pdht.h>

#define ASIZE 10

void pdht_poll(pdht_t *dht);

extern pdht_context_t *c;

int main(int argc, char **argv);

int main(int argc, char **argv) {
  pdht_t *ht;
  int pdhtc;
  uint64_t cval;
  
  setenv("PTL_IGNORE_UMMUNOTIFY", "1",1);

  // create hash table
  ht = pdht_create(sizeof(unsigned long), ASIZE * sizeof(double), PdhtModeStrict);

/*
  if (c->size != 2) {
    if (c->rank == 0) { 
      printf("requires two (and only two) processes to run\n");
    }
    goto done;
  }
*/

  pdhtc = pdht_counter_init(ht, 0);

  printf("counter: %d\n", pdhtc);

  pdht_barrier();

  for (int i=0; i<10; i++) {
    cval = pdht_counter_inc(ht, pdhtc, 1);
    printf("rank %d got %ld\n", c->rank, cval);
    fflush(stdout);
    pdht_barrier();
  }

  pdht_barrier();

done:
  pdht_free(ht);
}
