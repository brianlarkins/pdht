#include <sys/time.h>
#include <sys/resource.h>

#include <pdht.h>

#define ASIZE 10

void pdht_poll(pdht_t *dht);

extern pdht_context_t *c;

int main(int argc, char **argv);

int main(int argc, char **argv) {
  pdht_t *ht;
  pdht_status_t ret;
  unsigned long key = 10;
  double pbuf[ASIZE], gbuf[ASIZE];
  
  // create hash table
  ht = pdht_create(sizeof(unsigned long), ASIZE * sizeof(double), PdhtModeStrict);

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

    printf("%d: getting object\n", c->rank);
    printf("%d: gbuf: %p &gbuf: %p\n", c->rank, gbuf, &gbuf);
    ret = pdht_get(ht, &key, gbuf);
    if (ret == PdhtNotFound) {
      printf("hooray!\n");
    } else {
      printf("boo!\n");
    }

#if 0
    for (int i=0; i<ASIZE; i++) {
       printf("gbuf[%d] = %5.1f\n", i, gbuf[i]);
    }
#endif

  } else {

    printf("%d: calling poll()\n", c->rank); 
    pdht_poll(ht);
    pdht_barrier();

  }

  pdht_barrier();

done:
  pdht_free(ht);
}
