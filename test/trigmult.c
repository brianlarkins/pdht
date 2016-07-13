#include <sys/time.h>
#include <sys/resource.h>

#include <pdht.h>

#define ASIZE 10
#define ELEMS 12

void pdht_poll(pdht_t *dht);

extern pdht_context_t *c;

int main(int argc, char **argv);

int main(int argc, char **argv) {
  pdht_t *ht;
  unsigned long key = 10;
  double pbuf[ELEMS][ASIZE], gbuf[ELEMS][ASIZE];

  // create hash table
  ht = pdht_create(sizeof(unsigned long), ASIZE * sizeof(double), PdhtModeStrict);

  if (c->size != 2) {
    if (c->rank == 0) { 
      printf("requires two (and only two) processes to run\n");
    }
    goto done;
  }


  if (c->rank == 0) {

    for (int e=0; e<ELEMS; e++) {
      for (int i=0; i<ASIZE; i++) {
        pbuf[e][i] = (e*100.0) + i * 1.1;
      }
      key = 10 + e;

      printf("%d: putting object %d\n", c->rank, e);
      pdht_put(ht, &key, pbuf[e]);
    }
    pdht_barrier();


    for (int e=0; e<ELEMS; e++) {
      printf("%d: getting object %d\n", c->rank, e);

      key = 10 + e;
      //printf("%d: gbuf: %p &gbuf: %p\n", c->rank, gbuf, &gbuf);
      pdht_get(ht, &key, gbuf[e]);

      for (int i=0; i<ASIZE; i++) {
        printf("gbuf[%d][%d] = %5.1f\n", e, i, gbuf[e][i]);
      }
    }

  } else {
    print_count(ht, "before barrier:");
    pdht_barrier();
    print_count(ht, "after barrier:");
  }

  pdht_barrier();

done:
  pdht_free(ht);
}
