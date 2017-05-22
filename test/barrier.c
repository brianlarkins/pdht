#include <sys/time.h>
#include <sys/resource.h>

#include <pdht.h>

#define NELEMS 10

extern pdht_context_t *c;

int main(int argc, char **argv);

int main(int argc, char **argv) {
  pdht_t *ht;
  int val[NELEMS];
  int rval[NELEMS];
  //printf("pid: %d\n", getpid());


  ht = pdht_create(sizeof(int), sizeof(int), PdhtModeStrict);

#if 0
  printf("%d: before barrier\n", c->rank);
  fflush(stdout);
  pdht_barrier();
  printf("%d: after barrier\n", c->rank);
  printf("%d: before barrier\n", c->rank);
  fflush(stdout);
  pdht_barrier();
  printf("%d: after barrier\n", c->rank);
#endif


  for (int i=0; i < c->size; i++) {
    if (i == c->rank) {
      printf("hello from %d of %d!\n", c->rank, c->size);
      fflush(stdout);
    } 
    pdht_barrier();
  }


  for (int i=0; i < NELEMS; i++) {
    val[i] = 1+i;
  }

  for (int j=0; j < c->size; j++) {
    if (j == c->rank) {
      printf("rank %d values: ", c->rank);
      for (int i=0; i < NELEMS; i++) {
        printf("%d ", val[i]);
      }
      printf("\n");
      fflush(stdout);
    }
    pdht_barrier();
  }

  pdht_barrier();

  pdht_reduce(val, rval, PdhtReduceOpSum, IntType, NELEMS); 

  if (c->rank == 0) {
    printf("reduce values: ");
    for (int i=0; i < NELEMS; i++) {
      printf("%d ", rval[i]);
    }
    printf("\n");
  }

  pdht_barrier();
  
  if (c->rank == 0) {
    for (int i=0; i < NELEMS; i++) {
       val[i] = NELEMS * NELEMS * i;
    }
  } else {
    for (int i=0; i < NELEMS; i++) {
       val[i] = 0;
    }
  }
  pdht_broadcast(val, IntType, NELEMS);

  for (int j=0; j < c->size; j++) {
    if (j == c->rank) {
      printf("rank %d broadast: ", c->rank);
      for (int i=0; i < NELEMS; i++) {
        printf("%d ", val[i]);
      }
      printf("\n");
      fflush(stdout);
    }
    pdht_barrier();
  }
 
  pdht_barrier();

  pdht_allreduce(val, rval, PdhtReduceOpSum, IntType, NELEMS);

  for (int j=0; j < c->size; j++) {
    if (j == c->rank) {
      printf("rank %d allreduce: ", c->rank);
      for (int i=0; i < NELEMS; i++) {
        printf("%d ", rval[i]);
      }
      printf("\n");
      fflush(stdout);
    }
    pdht_barrier();
  }

  pdht_barrier();
done:
  pdht_free(ht);
}
