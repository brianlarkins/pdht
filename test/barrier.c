#include <sys/time.h>
#include <sys/resource.h>

#include <pdht.h>


extern pdht_context_t *c;

int main(int argc, char **argv);

int main(int argc, char **argv) {
  pdht_t *ht;
  
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
#if 1


  for (int i=0; i < c->size; i++) {
    if (i == c->rank) {
      printf("hello from %d of %d!\n", c->rank, c->size);
      fflush(stdout);
    } 
    pdht_barrier();
  }
#endif
  
  pdht_free(ht);
}
