#include <stdlib.h>
#include <pdht.h>


extern pdht_context_t *c;



int main(int argc, char *argv[]){
  pdht_t *ht;
  printf("failure in create?\n");
  fflush(stdout);
  ht = pdht_create(sizeof(int), sizeof(int),PdhtModeStrict);
  printf("c->rank : %d no failure in create\n",c->rank);
  fflush(stdout);
  pdht_barrier();
  
  int key1 = 0;
  int value = 10;
  if(c->rank == 0){
    for(int i = 0; i < 100;i++){
      value = i;
      pdht_put(ht,&key1,&value);
      key1++;
    }
  }
  if (c->rank == 2){
    key1 = 100;
    for(int i = 100; i < 200;i++){
      value = i;
      pdht_put(ht,&key1,&value);
      key1++;
    }

  }
  /*
  if (c->rank == 0){
    pdht_print_all(ht);

  }
  pdht_barrier();
  if(c->rank == 1){
    pdht_print_all(ht);
  }
  pdht_barrier();
  if(c->rank == 2){
    pdht_print_all(ht);

  }
  */
  key1 = 0;
  
  pdht_barrier();
  
  if(c->rank == 1){
    for(int i = 0;i < 200;i++){
  
      pdht_get(ht,&key1,&value);
      key1++;
    }

  }
  pdht_barrier();
  printf("c->rank : %d  fuckups : %d \n",c->rank,ht->fuckups); 
  fflush(stdout);
  pdht_free(ht);
  return 1;





}
