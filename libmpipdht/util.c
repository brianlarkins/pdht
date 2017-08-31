#include <pdht.h>

extern pdht_context_t *c;






int eprintf(const char *format, ...) {
  
  va_list ap;
  int ret;

  if (c->rank == 0) {
    va_start(ap, format);
    ret = vfprintf(stdout, format, ap);
    va_end(ap);
    fflush(stdout);
    return ret;
  }
  else
    return 0;
  
}


void pdht_print_all(pdht_t *dht){

  ht_t *s;
  char buf[PDHT_MAXKEYSIZE + dht->elemsize];
  int key;
  int value;
  for(s = dht->ht;s != NULL;s = (ht_t*)(s->hh.next)){
    memcpy(buf,s->value,PDHT_MAXKEYSIZE + dht->elemsize);
    printf("c->rank : %d key : %d value : %d \n",c->rank, *(int*)buf,*(int*)(buf+PDHT_MAXKEYSIZE));
  }



}


void pdht_print_stats(pdht_t *dht){
  printf("I dont do anything \n");


}
