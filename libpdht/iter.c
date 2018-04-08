/********************************************************/
/*                                                      */
/*  iter.c - PDHT hash function local iterators         */
/*                                                      */
/*  author: d. brian larkins                            */
/*  created: 2/13/18                                    */
/*                                                      */
/********************************************************/

#include <pdht_impl.h>

/**
 * @file
 * 
 * portals distributed hash function local iterator
 */

/** 
 * pdht_iterate() - construct an iterator over a pdht hash table
 *  @param dht hash table structure
 *  @param it a PDHT iterator structure
 *  @returns status of creation operation
 */
pdht_status_t pdht_iterate(pdht_t *dht, pdht_iter_t *it) {
  it->dht = dht;
  it->iterator = (char *)dht->ht;
  return PdhtStatusOK;
}



/** 
 * pdht_hasnext - checks to see if the next HT entry is legal to iterate over
 * @param it a PDHT iterator structure
 * @returns 1 if next entry is valid, 0 otherwise
 */
int pdht_hasnext(pdht_iter_t *it) {
  int ret = 0;
  _pdht_ht_entry_t     *phte;
  _pdht_ht_trigentry_t *thte;
  char *peek = it->iterator;

  switch (it->dht->pmode) {
    case PdhtPendingPoll:
      phte = (_pdht_ht_entry_t *)peek;
      if (phte->ame != PTL_INVALID_HANDLE)
        ret = 1;
      break;
    case PdhtPendingTrig:
      thte = (_pdht_ht_trigentry_t *)peek;
      if (thte->ame != PTL_INVALID_HANDLE)
        ret = 1;
      break;
  }
  //pdht_dprintf("returning %d: %p %p %d\n", ret, thte, &thte->ame, (int)thte->ame);
  return ret;
}



/**
 * pdht_getnext - get next entry of iteration over local part of an HT
 * @param it a PDHT iterator structure
 * @param key optional copy-out of matching key for entry
 * @returns pointer to HT entry
 */
void *pdht_getnext(pdht_iter_t *it, void **key)  {
  void *ret = NULL;
  _pdht_ht_entry_t     *phte;
  _pdht_ht_trigentry_t *thte;

  switch (it->dht->pmode) {
    case PdhtPendingPoll:
      phte = (_pdht_ht_entry_t *)it->iterator;
      if (phte->ame != PTL_INVALID_HANDLE) {
        ret = &phte->data;
        if (key)
          *key = &phte->key;
      }
      break;
    case PdhtPendingTrig:
      thte = (_pdht_ht_trigentry_t *)it->iterator;
      if (thte->ame != PTL_INVALID_HANDLE)
        ret =  &thte->data;
        if (key)
          *key = &thte->key;
      break;
  }
  it->iterator += it->dht->entrysize;
  return ret;
}
