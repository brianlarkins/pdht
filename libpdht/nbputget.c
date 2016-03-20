/********************************************************/
/*                                                      */
/*  nbputget.c - PDHT asynch non-blocking operations    */
/*                                                      */
/*  author: d. brian larkins                            */
/*  created: 3/20/16                                    */
/*                                                      */
/********************************************************/


/**
 * @file
 * 
 * portals distributed hash table non-blocking ops
 */


/**
 * pdht_nbput - asynchronously puts or overwrites an entry in the global hash table
 *   @param key - hash table key
 *   @param ksize - size of key
 *   @param value - value for table entry
 *   @returns handle for completion operations
 */
pdht_handle_t        pdht_nbput(void *key, int ksize, void *value) {
}



/**
 * pdht_nbget - asynchronously gets an entry from the global hash table
 *   @param key - hash table key
 *   @param ksize - size of key
 *   @param value - value of the global entry
 *   @returns handle for completion operations
 */
pdht_handle_t        pdht_nbget(void *key, int ksize, void **value) {
}
