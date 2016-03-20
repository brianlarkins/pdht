/********************************************************/
/*                                                      */
/*  putget.c - PDHT put / get operations                */
/*                                                      */
/*  author: d. brian larkins                            */
/*  created: 3/20/16                                    */
/*                                                      */
/********************************************************/


/**
 * @file
 * 
 * portals distributed hash table put/get ops
 */

/**
 * pdht_put - puts or overwrites an entry in the global hash table
 *   @param key - hash table key
 *   @param ksize - size of key
 *   @param value - value for table entry
 *   @returns status of operation
 */
pdht_status_t        pdht_put(void *key, int ksize, void *value) {
}



/**
 * pdht_get - gets an entry from the global hash table
 *   @param key - hash table key
 *   @param ksize - size of key
 *   @param value - value for table entry
 *   @returns status of operation
 */
pdht_status_t        pdht_get(void *key, int ksize, void **value) {
}
