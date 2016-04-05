/********************************************************/
/*                                                      */
/*  assoc.c - PDHT associative operations               */
/*                                                      */
/*  author: d. brian larkins                            */
/*  created: 3/20/16                                    */
/*                                                      */
/********************************************************/

#include <pdht.h>

/**
 * @file
 * 
 * portals distributed hash table associative operations
 */


/** 
 * pdht_acc() - associative accumulate operation into a hashed object
 *   @param key key of hash table entry to accumulate into
 *   @param ksize size of the key
 *   @param type data type of the object
 *   @param op operation to perform
 *   @param value to accumulate into entry
 */
pdht_status_t pdht_acc(pdht_t *dht, void *key, int ksize, pdht_datatype_t type, pdht_oper_t op, void *value) {
  // this is being ignored right now
  return PdhtStatusOK;
}



/** 
 * pdht_nbacc() - non-blocking associative accumulate operation into a hashed object
 *   @param key key of hash table entry to accumulate into
 *   @param ksize size of the key
 *   @param type data type of the object
 *   @param op operation to perform
 *   @param value to accumulate into entry
 */
pdht_handle_t pdht_nbacc(pdht_t *dht, void *key, int ksize, pdht_datatype_t type, pdht_oper_t op, void *value) {
  // this is being ignored right now
  return PdhtStatusOK;
}
