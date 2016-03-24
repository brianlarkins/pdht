/********************************************************/
/*                                                      */
/*  pdht_inline.h - PDHT inline functions               */
/*                                                      */
/*  author: d. brian larkins                            */
/*  created: 3/20/16                                    */
/*                                                      */
/********************************************************/


/**
 * @file
 * 
 * portals distributed hash table inline functions
 */

#pragma once


/* blocking shortcuts */

/**
 * pdht_puti - put an integer keyed value into the DHT
 *  @param k integer key
 *  @param v DHT entry to put
 *  @returns status of the operation
 */
static inline pdht_status_t pdht_puti(pdht_t *dht,int k, void *v) {
   return pdht_put(dht,k,sizeof(int), v);
}



/**
 * pdht_geti - get an integer keyed value from the DHT
 *  @param k integer key
 *  @param v DHT entry to get
 *  @returns status of the operation
 */
static inline pdht_status_t pdht_geti(pdht_t *dht,int k, void **v) {
   return pdht_get(dht,k,sizeof(int), v);
}



/**
 * pdht_putf - put an double keyed value into the DHT
 *  @param k double key
 *  @param v DHT entry to put
 *  @returns status of the operation
 */
static inline pdht_status_t pdht_putf(pdht_t *dht,double k, void *v) {
   return pdht_put(dht,k,sizeof(double), v);
}



/**
 * pdht_getf - get an double/float keyed value from the DHT
 *  @param k double key
 *  @param v DHT entry to get
 *  @returns status of the operation
 */
static inline pdht_status_t pdht_getf(pdht_t *dht,double k, void **v) {
   return pdht_get(dht,k,sizeof(double), v);
}


/* non-blocking shortcuts */

/**
 * pdht_nbputi - put an integer keyed value into the DHT
 *  @param k integer key
 *  @param v DHT entry to put
 *  @returns status of the operation
 */
static inline pdht_status_t pdht_nbputi(pdht_t *dht,int k, void *v) {
   return pdht_put(dht,k,sizeof(int), v);
}



/**
 * pdht_nbgeti - get an integer keyed value from the DHT
 *  @param k integer key
 *  @param v DHT entry to get
 *  @returns status of the operation
 */
static inline pdht_status_t pdht_nbgeti(pdht_t *dht,int k, void **v) {
   return pdht_get(dht,k,sizeof(int), v);
}



/**
 * pdht_nbputf - put an double keyed value into the DHT
 *  @param k double key
 *  @param v DHT entry to put
 *  @returns status of the operation
 */
static inline pdht_status_t pdht_putf(pdht_t *dht,double k, void *v) {
   return pdht_nbput(dht,k,sizeof(double), v);
}



/**
 * pdht_nbgetf - get an double/float keyed value from the DHT
 *  @param k double key
 *  @param v DHT entry to get
 *  @returns status of the operation
 */
static inline pdht_status_t pdht_nbgetf(pdht_t *dht,double k, void **v) {
   return pdht_nbget(dht,k,sizeof(double), v);
}
