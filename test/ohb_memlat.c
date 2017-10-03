/* Copyright (c) 2011-2016, The Ohio State University. All rights
 * reserved.
 *
 * This file is part of the OSU HiBD-Benchmarks software package
 * developed by the team members of The Ohio State University's
 * Network-Based Computing Laboratory (NBCL), headed by Professor
 * Dhabaleswar K. (DK) Panda.
 *
 * For detailed copyright and licensing information, please refer to
 * the copyright file COPYRIGHT in the top level directory.
 *
*/

#include <errno.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#ifdef HAVE_STRINGS_H
#include <strings.h>
#endif
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>

#include <pdht.h>
#include <city.h>


/*
 * program description and constants
 */
#define PROGRAM_NAME 	     "OHB Micro-benchmarks"
#define PROGRAM_DESCRIPTION  "Micro-benchmarks for Memcached"
#define VERSION		     "0.9.2"

#define MAX_SZ 		     (512*1024)
#define BENCH_PURE_SET       (1)
#define BENCH_PURE_GET       (2)
#define BENCH_MIX_SET_GET    (3)
#define BENCH_ALL	     (4)

/*
 * time function to measure latency 
 */
#define TIME() getMicrosecondTimeStamp()
static int64_t getMicrosecondTimeStamp()
{
    int64_t retval;
    struct timeval tv;
    if (gettimeofday(&tv, NULL)) {
        perror("gettimeofday");
        abort();
    }
    retval = ((int64_t)tv.tv_sec) * 1000000 + tv.tv_usec;
    return retval;
}


/* parse user-specified options */

/*
 * key/value pair to be used for micro-benchmark tests
 */
char *my_key1 = "my_key1";
char *my_key2 = "my_key2";
char my_value[MAX_SZ];

void ahash(pdht_t *dht, void *key, ptl_match_bits_t *mbits, uint32_t *ptindex, ptl_process_t *rank){
    *mbits = CityHash64((char *)key, dht->keysize);
    *ptindex = *mbits % dht->ptl.nptes;
    (*rank).rank = 1;
}





/* Main fuction for OHB Micro-benchmarks */
int main(int argc, char *argv[])
{
    int benchmark_type = BENCH_PURE_GET;
   
    int return_code= 0, i, sz, iters=10000, opt;
    int64_t begin=0, end=0;
    sz = 512; // default value size is 512    


    while((opt = getopt(argc,argv,"s:")) != -1){

        switch (opt)
        {
            case 's':
                sz = atoi(optarg);
               
        }
            

    }
    if (sz > MAX_SZ){
        goto done;
    }

    pdht_t *ht;


    
    
    int count = 0;


    if (benchmark_type ==BENCH_PURE_GET){
        size_t val_len;
        uint32_t flags;
        pdht_status_t rc;


        memset(my_value, 'v', MAX_SZ);
            

        pdht_config_t cfg;
        cfg.pendmode = PdhtPendingTrig;
        cfg.nptes = 1;
        cfg.maxentries = 25000;
        cfg.pendq_size = 10000;
        cfg.ptalloc_opts = PTL_PT_MATCH_UNORDERED;
        
        pdht_tune(PDHT_TUNE_ALL,&cfg);
            
        pdht_t *ht;
        ht = pdht_create(sizeof(my_key1),sz,PdhtModeStrict);
        pdht_sethash(ht,ahash);
        
        
        //printf("c->rank : %d calling fence\n",c->rank);
        pdht_fence(ht);
        //printf("c->rank out of the fence : %d\n",c->rank);

        if (c->rank == 0){
            rc= pdht_put(ht, my_key1, my_value);
        }

        pdht_fence(ht);    
        begin = TIME();
        if (c->rank == 0){
            for(i=0; i<iters; i++){
                pdht_get(ht, my_key1, my_value);
            }
        }
        end = TIME();
        if((rc == PdhtStatusOK) && (c->rank == 0)) {
         
            fprintf(stderr, "%10i byte:\t%11.2f usecs\n", sz, (1.0f * (end-begin)) / ((iters)));
        } 

        pdht_print_stats(ht);
        pdht_free(ht);
    }
done:
    return return_code;
}

