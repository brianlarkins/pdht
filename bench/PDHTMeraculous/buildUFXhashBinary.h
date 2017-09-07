#ifndef BUILD_UFX_HASH_H
#define BUILD_UFX_HASH_H

#include <stdio.h>
#include <sys/stat.h>
#include <assert.h>

#include "meraculous.h"
#include "packingDNAseq.h"
#include "kmer_hash.h"

hash_table_t* buildUFXhash(int64_t size, FILE *fd, memory_heap_t *memory_heap_res, int64_t myShare, int64_t dsize, int64_t dmin, int CHUNK_SIZE, int load_factor)
{
  //pdht_t *ht;
  hash_table_t *dist_hashtable;
  //   memory_heap_t memory_heap;
  //   shared list_t *lookup_res;
  int64_t chars_read, cur_chars_read , ptr, retval, i;
  int64_t buffer_size, offset;
  int64_t num_chars, chars_to_be_read;
  //   shared[1] int64_t *heap_sizes;
  int64_t *my_heap_sizes, my_heap_size, my_ufx_lines, idx;
  int *ufx_remote_thread, remote_thread;
  int64_t  hashval;


  //UPC_TICK_T start_read, end_read, start_storing, end_storing, start_setup, end_setup, start_calculation, end_calculation;
  pdht_timer_t read_timer, storing_timer, setup_timer, calculation_timer;

  //heap_sizes = (shared[1] int64_t*) upc_all_alloc(THREADS, sizeof(int64_t));
  //heap_sizes[MYTHREAD] = 0;


  //I need to find a way to figure out what lines in the data i am putting? How ??? need to know how many processes I have. Reduce? probably

  my_ufx_lines = (size / THREADS + size % THREADS);
  if (MYTHREAD == THREADS-1) {
    my_ufx_lines = (size / THREADS + size % THREADS);
  } else {
    my_ufx_lines = size / THREADS;
  }
  //this can probably get thrown away but idk
  //ufx_remote_thread = (int*) malloc(my_ufx_lines * sizeof(int));

  if (VERBOSE > 1) LOG("Thread %d: Preparing to read %ld UFX lines (%ld chars)\n", MYTHREAD, my_ufx_lines, chars_to_be_read);

#ifdef DEBUG
  int64_t kmer_count = 0;
#endif

  /* Initialize lookup-table --> necessary for packing routine */
  init_LookupTable();
  char **kmersarr;
  int *counts;
  char *lefts;
  char *rights;


  PDHT_START_ATIMER(read_timer);
  
  
  int64_t kmers_read = UFXRead(fd, dsize, &kmersarr, &counts, &lefts, &rights, my_ufx_lines, dmin, 0, MYTHREAD);

  if (MYTHREAD == 0){
    printf("Threads done with I/O\n");
  }
  PDHT_STOP_ATIMER(read_timer);



  PDHT_START_ATIMER(setup_timer);
  /* Create and initialize hashtable */
  dist_hashtable = pdht_create(KMER_LENGTH,sizeof(list_t),PdhtModeStrict);

  pdht_fence(dist_hashtable);

  PDHT_STOP_ATIMER(setup_timer);

  if (MYTHREAD == 0){
    printf("Threads done with setting-up\n");
  }

  ptr = 0;
  list_t new_entry;
  PDHT_START_ATIMER(store_timer);
  char exts[2];
  int my_offset;
  int num_per_process;
  num_per_process = kmers_read / c->size
  myoffset = num_per * c->rank;
  if (c->rank == (c->size - 1)){
    num_per_process += kmers_read % c->size;
  }


  while (ptr < num_per_process) {
    packSequence((unsigned char*) (kmersarr[my_offset + ptr]), new_entry.packed_key, KMER_LENGTH);

#ifdef MERACULOUS
    exts[0] = lefts[my_offset + ptr];
    exts[1] = rights[my_offset + ptr];
    new_entry.packed_extensions = convertExtensionsToPackedCode((unsigned char*) exts);
#endif

    new_entry.next = NULL;
    pdht_put(dist_hashtable,kmersarr[ptr],&new_entry);

    ptr++;
  }
  assert(idx == my_ufx_lines);



  pdht_fence(dist_hashtable);
  pdht_fence(dist_hashtable);
  
  PDHT_STOP_ATIMER(store_timer);
  #ifdef PROFILE 
  pdht_fence(dist_hashtable);
  if (c->rank == 0) {
    printf("\n************* SET - UP TIME *****************");
    printf("\nTime spent on setting up the distributed hash table is %f seconds\n", PDHT_READ_ATIMER(setup_timer));

  #ifdef DETAILED_IO_PROFILING
    printf("\n\n************* DETAILED TIMINGS FOR I/O AND STORING KMERS *****************\n");
  #endif
  }
  
  pdht_fence(dist_hashtable);
  //upc_barrier;

  #ifdef DETAILED_IO_PROFILING

  printf("Thread %d spent %f seconds on read I/O and %f seconds on storing %ld kmers\n", c->rank, PDHT_READ_ATIMER(read_timer), PDHT_READ_ATIMER(store_timer));

  #endif

  pdht_fence(dist_hashtable);
  //upc_barrier;

  DeAllocateAll(&kmersarr, &counts, &lefts, &rights, kmers_read);

  #ifdef BUCKET_BALANCE_PROFILING
  if (c->rank == 0) {
    printf("\n--------------------------------\n");
    printf("--------- Bucket balance -------\n");
    printf("--------------------------------\n");
  }
  if (c->rank == 0) {
    for (i=0; i<THREADS; i++)
      printf("Thread %ld has %ld elements in its buckets\n", i, memory_heap.heap_indices[i]);    
  }
  #endif

  #endif

  pdht_fence(dist_hashtable);
  return dist_hashtable;
}

#endif // BUILD_UFX_HASH_H

