/**
 * this is a re-write of the OpenSHMEM hash table benchmark.
 * the original was both awful and unreadable
 * this is only awful
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#include <shmem.h>

typedef unsigned long int numb;   
/* Should be 64 bit wide, to hold the square of: */
/* If you change this, also change "atol" in main */

#define ITERCOUNT 2

/* global vars */

int size, rank;
int hashlen, collisions, localhashlen, total_collisions;
int *hashtab, *hashcount;
long *pe_lock;

FILE *outfile;

int  p_wrk[_SHMEM_REDUCE_MIN_WRKDATA_SIZE];
long p_sync[_SHMEM_REDUCE_SYNC_SIZE];


numb N=NHASH;     /* Number of objects to hash */
size_t m=8;  /* Size of objects in bytes, rounded up to be a multiple of
                sizeof(numb) */

/*numb k;*/     /* Number of times each object is expected */



/* Some fast way to create funny data: */
#define modulus 1073741827
#define multipl 33554467

numb val = 1234567;



// initially scramble the object value based on my rank
void resetvalue()
{
  val = 1234567;

  for(int j=0;j < (rank*NHASH);j++){
    val = (val * multipl) % modulus;
  }
}

// obj == value
void fnew(numb *key, numb *obj) {
  numb hashlen = 2 * NHASH * size + 1;

  // val is a global "randomish" thing
  val = ((val * multipl) % modulus); // re-scramble

  *obj = val; 
  *key = *obj % hashlen + 1; // the key is the same as the value! brilliant.
}



double get_wtime()
{
  double t;
  struct timeval tv;

  gettimeofday(&tv, NULL);

  t = (tv.tv_sec*1000000LL + tv.tv_usec)/1000000.0;

  return t;
}


void inithash()
{
  size = shmem_n_pes();
  rank = shmem_my_pe();

  hashlen = 2*NHASH*size + 1;
  localhashlen = 2*NHASH + 1;
  collisions = 0;
  total_collisions = 0;

  hashtab   = shmem_malloc(localhashlen * sizeof *hashtab);
  hashcount = shmem_malloc(localhashlen * sizeof *hashcount);

  // one lock per process
  pe_lock   = shmem_malloc(size * sizeof(*pe_lock));
  memset(pe_lock, 0, size * sizeof(*pe_lock));
}



void finhash()
{
  shmem_free(hashtab);
  shmem_free(hashcount);
  shmem_free(pe_lock);
}



void hashlookup_with_amo(long k, long v)
{
  int dest_pe, dest_pos, local_count, local_hash;

  while (1) {
    dest_pe = k / localhashlen;

    if (dest_pe*localhashlen < k) {
      dest_pe += 1;
    }
    dest_pe -= 1;

    dest_pos = k - (dest_pe)*localhashlen;

    local_hash = shmem_int_cswap(&hashtab[dest_pos-1], 0, v, dest_pe);

    if (local_hash == 0 || local_hash == v) {
      /* update successful, so increment count and return */
      shmem_int_inc(&hashcount[dest_pos-1], dest_pe);

      return;
    } else {
      /* its a collision */
      collisions += 1;
      k += 1;
      if (k > hashlen) {
        k = 1;
      }
    }
  }
}



void hashlookup(long k, long v)
{
  int dest_pe, dest_pos, local_count, local_hash;

  while (1) {
    dest_pe = k / localhashlen;

    if (dest_pe*localhashlen < k) {
      dest_pe += 1;
    }
    dest_pe -= 1;

    dest_pos = k - (dest_pe)*localhashlen;

    /* lock the data */
    shmem_set_lock( &pe_lock[dest_pe] );

    shmem_int_get(&local_hash, &hashtab[dest_pos-1], 1, dest_pe);

    if (local_hash == 0) {
      /* insert the entry */
      int one = 1;
      shmem_int_put(&hashtab[dest_pos-1], (int *)&v, 1, dest_pe);
      shmem_int_put(&hashcount[dest_pos-1], &one, 1, dest_pe);

      /* unlock before return */
      shmem_clear_lock( &pe_lock[dest_pe] );
#ifdef _DEBUG
      fprintf(outfile, "%d writing %d to pe %d at position %d\n", rank, v, dest_pe, dest_pos);
#endif
      return;
    } else {
      /* check to see if it is a collision */
      if (local_hash == v) {
        /* its a repetition */
        shmem_int_inc(&hashcount[dest_pos-1], dest_pe);
        shmem_clear_lock( &pe_lock[dest_pe] );
#ifdef _DEBUG
        fprintf(outfile, "%d writing %d to pe %d at position %d\n", rank, v, dest_pe, dest_pos);
#endif
        return;
      } else {
        /* its a collision */
#ifdef _DEBUG
        fprintf(outfile, "%d writing %d to pe %d at position %d ***\n", rank, v, dest_pe, dest_pos);
        printf("%d:%d collision!\n", rank);
#endif
        collisions += 1;
        k += 1;
        if (k > hashlen) {
          k = 1;
        }
        shmem_clear_lock( &pe_lock[dest_pe] );
      }
    }

  }
}


int main()
{
  numb pb,key,value;
  int M;
  double t1, t2;
  int rank_id;

#ifdef _DEBUG
  char filename[20];
#endif

  shmem_init();

  for (int i = 0; i < _SHMEM_REDUCE_SYNC_SIZE; i++) {
    p_sync[i] = _SHMEM_SYNC_VALUE;
  }

  // initialize hash table
  inithash();

  shmem_barrier_all();

  if (rank == 0) {
     printf("starting OpenSHMEM hash table benchmark: NHASH is %d\n", NHASH);
  }

#ifdef _DEBUG
  sprintf(filename, "output.%d", rank);
  outfile = fopen(filename, "w+");
#endif

  t1 = get_wtime();

  M = NHASH / 100;
  if (M == 0) M = 1;

  rank_id = rank + 1;
  for (int iter=0; iter < ITERCOUNT;iter++) {
    resetvalue(size, rank_id);

#ifdef PROGRESS_BAR
    if (rank == 0) {
      printf("Pass %d: ", iter);
    }
#endif

    for (int i = 1; i <= NHASH; i++) {
      fnew(&key, &value);
      //printf("k: %lu v: %lu\n", key, value);
#ifndef USE_AMO
      hashlookup(key,value);
#else
      hashlookup_with_amo(key,value);
#endif
#ifdef PROGRESS_BAR
      if (rank == 0 && i%M == 1) {
        printf(".");
        fflush(stdout);
      }
#endif
     // printf("%d : %d\n", rank, i);
      //fflush(stdout);

    }

#ifdef PROGRESS_BAR
    if (rank == 0) {
      printf(" DONE!\n");
    }
#endif
    shmem_barrier_all();
  }

  shmem_barrier_all();

  t2 = get_wtime();

  if (rank == 0 && NHASH < 1000) {
    for (int j = 0; j < size; j++) {
      for (int i = 0; i < localhashlen; i++) {
        int count;
        shmem_int_get(&count, &hashcount[i], 1, j);
        if (count > 0) {
          int tab;
          shmem_int_get(&tab, &hashtab[i], 1, j);
          printf("%d %d %d %d\n", j+1, i+1, tab, count);
        }
      }
    }
  }

  shmem_int_sum_to_all(&total_collisions, &collisions, 1, 0, 0, size,
      p_wrk, p_sync);

  if (rank == 0) {
    printf("Avg # collisions: %12.2f\n", total_collisions/(1.0*size));
    printf("Total time is %8.3f sec\n", (t2-t1) );
  }

#ifdef _DEBUG
  fclose(outfile);
#endif
  finhash();

  shmem_finalize();

  return 0;
}
