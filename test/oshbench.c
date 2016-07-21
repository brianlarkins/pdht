#include <sys/time.h>
#include <sys/resource.h>
#include <pdht.h>

#define NHASH 10000
#define modulus 1073741827
#define multipl 33554467

//#define ASIZE 10

typedef unsigned long int numb;

extern pdht_context_t *c;

numb val;
int hashlen, collisions, localhashlen, total_collisions;

int eprintf(const char *format, ...);


void resetvalue() {
  val = 1234567;

  for(int j=0;j < (c->rank*NHASH);j++){
    val = (val * multipl) % modulus;
  }
}


// obj == value
void fnew(numb *key, numb *obj) {
  numb hashlen = 2 * NHASH * size + 1;

  // val is a global "randomish" thing
  val = ((val * multipl) % modulus); // re-scramble

  *obj = val; 
  *key = *obj % hashlen + 1; 
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
  hashlen = 2*NHASH*size + 1;
  localhashlen = 2*NHASH + 1;
  collisions = 0;
  total_collisions = 0;
}


void f_hash(pdht_t *dht, void *key, ptl_match_bits_t *mbits, ptl_process_t *rank) {
  int dest_pe;
  numb k;
  *mbits = *(ptl_match_bits_t *)key; // treat long key same as match bits

  k = *(numb *)key; // deal with as a long

  dest_pe = k / localhashlen;
  if (dest_pe * localhashlen < k)
    dest_pe += 1;
  dest_pe -= 1;
  (*rank).rank = dest_pe;
}

void hashlookup(pdht_t *ht, long k, long v) {
  int ret;
  numb vv;
  
  while (1) {

    ret = pdht_get(ht, &k, &vv);
    if (ret == PdhtStatusNotFound) {
     // insert the entry
     pdht_put(ht, &k, &v);
     return;
     
    } else {
      // check for repetition or collision
      if (vv == v) {
	// repetition -- do nothing
	;
#ifdef _DEBUG
	fprintf(outfile, "%d writing %d to pe %d at position %d\n", rank, v, dest_pe, dest_pos);
#endif

      } else {
	// collision
#ifdef _DEBUG
        fprintf(outfile, "%d writing %d to pe %d at position %d ***\n", rank, v, dest_pe, dest_pos);
        printf("%d:%d collision!\n", rank);
#endif
	collisions++;
	k = k > hashlen ? 1 : k + 1; // linear probe
      }
    }
  }
}



int main(int argc, char **argv) {
  pdht_t *ht;
  int M;
  numb key, value;
  unsigned long key = 10;
  double t1, t2;

  //double pbuf[ASIZE], gbuf[ASIZE];
  
  // create hash table
  ht = pdht_create(sizeof(unsigned long), NHASH * sizeof(double), PdhtModeStrict);

  inithash();

  if (c->size < 2) {
    eprintf("requires at least two processes to run\n");
    goto done;
  }

  pdht_sethash(ht, f_hash);

  pdht_barrier();

  eprintf("starting PDHT hash table benchmark: NHASH is %d\n", NHASH);


#ifdef _DEBUG
  sprintf(filename, "output.%d", rank);
  outfile = fopen(filename, "w+");
#endif

  t1 = get_wtime();

  M = NHASH / 100;
  if (M == 0) M = 1;

  rank_id = rank + 1; // seriously
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
      hashlookup(key,value);


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
    pdht_barrier();
  }

  pdht_barrier();

  t2 = get_wtime();

  shmem_int_sum_to_all(&total_collisions, &collisions, 1, 0, 0, size,
      p_wrk, p_sync);

  if (rank == 0) {
    printf("Avg # collisions: %12.2f\n", total_collisions/(1.0*size));
    printf("Total time is %8.3f sec\n", (t2-t1) );
  }

#ifdef _DEBUG
  fclose(outfile);
#endif

done:
  pdht_free(ht);
  return 0;
}
