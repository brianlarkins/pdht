#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <infiniband/verbs.h>



int main(void) {
  struct ibv_device **ibd;
  struct ibv_context *ctx;
  struct ibv_pd *pd;
  struct ibv_mr *mr;
  char *buf;
  u_int64_t size;
  u_int64_t inc = 1024 * 1204;
  inc *= 1024; // 1G each check
  int res = 1024 * 1024; // resolution of max (1MB)

  ibd = ibv_get_device_list(NULL);
  if (!ibd)  {
    perror("get devices");
    exit(1);
  }

  ctx = ibv_open_device(ibd[0]);
  if (!ctx)  {
    perror("open device");
    exit(1);
  }
  
  pd = ibv_alloc_pd(ctx);
  if (!pd) {
    fprintf(stderr, "alloc pd failed\n");
    exit(1);
  }

  char *rank = getenv("SLURM_PROCID");
  if (!rank) 
    rank = "seq";


  printf("running allocation test\n"); fflush(stdout);
  u_int64_t max = 64 * 1024;
  max *= 1024 * 1024; // 64 GB
  while (size < max) { 
    printf("%.1lf MB %s\n", (double)size/(double)(1024*1024), rank); fflush(stdout);
    buf = malloc(size);
    if (!buf) {
      perror("malloc");
      exit(1);
    }
    mr = ibv_reg_mr(pd, buf, size, IBV_ACCESS_LOCAL_WRITE | 
                                   IBV_ACCESS_REMOTE_WRITE | 
                                   IBV_ACCESS_REMOTE_READ);
    if (!mr) {
      if (inc > res) {
        size -= inc;
        inc /= 2; // halve the increment
        goto n;
      } else {
        printf("max allocation: %lu\n", size);
        exit(0);
      }
    }
    if (ibv_dereg_mr(mr)) {
       fprintf(stderr,"dereg failed: %lu\n", size);
    }
n:
    free(buf);
    size += inc;
  }

 ibv_dealloc_pd(pd);
 ibv_close_device(ctx);
}
