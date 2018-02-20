/********************************************************/
/*                                                      */
/*    diff3d.c                                          */
/*                                                      */
/*    author: d. brian larkins                          */
/*    id: $Id: diff_par.c 618 2007-03-06 23:53:38Z dinan $  */
/*                                                      */
/********************************************************/

#include <assert.h>
#include <execinfo.h>
#include <inttypes.h>
#include <math.h>
#include <signal.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

//#include "diffconst.h" 
//#include <pdht.h>

#include "tensor.h"
#include "diff3d.h"

#define DEFAULT_K TENSOR_DEFAULT_K

// small:    73
// med  :  2633
// large: 37449
#define THRESHOLD_TEST  .1
#define THRESHOLD_SMALL  1e-6
#define THRESHOLD_MEDIUM 1e-10
//#define THRESHOLD_LARGE  1e-14
#define THRESHOLD_LARGE 1e-16
#define INITIAL_LEVEL    2

#define PAR_LEVEL        3

#define TASK_AFFINITY_HIGH 0
#define TASK_AFFINITY_LOW  1
#define TASK_PRIORITY      0

// need to make timers global so we a. dont get destroyed by MPI
// and b. so we can actually time MPI
pdht_timer_t compress_timer, reconstruct_timer, initialization_timer, diff_timer;

//#define DIFF_UPDATES 

/********************************/
/* global shared variables      */
/********************************/
int ncounter = 0;
int pdhtcounter = 0;
int stcount = 0;
int totcount = 0;

func_t *f, *fprime;
int     defaultparlvl;
mtimers_t mtimers;
extern char *optarg;

int eprintf(const char *format, ...);

/********************************/
/* protos                       */
/********************************/
func_t   *init_function(int k, double thresh, double (* test)(double x, double y, double z), int initial_level);
void      fine_scale_projection(func_t *f, madkey_t *node, long initial_level);
void      fine_scale_project(func_t *f, madkey_t *node);
void      refine_fine_scale_projection(func_t *f, madkey_t *node, int initial_level);
void      refine_fine_scale_project(func_t *f, madkey_t *node);
tensor_t *gather_scaling_coeffs(func_t *f, madkey_t *node);
tensor_t *compress(func_t *f, madkey_t *nkey, int limit, int keep);
void      reconstruct(func_t *f, node_t *node, int limit);
void      recur_down(func_t *f, node_t *node, tensor_t *s);
void      par_compress(func_t *f, int limit);
void      par_reconstruct(func_t *f, int limit);
func_t   *par_diff(func_t *f, diffdim_t wrtdim, int thresh,  double (* test)(double x, double y, double z), int limit);
void      diff(func_t *f, diffdim_t wrtdim, node_t *node, func_t *fprime, node_t *dnode);
double    eval(func_t *f, madkey_t *node, double x, double y, double z);
void      summarize(func_t *f);
void      summarize_subtree(func_t *f, madkey_t *node, int depth, double *sums, double *diffs);
void      print_tree(func_t *f);
void      print_subtree(func_t *f, madkey_t *nkey, int indent, int childidx);

static void print_subtree_keys(madkey_t *keys, int len);


static double test1(double x, double y, double z);
//static double dtest1(double x, double y, double z);

void      usage(char **argv);
int       main(int argc, char **argv, char **envp);


void bthandler(int sig) {
  void *a[100];
  size_t size;
  size = backtrace(a, 100);
  fprintf(stderr, "c->rank : %d Error: signal: %d:\n", c->rank, sig);
  backtrace_symbols_fd(a,size, STDERR_FILENO);
  exit(1);
}

/*
 * initializes a new function octtree
 */
func_t *init_function(int k, double thresh, double (* test)(double x, double y, double z), int initial_level) {
  func_t *fun;
  int     i;
  int     parlvl = initial_level; // for now, parallel traversals are same as initial projection
  madkey_t root = { 0, 0, 0, 0 };
  uint64_t st;

  fun = malloc(sizeof(func_t));
  fun->ftree = create_tree(); // safe to call eprintf() after this

  fun->compressed = 0;
  fun->k         = k;
  fun->npt       = k;
  fun->thresh    = thresh;
  fun->f         = test;
  fun->max_level = 30;

  for (i=0; i<NDIM; i++) {
    fun->s0[i] = fun->s[0]; // XXX where is s[0] initialized?
    fun->vk[i] = fun->k;
    fun->vq[i] = fun->npt;
    fun->v2k[i] = 2*fun->k;
  }
  fun->work1 = tensor_create3d(fun->vk[0],fun->vk[1],fun->vk[2],TENSOR_ZERO);
  fun->work2 = tensor_create3d(fun->v2k[0], fun->vk[1], fun->vk[2],TENSOR_ZERO);
  fun->workq = tensor_create3d(fun->vq[0], fun->vq[1], fun->vq[2],TENSOR_ZERO);

  eprintf("   initializing twoscale, quadrature, dc_periodic");
  init_twoscale(fun);
  init_quadrature(fun);
  make_dc_periodic(fun);
  eprintf("...complete\n");

  /*
   * create an array for all nodes at level k
   * asize = 8^k, madkey_t cindices[8^k];
   * recursively fill out  idx -> { 2x+lx, 2y+ly, 2z+lz }
   * MD for counter in portals
   * do { 
   *     next = atomic_fetch_add(val, +1); 
   *     refine_fine_scale_projection(fun, next);
   * }
   */


  if (fun->f) {
    // cheating. set global func_t *f for || task
    f = fun;

    fun->stlen    = (int)pow(8,parlvl);
    fun->subtrees = calloc(fun->stlen, sizeof(madkey_t));
    pdhtcounter = pdht_counter_init(fun->ftree, 0);
    eprintf("  creating octree with initial projection depth of %d : %d nodes\n", parlvl, fun->stlen);
    // initial function projection @ parlvl 
    // - need to project to subtrees[] level + 1, because refine requires children
    

    PDHT_INIT_ATIMER(initialization_timer);
    fine_scale_projection(fun, &root, parlvl+1);
    //print_subtree_keys(fun->subtrees, fun->stlen);
    eprintf("   projection done.\n");
    pdht_fence(fun->ftree);


    pdht_barrier();
    //print_tree(fun);

    PDHT_START_ATIMER(initialization_timer);
    st = pdht_counter_inc(fun->ftree, fun->counter, 1);
    while (st < fun->stlen) {
      stcount = 0;
      refine_fine_scale_project(fun,  &fun->subtrees[st]);
#ifdef DIFF_UPDATES
      printf("%d: completed subtree %"PRIu64": <%ld,%ld,%ld> @ %ld - %d nodes\n",
          c->rank, st, fun->subtrees[st].x, fun->subtrees[st].y,
          fun->subtrees[st].z, fun->subtrees[st].level, stcount);
#endif
      st = pdht_counter_inc(fun->ftree, fun->counter, 1);
    }

    pdht_fence(fun->ftree);
    // timing!
#ifdef DIFF_UPDATES
    printf("   refinement done: %d nodes\n", totcount);
#endif
  }
  return fun;
}


/*
   - perhaps parallel tree construction
   - if mine, then pdht_insert(), instead of put()
   - keep list of nodes inserted
   - barrier
   - update child pointers on all nodes that were added
   - either mswap or cswap?
   */



/*
 * fine_scale_projection - creates octree top
 * projects function f->f to scaling basis on all nodes 
 * at initial_level
 * @param f the function
 * @param nkey address of the current node to recur from
 * @param initial_level, where to cutoff initial projection
 */
void fine_scale_projection(func_t *f, madkey_t *nkey, long initial_level) {

  madkey_t childkey;
  node_t node;
  long x, y, z;

  // create a new tree node at this level (always)
  node.a = *nkey; // copy current node key into node (coords)
  node.valid = madCoeffNone;
  node.children = 1; // true if have children

  // only rank 0 adds nodes to global space, everyone else is just
  // creating the work-sharing array
  if ((c->rank == 0) && (nkey->level != 0)) {
    //printf("%d: putting: <%ld,%ld,%ld> @ %ld\n", c->rank, nkey->x, nkey->y, nkey->z, nkey->level);
    if (pdht_put(f->ftree, nkey, &node) != PdhtStatusOK) { // store new node in PDHT
      printf("%d: fine_scale_projection: put error\n", c->rank);
      exit(1);
    }
  }


  // decide whether to create leaves, or keep recursing
  if (nkey->level == (initial_level-1)) {
    // stop, create leaves with scaling coefficients
    //this was going to far since ncounter was not being set back to 0
    f->subtrees[ncounter++] = *nkey; // save key at this level for parallel tree creation next

    // only rank 0 creates tree top
    if (c->rank == 0) fine_scale_project(f,nkey);

  } else {
    // keep recursing...

    // coords at next level
    x = nkey->x * 2;
    y = nkey->y * 2;
    z = nkey->z * 2;

    childkey.level = nkey->level + 1;

    for (int lx=0; lx<2; lx++) {
      for (int ly=0; ly<2; ly++) {
        for (int lz=0; lz<2; lz++) {
          childkey.x = x+lx;
          childkey.y = y+ly;
          childkey.z = z+lz;
          fine_scale_projection(f, &childkey, initial_level);
        }
      }
    }
  }
}



/* 
 * projects function f->f to all children of the specified node
 * @param f function to project
 * @param nkey key of node to project
 */
void fine_scale_project(func_t *f, madkey_t *nkey) {

  tensor_t *scoeffs = NULL, *tscoeffs = NULL;
  long      level = nkey->level;
  double    h       = 1.0/pow(2.0,level+1);
  double    scale   = sqrt(h);
  long x, y, z;
  madkey_t ckey;
  node_t cnode;
  double xlo, ylo, zlo;

  scale = scale*scale*scale;
  
  x = nkey->x * 2;
  y = nkey->y * 2;
  z = nkey->z * 2;


  //printf("%d: creating scaling: children of <%ld, %ld, %ld> @ %ld\n", c->rank, nkey->x,nkey->y,nkey->z,nkey->level);
  scoeffs = tensor_create3d(f->npt, f->npt, f->npt, TENSOR_NOZERO);
  
  for (int lx=0; lx<2; lx++) {
    xlo = (x+lx)*h;
    for (int ly=0; ly<2; ly++) {
      ylo = (y+ly)*h;
      for (int lz=0; lz<2; lz++) {
        zlo = (z+lz)*h;

        ckey.x = x+lx;
        ckey.y = y+ly;
        ckey.z = z+lz;
        ckey.level = nkey->level+1;

        //printf("%d: creating  <%ld,%ld,%ld> @ %ld\n", c->rank, ckey.x,ckey.y,ckey.z,ckey.level);
        


        fcube(f, f->npt, xlo, ylo, zlo, h, f->f, scoeffs);
        

        tensor_scale(scoeffs, scale);
        

        tscoeffs = transform3d(scoeffs, f->quad_phiw);
        

        cnode.a = ckey; // copy key
        cnode.valid = madCoeffScaling;
        cnode.children = 0;
             
        memcpy(&cnode.s, tscoeffs, sizeof(tensor3dk_t));
        

        // store child node
        if (pdht_put(f->ftree, &ckey, &cnode) != PdhtStatusOK) { // store new node in PDHT
          printf("%d: fine_scale_project: put error\n", c->rank);
          exit(1);
        }
        
        stcount++; totcount++;
        

        free(tscoeffs);



      }
    }
  }
}


void refine_fine_scale_project(func_t *f, madkey_t *nkey) {
  madkey_t ckey;
  node_t cnode;
  tensor_t *ss = NULL, *sf;
  double  dnorm;
  long    i,j,k;
  long    x,y,z;


  //printf("%d: refine_fine_scale_project: %ld : %ld %ld %ld\n", c->rank, nkey->level, nkey->x,nkey->y,nkey->z);

  // scaling coeffs _must_ exist at level n+1
  if (!(ss = gather_scaling_coeffs(f,nkey))) {
    printf("%d: fatal error getting scaling coeffs\n", c->rank);
    exit(1);
  }

  //ss = self.filter(ss);
  sf = filter(f,ss);

  // fill(ss, 0.0);
  for (i=0;i<f->k;i++) {
    for (j=0;j<f->k;j++) {
      for (k=0;k<f->k;k++) {
        tensor_set3d(sf,i,j,k,0.0); // fill scaling part of tensor with zeros
      }
    }
  }

  dnorm = normf(sf);

  tfree(ss); ss = NULL;
  tfree(sf); sf = NULL;

  // do we need to refine further?
  //printf("%d: %12.7f > %12.7f\n", c->rank, dnorm, f->thresh);
  if (dnorm > f->thresh) {
    // yes, for each child at n+1, project function to n+2, then refine at n+1
    for (int ix=0;ix<2;ix++) {
      for (int iy=0;iy<2;iy++) {
        for (int iz=0;iz<2;iz++) {
          ckey.x = 2 * nkey->x + ix;
          ckey.y = 2 * nkey->y + iy;
          ckey.z = 2 * nkey->z + iz;
          ckey.level = nkey->level + 1;

          //printf("%d: refining further: <%ld,%ld,%ld> @ %ld\n", c->rank, ckey.x, ckey.y, ckey.z, ckey.level);
          pdht_get(f->ftree, &ckey, &cnode);
          cnode.valid = (cnode.valid == madCoeffBoth) ? madCoeffWavelet : madCoeffNone;
          cnode.children = 1;
          pdht_update(f->ftree, &ckey, &cnode);
          fine_scale_project(f, &ckey);
          refine_fine_scale_project(f, &ckey);
        }
      }
    }
  }

}



/* 
 * gather scaling coeffs for node at level n from coeffs at n+1
 *
 * ss[0:9,0:9,0:9] = +0,+0,+0
 * ss[0:9,0:9,9:18] = +0,+0,+1
 * ss[0:9,9:18,0:9] = +0,+1,+0
 * ss[0:9,9:18,9:18] = +0,+1,+1
 * ss[9:18,0:9,0:9] = +1,+0,+0
 * ss[9:18,0:9,9:18] = +1,+0,+1
 * ss[9:18,9:18,0:9] = +1,+1,+0
 * ss[9:18,9:18,9:18] = +1,+1,+1
 *
 */
tensor_t *gather_scaling_coeffs(func_t *f, madkey_t *node) {
  //long level = get_level(f->ftree, node);
  tensor_t *ss = NULL;
  tensor_t *childsc = NULL;
  madkey_t ckey;
  node_t cnode;
  long ix,iy,iz, ixlo,iylo,izlo;
  long i,j,k;
  long count = 0;
  double t;
  
  ss = tensor_create3d(2*f->k,2*f->k,2*f->k, TENSOR_ZERO);
  
  // for each child
  for (ix=0;ix<2;ix++) {
    ixlo = ix*f->k;
    for (iy=0;iy<2;iy++) {
      iylo = iy*f->k;
      for (iz=0;iz<2;iz++) {
        izlo = iz*f->k;

        ckey.x = 2 * node->x + ix;
        ckey.y = 2 * node->y + iy;
        ckey.z = 2 * node->z + iz;
        ckey.level = node->level + 1;

        //printf("%d: asking for <%ld,%ld,%ld>@%ld\n", c->rank, ckey.x,ckey.y,ckey.z,ckey.level);

        if (pdht_get(f->ftree, &ckey, &cnode) != PdhtStatusOK) {
          printf("%d: gather: pdht_get error\n", c->rank);
          return NULL;
        }
        
        if ((cnode.valid == madCoeffScaling) || (cnode.valid == madCoeffBoth))
          childsc = tensor_copy((tensor_t *)&cnode.s);
        else
          childsc = NULL;

        assert(childsc);
        count++;
        //printf("count: %ld %ld %ld %ld\n", count, ixlo,iylo,izlo);
        // copy child scaling coeffs into ss
        for (i=0;i<f->k;i++) {
          for (j=0;j<f->k;j++) {
            //tensor_print(ss);
            for (k=0;k<f->k;k++) {
              t = tensor_get3d(childsc,i,j,k);
              //printf("%ld %ld,%ld,%ld **",(ixlo+i)*ss->h.stride[0]+(iylo+j)*ss->h.stride[1]+(izlo+k)*ss->h.stride[2],ixlo+i,iylo+j,izlo+k);
              tensor_set3d(ss, ixlo+i, iylo+j, izlo+k, t);
              //printf(":");
            }
            //printf("\n");
          }
        }
        tfree(childsc); childsc = NULL;
      }
    }
  }
  return ss;
}



/*
 * create compression tasks at parlvl
 *   - must call tc_process() collectively later
 *   - call compress(f,get_root(f->tree)); in thread 0 to finish
 */
void par_compress(func_t *f, int limit) {
  madkey_t rootkey = { 0, 0, 0, 0};
  uint64_t st;

  if (f->compressed)
    return;

  // bottom up parallel compression, seeding
  //   parallel tasks using subtrees setup in tree creation

  pdht_counter_reset(f->ftree, f->counter);


  st = pdht_counter_inc(f->ftree, f->counter, 1);
  
  while (st < f->stlen) {
    compress(f,  &f->subtrees[st], MAX_TREE_DEPTH+1, 1); // recur to leaf nodes
    //printf("%d: completed compress %ld: <%ld,%ld,%ld> @ %ld\n", 
    //   c->rank, st, f->subtrees[st].x, f->subtrees[st].y,
    //    f->subtrees[st].z, f->subtrees[st].level);
    st = pdht_counter_inc(f->ftree, f->counter, 1);
  }

  pdht_fence(f->ftree);

  // rank 0 finishes up tree top
  if (c->rank == 0) {
    compress(f, &rootkey, limit, 1);
  }
  pdht_fence(f->ftree);
  f->compressed = 1;
}



/**
 * compress - convert a MADNESS function tree to use difference coeffs
 * @param f function tree
 * @param nkey node address of subtree to start at 
 * @param limit recursion limit
 * @param keep true if tree/subtree root and want to keep scaling coeffs
 * @returns tensor containing scaling coeffs of current node
 */
tensor_t *compress(func_t *f, madkey_t *nkey, int limit, int keep) {
  tensor_t *ss = NULL, *sf = NULL, *sc = NULL;
  node_t node;
  node_t cnode;
  madkey_t ckey;
  long ix,iy,iz, ixlo,iylo,izlo;
  int i, j, k;
  double t;

  //printf("%d: compress of <%ld,%ld,%ld> @ %ld\n", c->rank, nkey->x,nkey->y,nkey->z,nkey->level);

  // get the current node
  if (pdht_get(f->ftree, nkey, &node) != PdhtStatusOK) {
    printf("%d: compress: pdht_get error\n", c->rank);
    exit(1);
  }


  // check to see if we're a leaf 
  // or interior node is at the limit depth
  if ((!node.children) || (nkey->level == limit)) {
    // base case, leaf node in octree 

    //print_madnode(&node); printf("\n");

    // copy our scaling coeffs for our caller
    sc = tensor_copy((tensor_t *)&node.s); // cast from tensor3dk_t
    node.valid = madCoeffNone; // mark scaling coeffs as invalid 

  } else {
    // landing spot for children scaling coeffs
    ss = tensor_create3d(2*f->k,2*f->k,2*f->k, TENSOR_ZERO);

    // interior node, recur down to bottom of octree

    // for each child of ours
    for (int ix=0;ix<2;ix++) {
      ixlo = ix*f->k;
      for (iy=0;iy<2;iy++) {
        iylo = iy*f->k;
        for (iz=0;iz<2;iz++) {
          izlo = iz*f->k;

          // compute child coords
          ckey.x = 2 * nkey->x + ix;
          ckey.y = 2 * nkey->y + iy;
          ckey.z = 2 * nkey->z + iz;
          ckey.level = nkey->level + 1;

          // compress each child subtree and collect the scaling coeffs
          // on the way back up
          sc = compress(f, &ckey, limit, 0); // nokeep scaling

          // put them in their rightful place
          for (i=0;i<f->k;i++) {
            for (j=0;j<f->k;j++) {
              for (k=0;k<f->k;k++) {
                t = tensor_get3d(sc,i,j,k);
                tensor_set3d(ss,ixlo+i,iylo+j,izlo+k, t);
              }
            }
          }
          tfree(sc);
        }
      }
    }

    // now we have all of the scaling coeffs, filter them
    sf = filter(f,ss); // two-scale op

    // copy our scaling coeffs into sc for our parent
    // fix up difference coeffs (no scaling coeffs in diff)
    // s = sf(0:k,0:k,0:k);    -- zero out these
    // d = sf(k:2k,k:2k,k:2k); -- these are valid
    sc = tensor_create3d(f->k,f->k,f->k,TENSOR_NOZERO);
    for (i=0;i<f->k;i++) {
      for (j=0;j<f->k;j++) {
        for (k=0;k<f->k;k++) {
          tensor_set3d(sc,i,j,k,tensor_get3d(sf,i,j,k)); // copy for parent
          tensor_set3d(sf,i,j,k,0.0); // fix up wavelet coeffs
        }
      }
    }
    // copy difference coeffs into current node
    memcpy(&node.d, sf, sizeof(tensor3d2k_t));

    // handle tree/subtree root nodes as special case
    if (keep) {
      // avoids cluttering loop above... sigh.
      memcpy(&node.s, sc, sizeof(tensor3dk_t)); // copy back scaling coeffs
      node.valid = madCoeffBoth; // both scaling and wavelet coeffs are at root
    } else {
      node.valid = madCoeffWavelet; // only wavelet coeffs are valid
    }

    // clean up our mess
    tfree(ss);
    tfree(sf);
  }

  // zero scaling coeffs at this node, unless root of tree or subtree
  if (!keep) {
    //memset(&node.s, 0, sizeof(tensor3dk_t));
    //only zero coeffs not the whole thing
    memset(node.s.array, 0, TENSOR_DEFAULT_K * TENSOR_DEFAULT_K * TENSOR_DEFAULT_K * sizeof(double));
  }

  // update global copy of our node
  pdht_update(f->ftree, nkey, &node);

  return sc;
}



void par_reconstruct(func_t *f, int limit) {
  uint64_t st;
  madkey_t rootkey = { 0, 0, 0, 0};
  node_t root, stnode;

  //   rank 0 does top of tree (limit = PAR_LEVEL)
  if (c->rank == 0) {
    
    // get root tree node
    if (pdht_get(f->ftree, &rootkey, &root) != PdhtStatusOK) {
      printf("%d: reconstruct root pdht_get error.\n", c->rank);
      exit(1);
    }
    //XXX check stuff out here
#ifdef DIFF_UPDATES
    tensor_print((tensor_t *)&root.d,1);
    tensor_print((tensor_t *)&root.s,0);
#endif
    // modifies root, but doesn't store into PDHT
    reconstruct(f, &root, limit);

    // update root tree node
    if (pdht_update(f->ftree, &rootkey, &root) != PdhtStatusOK) {
      printf("%d: reconstruct root pdht_update error.\n", c->rank);
      exit(1);
    }
  }
  pdht_fence(f->ftree);


  // parallel reconstruction starts at limit depth
  pdht_counter_reset(f->ftree, f->counter);
  st = pdht_counter_inc(f->ftree, f->counter, 1);
  while (st < f->stlen) {

    // fetch subtree node
    if (pdht_get(f->ftree, &f->subtrees[st], &stnode) != PdhtStatusOK) {
      printf("%d: reconstruct subtree pdht_get error.\n", c->rank);
      exit(1);
    }

    // reconstruct subtree 
    reconstruct(f,  &stnode, MAX_TREE_DEPTH+1); // recur to leaf nodes

    // store modifications to subtree root
    if (pdht_update(f->ftree, &f->subtrees[st], &stnode) != PdhtStatusOK) {
      printf("%d: reconstruct subtree pdht_update error.\n", c->rank);
      exit(1);
    }
#ifdef DIFF_UPDATES
    printf("%d: completed reconstruct %"PRIu64": <%ld,%ld,%ld> @ %ld\n",
        c->rank, st, f->subtrees[st].x, f->subtrees[st].y,
        f->subtrees[st].z, f->subtrees[st].level);
#endif
    st = pdht_counter_inc(f->ftree, f->counter, 1);
  }

  f->compressed = 0;
  pdht_fence(f->ftree);
}



/* reconstruct scaling basis tree from wavelet form
 * based on madness3/mra/mra.py
 */
void reconstruct(func_t *f, node_t *node, int limit) {
  tensor_t *du = NULL, *s = NULL, *d = NULL, *tmp = NULL;
  madkey_t ckey;
  node_t cnode;
  long   i,j,k;
  long   ix, iy, iz, iix, iiy, iiz;
  long   count = 0;


  // reconstruct assumes that scaling/difference coeffs exist at this level
  // and that difference coeffs exist at n+1
  // - upon completion scaling/diff coeffs exist at n+1

  // root has scaling+diff
  // interior at n has scaling+diff
  // interior at n+1 has diff
  // leaf nodes have no coeffs at all after compression

  // don't walk too far down the tree
  if (node->a.level >= limit) return;
#ifdef DIFF_UPDATES
  printf("%d reconstruct: %ld : %ld %ld %ld\n", c->rank, node->a.level, node->a.x,node->a.y,node->a.z);
#endif
  // if have both coeffs, then we must be an interior node
  if (node->valid == madCoeffBoth) {
    d = tensor_copy((tensor_t *)&node->d);
    s = tensor_copy((tensor_t *)&node->s); 

    // d(0:k,0:k,0:k) = s;
    for (i=0;i<f->k;i++) {
      for (j=0;j<f->k;j++) {
        for (k=0;k<f->k;k++) {
          tensor_set3d(d,i,j,k,tensor_get3d(s,i,j,k));
        }
      }
    }

    // two-scale gives us child scaling coeffs
    du = unfilter(f,d,0);

    // remove coeffs from this level of the tree
    memset(&node->d, 0, sizeof(tensor3d2k_t)); // caller must update to make permanent
    memset(&node->s, 0, sizeof(tensor3dk_t));  // caller must call pdht_update
    node->valid = madCoeffNone;                // caller must call pdht_update

    // tmp octant tensor to avoid dealing with slices
    // could save copying... set_wavelet(f,cnode,du[slice])
    tmp = tensor_create3d(f->k,f->k,f->k, TENSOR_NOZERO);

    // for each child of node
    for (ix=0;ix<2;ix++) {
      for (iy=0;iy<2;iy++) {
        for (iz=0;iz<2;iz++) {
          ckey.x = 2 * node->a.x + ix;
          ckey.y = 2 * node->a.y + iy;
          ckey.z = 2 * node->a.z + iz;
          ckey.level = node->a.level + 1;
          if (pdht_get(f->ftree, &ckey, &cnode) != PdhtStatusOK) {
            printf("%d: reconstruct: pdht_get error\n", c->rank);
            exit(1);
          }
          // for each octant of d[2*k,2*k,2*k]
          for (i=0;i<f->k;i++) {
            for (j=0;j<f->k;j++) {
              for (k=0;k<f->k;k++) {
                iix = ix*f->k; iiy = iy*f->k; iiz = iz*f->k;
                tensor_set3d(tmp,i,j,k,tensor_get3d(du,iix+i,iiy+j,iiz+k));
              }
            }
          }
          memcpy(&cnode.s, tmp, sizeof(tensor3dk_t));
          cnode.valid = (cnode.valid == madCoeffWavelet) ? madCoeffBoth : madCoeffScaling;

          // call reconstruct recursively XXX
          reconstruct(f, &cnode, limit);

          // this updates the tree for everything except the call from par_reconstruct
          if (pdht_update(f->ftree, &ckey, &cnode) != PdhtStatusOK) {
            printf("%d: reconstruct: pdht_update error\n", c->rank);
            exit(1);
          }
        }
      }
    }
  }

  if (d)
    tfree(d);
  if (s)
    tfree(s);
  if (du)
    tfree(du);
  if (tmp)
    tfree(tmp);
}



/*
 * what hath god wrought
 *
 */
tensor_t *get_coeffs(func_t *f, diffdim_t wrtdim, node_t *node, madkey_t neighbor) {
  long nextlvl, nextx, nexty, nextz;
  long twon;
  long flevel;
  madkey_t path[30];
  madkey_t fkey; // found key
  node_t fnode;  // found node
  char found = 0;

  if (!node)
    return NULL;

  //printf("get_coeffs: %ld %ld,%ld,%ld from %ld %ld,%ld,%ld\n",neighbor.level,neighbor.x,neighbor.y,neighbor.z, 
  //         node->a.level, node->a.x, node->a.y, node->a.z);

  // check left region boundaries
  if ((neighbor.x < 0) || (neighbor.y < 0) || (neighbor.z < 0))
    return tensor_create3d(f->k,f->k,f->k,TENSOR_ZERO); // empty box

  twon = pow(2.0, node->a.level);
  // check right region boundaries
  if ((neighbor.x >= twon) || (neighbor.y >= twon) || (neighbor.z >= twon))
    return tensor_create3d(f->k,f->k,f->k,TENSOR_ZERO); // empty box


  // if we found what we're looking for, return scaling coeffs
  if ((node->a.level == neighbor.level) &&
      (node->a.x == neighbor.x) &&
      (node->a.y == neighbor.y) &&
      (node->a.z == neighbor.z)) {

    // XXX syntax error tensor3dk_t * != tensor_t  ?
    return (tensor_t *)&node->s; // XXX this may not be safe
  }

  // locate neighbor, can either be above, below, or at same level
  nextlvl = neighbor.level;
  nextx = neighbor.x;
  nexty = neighbor.y;
  nextz = neighbor.z;

  // may need to walk neighbor node up to our level
  while (nextlvl > node->a.level) {
    nextlvl--;
    nextx /= 2; nexty /= 2; nextz /= 2;
  }

  flevel = nextlvl;              // remember final level
  path[nextlvl].x = nextx;       // remember final x coord
  path[nextlvl].y = nexty;       // remember final y coord
  path[nextlvl].z = nextz;       // remember final z coord
  path[nextlvl].level = nextlvl; // remember level


  // try to find 
  while (nextlvl > 0) {
    if (pdht_get(f->ftree, &path[nextlvl], &fnode) == PdhtStatusOK) {
      found = 1;
      break;
    }

    // not found, move on up the tree until we find success
    nextlvl--;
    nextx /= 2; nexty /= 2; nextz /= 2;
    path[nextlvl].x = nextx;       // remember final x coord
    path[nextlvl].y = nexty;       // remember final y coord
    path[nextlvl].z = nextz;       // remember final z coord
    path[nextlvl].level = nextlvl; // remember level
  }

  assert(found);

  while (nextlvl < flevel) {

    recur_down(f, &fnode, (tensor_t *)&fnode.s);

    // XXX - again, may have race condition
    if (pdht_get(f->ftree, &path[nextlvl], &fnode) != PdhtStatusOK) {
      printf("%d: get_coeffs pdht_get error- may be race condition (2).\n", c->rank);
      exit(1);
    }
    nextlvl++; // correct child coords should be in path[]
  }

  // XXX need to return tensor_t *
  return NULL;
}



/* 
 * see __look_out_below() in madness3/mra/mra.py
 */
void recur_down(func_t *f, node_t *node, tensor_t *s) {
  tensor_t *su = NULL, *tmp = NULL;
  madkey_t ckey;
  node_t cnode;
  long ix, iy, iz, iix, iiy, iiz;
  long i,j,k;
  long count = 0;

  //printf("recur_down: %ld %ld,%ld,%ld\n", level,x,y,z);

  // two-scale gives us child scaling coeffs at n+1
  su = unfilter(f,s,1); // sonly=1, we only have scaling coeffs

  // su = [0..k,0..k,0..k]

  // tmp octant tensor to avoid dealing with slices
  tmp = tensor_create3d(f->k,f->k,f->k, TENSOR_NOZERO);

  // for each child of node
  for (ix=0;ix<2;ix++) {
    for (iy=0;iy<2;iy++) {
      for (iz=0;iz<2;iz++) {
        ckey.x = 2 * node->a.x + ix;
        ckey.y = 2 * node->a.y + iy;
        ckey.z = 2 * node->a.z + iz;
        ckey.level = node->a.level + 1;

        // for each octant of su[2*k,2*k,2*k]
        for (i=0;i<f->k;i++) {
          for (j=0;j<f->k;j++) {
            for (k=0;k<f->k;k++) {
              iix = ix*f->k; iiy = iy*f->k; iiz = iz*f->k;
              tensor_set3d(tmp,i,j,k,tensor_get3d(su,iix+i,iiy+j,iiz+k));
            }
          }
        }
        // update child's scaling coeffs
        memcpy(&cnode.s, tmp, sizeof(tensor3dk_t));
        cnode.valid = (cnode.valid == madCoeffWavelet) ? madCoeffBoth : madCoeffScaling;

        // pdht_put() our updated child into f
        // XXX may cause race condition with recursive diff to follow
        if (pdht_put(f->ftree, &ckey, &cnode) != PdhtStatusOK) {
          printf("%d: diff f' pdht_put error.\n", c->rank);
          exit(1);
        }
      }
    }
  }

  if (su)
    free(su);
  if (tmp)
    free(tmp);
}




/*
 * run diff at tree-top and parallel execute diff() on subtrees
 */
func_t *par_diff(func_t *f, diffdim_t wrtdim, int thresh,  double (* test)(double x, double y, double z), int limit) {
  madkey_t rootkey = { 0, 0, 0, 0};
  node_t root, droot, stnode, stdnode;
  int parlvl = limit;
  uint64_t st;
  ncounter = 0;

  // create output tree
  fprime = (func_t *)malloc(sizeof(func_t));
  fprime->ftree = create_tree(); // safe to call eprintf() after this
  PDHT_START_ATIMER(diff_timer);
  fprime->compressed = 0;
  fprime->k         = DEFAULT_K;
  fprime->npt       = DEFAULT_K;
  fprime->thresh    = thresh;
  fprime->f         = test;
  fprime->max_level = 30;

  for (int i=0; i<NDIM; i++) {
    fprime->s0[i] = fprime->s[0]; // XXX check on this?
    fprime->vk[i] = fprime->k;
    fprime->vq[i] = fprime->npt;
    fprime->v2k[i] = 2*fprime->k;
  }
  fprime->work1 = tensor_create3d(fprime->vk[0],fprime->vk[1],fprime->vk[2],TENSOR_ZERO);
  fprime->work2 = tensor_create3d(fprime->v2k[0], fprime->vk[1], fprime->vk[2],TENSOR_ZERO);
  fprime->workq = tensor_create3d(fprime->vq[0], fprime->vq[1], fprime->vq[2],TENSOR_ZERO);

  eprintf("   initializing f' twoscale, quadrature, dc_periodic");
  init_twoscale(fprime);
  init_quadrature(fprime);
  make_dc_periodic(fprime);
  eprintf("...complete\n");

  fprime->stlen    = (int)pow(8,parlvl);
  fprime->subtrees = calloc(fprime->stlen, sizeof(madkey_t));
  pdhtcounter = pdht_counter_init(fprime->ftree, 0);
  pdht_barrier();
  pdht_fence(fprime->ftree);
  eprintf("  creating f' octree with initial projection depth of %d : %d nodes\n", parlvl, fprime->stlen);
  
  // initial function projection @ parlvl 
  // - need to project to subtrees[] level + 1, because refine requires children
  fine_scale_projection(fprime, &rootkey, parlvl+1);

  //print_subtree_keys(fun->subtrees, fun->stlen);
  eprintf("    projection done.\n",c->rank);
  
  pdht_fence(fprime->ftree);
  
  pdht_barrier();

  // parallel differentiation starts at limit depth
  pdht_counter_reset(f->ftree, f->counter);
  st = pdht_counter_inc(f->ftree, f->counter, 1);
  while (st < f->stlen) {

    // fetch f subtree node
    if (pdht_get(f->ftree, &f->subtrees[st], &stnode) != PdhtStatusOK) {
      printf("%d: diff subtree pdht_get error.\n", c->rank);
      exit(1);
    }

    //
    // fetch f' subtree node
    if (pdht_get(fprime->ftree, &fprime->subtrees[st], &stdnode) != PdhtStatusOK) {
      printf("%d: diff subtree pdht_get error.\n", c->rank);
      exit(1);
    }
    
    // diff subtree 
    diff(f, wrtdim, &stnode, fprime, &stdnode); // recur to leaf nodes
    
    // store modifications to f' subtree root (because it may have coeffs from diff)
    if (pdht_update(fprime->ftree, &fprime->subtrees[st], &stdnode) != PdhtStatusOK) {
      printf("%d: diff subtree pdht_update error.\n", c->rank);
      exit(1);
    }
#ifdef DIFF_UPDATES    
    printf("%d: completed difference %"PRIu64": <%ld,%ld,%ld> @ %ld\n",
        c->rank, st, f->subtrees[st].x, f->subtrees[st].y,
        f->subtrees[st].z, f->subtrees[st].level);
#endif
    st = pdht_counter_inc(f->ftree, f->counter, 1);
  }

  pdht_fence(f->ftree);
  return fprime;
}



void diff(func_t *f, diffdim_t wrtdim, node_t *node, func_t *fprime, node_t *dnode) {
  long level, x,y,z;
  double twon;
  long i,j,k, count = 0;
  int ix,iy,iz;
  long px, py, pz, mx, my, mz;
  tensor_t *r, *sm, *sp, *s0 = NULL;
  madkey_t mk, pk, ckey;
  node_t cnode, cdnode;
  // if no scaling coeffs, keep recursing
  if ((node->valid != madCoeffScaling) && (node->valid != madCoeffBoth)) {
    // for each child of node
    for (ix=0;ix<2;ix++) {
      for (iy=0;iy<2;iy++) {
        for (iz=0;iz<2;iz++) {
          ckey.x = 2 * node->a.x + ix;
          ckey.y = 2 * node->a.y + iy;
          ckey.z = 2 * node->a.z + iz;
          ckey.level = node->a.level + 1;

          // fetch f subtree node
          if (pdht_persistent_get(f->ftree, &ckey, &cnode) != PdhtStatusOK) {
            printf("%d: diff f pdht_get error.\n", c->rank);
            exit(1);
          }

          // recursively call diff(f, wrtdim, cnode, fprime, cdnode
          diff(f, wrtdim, &cnode, fprime, &cdnode);

          // pdht_put() a copy of cdnode out to PGLAS
          if (pdht_put(fprime->ftree, &ckey, &cdnode) != PdhtStatusOK) {
            printf("%d: diff f' pdht_put error.\n", c->rank);
            exit(1);
          }
        }
      }
    }
  } else {

    // determine axis of differentiaton
    pk = mk = node->a;
    switch (wrtdim) {
      case Diff_wrtX:
        pk.x += 1;
        mk.x -= 1;
        break;
      case Diff_wrtY:
        pk.y += 1;
        mk.y -= 1;
        break;
      case Diff_wrtZ:
        pk.z += 1;
        mk.z -= 1;
        break;
    }


    // chase down coeffs

    sm = get_coeffs(f, wrtdim, node, mk);
    sp = get_coeffs(f, wrtdim, node, pk);

    s0 = (tensor_t *)&node->s; 

    // find anything?
    if (sm && s0 && sp) {
      twon = pow(2.0,level);

      // rm,rp [0..k,0..k]
      // sm,sp,s0 [0..k,0..k,0..k]

      // s  [9,9,9]
      // ss [18,18,18]

      //MSTART_TIMER(mvmult);
      // rm,r0,rp from make_dc_periodic()
      // ind holds indicies to contract in inner products
      r = inner(f->rp,sm,1,0,NULL);
      inner(f->r0,s0,1,0,r); // work on r in place
      inner(f->rm,sp,1,0,r); // work on r in place
      //MSTOP_TIMER(tensor);

      //mvmult_scale(r,twon); // XXX THIS NEEDS FIXED.

      // XXX this needs fixed too
      void *tmp;
      memcpy(&dnode->s, r, sizeof(tensor3dk_t));
      dnode->valid = (dnode->valid == madCoeffWavelet) ? madCoeffBoth : madCoeffScaling;

      // pdht_put/update of dnode scaling coeffs is done by parent

    } else {
      // these are not the coeffs you are looking for.  

      // project ourselves down to the next level
      // -- may modify f tree, race condition with remote child progress thread
      //    XXX - maybe have a "persistent pdht_get()"
      recur_down(f, node, s0);

      // for each child of node
      for (ix=0;ix<2;ix++) {
        for (iy=0;iy<2;iy++) {
          for (iz=0;iz<2;iz++) {
            ckey.x = 2 * node->a.x + ix;
            ckey.y = 2 * node->a.y + iy;
            ckey.z = 2 * node->a.z + iz;
            ckey.level = node->a.level + 1;

            // fetch f subtree node
            if (pdht_persistent_get(f->ftree, &ckey, &cnode) != PdhtStatusOK) {
              printf("%d: diff f pdht_get error- may be race condition.\n", c->rank);
              exit(1);
            }

            // recursively call diff(f, wrtdim, cnode, fprime, cdnode
            diff(f, wrtdim, &cnode, fprime, &cdnode);

            // pdht_put() a copy of cdnode out to PGLAS
            if (pdht_put(fprime->ftree, &ckey, &cdnode) != PdhtStatusOK) {
              printf("%d: diff f' pdht_put error.\n", c->rank);
              exit(1);
            }
          }
        }
      }
    }
  }
}



#if 0

// ripped off from madness3/mra/mra.py
static double test1(double x, double y, double z) {
  double alpha = 5.0;
  x = x-0.5;
  y = y-0.5;
  z = z-0.5;
  return exp(-alpha*(x*x+y*y+z*z))*5.0;
}



#if 0
static double dtest1(double x, double y, double z) {
  return 0.0;
}
#endif



void summarize(func_t *f) {
  double sums[60], diffs[60];
  int count = 0;
  int printscale = 0;
  int printdiff = 0;


  for (int i=0;i<60;i++)
    sums[i] = diffs[i] = 0.0;

  summarize_subtree(f, get_root(f->ftree), 0, sums, diffs);

  for (int i=0;i<30;i++) {
    count = sums[30+i];
    if (count > 0) {
      if (printscale == 0) {
        printf("sum coefficients\n");
        printscale++;
      }
      printf("  level %3d #boxes=%6d norm=%.8e\n", i, count, sqrt(sums[i]));
    }
  }

  for (int i=0;i<30;i++) {
    count = diffs[30+i];
    if (count > 0) {
      if (printdiff == 0) {
        printf("difference coefficients\n");
        printdiff++;
      }
      printf("  level %3d #boxes=%6d norm=%.8e\n", i, count, sqrt(diffs[i]));
    }
  }
}


void summarize_subtree(func_t *f, gt_cnp_t *node, int depth, double *sums, double *diffs) {
  tensor_t *s, *d;
  gt_cnp_t *cnode = NULL;
  double t;
  long i;
  //  long x,y,z;

  if (depth > 30)
    return;

  if ((s = get_scaling(f,node))) {
    t = normf(s);
    //    get_xyzindex(f->ftree,node,&x,&y,&z);
    //printf("%ld %ld,%ld,%ld : %.0e\n",get_level(f->ftree,node),x,y,z, t);
    sums[depth] += t*t;
    sums[depth+30] += 1.0; // counter
    tfree(s);
  }
  if ((d = get_wavelet(f,node))) {
    t = normf(d);
    diffs[depth] += t*t;
    diffs[depth+30] += 1.0; // counter
    tfree(d);
  }

  for (i=0;i<8;i++) {
    cnode = get_child(f->ftree, node, i);
    if (cnode) {
      summarize_subtree(f,cnode,depth+1,sums,diffs);
      tfree(cnode);
    }
  }
}

#endif

void print_tree(func_t *f) {
  madkey_t k = { 0, 0, 0, 0 };
  if (c->rank == 0) {
    printf("printing octree\n");
    print_subtree(f, &k, 2, 0);
  }
  pdht_barrier();
}



void print_subtree(func_t *f, madkey_t *nkey, int indent, int childidx) {
  pdht_status_t ret;
  char *spaces = NULL;
  node_t node;
  madkey_t ckey;
  long x, y, z;
  tensor_t *s, *d;
  int c=0;

  spaces = malloc(indent+1);
  for (int i=0;i<indent;i++)
    spaces[i] = '.';
  spaces[indent] = '\0';

  ret = pdht_get(f->ftree, nkey, &node);
  if (ret != PdhtStatusOK) {
    printf("%snode: %ld, %ld,%ld,%ld not found\n", spaces, nkey->level, nkey->x, nkey->y, nkey->z);
    return;
  }

  char flag = '-';
  switch(node.valid) {
    case madCoeffBoth:
      flag = 'b';
      break;
    case madCoeffScaling:
      flag = 's';
      break;
    case madCoeffWavelet:
      flag = 'd';
      break;
    case madCoeffNone:
      break;
  }
  printf("%snode: %ld  %ld,%ld,%ld f:%c", spaces, node.a.level, node.a.x, node.a.y,node.a.z, flag);
  printf(" (%d) %c",childidx, node.children ? 'c' : '-');

  s = get_scaling(f, &node);
  d = get_wavelet(f, &node);

  if (s) {
    printf(" s");
    free(s);
  }

  if (d) {
    printf(" d");
    free(d);
  }

  printf("\n");

  x = nkey->x * 2;
  y = nkey->y * 2;
  z = nkey->z * 2;

  if (node.children) {
    for (int lx=0; lx<2; lx++) {
      for (int ly=0; ly<2; ly++) {
        for (int lz=0; lz<2; lz++) {
          ckey.x = x+lx;
          ckey.y = y+ly;
          ckey.z = z+lz;
          ckey.level = nkey->level+1;
#ifdef DIFF_UPDATES
          print_subtree(f, &ckey, indent+2, c);
#endif          
          c++;
        }
      }
    }
  }

  free(spaces);
  return;

}


static void print_subtree_keys(madkey_t *keys, int len) {
  if (c->rank == 0) {
    for (int i=0; i<len; i++) {
      printf("keys[%d]: x: %lu y: %lu z: %lu lvl: %lu\n", i, keys[i].x, keys[i].y, keys[i].z, keys[i].level);
    }
    printf("subtrees contains %d elements.\n", len);
  }
}


// ripped off from madness3/mra/mra.py
static double test1(double x, double y, double z) {
  double alpha = 5.0;
  x = x-0.5;
  y = y-0.5;
  z = z-0.5;
  return exp(-alpha*(x*x+y*y+z*z))*5.0;
}



void usage(char **argv) {
  printf("PDHT Parallel 3-D Madness -- Differentiation\n");
  printf("  Usage: %s [args]\n\n", argv[0]);
  printf("Options:\n");
  printf("  -[s,m,l]        Problem Size: Small, Medium, or Large\n");
}



int main(int argc, char **argv, char **envp) {
  int arg;
  int k = DEFAULT_K; // == 
  double (* test)(double x, double y, double z);
  int chunksize, caching, evaluateme = 0;
  double threshold = THRESHOLD_TEST;
  double t_compress = 0.0, t_reconstruct = 0.0, t_diff = 0.0;
  double expected, actual;
  double local[4];
  double min[4], max[4], avg[4];
  

  pdht_config_t cfg;

  PDHT_INIT_ATIMER(compress_timer);
  PDHT_INIT_ATIMER(reconstruct_timer);

  PDHT_INIT_ATIMER(diff_timer);

  chunksize = DEFAULT_CHUNKSIZE;
  defaultparlvl = INITIAL_LEVEL;

  // deal with cli args
  while ((arg = getopt(argc, argv, "ehsmlt:C:cp:q")) != -1) {
    switch (arg) {
      case 'c':
        caching = 1;
        break;
      case 'C':
        chunksize = atoi(optarg);
        break;
      case 'e':
        evaluateme = 1;
        break;
      case 's':
        threshold = THRESHOLD_SMALL;
        break;
      case 'm':
        threshold = THRESHOLD_MEDIUM;
        break;
      case 'l':
        threshold = THRESHOLD_LARGE;
        //defaultparlvl = 3;
        break;
      case 'h':
        usage(argv);
      case 'q':
        cfg.quiet = 1;
        break;
      case 'p':
        cfg.pendq_size = atoi(optarg);
        break;
      default:
        printf("%s: unknown option: -%c\n\n", argv[0], arg);
        usage(argv);
        exit(1);
    }
  }
  signal(SIGSEGV, bthandler);

  test  = test1;
  //
  // parallel tree creation
  //
  cfg.nptes = 1;
  cfg.pendmode = PdhtPendingTrig;
  cfg.maxentries = 50000;
  cfg.pendq_size = 20000;
  cfg.ptalloc_opts = PTL_PT_MATCH_UNORDERED;
  cfg.rank = PDHT_DEFAULT_RANK_HINT;
  
  pdht_tune(PDHT_TUNE_ALL, &cfg);
  //init timer gets started in init_function
  f = init_function(k, threshold, test1, defaultparlvl);
  PDHT_STOP_ATIMER(initialization_timer);
  eprintf("function tree initialization complete.\n");
  //print_tree(f);
  
//  printf("c->rank : %d pid : %d \n", c->rank, getpid());
  pdht_barrier();
  
  eprintf("compress.\n");
  PDHT_START_ATIMER(compress_timer);
  par_compress(f, defaultparlvl);
  PDHT_STOP_ATIMER(compress_timer);
  
  pdht_barrier();

  eprintf("reconstruct.\n");

  PDHT_START_ATIMER(reconstruct_timer);
  par_reconstruct(f, defaultparlvl);
  PDHT_STOP_ATIMER(reconstruct_timer);
  // have to initialize fprime for differentiation
  //   - don't need to init_function entire tree
  //   - just create tree-top and add subtrees
  
  pdht_barrier();
  
  eprintf("diff.\n");
  //diff timer gets started in diff
//  fprime = par_diff(f,Diff_wrtX, threshold,  test1, defaultparlvl);
  PDHT_STOP_ATIMER(diff_timer);
  
  pdht_barrier();

  local[0] = PDHT_READ_ATIMER_SEC(initialization_timer);
  local[1] = PDHT_READ_ATIMER_SEC(compress_timer);
  local[2] = PDHT_READ_ATIMER_SEC(reconstruct_timer);
  local[3] = PDHT_READ_ATIMER_SEC(diff_timer);
  local[3] = 0;

  pdht_allreduce(local, &avg, PdhtReduceOpSum, DoubleType, 4);
  pdht_allreduce(local, &min, PdhtReduceOpMin, DoubleType, 4);
  pdht_allreduce(local, &max, PdhtReduceOpMax, DoubleType, 4);

  
  if(c->rank == 0){

    printf("MADNESS TIMING size : %d \n", c->size);
#ifndef MPI
    printf("initialiation   min : %12.7f avg : %12.7f max : %12.7f \n", min[0], avg[0] / c->size, max[0]);
    printf("compress        min : %12.7f avg : %12.7f max : %12.7f \n", min[1], avg[1] / c->size, max[1]);
    printf("reconstruct     min : %12.7f avg : %12.7f max : %12.7f \n", min[2], avg[2] / c->size, max[2]);
    printf("diff            min : %12.7f avg : %12.7f max : %12.7f \n", min[3], avg[3] / c->size, max[3]);
#else
    printf("initialiation   min : %12.7f avg : %12.7f max : %12.7f \n", min[0], avg[0] * 2 / c->size, max[0]);
    printf("compress        min : %12.7f avg : %12.7f max : %12.7f \n", min[1], avg[1] * 2 / c->size, max[1]);
    printf("reconstruct     min : %12.7f avg : %12.7f max : %12.7f \n", min[2], avg[2] * 2 / c->size, max[2]);
    printf("diff            min : %12.7f avg : %12.7f max : %12.7f \n", min[3], avg[3] * 2 / c->size, max[3]);
#endif
  }
  pdht_print_stats(f->ftree);
  eprintf("complete.\n");
  fflush(stdout);
  pdht_print_stats(f->ftree);
  return 1;
  exit(0);
}



void print_madnode(void *node) {
  node_t *n = (node_t *)node;
  tensor_t *s, *d;

  printf(" node: <%ld,%ld,%ld>@ %ld %2x", n->a.x, n->a.y,n->a.z, n->a.level, n->children);

  s = get_scaling(f, n);
  d = get_wavelet(f, n);

  if (s) {
    printf(" s");
    free(s);
  }

  if (d) {
    printf(" d");
    free(d);
  }
}


// XXX: call eval() to check results


// compress

// reconstruct

// diff

// 
// check our results
//
#if 0
//print_tree(f->ftree);
if (evaluateme) {
  if (gt_context->mythread == 0) {
    print_tree(f->ftree);
    summarize(f);
    actual = eval(f,get_root(f->ftree), 0.33, 0.44, 0.55);
    expected = test(0.33,0.44,0.55);
    printf("expected: %f actual: %f\n", expected, actual);

    //nodecount = 0;
    //for (i=0;i<THREADS;i++) {
    //nodecount += nodecounter[i];
    //printf("   nodes[%ld]: %ld\n", i, nodecounter[i]);
    //}
    //printf("nodes in octtree: %ld\n", nodecount);
  }
}

gt_abort(gt_context);

//
// compress 
//

if (caching)
  gt_enable_caching(f->ftree);

  gtc_reset(madtc);
  if (gt_context->mythread == 0) {
    eprintf("compress.\n");
    // recur down to defaultparlvl and seed task pool
    //    with 2^parlvl subtrees to compress (64 or 512)
    par_compress(f,get_root(f->ftree),defaultparlvl);
  }
gtc_process(madtc);

// let thread 0 finish up top of tree
if (gt_context->mythread == 0)
  compress(f,get_root(f->ftree));

if (caching)
  gt_disable_caching(f->ftree);

  f->compressed = 1;
  gcl_barrier();

  eprintf("compress complete.\n");
  t_compress = READ_TIMER(tc_process);

  //
  // parallel reconstruction
  //  

if (caching)
  gt_enable_caching(f->ftree);
  gtc_reset(madtc);
  eprintf("reconstruct.\n");
  if (gt_context->mythread == 0) {
    create_reconstruct_task(madtc, get_root(f->ftree));
    f->compressed = 0;
  }
gtc_process(madtc);
t_reconstruct = READ_TIMER(tc_process);
gcl_barrier();
eprintf("reconstruct complete.\n");

if (caching)
  gt_disable_caching(f->ftree);

#if 0
  //
  // pretend to care about the results again
  // 
  if (MYTHREAD == 0) {
    actual = eval(f,get_root(f->ftree), 0.33, 0.44, 0.55);
    expected = test(0.33,0.44,0.55);
    printf("expected: %f actual: %f\n", expected, actual);
  }
gcl_barrier();
#endif

//
// differentiation
// 

exit(0);

if (gt_context->mythread == 0)
  eprintf("creating output tree.\n");
  fprime = init_function(k, threshold, NULL, INITIAL_LEVEL, chunksize);

  gtc_reset(madtc);
  eprintf("diff.\n");
  if (gt_context->mythread == 0) {
    create_diff_task(madtc, Diff_wrtX, get_root(f->ftree), get_root(fprime->ftree));
  }
gtc_process(madtc);
gcl_barrier();
eprintf("diff complete.\n");
t_diff = READ_TIMER(tc_process);

//print_tree(fprime->ftree);
return 0;
}
#endif
