/********************************************************/
/*                                                      */
/*    diff3d.c                                          */
/*                                                      */
/*    author: d. brian larkins                          */
/*    id: $Id: diff_par.c 618 2007-03-06 23:53:38Z dinan $  */
/*                                                      */
/********************************************************/

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

//#include "diffconst.h" 
#include <pdht.h>

#include "tensor.h"
#include "diff3d.h"

#define DEFAULT_K TENSOR_DEFAULT_K

// small:    73
// med  :  2633
// large: 37449
#define THRESHOLD_SMALL  1e-6
#define THRESHOLD_MEDIUM 1e-10
#define THRESHOLD_LARGE  1e-14
#define INITIAL_LEVEL    3

#define PAR_LEVEL        3

#define TASK_AFFINITY_HIGH 0
#define TASK_AFFINITY_LOW  1
#define TASK_PRIORITY      0

/********************************/
/* global shared variables      */
/********************************/
int ncounter = 0;

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
void      compress(func_t *f, madkey_t *node);
void      reconstruct(func_t *f, madkey_t *node);
void      diff(func_t *f, diffdim_t wrtdim, madkey_t *node, func_t *fprime, madkey_t *dnode);
double    eval(func_t *f, madkey_t *node, double x, double y, double z);
void      summarize(func_t *f);
void      summarize_subtree(func_t *f, madkey_t *node, int depth, double *sums, double *diffs);
#if 0
void      print_tree(gt_tree_t tree);
void      print_subtree(gt_tree_t tree, gt_cnp_t *t, int indent, int childidx);
#endif

static double test1(double x, double y, double z);
//static double dtest1(double x, double y, double z);

void      usage(char **argv);
int       main(int argc, char **argv, char **envp);

/*
 * initializes a new function octtree
 */
func_t *init_function(int k, double thresh, double (* test)(double x, double y, double z), int initial_level) {
  func_t *fun;
  int     i;
  int     parlvl = defaultparlvl;
  madkey_t root = { 0, 0, 0, 0 };

  fun = malloc(sizeof(func_t));
  fun->k         = k;
  fun->npt       = k;
  fun->thresh    = thresh;
  fun->f         = test;
  fun->max_level = 30;

  for (i=0; i<NDIM; i++) {
    fun->s0[i] = fun->s[0];
    fun->vk[i] = fun->k;
    fun->vq[i] = fun->npt;
    fun->v2k[i] = 2*fun->k;
  }
  fun->work1 = tensor_create3d(fun->vk[0],fun->vk[1],fun->vk[2],TENSOR_ZERO);
  fun->work2 = tensor_create3d(fun->v2k[0], fun->vk[1], fun->vk[2],TENSOR_ZERO);
  fun->workq = tensor_create3d(fun->vq[0], fun->vq[1], fun->vq[2],TENSOR_ZERO);

  eprintf("   initializing twoscale, quadrature, dc_periodic\n");
  init_twoscale(fun);
  init_quadrature(fun);
  make_dc_periodic(fun);

  fun->compressed = 0;
  fun->ftree = create_tree();

  pdht_barrier();

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

    if (c->rank == 0) {
      // initial function projection @ parlvl 
      fine_scale_projection(fun, &root, parlvl+1);  // xxx

      printf("   projection done.\n");
      //print_tree(fun->ftree);
      // XXX this!
      //refine_fine_scale_projection(fun, &root,parlvl);
    }

    // timing!
    //gtc_process(madtc);
    printf("   refinement done.\n");
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
 * fine_scale_projection - recursive wrapper for fine_scale_project 
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
  node.children = 0xff; // turn on all 8-bits == 8 children.
  pdht_put(f->ftree, nkey, &node); // store new node in PDHT


  // decide whether to create leaves, or keep recursing
  if (nkey->level == (initial_level-1)) {
    // stop, create leaves with scaling coefficients
    fine_scale_project(f,nkey);

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

  printf("creating scaling: children of <%ld, %ld, %ld> @ %ld\n", nkey->x,nkey->y,nkey->z,nkey->level);

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

        fcube(f, f->npt, xlo, ylo, zlo, h, f->f, scoeffs);
        tensor_scale(scoeffs, scale);
        tscoeffs = transform3d(scoeffs, f->quad_phiw);
        cnode.a = ckey; // copy key
        cnode.valid = madCoeffScaling;
        memcpy(&cnode.s, tscoeffs, sizeof(tensor3dk_t));

        // store child node
        ncounter++;
        pdht_put(f->ftree, &ckey, &cnode);

        free(tscoeffs);
      }
    }
  }
}

#if 0
/*
 * projects function f->f to all children of node
 */
void fine_scale_project(func_t *f, node_t *node) {
  long      level   = get_level(f->ftree, node);
  tensor_t *scoeffs = NULL, *tscoeffs = NULL;
  double    h       = 1.0/pow(2.0,level+1);
  double    scale   = sqrt(h);
  gt_cnp_t *cnode;
  long      lx,ly,lz,ix,iy,iz;
  double    xlo,ylo,zlo;
  //long      twoton = pow(2.0,level);
  long      count = 0;
  int       freeflag = 1;

  scale = scale*scale*scale;

  get_xyzindex(f->ftree,node,&lx,&ly,&lz);

  printf("creating scaling: %ld %ld,%ld,%ld\n", level,lx,ly,lz);

  lx *= 2; ly *= 2; lz *= 2;

  scoeffs = tensor_create3d(f->npt,f->npt,f->npt,TENSOR_NOZERO);

  // for each child of node
  for (ix=0;ix<2;ix++) {
    xlo = (lx+ix)*h;
    for (iy=0;iy<2;iy++) {
      ylo = (ly+iy)*h;
      for (iz=0;iz<2;iz++) {
        zlo = (lz+iz)*h;
        cnode = get_child(f->ftree,node,count); 

        // create child if needed
        if (!cnode) {
          cnode = set_child(f->ftree, node, level+1, lx+ix, ly+iy, lz+iz, count);
          freeflag = 0;
        }

        fcube(f, f->npt, xlo, ylo, zlo, h, f->f, scoeffs); // f->quad_x read through f
        tensor_scale(scoeffs, scale);
        tscoeffs = transform3d(scoeffs, f->quad_phiw);
        set_scaling(f, cnode, tscoeffs);
        tfree(tscoeffs); 
        if (freeflag)
          tfree(cnode);
        count++;
      }
    }
  }
  tfree(scoeffs);
}

#endif

#if 0
// PDHT VERSION
void refine_fine_scale_projection(funct_t *f, madkey_t *nkey, int initial_level) {
  long level;
  madkey_t ckey;
  node_t cnode;
  pdht_status_t ret;

  if (nkey.level > f->max_level)
    return;
  if (nkey.level < initial_level) {
    // recur down to initial level
    for (int i=0;i<8;i++) {
      // get child
      ckey.x = nkey.x+lx;
      ckey.y = nkey.y+ly;
      ckey.z = nkey.z+lz;
      ckey.level = nkey.level+1;	  

      ret = pdht_get(f->ftree, &ckey, &cnode);
      if (ret == PdhtStatusOK) {
        refine_fine_scale_projection(f, &ckey, initial_level);
      }
    }
    return;
  }

  // create tasks

}


void refine_fine_scale_project(func_t *f, madkey_t *nkey) {
  gt_cnp_t *cnode = NULL;
  tensor_t *ss = NULL, *sf;
  double  dnorm;
  long    i,j,k;
  long    x,y,z;


  //printf("refine_fine_scale_project: %ld : %ld %ld %ld\n", get_level(f->ftree, node), nkey.x,nkey.y,nkey.z);

  // scaling coeffs _must_ exist at level n+1
  ss = gather_scaling_coeffs(f,node);

  sf = filter(f,ss);

  // fill(ss, 0.0);
  for (i=0;i<f->k;i++) {
    for (j=0;j<f->k;j++) {
      for (k=0;k<f->k;k++) {
        tensor_set3d(sf,i,j,k,0.0);
      }
    }
  }

  dnorm = normf(sf);

  tfree(ss); ss = NULL;
  tfree(sf); sf = NULL;

  // do we need to refine further?
  if (dnorm > f->thresh) {
    // yes, for each child at n+1, project function to n+2, then refine at n+1
    for (i=0;i<8;i++) {
      cnode = get_child(f->ftree, node, i);
      set_scaling(f,cnode,NULL); // delete scaling coeffs from children
      fine_scale_project(f, cnode); // creates nodes at n+2

      // don't create tasks here, should have all been done higher up

      refine_fine_scale_project(f, cnode);
      tfree(cnode); cnode = NULL;
    }
  }
}

#endif


#if 0
void refine_fine_scale_projection(func_t *f, gt_cnp_t *node, int initial_level) {
  long level = get_level(f->ftree, node);
  gt_cnp_t *cnode = NULL;
  task_t *task;
  mad_task_t *madtask;
  long    i;
  long x,y,z;

  get_xyzindex(f->ftree,node,&x,&y,&z);
  //printf("%d refine_fine_scale_projection %ld %ld %ld %ld \n", MYTHREAD, level, x,y,z);

  if (level > f->max_level) 
    return;
  if (level < initial_level) {
    // simply recur down to where the initial coeffs live
    for (i=0;i<8;i++) {
      cnode = get_child(f->ftree, node, i);
      if (cnode) {
        refine_fine_scale_projection(f, cnode, initial_level);
        tfree(cnode);
      }
    }
    return;
  }

  get_xyzindex(f->ftree,node,&x,&y,&z);
  //printf("added tree creation task @ %ld %ld,%ld,%ld\n", get_level(f->ftree,node),x,y,z);

  task = gtc_task_create(sizeof(mad_task_t), refine_project_handle);
  madtask = (mad_task_t *)gtc_task_body(task);
  // create tasks at initial_level
  gt_cnp_copy(f->ftree, node, &madtask->node);
  gtc_add(madtc, madtc->mythread, task, TASK_PRIORITY, TASK_AFFINITY_HIGH);
  gtc_task_destroy(task);
}



void refine_fine_scale_project_wrapper(tc_t *tc, task_t *closure) {
  mad_task_t *madtask;

  madtask = (mad_task_t *)gtc_task_body(closure);
  refine_fine_scale_project(f, &madtask->node);
}
#endif

#if 0 
// PDHT
void refine_fine_scale_project(func_t *f, madkey_t *node) {
  madkey_t *cnode = NULL;
  tensor_t *ss = NULL, *sf;
  double  dnorm;
  long    i,j,k;
  long    x,y,z;


  get_xyzindex(f->ftree, node, &x,&y,&z);
  //printf("refine_fine_scale_project: %ld : %ld %ld %ld\n", get_level(f->ftree, node), x,y,z);

  // scaling coeffs _must_ exist at level n+1
  ss = gather_scaling_coeffs(f,node);

  sf = filter(f,ss);

  // fill(ss, 0.0);
  for (i=0;i<f->k;i++) {
    for (j=0;j<f->k;j++) {
      for (k=0;k<f->k;k++) {
        tensor_set3d(sf,i,j,k,0.0);
      }
    }
  }

  dnorm = normf(sf);

  tfree(ss); ss = NULL;
  tfree(sf); sf = NULL;

  // do we need to refine further?
  if (dnorm > f->thresh) {
    // yes, for each child at n+1, project function to n+2, then refine at n+1
    for (i=0;i<8;i++) {
      cnode = get_child(f->ftree, node, i);
      set_scaling(f,cnode,NULL); // delete scaling coeffs from children
      fine_scale_project(f, cnode); // creates nodes at n+2

      // don't create tasks here, should have all been done higher up

      refine_fine_scale_project(f, cnode);
      tfree(cnode); cnode = NULL;
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
  gt_cnp_t *cnode = NULL;
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
        cnode = get_child(f->ftree,node,count);
        childsc = get_scaling(f,cnode);
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
        tfree(childsc); childsc = NULL; tfree(cnode); cnode = NULL;
      }
    }
  }
  return ss;
}
#endif


#if 0
void create_compress_task(tc_t *tc, gt_cnp_t *node) {
  mad_task_t *madtask;
  task_t     *task;

  task = gtc_task_create(sizeof(mad_task_t), compress_handle);
  madtask = (mad_task_t *)gtc_task_body(task);
  gt_cnp_copy(f->ftree, node, &madtask->node);
  gtc_add(madtc, madtc->mythread, task, TASK_PRIORITY, TASK_AFFINITY_HIGH);
  gtc_task_destroy(task);
}



void compress_wrapper(tc_t *tc, task_t *closure) {
  mad_task_t *madtask;

  madtask = (mad_task_t *)gtc_task_body(closure);
  compress(f, &madtask->node);
}


/*
 * create compression tasks at parlvl
 *   - must call tc_process() collectively later
 *   - call compress(f,get_root(f->tree)); in thread 0 to finish
 */
void par_compress(func_t *f, gt_cnp_t *node, int parlvl) {
  long level = get_level(f->ftree, node);
  gt_cnp_t *cnode = NULL;
  long i;

  if (f->compressed)
    return;

  if (level >= parlvl)
    create_compress_task(madtc, node);
  else {
    for (i=0;i<8;i++) {
      cnode = get_child(f->ftree,node,i);
      if (cnode) {
        par_compress(f,cnode,parlvl);
        tfree(cnode); 
      }
    }
  }
}




void compress(func_t *f, gt_cnp_t *node) {
  tensor_t *ss = NULL, *sf = NULL, *s = NULL;
  gt_cnp_t *cnode = NULL;
  long    i,j,k;
  //long    x,y,z;


  //get_xyzindex(f->ftree, node, &x,&y,&z);
  //printf("compress: %ld : %ld %ld %ld\n", get_level(f->ftree, node), x,y,z);

  if (f->compressed)
    return;

  MSTART_TIMER(mvmult);

  // recursive case: (parallelize here)
  // find leaf nodes with scaling coeffs

  // check first child, all or no children will have scaling coeffs
  cnode = get_child(f->ftree,node,0);

  // sometimes we seg fault here... 
  //if (!cnode) 
  //print_tree(f->ftree);

  MSTOP_TIMER(mvmult);
  if (!has_scaling(f->ftree, cnode)) {
    tfree(cnode);
    for (i=0;i<8;i++) {
      cnode = get_child(f->ftree,node,i);
      compress(f, cnode);
      tfree(cnode);
    }
  } else if (cnode)
    tfree(cnode);

  // base case: (and bottom-up continuation)

  // at parent of leaf nodes - gather scaling @n+1 to n
  ss = gather_scaling_coeffs(f,node);

  MSTART_TIMER(mvmult);

  sf = filter(f,ss); // two-scale
  tfree(ss); ss = NULL;

  s = tensor_create3d(f->k,f->k,f->k,TENSOR_NOZERO);

  // s = sf(0:k,0:k,0:k);
  // d = sf(k:2k,k:2k,k:2k);
  for (i=0;i<f->k;i++) {
    for (j=0;j<f->k;j++) {
      for (k=0;k<f->k;k++) {
        tensor_set3d(s,i,j,k,tensor_get3d(sf,i,j,k));
      }
    }
  }

  MSTOP_TIMER(mvmult);

  // check mra.py (existing scaling coeffs are accumulated into via gaxpy)
  set_scaling(f,node,s);

  MSTART_TIMER(mvmult);
  // fill(ss, 0.0)
  for (i=0;i<f->k;i++) {
    for (j=0;j<f->k;j++) {
      for (k=0;k<f->k;k++) {
        tensor_set3d(sf,i,j,k,0.0);
      }
    }
  }

  MSTOP_TIMER(mvmult);

  // check mra.py (existing wavelet coeffs are accumulated into via gaxpy)
  set_wavelet(f,node,sf);

  for (i=0;i<8;i++)
    set_scaling(f, get_child(f->ftree,node,i), NULL);

  tfree(sf);
  tfree(s);
}


//gt_cnp_t *testo = NULL;
void create_reconstruct_task( tc_t *tc, gt_cnp_t *node) {
  mad_task_t *madtask;
  task_t     *task;
  //  long x,y,z;

  //if (!node)
  //print_tree(f->ftree);

  //get_xyzindex(f->ftree,node, &x,&y,&z);
  //printf("%d: reconstruct task created: %ld %ld,%ld,%ld : %p\n", MYTHREAD, get_level(f->ftree, node),x,y,z, node);

  task = gtc_task_create(sizeof(mad_task_t), reconstruct_handle);
  madtask = (mad_task_t *)gtc_task_body(task);
  gt_cnp_copy(f->ftree, node, &madtask->node);
  gtc_add(madtc, madtc->mythread, task, TASK_PRIORITY, TASK_AFFINITY_HIGH);
  gtc_task_destroy(task);
}



void reconstruct_wrapper(tc_t *tc, task_t *closure) {
  mad_task_t *madtask;
  //long x,y,z; 
  //get_xyzindex(f->ftree, &node, &x,&y,&z);
  //printf("%d reconstruct task: %ld : %ld %ld %ld\n", MYTHREAD, get_level(f->ftree, &node), x,y,z);

  madtask = (mad_task_t *)gtc_task_body(closure);
  reconstruct(f, &madtask->node);
}



/* reconstruct scaling basis tree from wavelet form
 * based on madness3/mra/mra.py
 */
void reconstruct(func_t *f, gt_cnp_t *node) {
  tensor_t *du = NULL, *s = NULL, *d = NULL, *tmp = NULL;
  gt_cnp_t *cnode = NULL;
  long    i,j,k;
  long    ix, iy, iz, iix, iiy, iiz;
  long    count = 0;
  long    level = get_level(f->ftree, node);
  long    x,y,z;

  get_xyzindex(f->ftree, node, &x,&y,&z);
  //printf("%d reconstruct: %ld : %ld %ld %ld\n", MYTHREAD, level, x,y,z);

  if ((d = get_wavelet(f, node))) {
    s = get_scaling(f, node);

    assert(s);

    MSTART_TIMER(mvmult);
    // d(0:k,0:k,0:k) = s;
    for (i=0;i<f->k;i++) {
      for (j=0;j<f->k;j++) {
        for (k=0;k<f->k;k++) {
          tensor_set3d(d,i,j,k,tensor_get3d(s,i,j,k));
        }
      }
    }

    //printf(" %d unfilter: %ld : %ld %ld %ld\n", MYTHREAD, level, x,y,z);

    // two-scale gives us child scaling coeffs
    du = unfilter(f,d,0);

    MSTOP_TIMER(mvmult);

    // remove coeffs from this level of the tree
    set_wavelet(f, node, NULL);
    set_scaling(f, node, NULL);

    //printf(" %d clear coeffs\n", MYTHREAD);

    MSTART_TIMER(mvmult);
    // tmp octant tensor to avoid dealing with slices
    // could save copying... set_wavelet(f,cnode,du[slice])
    tmp = tensor_create3d(f->k,f->k,f->k, TENSOR_NOZERO);
    MSTOP_TIMER(mvmult);

    // for each child of node
    for (ix=0;ix<2;ix++) {
      for (iy=0;iy<2;iy++) {
        for (iz=0;iz<2;iz++) {
          cnode = get_child(f->ftree, node, count);
          count++;
          MSTART_TIMER(mvmult);
          // for each octant of d[2*k,2*k,2*k]
          for (i=0;i<f->k;i++) {
            for (j=0;j<f->k;j++) {
              for (k=0;k<f->k;k++) {
                iix = ix*f->k; iiy = iy*f->k; iiz = iz*f->k;
                tensor_set3d(tmp,i,j,k,tensor_get3d(du,iix+i,iiy+j,iiz+k));
              }
            }
          }
          MSTOP_TIMER(mvmult);
          // update child's scaling coeffs
          set_scaling(f, cnode, tmp); 
          tfree(cnode); cnode = NULL;
        }
      }
    }
    // parallelize me. recursive call to next lower level

    // DO NOT spawn tasks that may overlap lower in the tree.
    for (i=0;i<8;i++) {
      cnode = get_child(f->ftree, node, i);

      //if (cnode) {
      get_xyzindex(f->ftree,cnode,&x,&y,&z);
      //printf("  %d: get_child: %ld :: %ld %ld,%ld,%ld\n", MYTHREAD, i, get_level(f->ftree,cnode),x,y,z);
      //} else {
      //printf("node: %p cnode: %p : %ld,%ld :: %ld,%ld\n", node, cnode, node->ci, node->ni, cnode->ci, cnode->ni);
      //}

      if (level < PAR_LEVEL) {
        MSTART_TIMER(tcreate);
        create_reconstruct_task(madtc,cnode);
        MSTOP_TIMER(tcreate);
      } else 
        reconstruct(f,cnode);
      tfree(cnode); cnode = NULL;
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



void print_tree(gt_tree_t tree) {
  print_subtree(tree, get_root(tree),2,0);
}



void print_subtree(gt_tree_t tree, gt_cnp_t *t, int indent, int childidx) {
  char *spaces = NULL;
  gt_cnp_t *child = NULL;
  long x, y, z;
  tensor_t *s, *d;
  int i;

  spaces = malloc(indent+1);
  for (i=0;i<indent;i++)
    spaces[i] = '.';
  spaces[indent] = '\0';

  get_xyzindex(tree,t,&x,&y,&z);
  printf("%snode: %ld  %ld,%ld,%ld", spaces, get_level(tree,t),x,y,z);

  s = get_scaling(f, t);
  d = get_wavelet(f, t);

  if (s) {
    printf(" s");
    free(s);
  }

  if (d) {
    printf(" d");
    free(d);
  }

  printf(" (%d) \n",childidx);
  for (i=0;i<8;i++) {
    child = get_child(tree,t,i);
    if (child) {
      print_subtree(tree,child,indent+2,i);
      tfree(child);
    }
  }
  free(spaces);
  return;
}

#endif 

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
  printf("  -C <chunksize>  GT Chunksize\n");
}



int main(int argc, char **argv, char **envp) {
  int arg;
  int k = DEFAULT_K; // == 
  double (* test)(double x, double y, double z);
  int chunksize, caching, evaluateme = 0;
  double threshold = THRESHOLD_SMALL;
  double t_compress = 0.0, t_reconstruct = 0.0, t_diff = 0.0;
  double expected, actual;

  chunksize = DEFAULT_CHUNKSIZE;
  defaultparlvl = 2;

  eprintf("initializing global trees.\n");


  // deal with cli args
  while ((arg = getopt(argc, argv, "ehsmlt:C:c")) != -1) {
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
      default:
        printf("%s: unknown option: -%c\n\n", argv[0], arg);
        usage(argv);
        exit(1);
    }
  }

  test  = test1;
  //
  // parallel tree creation
  //

  eprintf("initializing function tree.\n");
  f = init_function(k, threshold, test1, INITIAL_LEVEL);
  eprintf("initializing function tree complete.\n");

  exit(0);
}
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
