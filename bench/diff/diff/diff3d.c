/********************************************************/
/*                                                      */
/*    diff3d.c                                          */
/*                                                      */
/*    author: d. brian larkins                          */
/*    id: $Id: diff_par.c 618 2007-03-06 23:53:38Z dinan $  */
/*                                                      */
/********************************************************/

// TODO:
//
//   - determine caching semantics and fix wrt reconstruct/diff
//   - could cache chunk home directory in local array
//   - update task processing to have ready/waiting queues
//   - update get/put operations to use non-blocking i/o
//     - tasks waiting for _remote_ i/o should move themselves to waiting queue
//
//   - deal with recur_down writing to global tree
//   - find numerical error in either compress/reconstruct
//     - check normf() on tree coeffs in scaling basis (print_subtree())
//   - figure out if we can avoid any locking for updates

//   - need runs for 1,2,4,8,16,32 procs
//     - par tree creation
//     - reconstruct
//     - compression
//     - differentiation


// indep parallel
//   reads  : read from cache only except for miss
//   writes : update cached version and master copy
//            async if possible
//   rlock  : no-op if read is from cache
//   wlock  : upc_lock()/upc_unlock() chunk

// overlapping parallel
//   reads  : read current from master only ?
//   writes : update master copy
//   rlock  : upc_lock()/upc_unlock() chunk
//   wlock  : upc_lock()/upc_unlock() chunk

// qsub -I -lsize=4 -lwalltime=00:30:00 -A CHM022

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

//#include "diffconst.h" 
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
gt_context_t *gt_context;

static task_handle_t refine_project_handle;
static task_handle_t compress_handle;
static task_handle_t reconstruct_handle;
static task_handle_t diff_handle;
tc_t                *madtc; // task collection.

func_t *f, *fprime;
int     defaultparlvl;
mtimers_t mtimers;
extern char *optarg;


/********************************/
/* protos                       */
/********************************/
func_t   *init_function(int k, double thresh, double (* test)(double x, double y, double z), 
                        int initial_level, int chunksize);
void      fine_scale_projection(func_t *f, gt_cnp_t *node, long initial_level);
void      fine_scale_project(func_t *f, gt_cnp_t *node);
void      refine_fine_scale_projection(func_t *f, gt_cnp_t *node, int initial_level);
void      refine_fine_scale_project(func_t *f, gt_cnp_t *node);
tensor_t *gather_scaling_coeffs(func_t *f, gt_cnp_t *node);
void      compress(func_t *f, gt_cnp_t *node);
void      reconstruct(func_t *f, gt_cnp_t *node);
void      diff(func_t *f, diffdim_t wrtdim, gt_cnp_t *node, func_t *fprime, gt_cnp_t *dnode);
double    eval(func_t *f, gt_cnp_t *node, double x, double y, double z);
void      summarize(func_t *f);
void      summarize_subtree(func_t *f, gt_cnp_t *node, int depth, double *sums, double *diffs);
void      print_tree(gt_tree_t tree);
void      print_subtree(gt_tree_t tree, gt_cnp_t *t, int indent, int childidx);

void      refine_fine_scale_project_wrapper(tc_t *tc, task_t *closure);
void      create_compress_task(tc_t *tc, gt_cnp_t *node);
void      create_reconstruct_task(tc_t *tc, gt_cnp_t *node);
void      create_diff_task(tc_t *tc, diffdim_t wrtdim, gt_cnp_t *node, gt_cnp_t *dnode);
void      compress_wrapper(tc_t *tc, task_t *closure);
void      reconstruct_wrapper(tc_t *tc, task_t *closure);
void      diff_wrapper(tc_t *tc, task_t *closure);


static double test1(double x, double y, double z);
//static double dtest1(double x, double y, double z);

void      usage(char **argv);
int       main(int argc, char **argv, char **envp);

func_t *init_function(int k, double thresh, double (* test)(double x, double y, double z), 
                      int initial_level, int chunksize) {
  func_t *fun;
  int     i;
  int     parlvl = defaultparlvl;

  fun = talloc(sizeof(func_t));
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
  fun->ftree = create_tree(gt_context, chunksize);
  gcl_barrier();
  if (fun->f) {
    // cheating. set global func_t *f for || task
    f = fun;

    if (gt_context->mythread == 0) {
      // initial function projection @ parlvl 
      fine_scale_projection(fun, get_root(fun->ftree),parlvl+1);
      printf("   projection done.\n");
      print_tree(fun->ftree);
      refine_fine_scale_projection(fun, get_root(fun->ftree),parlvl);
    }

    // timing!
    gtc_process(madtc);
    printf("   refinement done.\n");
  }
  return fun;
}




/*
 * recursive wrapper for fine_scale_project 
 * projects function f->f to scaling basis on all nodes 
 * at initial_level
 */
void fine_scale_projection(func_t *f, gt_cnp_t *node, long initial_level) {
  long level = get_level(f->ftree, node);
  gt_cnp_t *child = NULL;
  long i=0, x, y, z, lx, ly, lz;
  
  if (level == initial_level - 1) {
    // base case
    fine_scale_project(f, node);
  } else {
    // recur down along each child
    get_xyzindex(f->ftree, node, &x, &y, &z);

    x *= 2; y *= 2; z *= 2;

    for (lx=0;lx<2;lx++) {
      for (ly=0;ly<2;ly++) {
	for (lz=0;lz<2;lz++) {
	  child = set_child(f->ftree, node, level+1, x+lx, y+ly, z+lz, i);
	  fine_scale_projection(f, child, initial_level);
	  i++;
	}
      }
    }
  }
}



/*
 * projects function f->f to all children of node
 */
void fine_scale_project(func_t *f, gt_cnp_t *node) {
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
	if (freeflag) tfree(cnode);
	count++;
      }
    }
  }
  tfree(scoeffs);
}



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



void refine_fine_scale_project(func_t *f, gt_cnp_t *node) {
  gt_cnp_t *cnode = NULL;
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
tensor_t *gather_scaling_coeffs(func_t *f, gt_cnp_t *node) {
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


gt_cnp_t *testo = NULL;
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



/* 
 * see __look_out_below() in madness3/mra/mra.py
 */
void recur_down(func_t *f, gt_cnp_t *node, long level, long x, long y, long z, tensor_t *s) {
  tensor_t *su = NULL, *tmp = NULL;
  gt_cnp_t *cnode = NULL;
  long ix, iy, iz, iix, iiy, iiz;
  long i,j,k;
  long count = 0;

  //printf("recur_down: %ld %ld,%ld,%ld\n", level,x,y,z);
  
  
  // two-scale gives us child scaling coeffs at n+1
  STOP_TIMER(gcoeff);
  MSTART_TIMER(mvmult);
  su = unfilter(f,s,1); // sonly=1, we only have scaling coeffs
  MSTOP_TIMER(mvmult);
  START_TIMER(gcoeff);

  // su = [0..k,0..k,0..k]

  // tmp octant tensor to avoid dealing with slices
  tmp = tensor_create3d(f->k,f->k,f->k, TENSOR_NOZERO);

  // for each child of node
  for (ix=0;ix<2;ix++) {
    for (iy=0;iy<2;iy++) {
      for (iz=0;iz<2;iz++) {

	START_TIMER(alloc);
	cnode = set_child(f->ftree, node, level+1, 2*x+ix, 2*y+iy, 2*z+iz, count);
	STOP_TIMER(alloc);
	count++;

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
	set_scaling(f, cnode, tmp); 
      }
    }
  }

  if (su)
    tfree(su);
  if (tmp)
    tfree(tmp);
}



tensor_t *get_coeffs(func_t *f, diffdim_t wrtdim, gt_cnp_t *node, long level, long x, long y, long z) {	
  gt_cnp_t *pnode = NULL, *cnode = NULL, *tmp = NULL;
  long mylevel, myx, myy, myz;
  long nextlvl, nextx, nexty, nextz;
  long nx[30], ny[30], nz[30]; // path trace holders
  tensor_t *s = NULL, *ts;
  long twon;
  int  whichchild = 0;
  int  xdir = 1, ydir = 1, zdir = 1, xdiff = 0, ydiff = 0, zdiff = 0;
  int  i;

  if (!node)
    return NULL;

  pnode = talloc(sizeof(gt_cnp_t));
  gt_cnp_copy(f->ftree, node, pnode);


  // check left region boundaries
  if ((x < 0) || (y < 0) || (z < 0))
    return tensor_create3d(f->k,f->k,f->k,TENSOR_ZERO); // empty box

  mylevel = get_level(f->ftree, node);
  get_xyzindex(f->ftree, node, &myx, &myy, &myz);

  //printf("get_coeffs: %ld %ld,%ld,%ld from %ld %ld,%ld,%ld\n",level,x,y,z, mylevel, myx, myy,myz);
  
  
  // check right region boundaries
  twon = pow(2.0,level);
  if ((x >= twon) || (y >= twon) || (z >= twon))
    return tensor_create3d(f->k,f->k,f->k,TENSOR_ZERO); // empty box


  // if we've found what we're looking for, return
  if (node && (mylevel == level) && (myx == x) && (myy == y) && (myz == z))
    return get_scaling(f, node);

  // find unifying parent

  nextlvl = level; // tracks path of desired node
  nextx = x; nexty = y; nextz = z;
  nx[nextlvl] = nextx; ny[nextlvl] = nexty; nz[nextlvl] = nextz;

  // may need to walk desired node to our level
  while (nextlvl > mylevel) {
    nextlvl--;
    nextx /= 2; nexty /= 2; nextz /= 2;
    nx[nextlvl] = nextx; ny[nextlvl] = nexty; nz[nextlvl] = nextz;
  } 

  // index difference at common ancestor
  
  // ONLY HANDLING DIFF WRT X RIGHT NOW. Y&Z MAY NOT BE CORRECT
  switch (wrtdim) {
  case Diff_wrtX:
    xdiff = 4;
    break;
  case Diff_wrtY:
    ydiff = 2;
    break;
  case Diff_wrtZ:
    zdiff = 1;
    break;
  }

  // direction difference at common ancestor
  if (nextx < myx)
    xdir = -1;
  
  if (nexty < myy)
    ydir = -1;

  if (nextz < myz)
    zdir = -1;

  // walk up tree
  while (nextlvl > 0) {
    if ((nextx == myx) && (nexty == myy) && (nextz == myz))
      break;
    nextlvl--;
    nextx /= 2; nexty /= 2; nextz /= 2;
    nx[nextlvl] = nextx; ny[nextlvl] = nexty; nz[nextlvl] = nextz;
    myx /= 2; myy /= 2; myz /= 2;
    tmp = pnode;
    pnode = get_parent(f->ftree, tmp);
    whichchild = child_index(f->ftree, pnode, tmp); // XXX - gt_child_index()
    tfree(tmp);
  }

  //  gt_node_lock(pnode); 
  // should be whichchild + (4*[1/-1]) + 0 + 0
  whichchild += (xdiff*xdir)+(ydiff*ydir)+(zdiff*zdir);
  cnode = get_child(f->ftree,pnode, whichchild);
  //  gt_node_unlock(pnode);  

  tfree(pnode);

  //printf(" common ancestor: %ld %ld,%ld,%ld\n", nextlvl, nextx,nexty,nextz);

  // follow stored path down to goal
  for (i=nextlvl;i<level;i++) {
    whichchild = 0;

    if (2*nx[i] != nx[i+1])
      whichchild += 4;
    
    if (2*ny[i] != ny[i+1])
      whichchild += 2;

    if (2*nz[i] != nz[i+1])
      whichchild += 1;

    // printf("    %d %ld,%ld,%ld : %d\n",i,nx[i],ny[i],nz[i], whichchild);
 
    // path may not actually exist. recur down where needed.
    if (!(has_child(f->ftree,cnode,whichchild))) {
      s = get_scaling(f,cnode);
      assert(s);
      recur_down(f,cnode,i,nx[i],ny[i],nz[i],s);
      tfree(s);
    }
    
    // child must now exist
    tmp = cnode;
    cnode = get_child(f->ftree, tmp, whichchild);
    tfree(tmp);
  }
  ts = get_scaling(f,cnode);
  tfree(cnode);
  return ts;
}



void create_diff_task(tc_t *tc, diffdim_t wrtdim, gt_cnp_t *node, gt_cnp_t *dnode) {
  mad_task_t *madtask;
  task_t     *task;

  task = gtc_task_create(sizeof(mad_task_t), diff_handle);
  madtask = (mad_task_t *)gtc_task_body(task);
  madtask->wrtdim = wrtdim;
  gt_cnp_copy(f->ftree, node, &madtask->node);
  gt_cnp_copy(f->ftree, dnode, &madtask->dnode);
  gtc_add(madtc, madtc->mythread, task, TASK_PRIORITY, TASK_AFFINITY_HIGH);
  gtc_task_destroy(task);
}



void diff_wrapper(tc_t *tc, task_t *closure) {
  mad_task_t *madtask;

  madtask = (mad_task_t *)gtc_task_body(closure);
  diff(f, (diffdim_t)madtask->wrtdim, &madtask->node, fprime, &madtask->dnode);
}



void diff(func_t *f, diffdim_t wrtdim, gt_cnp_t *node, func_t *fprime, gt_cnp_t *dnode) {
  gt_cnp_t *cnode, *cdnode;
  long level, x,y,z;
  double twon;
  long i,j,k, count = 0;
  long px, py, pz, mx, my, mz;
  tensor_t *r, *sm, *sp, *s0;

  // get level and index
  level = get_level(f->ftree,node);
  get_xyzindex(f->ftree,node,&x,&y,&z);

  // if no scaling coeffs, recur down each child in || 
  if (!has_scaling(f->ftree, node)) {
    for (i=0;i<2;i++) {
      for (j=0;j<2;j++) {
	for (k=0;k<2;k++) {
	  cnode = get_child(f->ftree, node, count);
	  if (!cnode)
	    print_tree(f->ftree);
	  START_TIMER(alloc);
	  cdnode = set_child(fprime->ftree, dnode, level+1,2*x+i,2*y+j,2*z+k,count);
	  STOP_TIMER(alloc);
	  
	  if (level < PAR_LEVEL) {
	    create_diff_task(madtc, wrtdim, cnode, cdnode);
	  } else
	    diff(f,wrtdim,cnode,fprime,cdnode);

	  tfree(cnode);
	  count++;
	}
      }

    }
  } else {
    // determine axis of differentiaton
    px = mx = x; 
    py = my = y;
    pz = mz = z;
    switch (wrtdim) {
    case Diff_wrtX:
      px += 1;
      mx -= 1;
      break;
    case Diff_wrtY:
      py += 1;
      my -= 1;
      break;
    case Diff_wrtZ:
      pz += 1;
      mz -= 1;
      break;
    }

    // chase down coeffs
    START_TIMER(gcoeff);
    sm = get_coeffs(f,wrtdim,node,level,mx,my,mz);
    sp = get_coeffs(f,wrtdim,node,level,px,py,pz);
    s0 = get_scaling(f,node);
    STOP_TIMER(gcoeff);

    // find anything?
    if (sm && s0 && sp) {
      twon = pow(2.0,level);

      // rm,rp [0..k,0..k]
      // sm,sp,s0 [0..k,0..k,0..k]

      // s  [9,9,9]
      // ss [18,18,18]

      MSTART_TIMER(mvmult);
      // rm,r0,rp from make_dc_periodic()
      // ind holds indicies to contract in inner products
      r = inner(f->rp,sm,1,0,NULL);
      inner(f->r0,s0,1,0,r); // work on r in place
      inner(f->rm,sp,1,0,r); // work on r in place
      MSTOP_TIMER(tensor);

      //mvmult_scale(r,twon); // XXX THIS NEEDS FIXED.
      set_scaling(f,dnode,r);
      
    // these are not the coeffs you are looking for.  
    } else {

      // project ourselves down to the next level...
      recur_down(f,node,level,x,y,z,s0);

      count = 0;
      // and then differentiate from there.
      for (i=0;i<2;i++) {
	for (j=0;j<2;j++) {
	  for (k=0;k<2;k++) {
	    cnode = get_child(f->ftree, node, count);
	    START_TIMER(alloc);
	    cdnode = set_child(fprime->ftree, dnode, level+1,2*x+i,2*y+j,2*z+k,count);
	    STOP_TIMER(alloc);
	    diff(f,wrtdim,cnode,fprime,cdnode); // parallelize here?
	    tfree(cnode);
	    count++;
	  }
	}
      }
    }
  }
}



/* 
 * evaluate f(x,y,z) madness3/mra/cfcube.c
 */
double eval(func_t *f, gt_cnp_t *node, double x, double y, double z) {
  long level;
  double px[100],py[100],pz[100], sum = 0.0;
  gt_cnp_t *curnode = NULL, *tnode = NULL;
  tensor_t *s = NULL;
  long ix, iy, iz;
  long p,q,r;
  long index;
  double xx,yy,zz;
  double *ptr;
  double twon, twoinv;
  double aa,bb,cc;
  
  // only work on uncompressed functions
  // code for compressed is in mra.py:1172
  if (f->compressed)
    return 0.0;

  // no aliases...
  curnode = talloc(sizeof(gt_cnp_t));
  gt_cnp_copy(f->ftree, node,curnode);

  while (!has_scaling(f->ftree,curnode)) {
    // recur down tree from root until we find scaling coeffs
    // find the child containing the spatial box with x,y,z
    get_xyzindex(f->ftree,curnode,&ix,&iy,&iz);
    level = get_level(f->ftree, curnode);
    
    if (level > f->max_level)
      return 0.0;

    index = 0;

    //printf("eval %ld:%ld,%ld,%ld for %f,%f,%f \n",level,ix,iy,iz,x,y,z);

    // #nodes on next level down
    twon = pow(2.0,level+1);
    twoinv = 1.0/twon;

    ix *= 2; iy *= 2; iz *= 2;
    
    aa = ((double) ix)*twoinv;
    bb = ((double) iy)*twoinv;
    cc = ((double) iz)*twoinv;

    if (z > (cc+twoinv)) {
      index += 1; 
      iz +=1;
    }

    if (y > (bb+twoinv)) {
      index += 2; 
      iy +=1;
    }

    if (x > (aa+twoinv)) {
      index += 4; 
      ix +=1;
    }
    
    tnode = get_child(f->ftree, curnode, index);
    tfree(curnode);
    curnode = tnode;
  }

  level = get_level(f->ftree, curnode);

  // found scaling coeffs, compute function value.
  s = get_scaling(f, curnode);
  tfree(curnode);
  ptr = s->array; // tensor abstraction violation... 

  // hoodoo magic
  xx = x * twon; yy = y * twon; zz = z * twon;
  xx = xx - ix; yy = yy - iy; zz = zz - iz;
   
  phi(xx,f->k,px);
  phi(yy,f->k,py);
  phi(zz,f->k,pz);

  for (p=0;p<f->k;p++) {
    for (q=0;q<f->k;q++) {
      for (r=0;r<f->k;r++) {
	sum = sum + *ptr++ * px[p]*py[q]*pz[r];
      }
    }
  }
  tfree(s);
  return sum*pow(2.0,1.5*level);
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



void usage(char **argv) {
  printf("GT Parallel 3-D Madness -- Differentiation\n");
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

  gt_context = gt_init(&argc, &argv);

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
  madtc = gtc_create(sizeof(mad_task_t), 1, 1024, MPI_COMM_WORLD);
  refine_project_handle = gtc_register(madtc, refine_fine_scale_project_wrapper);
  compress_handle       = gtc_register(madtc, compress_wrapper);
  reconstruct_handle    = gtc_register(madtc, reconstruct_wrapper);
  diff_handle           = gtc_register(madtc, diff_wrapper);

  //
  // parallel tree creation
  //

  eprintf("initializing function tree.\n");
  f = init_function(k, threshold, test, INITIAL_LEVEL, chunksize);
  eprintf("initializing function tree complete.\n");

  // 
  // check our results
  //

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
