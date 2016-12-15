/********************************************************/
/*                                                      */
/*    operators.c                                      */
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

// function operators
// - square() - squaring of a function // do-all ||ism
// - diff() - differentiate on an axis - top-down, stencil access
// - mul()  - top-down recursive descent of two functions
// - transform() - complicated
// - gaxpy_oop() - alpha*left + beta*right - do-all ||ism
// - sub() - gaxpy(1, left, -1, right)
// - scale() - *= - scaling by a constant - do-all ||ism
// - conj() - complex conjugate of input func - do-all ||ism
// - inner() - element-wise inner product of two function vectors
//           - finds matching coeffs in two trees, order not important?
//           - depends on hash table lookup semantics

// these iterate over all coeffs or values in the hash table
// - binary_op() 
// - unary_op() - inplace unary operation on function values
// - unary_op_coeffs() - inplace unary op applied to coeffs


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
