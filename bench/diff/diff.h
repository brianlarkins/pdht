/********************************************************/
/*                                                      */
/*    diff.h                                            */
/*                                                      */
/*    author: d. brian larkins                          */
/*    id: $Id: diff.h 611 2007-03-06 15:28:37Z brian $  */
/*                                                      */
/********************************************************/

#ifndef _diff_h
#define _diff_h

#include "diffconst.h"
#include "difftypes.h"

#define MAD_CHILD_LEFT  0
#define MAD_CHILD_RIGHT 1

#define MAX_NODE_LIMIT 20
#define MAX_TREE_DEPTH 10

#define DEFAULT_CHUNKSIZE 1024

/* function structure */
struct func_s {
  int        k;           // first k legendre polynomials
  double     thresh;      // truncation threshold for wavelet coeffs
  double   (*f)(double x);// projection function
  int        max_level;   // don't refine below this level
  int        compressed;  // tree compressed? or not
  double     hg [10][10];    // twoscale relation
  double     hgT[10][10];    // transpose
  double     quad_w[5];      // quadrature weights
  double     quad_x[5];      // quadrature points
  int        quad_npt;       // number of quadrature points (5, really)
  double     quad_phi[5][5]; //
  double     quad_phiT[5][5];//
  double     quad_phiw[5][5];// 
  double     rm[5][5];       // block for derivative operator
  double     r0[5][5];       // block for derivative operator
  double     rp[5][5];       // block for derivative operator
  gt_tree_t *ftree;          // global function tree
}; 
typedef struct func_s func_t;


// from diff_par.h
struct worker_timers_s {
  double total_last,   total;
  double idle_last,    idle;
  double work_last,    work;
  double mvmult_last,  mvmult;
  double rdown_last,   rdown;
  double gcoeff_last,  gcoeff;
  double tcreate_last, tcreate;
  double tcadd_last,   tcadd;
  double alloc_last,   alloc;
  double alloc1_last,   alloc1;
  double alloc2_last,   alloc2;
  double alloc3_last,   alloc3;
};

gt_tree_t     *create_tree(shared gt_context_t *context, int chunksize);
gt_cnp_t      *get_root(gt_tree_t *ftree);
gt_cnp_t      *get_parent(gt_tree_t *ftree, gt_cnp_t *node, gt_cnp_t *buf);
gt_cnp_t      *get_child(gt_tree_t *ftree, gt_cnp_t *node, gt_cnp_t *buf, int childidx);
treedata_t    *get_data(gt_tree_t *ftree, gt_cnp_t *node);
int            get_level(gt_tree_t *ftree, gt_cnp_t *node);
int            get_index(gt_tree_t *ftree, gt_cnp_t *node);
int            has_scaling(gt_tree_t *ftree, gt_cnp_t *node);
int            has_wavelet(gt_tree_t *ftree, gt_cnp_t *node);
double        *get_scaling(func_t *f, gt_tree_t *ftree, gt_cnp_t *node);
double        *get_wavelet(func_t *f, gt_tree_t *ftree, gt_cnp_t *node);
#define        get_left(t,n,b)	    get_child(t,n,b, MAD_CHILD_LEFT)
#define	       get_right(t,n,b)     get_child(t,n,b, MAD_CHILD_RIGHT)

void           set_root(gt_tree_t *ftree, gt_cnp_t *root);
gt_cnp_t      *set_child(gt_tree_t *ftree, gt_cnp_t *parent, int level, int index, int childidx);
void           free_node(gt_tree_t *ftree, gt_cnp_t *node);
void           set_level(gt_tree_t *ftree, gt_cnp_t *node, int level);
void           set_index(gt_tree_t *ftree, gt_cnp_t *node, int index);
void           set_scaling(func_t *f, gt_tree_t *ftree, gt_cnp_t *node, double *scoeffs);
void           set_wavelet(func_t *f, gt_tree_t *ftree, gt_cnp_t *node, double *dcoeffs);
#define        set_left( t,p,l,i)	set_child(t, p, l, i, MAD_CHILD_LEFT)
#define        set_right(t,p,l,i)	set_child(t, p, l, i, MAD_CHILD_RIGHT)
void           print_tree(gt_tree_t *tree);
void           print_subtree(gt_tree_t *tree, gt_cnp_t *t, int indent);
void	       print_node(gt_cnp_t *t);

double dp_last;
double dp_total = 0.0;

/******************************/
/* diff tree helper functions */
/******************************/
gt_tree_t *create_tree(shared gt_context_t *context, int chunksize) { 
  gt_tree_t *gtree;

  gtree =  gt_tree_init(context, chunksize);
  // use local open strategy
  gt_placement_localopen(gtree);
  if (MYTHREAD == 0) {
    gtree->root = gt_node_alloc(gtree, NULL, NULL);
    
  }
  return gtree;
}


gt_cnp_t *get_root(gt_tree_t *ftree) {
  return ftree->root;
}

void set_root(gt_tree_t *ftree, gt_cnp_t *root) {
  ftree->root = root;
}

gt_cnp_t *set_child(gt_tree_t *ftree, gt_cnp_t *parent, int level, int index, int childidx) {
  gt_cnp_t *child;
  tree_t *pnode, *cnode;

  pnode = gt_get_node(ftree, parent);  

  child = gt_node_alloc(ftree, &pnode->children[childidx], parent);
  cnode = gt_get_node(ftree, child);

  cnode->parent = *parent; // struct copy

  for (int i=0;i<GT_NUM_CHILDREN;i++)
    cnode->children[i].ptype = GclCnpPtypeInactive;
  cnode->data.level = level;
  cnode->data.index = index;
  cnode->data.valid = madCoeffNone;

  gt_put_node(ftree, child, cnode);
  gt_put_node(ftree, parent, pnode);
  return child;
}

void free_node(gt_tree_t *ftree, gt_cnp_t *node) {
  gt_cnp_t buf;
  tree_t *p = gt_get_node(ftree, get_parent(ftree, node, &buf));

  for (int i=0;i<GT_NUM_CHILDREN;i++) {
    if (&p->children[i] == node)
      p->children[i].ptype = GclCnpPtypeInactive;
  }
  gt_put_node(ftree, node, p);
  gt_node_free(ftree, node);
}

gt_cnp_t *get_parent(gt_tree_t *ftree, gt_cnp_t *node, gt_cnp_t *buf) { 
  tree_t *t = gt_get_node(ftree, node);
  
  if (!buf)
    buf = talloc(sizeof(gt_cnp_t));
  buf->ptype = t->parent.ptype;
  buf->ci    = t->parent.ci;
  buf->ni    = t->parent.ni;
  gt_finish_node(ftree, t);
  return buf;
}

gt_cnp_t *get_child(gt_tree_t *ftree, gt_cnp_t *node, gt_cnp_t *buf, int childidx) { 
  tree_t *t = gt_get_node(ftree, node);

  if (t->children[childidx].ptype != GclCnpPtypeInactive) {
    if (!buf)
      buf = talloc(sizeof(gt_cnp_t));
    buf->ptype = t->children[childidx].ptype;
    buf->ci    = t->children[childidx].ci;
    buf->ni    = t->children[childidx].ni;
    gt_finish_node(ftree, t);
    return buf;

  } else {
    return NULL;
  }
}


int get_level(gt_tree_t *ftree, gt_cnp_t *node) { 
  tree_t *t =  gt_get_node(ftree, node); 
  int lvl = t->data.level; 
  gt_finish_node(ftree, t);
  return lvl;
}

int get_index(gt_tree_t *ftree, gt_cnp_t *node) { 
  tree_t *t =  gt_get_node(ftree, node); 
  int idx = t->data.index;
  gt_finish_node(ftree, t);
  return idx;
}

int has_scaling(gt_tree_t *ftree, gt_cnp_t *node) {
  tree_t *t =  gt_get_node(ftree, node); 
  int ret = ((t->data.valid == madCoeffScaling) || (t->data.valid == madCoeffBoth));
  gt_finish_node(ftree, t);
  return ret;
}

int has_wavelet(gt_tree_t *ftree, gt_cnp_t *node) {
  tree_t *t =  gt_get_node(ftree, node); 
  int ret = ((t->data.valid == madCoeffWavelet) || (t->data.valid == madCoeffBoth));
  gt_finish_node(ftree, t);
  return ret;
}

double *get_scaling(func_t *f, gt_tree_t *ftree, gt_cnp_t *node) { 
  tree_t *t =  gt_get_node(ftree, node);
  double *copy = NULL;
  int i;

  if ((t->data.valid == madCoeffScaling) || (t->data.valid == madCoeffBoth)) {
    copy = (double *)calloc(f->k, sizeof(double));
    for (i=0;i<f->k;i++)
      copy[i] = t->data.s[i];
    gt_finish_node(ftree, t);
    return copy;
  } else 
    return NULL;
}


double *get_wavelet(func_t *f, gt_tree_t *ftree, gt_cnp_t *node) { 
  tree_t *t =  gt_get_node(ftree, node);
  double *copy = NULL;
  int i;

  if ((t->data.valid == madCoeffWavelet) || (t->data.valid == madCoeffBoth)) {
    copy = (double *)calloc(f->k, sizeof(double));
    for (i=0;i<f->k;i++)
      copy[i] = t->data.d[i];
    gt_finish_node(ftree, t);
    return copy;
  } else 
    return NULL;
}


void set_level(gt_tree_t *ftree, gt_cnp_t *node, int level) {
  tree_t *t =  gt_get_node(ftree, node);
  t->data.level = level;
  gt_put_node(ftree, node, t);
}

void set_index(gt_tree_t *ftree, gt_cnp_t *node, int index) {
  tree_t *t =  gt_get_node(ftree, node);
  t->data.index = index;
  gt_put_node(ftree, node, t);
}

void set_scaling(func_t *f, gt_tree_t *ftree, gt_cnp_t *node, double *scoeffs) {
  tree_t *t =  gt_get_node(ftree, node);
  int i;

  // removing or adding?
  if (!scoeffs) {
    t->data.valid = (t->data.valid == madCoeffBoth) ? madCoeffWavelet : madCoeffNone;

  } else {
    for (i=0;i<f->k;i++)
      t->data.s[i] = scoeffs[i];

    // technically need to check for both && wavelet
    t->data.valid = (t->data.valid == madCoeffWavelet) ? madCoeffBoth : madCoeffScaling;
  }
  gt_put_node(ftree, node, t);
}

void set_wavelet(func_t *f, gt_tree_t *ftree, gt_cnp_t *node, double *dcoeffs) {
  tree_t *t =  gt_get_node(ftree, node);
  int i;

  // removing or adding?
  if (!dcoeffs) {
    t->data.valid = (t->data.valid == madCoeffBoth) ? madCoeffScaling : madCoeffNone;

  } else {
    for (i=0;i<f->k;i++)
      t->data.d[i] = dcoeffs[i];
    t->data.valid = (t->data.valid == madCoeffScaling) ? madCoeffBoth : madCoeffWavelet;
  }
  gt_put_node(ftree, node, t);
}

#if 0
void print_node(gt_cnp_t *t) {
  if (t) {
    //printf("global pointer: %p chunk: %p idx: %d\n", t, t->cp, t->ni);
    tree_t *node = gt_get_node(t);
    if (node)
      printf("   lvl: %d idx: %d \n", 
         get_level(t), get_index(t));
    else
      printf("  node: (null)\n");
  } else
    printf("(null)\n");

}
#endif

void reset_stats(void) {
  dp_total = 0.0;
}

void print_stats(void) {
  printf("TOTAL TIME: %f\n", dp_total);
}

void      filter_inplace(func_t *f, tensor_t *s);
void      unfilter_inplace(func_t *f, tensor_t *s);
tensor_t *transform3d_inplace(tensor_t *s, tensor_t *c, tensor_t *work);
double    truncate_tol(func_t *f, double tol, long level);
void      mTxm(long dimi, long dimj, long dimk, double *c, double *a, double *b);

void      fcube(func_t *f, long n, long lx, long ly, long lz, double (*f)(double,double,double), tensor_t *fcube);
void      vfcube(func_t *f, long n, long lx, long ly, long lz, double (*vf)(long l, double *,double *,double *), tensor_t *vfcube);



#endif // _diff_h
