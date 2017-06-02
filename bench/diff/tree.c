/**********************************************************/
/*                                                        */
/*    tree.c                                              */
/*                                                        */
/*    author: d. brian larkins                            */
/*    id: $Id: tensor.c 618 2007-03-06 23:53:38Z dinan $  */
/*                                                        */
/**********************************************************/
#define _XOPEN_SOURCE 600
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <tensor.h>
#include <diff3d.h>


#define MASK5  0x000000000000001f
#define MASK19 0x000000000007ffff
#define MASK20 0x00000000000fffff

/*
 *  hash function for madness tree coordinates
 */
void treehash(pdht_t *dht, void *key, ptl_match_bits_t *mbits, uint32_t *ptindex, ptl_process_t *rank) {
  madkey_t *kp = (madkey_t *)key;
  uint64_t temp = 0;

  /*
   * MADNESS keys are 4 long values: x, y, z, and level in tree (256 bits)
   * PDHT limits key size to 64-bits, so we grab low-order bits from each key field
   *   x     = 20 bits from 64-bit long                     [44..63]
   *   y     = 20 bits from 64-bit long                     [24..43]
   *   z     = 19 bits from 64-bit long                     [5..23] 
   *   level =  5 bits from 64-bit long  (max depth is 2^5) [0..4]
   */ 

  *mbits = 0;
  temp = (kp->x & MASK20) << 44;
  *mbits |= temp;
  temp = (kp->y & MASK20) << 24;
  *mbits |= temp;
  temp = (kp->z & MASK19) << 5;
  *mbits |= temp;
  temp = (kp->level & MASK5);
  *mbits |= temp;


  (*rank).rank = 1;
  //(*rank).rank = (*mbits) % c->size;
  *ptindex = (*mbits) % dht->ptl.nptes;
  //printf("%d: hashing <%ld,%ld,%ld> @%ld  (%lx @ %d pte: %d)\n", c->rank, kp->x,kp->y,kp->z,kp->level, *mbits, (*rank).rank, *ptindex);
}



/*
 * create_tree - collective call to create function tree
 */
pdht_t *create_tree(void) {
  pdht_t *gtree;
  madkey_t k;
  node_t n;
  gtree = pdht_create(sizeof(madkey_t), sizeof(node_t), PdhtModeStrict);

  pdht_sethash(gtree, treehash);
  // maybe set tunable defaults here

  if (c->rank == 0) {
    // create root node
    k.x = 0; k.y = 0; k.z = 0; k.level = 0;

    memset(&n, 0, sizeof(node_t)); // zero everything
    n.a = k; // set logical address
    n.children = 0xff;

    pdht_put(gtree, &k, &n);
  }

  return gtree;
}


node_t *get_root(pdht_t *gtree) {
  madkey_t k = { 0, 0, 0, 0 };
  node_t *root = malloc(sizeof(node_t));
  pdht_get(gtree, &k, root);
  return root;
}


tensor_t *get_scaling(func_t *f, node_t *node) {
  assert(node);
  if ((node->valid == madCoeffScaling) || (node->valid == madCoeffBoth)) {
    return tensor_copy((tensor_t *)&node->s);
  } else
    return NULL;
}


tensor_t *get_wavelet(func_t *f, node_t *node) {
  assert(node);
  if ((node->valid == madCoeffWavelet) || (node->valid == madCoeffBoth)) {
    return tensor_copy((tensor_t *)&node->d);
  } else
    return NULL;
}


#if 0

gt_cnp_t *get_parent(gt_tree_t ftree, gt_cnp_t *node) {
  return gt_get_parent(ftree, node);
}


gt_cnp_t *get_child(gt_tree_t ftree, gt_cnp_t *node, int childidx) {
  return gt_get_child_cnp(ftree, node, childidx);
  //return gt_get_child(ftree, node, childidx);
}



int child_index(gt_tree_t ftree, gt_cnp_t *pnode, gt_cnp_t *cnode) {
  return gt_child_index(ftree, pnode, cnode);
}



long get_level(node_t *node) {
  return node->a.lvl;
}



void get_xyzindex(node_t *node, long *x, long *y, long *z) {
  *x = node->a.x;
  *y = node->a.y;
  *z = node->a.z;
}



long get_xindex(node_t *node) {
  return node->a.x;
}



long get_yindex(node_t *node) {
  return node->a.y;
}



long get_zindex(node_t *node) {
  return node->a.z;
}



int has_child(gt_tree_t ftree, gt_cnp_t *node, int whichchild) {
  tree_t *t = gt_get_node(ftree,node);
  int ret = (gt_is_active(&(t->children[whichchild])));
  gt_finish_node(ftree, t);
  return ret;
}


int has_scaling(gt_tree_t ftree, gt_cnp_t *node) {
  tree_t *t = gt_get_node(ftree,node);
  int ret = ((t->data.valid == madCoeffScaling) || (t->data.valid == madCoeffBoth));
  gt_finish_node(ftree, t);
  return ret;
}



int has_wavelet(gt_tree_t ftree, gt_cnp_t *node) {
  tree_t *t = gt_get_node(ftree,node);
  int ret = ((t->data.valid == madCoeffWavelet) || (t->data.valid == madCoeffBoth));
  gt_finish_node(ftree, t);
  return ret;
}



tensor_t *get_scaling(func_t *f, gt_cnp_t *node) {
  assert(node);
  tree_t   *t = gt_get_node(f->ftree,node);
  assert(t);
  //printf(" %d: scaling\n", context->mythread);
  if ((t->data.valid == madCoeffScaling) || (t->data.valid == madCoeffBoth)) {
    return tensor_copy((tensor_t *)&t->data.s);
  } else
    return NULL;
}



tensor_t *get_wavelet(func_t *f, gt_cnp_t *node) {
  assert(node);
  tree_t   *t = gt_get_node(f->ftree,node);
  assert(t);
  //printf(" %d: wavelet\n", context->mythread);
  if ((t->data.valid == madCoeffWavelet) || (t->data.valid == madCoeffBoth)) {
    return tensor_copy((tensor_t *)&t->data.d);
  } else
    return NULL;
}



#if 0
void set_root(gt_tree_t ftree, gt_cnp_t *root) {
  gtrees[ftree].root = root;
}
#endif


gt_cnp_t *set_child(gt_tree_t ftree, gt_cnp_t *parent, long level, long x, long y, long z, int childidx) {
  gt_cnp_t *child;
  tree_t *pnode, *cnode;

  pnode = gt_get_node(ftree, parent);

  child = gt_node_alloc(ftree, &pnode->children[childidx], parent);
  cnode = gt_get_node(ftree, child);

  gt_cnp_copy(ftree, parent,&(cnode->parent)); // XXX - can be fixed up

  gt_clear_children(ftree, cnode);
  cnode->data.level = level;
  cnode->data.x = x;
  cnode->data.y = y;
  cnode->data.z = z;
  cnode->data.valid = madCoeffNone;

  gt_put_node(ftree, child, cnode);
  gt_put_node(ftree, parent, pnode);

  return child;
}



void free_node(gt_tree_t ftree, gt_cnp_t *node) {
  gt_cnp_t *parent = get_parent(ftree, node);
  long i;
  tree_t *p = gt_get_node(ftree, parent);

  i = gt_child_index(ftree, parent, node);

  gt_set_inactive(&(p->children[i]));
  gt_put_node(ftree,node,p);
  gt_node_free(ftree,node);
}



void set_level(gt_tree_t ftree, gt_cnp_t *node, long level) {
  tree_t *t = gt_get_node(ftree,node);
  t->data.level = level;
  gt_put_node(ftree, node,t);
}



void set_xyzindex(gt_tree_t ftree, gt_cnp_t *node, long x, long y, long z) {
  tree_t *t = gt_get_node(ftree,node);
  t->data.x = x;
  t->data.y = y;
  t->data.z = z;
  gt_put_node(ftree, node,t);
}



void set_xindex(gt_tree_t ftree, gt_cnp_t *node, long x) {
  tree_t *t = gt_get_node(ftree,node);
  t->data.x = x;
  gt_put_node(ftree, node,t);
}



void set_yindex(gt_tree_t ftree, gt_cnp_t *node, long y) {
  tree_t *t = gt_get_node(ftree,node);
  t->data.y = y;
  gt_put_node(ftree, node,t);
}



void set_zindex(gt_tree_t ftree, gt_cnp_t *node, long z) {
  tree_t *t = gt_get_node(ftree,node);
  t->data.z = z;
  gt_put_node(ftree, node,t);
}



void set_scaling(func_t *f, gt_cnp_t *node, tensor_t *scoeffs) {
  assert(node);
  tree_t *t = gt_get_node(f->ftree, node);
  assert(t);

  if (!scoeffs) {
    t->data.valid = (t->data.valid == madCoeffBoth) ? madCoeffWavelet : madCoeffNone;

  } else {
    memcpy(&t->data.s, scoeffs, sizeof(tensor3dk_t));

    // technically need to check for both && wavelet
    t->data.valid = (t->data.valid == madCoeffWavelet) ? madCoeffBoth : madCoeffScaling;
  }
  gt_put_node(f->ftree,node,t);
}



void set_wavelet(func_t *f, gt_cnp_t *node, tensor_t *dcoeffs) {
  assert(node);
  tree_t *t = gt_get_node(f->ftree, node);
  assert(t);
  if (!dcoeffs) {
    t->data.valid = (t->data.valid == madCoeffBoth) ? madCoeffScaling : madCoeffNone;

  } else {
    // fucking huge!
    memcpy(&t->data.d, dcoeffs, sizeof(tensor3d2k_t));

    // technically need to check for both && wavelet
    t->data.valid = (t->data.valid == madCoeffScaling) ? madCoeffBoth : madCoeffWavelet;
  }
  gt_put_node(f->ftree,node,t);
}


void print_node(gt_tree_t ftree, gt_cnp_t *t) {
  gt_cnp_t *p;
  tree_t   *n;
  int i;

  n = gt_get_node(ftree,t);
  printf("node: %ld %ld,%ld,%ld ", n->data.level, n->data.x, n->data.y, n->data.z);
  if ((n->data.valid == madCoeffScaling) || (n->data.valid == madCoeffBoth))
    printf("s ");
  if ((n->data.valid == madCoeffWavelet) || (n->data.valid == madCoeffBoth))
    printf("d ");
  printf("\n");

  printf("  node pointer: ci: %ld ni: %ld : %p\n", t->ci, t->ni, t);
  gt_finish_node(ftree,n);

  p = get_parent(ftree, t);
  printf("    parent: ci: %ld ni: %ld\n", p->ci, p->ni);
  tfree(p);
  for (i=0;i<8;i++) {
    p = get_child(ftree, t, i);
    if (p)
      printf("    child : ci: %ld ni: %ld\n", p->ci, p->ni);
    else
      printf("    child : (nil)\n");
    tfree(p);
  }

}




void set_children(gt_tree_t ftree, gt_cnp_t *node) {
  long level = get_level(ftree, node);
  long i,j,k;
  long x,y,z;
  long child = 0;

  x = 2*get_xindex(ftree, node);
  y = 2*get_yindex(ftree, node);
  z = 2*get_zindex(ftree, node);

  for (i=0;i<2;i++) {
    for (j=0;j<2;j++) {
      for (k=0;k<2;k++) {
        set_child(ftree, node, level+1, x+i, y+j, z+k, child++);
      }
    }
  }
}
#endif
