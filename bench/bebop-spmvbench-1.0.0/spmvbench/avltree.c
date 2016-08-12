#include "avltree.h"
#include <stdlib.h>
#include <smvm_malloc.h>
#include <stdio.h>

/**
 * Lexicographical comparison function for a pair of ints.  
 * Follows the name return convention as (e.g.) strcmp.
 *
 * @return  -1 if p1 < p2, 0 if p1 == p2, and +1 if p1 > p2.
 */
static int
compare (struct int_pair p1, struct int_pair p2)
{
  if (p1.first < p2.first)
    return -1;
  else if (p1.first == p2.first)
    {
      if (p1.second < p2.second)
	return -1;
      else if (p1.second == p2.second)
	return 0;
      else			/* p1.second > p2.second */
	return 1;
    }
  /* else *//* Compilers might whine about "no return statement" */
  return 1;
}


/**
 * Uses the above lexicographical comparison to find the max of two int pairs.
 */
/*
static struct int_pair
max (struct int_pair lhs, struct int_pair rhs)
{
  if (compare (lhs,rhs) < 0)
    return rhs;
  else
    return lhs;
}
*/
/* The standard int max. */
static int
max (int lhs, int rhs)
{
  return lhs < rhs ? rhs : lhs;
}

/**
 * Convenient shorthand. 
 */
static element_type
element (avl_position p)
{
  return p->element;
}




/*******************************************************
 * Now follows the actual implementation of the tree.  *
 *******************************************************/

avl_tree
make_empty (avl_tree T)
{
  if (T != NULL)
    {
      make_empty (T->left);
      make_empty (T->right);
      free (T);
    }
  return NULL;
}

avl_position
find (element_type X, avl_tree T)
{
  if (T == NULL)
    return NULL;

  /* if (X < element (T)) */
  if (compare (X, element (T)) < 0)
    return find (X, T->left);
  else
    {
      /* if (X > element (T)) */
      if (compare (X, element (T)) > 0)
	return find (X, T->right);
      else
	return T;
    }
}

avl_position
find_min (avl_tree T)
{
  if (T == NULL)
    return NULL;
  else
    {
      if (T->left == NULL)
	return T;
      else
	return find_min (T->left);
    }
}

avl_position
find_max (avl_tree T)
{
  if (T != NULL)
    while (T->right != NULL)
      T = T->right;

  return T;
}

static int
height (avl_position P)
{
  if (P == NULL)
    return -1;
  else
    return P->height;
}


/**
 * Should only be called if K2 has a left child.  Do a rotation
 * between K2 and its left child.  Update heights, and return 
 * the new root.
 */
static avl_position
single_rotate_with_left (avl_position K2)
{
  avl_position K1;

  K1 = K2->left;
  K2->left = K1->right;
  K1->right = K2;

  K2->height = max (height (K2->left), height (K2->right)) + 1;
  K1->height = max (height (K1->left), height (K2)) + 1;

  return K1;			/* New root */
}


/**
 * Should be called only if K1 has a right child.  Do a rotate
 * between K1 and its right child.  Update heights and return 
 * new root.
 */
static avl_position
single_rotate_with_right (avl_position K1)
{
  avl_position K2;

  K2 = K1->right;
  K1->right = K2->left;
  K2->left = K1;

  K1->height = max (height (K1->left), height (K1->right)) + 1;
  K2->height = max (height (K2->right), height (K1)) + 1;

  return K2;			/* New root */
}



/**
 * Should be called only if K3 has a left child and K3's left child
 * has a right child.  Do the left-right double rotation.  Update 
 * heights, and return new root.
 */
static avl_position
double_rotate_with_left (avl_position K3)
{
  /* Rotate between K1 and K2 */
  K3->left = single_rotate_with_right (K3->left);

  /* Rotate between K3 and K2 */
  return single_rotate_with_left (K3);
}


/**
 * Should be called only if K1 has a right child, and K1's right child
 * has a left child.  Do the right-left double rotation.  Update heights,
 * and return new root.
 */
static avl_position
double_rotate_with_right (avl_position K1)
{
  /* Rotate between K3 and K2 */
  K1->right = single_rotate_with_left (K1->right);

  /* Rotate between K1 and K2 */
  return single_rotate_with_right (K1);
}


avl_tree
insert (element_type X, avl_tree T)
{
  if (T == NULL)
    {
      /* Create and return a one-node tree */
      T = smvm_malloc (sizeof (struct avl_node));
      T->element = X;
      T->height = 0;
      T->left = T->right = NULL;
    }
  else
    /* if (X < T->element) */
  if (compare (X, element (T)) < 0)
    {
      T->left = insert (X, T->left);
      if (height (T->left) - height (T->right) == 2)
	/* if (X < element (T->left)) */
	if (compare (X, element (T->left)) < 0)
	  T = single_rotate_with_left (T);
	else
	  T = double_rotate_with_left (T);
    }
  else
    /* if (X > element (T)) */
  if (compare (X, element (T)) > 0)
    {
      T->right = insert (X, T->right);
      if (height (T->right) - height (T->left) == 2)
	/* if (X > element (T->right)) */
	if (compare (X, element (T->right)) > 0)
	  T = single_rotate_with_right (T);
	else
	  T = double_rotate_with_right (T);
    }
  /* Else X is in the tree already; we'll do nothing */

  T->height = max (height (T->left), height (T->right)) + 1;
  return T;
}



element_type
retrieve (avl_position P)
{
  return P->element;
}


void
test_visit (avl_tree T, void *data)
{
  if (T == NULL)
    return;			/* Sanity check (traversal fn should not visit NULL nodes) */

  printf ("(%d,%d)\n", (T->element).first, (T->element).second);
}



void
traverse_in_order (avl_tree T, visit_function visit, void *data)
{
  if (T == NULL)
    return;			/* base case */

  traverse_in_order (T->left, visit, data);
  visit (T, data);
  traverse_in_order (T->right, visit, data);
}
