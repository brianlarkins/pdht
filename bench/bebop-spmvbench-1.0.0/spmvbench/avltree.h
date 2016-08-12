#ifndef _avl_tree_h
#define _avl_tree_h
/**
 * @file avltree.h
 * @author mfh
 * @date Time-stamp: <2004-03-23 09:58:38 mhoemmen>
 *
 * @brief An AVL tree holding pairs of integers.
 *
 * An AVL tree holding pairs of integers.  Supplied with a dictionary
 * interface.  Borrowed from the Weiss data structures book web site (which
 * see) and modified to store pairs of integers with a lexicographic
 * comparison.
 ****************************************************************************/

/** A pair of integers. */
struct int_pair
{
  /** First element in the pair. */
  int first;  
  /** second element in the pair. */
  int second; 
};

/** The element in each avl_node is an integer pair.  @see avl_node */
typedef struct int_pair element_type;

/**
 * @struct avl_node
 * Node of an AVL tree.
 * 
 * @note In the original version of avl_node, we were able to hide the
 * internals.  Now we can't, since we need them to write visit functions for
 * tree traversal.  But there isn't anything unexpected in the implementation.
 */
struct avl_node
{
  /** The left subtree. */
  struct avl_node* left;  

  /** The right subtree. */
  struct avl_node* right; 

  /** The datum stored in a node. */
  element_type element;   

  /** Tree height of this node.   */
  int height;             
};


typedef struct avl_node *avl_position;
typedef struct avl_node *avl_tree;

/**
 * Clears the contents of T (if there are any), and returns NULL, 
 * the empty tree.
 */
avl_tree 
make_empty (avl_tree T);

/**
 * Returns NULL if X is not in the tree T.  Otherwise, returns
 * the position of X in T.
 */
avl_position 
find (element_type X, avl_tree T);

/**
 * Returns NULL if T is empty, otherwise returns the position of the
 * minimum element in the tree.
 */
avl_position 
find_min (avl_tree T);

/**
 * Returns NULL if T is empty, otherwise returns the position of the
 * maximum element in the tree.
 */
avl_position 
find_max (avl_tree T);

/**
 * `T = Insert (X,T)' is the correct usage.  Inserts X into T.
 */
avl_tree 
insert (element_type X, avl_tree T);

/**
 * Returns the value stored at position p in an avl_tree.
 */
element_type 
retrieve (avl_position p);

/**
 * Type of a function used to implement visiting each node, and doing
 * something at that node.  The function takes two arguments: first, the
 * current subtree, and second, a void pointer to a working set of data.  Your
 * function is responsible for interpreting the data correctly (this is an
 * unfortunate consequence of using a language like C, which is not
 * type-safe).  The function returns nothing.
 * 
 * @note For programmers not used to function pointers: A typedef for a
 * function pointer looks different than a usual typedef: the type's name is
 * `visit_function'.
 * 
 * @warn The structure of the tree must not be modified during traversal, or 
 * results are undefined.
 */
typedef  void(*visit_function)  (avl_tree, void*);


/**
 * A simple visiting function that just prints out the value at the node.
 * Useful for testing.  Pass in NULL (or anything) for data (no reads or
 * writes to data).
 */
void 
test_visit (avl_tree T, void* data);


/**
 * Traverses the tree in order, visiting each node using the given visit
 * function and using the given data during the visit.  Your visit function is
 * responsible for interpreting the data correctly.
 * 
 * @warn The structure of the tree must not be modified during traversal, or 
 * results are undefined.
 *
 * @param T       Current subtree.
 * @param visit   Function that `visits' (does operations on) the node.
 * @param data    Data to pass into the visit function.
 * 
 * @see visit_function
 */
void
traverse_in_order (avl_tree T, visit_function visit, void* data);


#endif  /* NOT _avl_tree_h */

