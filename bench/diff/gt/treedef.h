/********************************************************/
/*                                                      */
/*    treedef.h                                         */
/*                                                      */
/*    author: d. brian larkins                          */
/*    id: $Id: diff.h 611 2007-03-06 15:28:37Z brian $  */
/*                                                      */
/********************************************************/

#ifndef _treedef_h
#define _treedef_h

#include <unistd.h>
#include <sys/types.h>

#include <gt.h>
#include <tensor.h>

#define NUM_CHILDREN 8

enum coeff_e { madCoeffNone, madCoeffScaling, madCoeffWavelet, madCoeffBoth };
typedef enum coeff_e coeff_t;

struct treedata_s {
  long level;
  long x;
  long y;
  long z;
  coeff_t valid;
  tensor3dk_t  s;  // scaling (sum) coefficients
  tensor3d2k_t d;  // wavelet (difference) coefficients
};
typedef struct treedata_s treedata_t;

struct tree_s {
  u_int8_t     flags;
  gt_cnp_t     parent;
  gt_cnp_t     children[NUM_CHILDREN];
  treedata_t   data;
};
typedef struct tree_s tree_t;


#endif // _treedef_h
