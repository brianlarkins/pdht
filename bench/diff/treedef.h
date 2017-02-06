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

struct madkey_s {
  long x;
  long y;
  long z;
  long level;
};
typedef struct madkey_s madkey_t;

struct node_s {
  madkey_t     a;  // logical address of node
  char         children; // boolean flags for the existence of children
  coeff_t      valid;
  tensor3dk_t  s;  // scaling (sum) coefficients
  tensor3d2k_t d;  // wavelet (difference) coefficients
};
typedef struct node_s node_t;



#endif // _treedef_h
