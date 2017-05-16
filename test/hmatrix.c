/*
 * 
 */
#include <pdht.h>

#define BLOCK_SIZE 200
#define N 800

struct cmatrix_s {
  int n;
  std::vector<double> elements;
  struct cmatrix_s *children[4];
}

create_matrix(int n, int base1, int base2, int type);


matrix *create_matrix(int n, int base1, int base2, int type) {
  matrix children[4];

  if (n <= BLOCK_SIZE) {
    // lowest level

    // A = allocate new CMatrix
    // set size to N*N
    for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
         idx1 = base1 + i;
         idx2 = base2 + j;
         A[i*n+j] = matElementFunc(type, idx1, idx2);
      }
    }
    return A;

  } else {
    // not lowest level
    assert((n % 2) != 0);
    half = n / 2;
    
    for (int i1 = 0; i1 < 2; i1++) {
      cid_base_i1 = base1 + i1 * half;
      for (int i2 = 0; i2 < 2; i2++) {
        cid_base_i2 = base2 + i2 * half;
        children[i1*2+i2] = create_matrix(half, cid_base_i1, cid_base_i2, type);
      }
    }
    return create_matrix_from_ids(n, children[0], children[1], children[2], children[3]);
  }
}

int main(int argc, char **argv) {

   // initialize
   //
   // create matrix A
   // createMatrix(N,0,0,1);
   //
   // create matrix B
   // createMatrix(N,0,0,2);
   //
   // C = A * B
   //
   //
   // verify 20 elements
   // pick random i,j. compute C[i][j] += A[i][k] * B[k][j];

}
