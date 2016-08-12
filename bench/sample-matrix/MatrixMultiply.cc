/* Sample code using the Chunks and Tasks programming model.
 * See <http://www.chunks-and-tasks.org>.
 */
#include "MatrixMultiply.h"
#include "MatrixAdd.h"
#include "CreateMatrixFromIds.h"
#include <cstring>

CHT_TASK_TYPE_IMPLEMENTATION((MatrixMultiply));
cht::ID MatrixMultiply::execute(CMatrix const & A, CMatrix const & B) {
  int nA = A.n;
  int nB = B.n;
  if(nA != nB)
    throw std::runtime_error("Error in MatrixMultiply::execute: (nA != nB).");
  int n = nA;
  if(n <= CMatrix::BLOCK_SIZE) {
    // Lowest level
    CMatrix* C = new CMatrix();
    C->n = n;
    C->elements.resize(n*n);
    memset(&C->elements[0], 0, n*n*sizeof(double));
    for(int i = 0; i < n; i++)
      for(int k = 0; k < n; k++)
	for(int j = 0; j < n; j++) {
	  double Aik = A.elements[i*n+k];
	  double Bkj = B.elements[k*n+j];
	  C->elements[i*n+j] += Aik * Bkj;
	}
    return registerChunk(C, cht::persistent);
  }
  else {
    for(int i = 0; i < 4; i++) {
      if(A.children[i] == cht::CHUNK_ID_NULL)
	throw std::runtime_error("Error in MatrixMultiply::execute: CHUNK_ID_NULL found for A.");
      if(B.children[i] == cht::CHUNK_ID_NULL)
	throw std::runtime_error("Error in MatrixMultiply::execute: CHUNK_ID_NULL found for B.");
    }
    // Not lowest level
    cht::ID childTaskIDs[4];
    for(int i = 0; i < 2; i++)
      for(int j = 0; j < 2; j++) {
	cht::ID childTaskIDsForSum[2];
	for(int k = 0; k < 2; k++)
	  childTaskIDsForSum[k] = registerTask<MatrixMultiply>(A.children[i*2+k], B.children[k*2+j]);
	childTaskIDs[i*2+j] = registerTask<MatrixAdd>(childTaskIDsForSum[0], childTaskIDsForSum[1]);
      }
    cht::ChunkID cid_n = registerChunk( new CInt(n) );
    return registerTask<CreateMatrixFromIds>(cid_n, childTaskIDs[0], childTaskIDs[1], childTaskIDs[2], childTaskIDs[3], cht::persistent);
  }
} // end execute
