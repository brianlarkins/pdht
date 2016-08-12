/* Sample code using the Chunks and Tasks programming model.
 * See <http://www.chunks-and-tasks.org>.
 */
#include "CreateMatrix.h"
#include "CreateMatrixFromIds.h"
#include "MatrixElementValues.h"
#include <cmath>

CHT_TASK_TYPE_IMPLEMENTATION((CreateMatrix));
cht::ID CreateMatrix::execute(CInt const & matSize,
			      CInt const & baseIdx1,
			      CInt const & baseIdx2,
			      CInt const & matType) {
  int n = matSize;
  if(n <= CMatrix::BLOCK_SIZE) {
    // Lowest level
    CMatrix* A = new CMatrix();
    A->n = n;
    A->elements.resize(n*n);
    for(int i = 0; i < n; i++)
      for(int j = 0; j < n; j++) {
	int idx1 = baseIdx1 + i;
	int idx2 = baseIdx2 + j;
	A->elements[i*matSize+j] = matElementFunc(matType, idx1, idx2);
      }
    return registerChunk(A, cht::persistent);
  }
  else {
    // Not lowest level
    if(matSize % 2 != 0)
      throw std::runtime_error("Error in CreateMatrix::execute: matSize not divisible by 2.");
    int nHalf = matSize / 2;
    cht::ChunkID cid_nHalf = registerChunk( new CInt(nHalf) );
    cht::ID childTaskIDs[4];
    for(int i1 = 0; i1 < 2; i1++) {
      cht::ChunkID cid_baseIdx_i1 = registerChunk( new CInt(baseIdx1+i1*nHalf) );
      for(int i2 = 0; i2 < 2; i2++) {
	cht::ChunkID cid_baseIdx_i2 = registerChunk( new CInt(baseIdx2+i2*nHalf) );
	childTaskIDs[i1*2+i2] = registerTask<CreateMatrix>(cid_nHalf, cid_baseIdx_i1, cid_baseIdx_i2, getInputChunkID(matType));
      }
    }
    return registerTask<CreateMatrixFromIds>(getInputChunkID(matSize), childTaskIDs[0], childTaskIDs[1], childTaskIDs[2], childTaskIDs[3], cht::persistent);
  }
} // end execute
