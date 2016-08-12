/* Sample code using the Chunks and Tasks programming model.
 * See <http://www.chunks-and-tasks.org>.
 */
#include "MatrixAdd.h"
#include "CreateMatrixFromIds.h"

CHT_TASK_TYPE_IMPLEMENTATION((MatrixAdd));
cht::ID MatrixAdd::execute(CMatrix const & A, CMatrix const & B) {
  int nA = A.n;
  int nB = B.n;
  if(nA != nB)
    throw std::runtime_error("Error in MatrixAdd::execute: (nA != nB).");
  int n = nA;
  if(n <= CMatrix::BLOCK_SIZE) {
    // Lowest level
    CMatrix* C = new CMatrix();
    C->n = n;
    C->elements.resize(n*n);
    assert(A.elements.size() == n*n);
    assert(B.elements.size() == n*n);
    for(int i = 0; i < n; i++)
      for(int j = 0; j < n; j++)
	C->elements[i*n+j] = A.elements[i*n+j] + B.elements[i*n+j];
    return registerChunk(C, cht::persistent);
  }
  else {
    // Not lowest level
    cht::ID childTaskIDs[4];
    for(int i = 0; i < 2; i++)
      for(int j = 0; j < 2; j++)
	childTaskIDs[i*2+j] = registerTask<MatrixAdd>(A.children[i*2+j], B.children[i*2+j]);
    cht::ChunkID cid_n = registerChunk( new CInt(n) );
    return registerTask<CreateMatrixFromIds>(cid_n, childTaskIDs[0], childTaskIDs[1], childTaskIDs[2], childTaskIDs[3], cht::persistent);
  }
} // end execute
