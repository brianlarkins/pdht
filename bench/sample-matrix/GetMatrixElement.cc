/* Sample code using the Chunks and Tasks programming model.
 * See <http://www.chunks-and-tasks.org>.
 */
#include "GetMatrixElement.h"

CHT_TASK_TYPE_IMPLEMENTATION((GetMatrixElement));
cht::ID GetMatrixElement::execute(CMatrix const & A, CInt const & idx1, CInt const & idx2) {
  int n = A.n;
  if(n <= CMatrix::BLOCK_SIZE) {
    // Lowest level
    return registerChunk( new CDouble(A.elements[idx1*n+idx2]), cht::persistent);
  }
  else {
    // Not lowest level
    // Check which child contains the requested matrix element
    int nHalf = n/2;
    int childIdx1 = 0;
    if(idx1 >= nHalf)
      childIdx1 = 1;
    int childIdx2 = 0;
    if(idx2 >= nHalf)
      childIdx2 = 1;
    int idx1_child = idx1 - childIdx1*nHalf;
    int idx2_child = idx2 - childIdx2*nHalf;
    cht::ChunkID cid_idx1_child = registerChunk( new CInt(idx1_child) );
    cht::ChunkID cid_idx2_child = registerChunk( new CInt(idx2_child) );
    cht::ChunkID cid_child = A.children[childIdx1*2+childIdx2];
    return registerTask<GetMatrixElement>(cid_child, cid_idx1_child, cid_idx2_child, cht::persistent);
  }
} // end execute
