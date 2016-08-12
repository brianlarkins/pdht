/* Sample code using the Chunks and Tasks programming model.
 * See <http://www.chunks-and-tasks.org>.
 */
#include "CreateMatrixFromIds.h"

CHT_TASK_TYPE_IMPLEMENTATION((CreateMatrixFromIds));
cht::ID CreateMatrixFromIds::execute(CInt const & n, cht::ChunkID const & id1, cht::ChunkID const & id2, cht::ChunkID const & id3, cht::ChunkID const & id4) {
  CMatrix* A = new CMatrix();
  A->n = n;
  A->children[0] = copyChunk(id1);
  A->children[1] = copyChunk(id2);
  A->children[2] = copyChunk(id3);
  A->children[3] = copyChunk(id4);
  return registerChunk(A, cht::persistent);
}
