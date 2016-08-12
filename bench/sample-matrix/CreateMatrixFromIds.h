/* Sample code using the Chunks and Tasks programming model.
 * See <http://www.chunks-and-tasks.org>.
 */
#include "chunks_and_tasks.h"
#include "CMatrix.h"
#include "CInt.h"

struct CreateMatrixFromIds: public cht::Task {
  cht::ID execute(CInt const &, cht::ChunkID const &, cht::ChunkID const &, cht::ChunkID const &, cht::ChunkID const &);
  CHT_TASK_INPUT((CInt, cht::ChunkID, cht::ChunkID, cht::ChunkID, cht::ChunkID));
  CHT_TASK_OUTPUT((CMatrix));
  CHT_TASK_TYPE_DECLARATION;
};
