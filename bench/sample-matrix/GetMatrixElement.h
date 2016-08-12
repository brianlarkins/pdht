/* Sample code using the Chunks and Tasks programming model.
 * See <http://www.chunks-and-tasks.org>.
 */
#include "chunks_and_tasks.h"
#include "CInt.h"
#include "CDouble.h"
#include "CMatrix.h"

struct GetMatrixElement: public cht::Task {
  cht::ID execute(CMatrix const &, CInt const &, CInt const &);
  CHT_TASK_INPUT((CMatrix, CInt, CInt));
  CHT_TASK_OUTPUT((CDouble));
  CHT_TASK_TYPE_DECLARATION;
};
