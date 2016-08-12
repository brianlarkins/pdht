/* Sample code using the Chunks and Tasks programming model.
 * See <http://www.chunks-and-tasks.org>.
 */
#include "chunks_and_tasks.h"
#include "CInt.h"
#include "CMatrix.h"

struct CreateMatrix: public cht::Task {
  cht::ID execute(CInt const &, CInt const &, CInt const &, CInt const &);
  CHT_TASK_INPUT((CInt, CInt, CInt, CInt));
  CHT_TASK_OUTPUT((CMatrix));
  CHT_TASK_TYPE_DECLARATION;
};
