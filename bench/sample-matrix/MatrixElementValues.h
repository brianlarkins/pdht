/* Sample code using the Chunks and Tasks programming model.
 * See <http://www.chunks-and-tasks.org>.
 */
#include <cmath>
const int MATRIX_TYPE_A = 1;
const int MATRIX_TYPE_B = 2;
static double matElementFunc(int matType, int i, int j) {
  if(matType == MATRIX_TYPE_A)
    return sin(0.3 + 0.01*i + 0.123*j) + cos(0.4*i) + 0.01 * (i % 2) + 0.02 * (j % 3);
  else if(matType == MATRIX_TYPE_B)
    return cos(0.1 + 0.07*i + 0.432*j) + sin(0.3*i) + 0.05 * (i % 4) + 0.01 * (j % 3);
  else
    return 0;
}
