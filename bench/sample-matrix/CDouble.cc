/* Sample code using the Chunks and Tasks programming model.
 * See <http://www.chunks-and-tasks.org>.
 */
#include <cstring>
#include "CDouble.h"

CHT_CHUNK_TYPE_IMPLEMENTATION((CDouble));
void CDouble::writeToBuffer(char * dataBuffer, size_t const bufferSize) const {
  if (bufferSize != getSize())
    throw std::runtime_error("Wrong buffer size to CDouble::writeToBuffer.");
  memcpy(dataBuffer, &x, sizeof(double));
}
size_t CDouble::getSize() const {
  return sizeof(double);
}
void CDouble::assignFromBuffer(char const * dataBuffer, size_t const bufferSize) {
  if (bufferSize != getSize())
    throw std::runtime_error("Wrong buffer size to CDouble::assign_from_buffer.");
  memcpy(&x, dataBuffer, sizeof(double));
}
size_t CDouble::memoryUsage() const {
  return getSize();
}
