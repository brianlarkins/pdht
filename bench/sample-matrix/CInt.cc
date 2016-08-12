/* Sample code using the Chunks and Tasks programming model.
 * See <http://www.chunks-and-tasks.org>.
 */
#include <cstring>
#include "CInt.h"

CHT_CHUNK_TYPE_IMPLEMENTATION((CInt));
void CInt::writeToBuffer(char * dataBuffer, size_t const bufferSize) const {
  if (bufferSize != getSize())
    throw std::runtime_error("Wrong buffer size to CInt::writeToBuffer.");
  memcpy(dataBuffer, &x, sizeof(int));
}
size_t CInt::getSize() const {
  return sizeof(int);
}
void CInt::assignFromBuffer(char const * dataBuffer, size_t const bufferSize) {
  if (bufferSize != getSize())
    throw std::runtime_error("Wrong buffer size to CInt::assign_from_buffer.");
  memcpy(&x, dataBuffer, sizeof(int));
}
size_t CInt::memoryUsage() const {
  return getSize();
}
