/* Sample code using the Chunks and Tasks programming model.
 * See <http://www.chunks-and-tasks.org>.
 */
#include <cstring>
#include "CMatrix.h"

CHT_CHUNK_TYPE_IMPLEMENTATION((CMatrix));
void CMatrix::writeToBuffer(char * dataBuffer, size_t const bufferSize) const {
  if (bufferSize != getSize())
    throw std::runtime_error("Wrong buffer size to CMatrix::writeToBuffer.");
  if(n <= BLOCK_SIZE) {
    // Lowest level
    memcpy(dataBuffer, &n, sizeof(int));
    memcpy(dataBuffer+sizeof(int), &elements[0], n*n*sizeof(double));
  }
  else {
    // Not lowest level
    memcpy(dataBuffer, &n, sizeof(int));
    memcpy(dataBuffer+sizeof(int), &children[0], 4*sizeof(cht::ChunkID));
  }
}
size_t CMatrix::getSize() const {
  if(n <= BLOCK_SIZE) {
    // Lowest level
    return sizeof(int) + n*n*sizeof(double);
  }
  else {
    // Not lowest level
    return sizeof(int) + 4*sizeof(cht::ChunkID);
  }
}
void CMatrix::assignFromBuffer(char const * dataBuffer, size_t const bufferSize) {
  if (bufferSize < sizeof(int))
    throw std::runtime_error("Wrong buffer size to CMatrix::assign_from_buffer.");
  memcpy(&n, dataBuffer, sizeof(int));
  if(n <= BLOCK_SIZE) {
    // Lowest level
    if(bufferSize != sizeof(int) + n*n*sizeof(double))
      throw std::runtime_error("Wrong buffer size to CMatrix::assign_from_buffer.");
    elements.resize(n*n);
    memcpy(&elements[0], dataBuffer+sizeof(int), n*n*sizeof(double));
  }
  else {
    // Not lowest level
    if(bufferSize != sizeof(int) + 4*sizeof(cht::ChunkID))
      throw std::runtime_error("Wrong buffer size to CMatrix::assign_from_buffer.");
    memcpy(&children, dataBuffer+sizeof(int), 4*sizeof(cht::ChunkID));
  }
}
size_t CMatrix::memoryUsage() const {
  return getSize();
}
void CMatrix::getChildChunks(std::list<cht::ChunkID> & childChunkIDs) const {
  if(elements.size() != 0) {
    // Lowest level. Nothing to do in this case.
  }
  else {
    // Not lowest level
    for(int i = 0; i < 4; i++)
      childChunkIDs.push_back(children[i]);
  }
}
