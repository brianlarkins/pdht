/* Sample code using the Chunks and Tasks programming model.
 * See <http://www.chunks-and-tasks.org>.
 */
#ifndef CMATRIX_HEADER
#define CMATRIX_HEADER

#include "chunks_and_tasks.h"
struct CMatrix: public cht::Chunk {
  // Functions required for a Chunk
  void writeToBuffer(char * dataBuffer, size_t const bufferSize) const;
  size_t getSize() const;
  void assignFromBuffer(char const * dataBuffer, size_t const bufferSize);
  void getChildChunks(std::list<cht::ChunkID> & childChunkIDs) const;
  size_t memoryUsage() const;
  // CMatrix specific functionality
  static const int BLOCK_SIZE = 200;
  CMatrix() { }
  int n; // matrix dimension
  std::vector<double> elements; // matrix elements, if lowest level
  cht::ChunkID children[4]; // 2x2 matrix of ids for child matrices, if not lowest level
  CHT_CHUNK_TYPE_DECLARATION;
};

#endif
