/* Sample code using the Chunks and Tasks programming model.
 * See <http://www.chunks-and-tasks.org>.
 */
#ifndef CINT_HEADER
#define CINT_HEADER

#include "chunks_and_tasks.h"
struct CInt: public cht::Chunk {
  // Functions required for a Chunk
  void writeToBuffer(char * dataBuffer, size_t const bufferSize) const;
  size_t getSize() const;
  void assignFromBuffer(char const * dataBuffer, size_t const bufferSize);
  size_t memoryUsage() const;
  // CInt specific functionality
 CInt(int x_) : x(x_) { }
  CInt() { }
  operator int() const { return x; }
 private:
  int x; // The number itself
  CHT_CHUNK_TYPE_DECLARATION;
};

#endif
