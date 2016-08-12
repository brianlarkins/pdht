/* Sample code using the Chunks and Tasks programming model.
 * See <http://www.chunks-and-tasks.org>.
 */
#ifndef CDOUBLE_HEADER
#define CDOUBLE_HEADER

#include "chunks_and_tasks.h"
struct CDouble: public cht::Chunk {
  // Functions required for a Chunk
  void writeToBuffer(char * dataBuffer, size_t const bufferSize) const;
  size_t getSize() const;
  void assignFromBuffer(char const * dataBuffer, size_t const bufferSize);
  size_t memoryUsage() const;
  // CDouble specific functionality
 CDouble(double x_) : x(x_) { }
  CDouble() { }
  operator double() const { return x; }
 private:
  double x; // The number itself
  CHT_CHUNK_TYPE_DECLARATION;
};

#endif
