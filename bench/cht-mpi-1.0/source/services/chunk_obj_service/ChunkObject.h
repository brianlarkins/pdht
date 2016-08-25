/* CHT-MPI, version 1.0. An MPI-based Chunks and Tasks library
 *                       implementation.
 * For copyright and license information, see below under "Copyright and
 * license".
 * 
 * Primary academic reference: 
 * Chunks and Tasks: A programming model for parallelization of dynamic
 * algorithms,
 * Emanuel H. Rubensson and Elias Rudberg,
 * Parallel Computing 00, 00 (2013),
 * <http://dx.doi.org/10.1016/j.parco.2013.09.006>
 * 
 * For further information about Chunks and Tasks, see
 * <http://www.chunks-and-tasks.org>.
 * 
 * === Copyright and license ===
 * 
 * Copyright (c) 2009-2014 Emanuel H. Rubensson and Elias Rudberg. All
 *                         rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 * - Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 * 
 * - Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclaimer listed
 *   in this license in the documentation and/or other materials provided
 *   with the distribution.
 * 
 * - Neither the name of the copyright holders nor the names of its
 *   contributors may be used to endorse or promote products derived from
 *   this software without specific prior written permission.
 * 
 * The copyright holders provide no reassurances that the source code
 * provided does not infringe any patent, copyright, or any other
 * intellectual property rights of third parties. The copyright holders
 * disclaim any liability to any recipient for claims brought against
 * recipient by any third party for infringement of that parties
 * intellectual property rights.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#ifndef CHUNK_OBJECT_HEADER
#define CHUNK_OBJECT_HEADER

#include <string>
#include <list>

namespace cht {
  
  struct BaseObj {
    // Perhaps this struct should lie outside chunk stuff directory
    static std::string get_class_id() {return "BaseObj";}
    //    virtual void dummy(){}
  };

// Note: the ChunkID struct needs to define operator< so that it can be
// used as key for std::map.
  struct ChunkID : public BaseObj {
    struct MajorOwnerInfo {
      int rank;
      size_t size;
    MajorOwnerInfo() : rank(-1), size(0) { }
    };
  int ownerRank;
  int creatorRank;
  int creatorLocalID;
  size_t size;
  int chunkTypeID;
  int childDepth; // 0 means no child chunks, 1 means one level of child chunks, etc.
  size_t totalChunkCount; // total number of chunks including children recursively, 1 if no children.
  size_t totalDeepSize; // size including all child chunks recursively
  static const int MAX_NO_OF_OWNERS_IN_LIST = 3;
  MajorOwnerInfo majorOwners[MAX_NO_OF_OWNERS_IN_LIST];
  ChunkID();
  bool operator<  ( const ChunkID & x ) const;
  bool operator==  ( const ChunkID & x ) const;
  bool operator!=  ( const ChunkID & x ) const;
  std::string str() const;
  static std::string get_class_id();
  };
  static const ChunkID CHUNK_ID_NULL;


 class Chunk : public BaseObj {
 protected:
  void deleteChunk(ChunkID cid) const;
 public:
  static std::string object_base_type_id() {
    return "Chunk";
  }
  // FIXME? Pass cht::vector to avoid communicating length?
  virtual void writeToBuffer( char * dataBuffer, 
			      size_t const bufferSize ) const = 0;
  virtual size_t getSize() const = 0;
  virtual void assignFromBuffer ( char const * dataBuffer, 
				  size_t const bufferSize) = 0;
  virtual size_t memoryUsage() const {return getSize();}
  // Subchunks: implemented via a (virtual) function getChildChunks() returning 
  // list of ID:s of subchunks
  virtual void getChildChunks(std::list<ChunkID> & childChunkIDs) const {}
  // Inherited classes are required to define
  // static std::string get_class_id();
  static std::string get_class_id() {return object_base_type_id();}
  
  // Note: a virtual destructor is needed here so that instances of 
  // chunk object classes derived from this base class can be deleted 
  // using a pointer to the base class.
  virtual ~Chunk() { }
};


} // end namespace cht

#endif
